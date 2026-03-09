import {
  BatchIC50Request,
  BatchIC50Response,
  BatchMetabolomicsSnapshotResponse,
  BatchResponse,
  BatchSnapshotRequest,
  BatchSnapshotResponse,
  BatchSweepRequest,
  EnrichmentRequest,
  EnrichmentResponse,
  GeneExpressionResult,
  IC50Prediction,
  IC50PredictionRequest,
  KeyGeneRequest,
  KeyGeneResponse,
  MetabolomicsAnalyzeParams,
  MetabolomicsDiffPlotResponse,
  MetabolomicsTaskResponse,
  MetabolomicsTaskStatusResponse,
  PredictionRequest,
  PredictionResponse,
  SingleMetabolomicsRequest,
  SingleMetabolomicsResponse,
  StringNetworkRequest,
  StringNetworkResponse,
  SweepIC50Request,
  SweepIC50Response,
} from '../types';
import { buildApiUrl, mapNetworkError } from './apiConfig';

const apiFetch = async (path: string, init?: RequestInit): Promise<Response> => {
  try {
    return await fetch(buildApiUrl(path), init);
  } catch (error) {
    throw mapNetworkError(error);
  }
};

const readErrorDetail = async (response: Response): Promise<string> => {
  let detail = `HTTP error! status: ${response.status}`;

  try {
    const payload = await response.json();
    if (payload?.detail) {
      return `${detail} - ${payload.detail}`;
    }
    if (payload?.error) {
      return `${detail} - ${payload.error}`;
    }
  } catch {
    // Fall back to plain text below.
  }

  try {
    const text = await response.text();
    if (text) {
      return `${detail} - ${text}`;
    }
  } catch {
    // Ignore body parsing errors.
  }

  return detail;
};

const ensureOk = async (response: Response): Promise<Response> => {
  if (!response.ok) {
    throw new Error(await readErrorDetail(response));
  }
  return response;
};

const requestJson = async <T>(path: string, init?: RequestInit): Promise<T> => {
  const response = await ensureOk(await apiFetch(path, init));
  return await response.json() as T;
};

const postJson = async <T>(path: string, body: unknown): Promise<T> => {
  return await requestJson<T>(path, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
    },
    body: JSON.stringify(body),
  });
};

const requestBlob = async (path: string, init?: RequestInit): Promise<Blob> => {
  const response = await ensureOk(await apiFetch(path, init));
  return await response.blob();
};

const transformBackendData = (backendData: unknown): GeneExpressionResult[] => {
  if (!Array.isArray(backendData)) {
    return [];
  }

  return backendData.flatMap((item) => {
    if (!Array.isArray(item) || item.length < 2) {
      return [];
    }

    const [geneId, expressionLevel] = item;
    return [{
      geneId: String(geneId),
      expressionLevel: Number(expressionLevel),
    }];
  });
};

export const ASCENDService = {
  predictExpression: async (request: PredictionRequest): Promise<PredictionResponse> => {
    const backendRequest = {
      AnalysisMode: 'single',
      smiles: request.smiles,
      time: request.time.toString(),
      dose: request.dose.toString(),
      cellLineMode: request.cellLineMode,
      cellLineId: request.cellLineId,
      customExpression: request.customExpression,
    };

    const response = await postJson<{
      success: boolean;
      data: unknown;
      error?: string;
      metadata: {
        inferenceTime: string;
        modelVersion: string;
      };
    }>('/api/Analysis', backendRequest);

    if (!response.success) {
      throw new Error(response.error || 'Transcriptomics prediction failed.');
    }

    return {
      success: true,
      data: transformBackendData(response.data),
      metadata: response.metadata,
    };
  },

  predictIC50: async (request: IC50PredictionRequest): Promise<{ success: boolean; data: IC50Prediction; }> => {
    return await postJson<{ success: boolean; data: IC50Prediction; }>('/api/ic50/predict', request);
  },

  predictBatchIC50: async (request: BatchIC50Request): Promise<BatchIC50Response> => {
    return await postJson<BatchIC50Response>('/api/ic50/batch/predict', request);
  },

  predictSweepIC50: async (request: SweepIC50Request): Promise<SweepIC50Response> => {
    return await postJson<SweepIC50Response>('/api/ic50/sweep/predict', request);
  },

  predictBatchSweep: async (request: BatchSweepRequest): Promise<BatchResponse> => {
    const backendRequest = {
      AnalysisMode: 'batch',
      BatchMethod: 'sweep',
      SweepVariable: request.sweepVariable,
      range: {
        start: request.range.start,
        end: request.range.end,
        steps: request.range.steps,
      },
      fixedParamValue: request.fixedParamValue,
      smiles: request.smiles,
      cellLineId: request.cellLineId,
    };

    const response = await postJson<BatchResponse & { error?: string; }>('/api/Analysis', backendRequest);
    if (!response.success) {
      throw new Error(response.error || 'Batch sweep failed.');
    }
    return response;
  },

  processBatchUpload: async (file: File, cellLineId?: string): Promise<BatchResponse> => {
    const formData = new FormData();
    formData.append('file', file);
    if (cellLineId) {
      formData.append('cellLineId', cellLineId);
    }

    const result = await requestJson<BatchResponse & { error?: string; }>('/api/analysis/batch/upload', {
      method: 'POST',
      body: formData,
    });

    if (!result.success) {
      throw new Error(result.error || 'File upload processing failed.');
    }

    return result;
  },

  startMetabolomicsAnalysis: async (params: MetabolomicsAnalyzeParams): Promise<MetabolomicsTaskResponse> => {
    const formData = new FormData();
    formData.append('file', params.file);
    formData.append('epochs', params.epochs.toString());
    formData.append('imputation', params.imputation ? 'true' : 'false');
    formData.append('n_clusters', params.nClusters.toString());

    return await requestJson<MetabolomicsTaskResponse>('/api/metabolomics/analyze', {
      method: 'POST',
      body: formData,
    });
  },

  getMetabolomicsStatus: async (taskId: string): Promise<MetabolomicsTaskStatusResponse> => {
    return await requestJson<MetabolomicsTaskStatusResponse>(`/api/metabolomics/status/${taskId}`);
  },

  generateMetabolomicsDiffPlot: async (
    jobId: string,
    c1: number,
    c2: number,
  ): Promise<MetabolomicsDiffPlotResponse> => {
    return await requestJson<MetabolomicsDiffPlotResponse>(
      `/api/metabolomics/diff_plot?job_id=${encodeURIComponent(jobId)}&c1=${c1}&c2=${c2}`,
    );
  },

  getMetabolomicsResultUrl: (jobId: string, filename: string) => {
    return buildApiUrl(`/api/metabolomics/results/${jobId}/${encodeURIComponent(filename)}`);
  },

  fetchBatchMetabolomicsCsv: async (fileId: string): Promise<Blob> => {
    return await requestBlob(`/api/analysis/batch/export/${fileId}`);
  },

  fetchSweepMetabolomicsCsv: async (request: BatchSweepRequest): Promise<Blob> => {
    const payload = {
      sweepVariable: request.sweepVariable,
      range: request.range,
      fixedParamValue: request.fixedParamValue,
      smiles: request.smiles,
      cellLineId: request.cellLineId,
    };

    return await requestBlob('/api/analysis/batch/sweep/export', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify(payload),
    });
  },

  fetchBatchEnrichmentSummary: async (fileId: string): Promise<GeneExpressionResult[]> => {
    const result = await requestJson<{ success: boolean; data: unknown; error?: string; }>(
      `/api/enrichment/batch/summary/${fileId}`,
    );

    if (!result.success) {
      throw new Error(result.error || 'Failed to load batch enrichment summary.');
    }

    return transformBackendData(result.data);
  },

  fetchSweepEnrichmentSummary: async (request: BatchSweepRequest): Promise<GeneExpressionResult[]> => {
    const payload = {
      sweepVariable: request.sweepVariable,
      range: request.range,
      fixedParamValue: request.fixedParamValue,
      smiles: request.smiles,
      cellLineId: request.cellLineId,
    };

    const result = await postJson<{ success: boolean; data: unknown; error?: string; }>(
      '/api/enrichment/sweep/summary',
      payload,
    );

    if (!result.success) {
      throw new Error(result.error || 'Failed to load sweep enrichment summary.');
    }

    return transformBackendData(result.data);
  },

  fetchStringNetwork: async (request: StringNetworkRequest): Promise<StringNetworkResponse> => {
    const response = await postJson<StringNetworkResponse & { error?: string; }>(
      '/api/enrichment/string_network',
      request,
    );

    if (!response.success) {
      throw new Error(response.error || 'STRING network failed.');
    }

    return response;
  },

  getBatchSnapshot: async (request: BatchSnapshotRequest): Promise<BatchSnapshotResponse> => {
    const response = await postJson<BatchSnapshotResponse & { error?: string; }>(
      '/api/workflow/batch/snapshot',
      request,
    );

    if (!response.success) {
      throw new Error(response.error || 'Batch snapshot failed.');
    }

    return response;
  },

  getBatchMetabolomicsSnapshot: async (jobId: string, pointIndex: number): Promise<BatchMetabolomicsSnapshotResponse> => {
    const result = await requestJson<BatchMetabolomicsSnapshotResponse & { error?: string; }>(
      `/api/workflow/batch/metabolomics_snapshot?job_id=${encodeURIComponent(jobId)}&point_index=${pointIndex}`,
    );

    if (!result.success) {
      throw new Error(result.error || 'Batch metabolomics snapshot failed.');
    }

    return result;
  },

  runSingleWorkflowMetabolomics: async (request: SingleMetabolomicsRequest): Promise<SingleMetabolomicsResponse> => {
    const payload = {
      expression: request.expression.map(({ geneId, expressionLevel }) => ({ geneId, expressionLevel })),
      epochs: request.epochs ?? 100,
      imputation: request.imputation ?? false,
    };

    const response = await postJson<SingleMetabolomicsResponse & { error?: string; }>(
      '/api/workflow/single/metabolomics',
      payload,
    );

    if (!response.success) {
      throw new Error(response.error || 'Single metabolomics workflow failed.');
    }

    return response;
  },

  identifyKeyGenes: async (request: KeyGeneRequest): Promise<KeyGeneResponse> => {
    const response = await postJson<KeyGeneResponse & { error?: string; }>(
      '/api/workflow/single/key_genes',
      request,
    );

    if (!response.success) {
      throw new Error(response.error || 'Key gene analysis failed.');
    }

    return response;
  },

  runEnrichment: async (request: EnrichmentRequest): Promise<EnrichmentResponse> => {
    const response = await postJson<EnrichmentResponse & { error?: string; }>(
      '/api/enrichment/run',
      request,
    );

    if (!response.success) {
      throw new Error(response.error || 'Enrichment failed.');
    }

    return response;
  },
};
