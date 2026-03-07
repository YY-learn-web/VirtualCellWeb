
import { PredictionRequest, PredictionResponse, GeneExpressionResult, BatchSweepRequest, BatchResponse, SweepPoint, IC50PredictionRequest, IC50Prediction, BatchIC50Request, BatchIC50Response, SweepIC50Request, SweepIC50Response, MetabolomicsAnalyzeParams, MetabolomicsTaskResponse, MetabolomicsTaskStatusResponse, MetabolomicsDiffPlotResponse, EnrichmentRequest, EnrichmentResponse } from '../types';

/**
 * API SERVICE - Connects to FastAPI Backend
 */

// FastAPI服务器地址
const API_BASE_URL = 'http://localhost:8080';

// 通用请求函数
const apiRequest = async (endpoint: string, data: any): Promise<any> => {
  try {
    const response = await fetch(`${API_BASE_URL}${endpoint}`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify(data)
    });

    if (!response.ok) {
      let detail = `HTTP error! status: ${response.status}`;
      try {
        const err = await response.json();
        if (err?.detail) {
          detail = `${detail} - ${err.detail}`;
        }
      } catch (e) {
        try {
          const text = await response.text();
          if (text) {
            detail = `${detail} - ${text}`;
          }
        } catch (inner) {
          // ignore
        }
      }
      throw new Error(detail);
    }

    return await response.json();
  } catch (error) {
    console.error('API请求失败:', error);
    throw error;
  }
};

// 转换后端数据格式到前端格�?
const transformBackendData = (backendData: any): GeneExpressionResult[] => {
  if (Array.isArray(backendData)) {
    return backendData.map(([geneId, expressionLevel]) => ({
      geneId,
      expressionLevel: parseFloat(expressionLevel)
    }));
  }
  return [];
};

export const ASCENDService = {
  predictExpression: async (request: PredictionRequest): Promise<PredictionResponse> => {
    console.log("Sending request to ASCEND Core...", request);

    // 构建后端请求格式
    const backendRequest = {
      AnalysisMode: 'single',
      smiles: request.smiles,
      time: request.time.toString(),
      dose: request.dose.toString(),
      cellLineMode: request.cellLineMode,
      cellLineId: request.cellLineId,
      customExpression: request.customExpression
    };

    // 发送请求到后端
    const response = await apiRequest('/api/Analysis', backendRequest);
    
    if (!response.success) {
      throw new Error(response.error || '预测失败');
    }

    // 转换数据格式
    const transformedData = transformBackendData(response.data);

    return {
      success: response.success,
      data: transformedData,
      metadata: {
        inferenceTime: response.metadata.inferenceTime,
        modelVersion: response.metadata.modelVersion
      }
    };
  },
  predictIC50: async (request: IC50PredictionRequest): Promise<{ success: boolean; data: IC50Prediction; }> => {
    const response = await apiRequest('/api/ic50/predict', request);
    return response;
  },
  predictBatchIC50: async (request: BatchIC50Request): Promise<BatchIC50Response> => {
    const response = await apiRequest('/api/ic50/batch/predict', request);
    return response;
  },
  predictSweepIC50: async (request: SweepIC50Request): Promise<SweepIC50Response> => {
    const response = await apiRequest('/api/ic50/sweep/predict', request);
    return response;
  },

  predictBatchSweep: async (request: BatchSweepRequest): Promise<BatchResponse> => {
    console.log("Sending Batch Sweep Request...", request);

    // 构建后端请求格式
    const backendRequest = {
      AnalysisMode: 'batch',
      BatchMethod: 'sweep',
      SweepVariable: request.sweepVariable,
      range: {
        start: request.range.start,
        end: request.range.end,
        steps: request.range.steps
      },
      fixedParamValue: request.fixedParamValue,
      smiles: request.smiles,
      cellLineId: request.cellLineId
    };

    // 发送请求到后端
    const response = await apiRequest('/api/Analysis', backendRequest);
    
    if (!response.success) {
      throw new Error(response.error || 'Batch scanning failed.');
    }

    return {
      success: response.success,
      type: response.type,
      data: response.data 
    };
  },

  processBatchUpload: async (file: File, cellLineId?: string): Promise<BatchResponse> => {
    console.log("Processing Batch Upload...", file.name);

    const formData = new FormData();
    formData.append('file', file);
    if (cellLineId) {
      formData.append('cellLineId', cellLineId);
    }

    const response = await fetch(`${API_BASE_URL}/api/analysis/batch/upload`, {
      method: 'POST',
      body: formData
    });

    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }

    const result = await response.json();
    
    if (!result.success) {
      throw new Error(result.error || 'File upload processing failed.');
    }
    
    return {
      success: result.success,
      type: result.type,
      data: result.data
    };
  },

  startMetabolomicsAnalysis: async (params: MetabolomicsAnalyzeParams): Promise<MetabolomicsTaskResponse> => {
    const formData = new FormData();
    formData.append('file', params.file);
    formData.append('epochs', params.epochs.toString());
    formData.append('imputation', params.imputation ? 'true' : 'false');
    formData.append('n_clusters', params.nClusters.toString());

    const response = await fetch(`${API_BASE_URL}/api/metabolomics/analyze`, {
      method: 'POST',
      body: formData
    });

    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }

    return await response.json();
  },

  getMetabolomicsStatus: async (taskId: string): Promise<MetabolomicsTaskStatusResponse> => {
    const response = await fetch(`${API_BASE_URL}/api/metabolomics/status/${taskId}`);
    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }
    return await response.json();
  },

  generateMetabolomicsDiffPlot: async (
    jobId: string,
    c1: number,
    c2: number
  ): Promise<MetabolomicsDiffPlotResponse> => {
    const response = await fetch(
      `${API_BASE_URL}/api/metabolomics/diff_plot?job_id=${encodeURIComponent(jobId)}&c1=${c1}&c2=${c2}`
    );
    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }
    return await response.json();
  },

  getMetabolomicsResultUrl: (jobId: string, filename: string) => {
    return `${API_BASE_URL}/api/metabolomics/results/${jobId}/${filename}`;
  },

  fetchBatchMetabolomicsCsv: async (fileId: string): Promise<Blob> => {
    const response = await fetch(`${API_BASE_URL}/api/analysis/batch/export/${fileId}`);
    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }
    return await response.blob();
  },

  fetchSweepMetabolomicsCsv: async (request: BatchSweepRequest): Promise<Blob> => {
    const payload = {
      sweepVariable: request.sweepVariable,
      range: request.range,
      fixedParamValue: request.fixedParamValue,
      smiles: request.smiles,
      cellLineId: request.cellLineId
    };
    const response = await fetch(`${API_BASE_URL}/api/analysis/batch/sweep/export`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json'
      },
      body: JSON.stringify(payload)
    });
    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }
    return await response.blob();
  },

  fetchBatchEnrichmentSummary: async (fileId: string): Promise<GeneExpressionResult[]> => {
    const response = await fetch(`${API_BASE_URL}/api/enrichment/batch/summary/${fileId}`);
    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }
    const result = await response.json();
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
      cellLineId: request.cellLineId
    };
    const response = await fetch(`${API_BASE_URL}/api/enrichment/sweep/summary`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json'
      },
      body: JSON.stringify(payload)
    });
    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }
    const result = await response.json();
    if (!result.success) {
      throw new Error(result.error || 'Failed to load sweep enrichment summary.');
    }
    return transformBackendData(result.data);
  },

  runEnrichment: async (request: EnrichmentRequest): Promise<EnrichmentResponse> => {
    const response = await apiRequest('/api/enrichment/run', request);
    if (!response.success) {
      throw new Error(response.error || 'Enrichment failed.');
    }
    return response;
  }
};
