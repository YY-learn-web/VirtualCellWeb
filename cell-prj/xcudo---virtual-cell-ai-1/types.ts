
// Defines the input structure for the prediction model
export interface GeneExpressionPayload {
  geneId: string;
  expressionLevel: number;
}

export interface IC50PredictionRequest {
  cellLineId: string;
  smiles?: string;
  expression: GeneExpressionPayload[];
  metadata?: Record<string, unknown>;
}

export interface IC50Prediction {
  lnIc50: number;
  ic50: number;
  unit: string;
  confidence?: number;
  toolVersion?: string;
  extra?: Record<string, unknown>;
}


export interface PredictionRequest {
  smiles: string;
  time: number; // 0-360 hours
  dose: number; // 0-30000 uM
  cellLineMode: 'preset' | 'custom';
  cellLineId?: string; // For preset
  customExpression?: number[]; // For custom upload (simplified as array of floats)
}

// Defines a single gene's output data for visualization
export interface GeneExpressionResult {
  geneId: string;
  expressionLevel: number; // The predicted value
  controlLevel?: number;   // Optional: comparison to baseline
}

// Response from the backend
export interface PredictionResponse {
  success: boolean;
  data: GeneExpressionResult[];
  metadata: {
    inferenceTime: string;
    modelVersion: string;
  };
}

// --- BATCH MODE TYPES ---

export type AnalysisMode = 'single' | 'batch';
export type BatchMethod = 'sweep' | 'upload';
export type SweepVariable = 'time' | 'dose';

export interface BatchSweepRequest {
  smiles: string;
  cellLineId: string;
  fixedParamValue: number; // The value of the non-sweeping parameter
  sweepVariable: SweepVariable;
  range: {
    start: number;
    end: number;
    steps: number;
  };
}

export interface SweepPoint {
  xValue: number; // The time or dose value
  genes: { [geneId: string]: number }; // Key is gene ID, value is expression
}

export interface BatchResponse {
  success: boolean;
  type: BatchMethod;
  data: SweepPoint[] | BatchUploadData;
}

export interface BatchUploadData {
  processedRows: number;
  downloadUrl: string;
  fileId: string;
}

export interface BatchIC50Request {
  fileId: string;
  cellLineId?: string;
  previewCount?: number;
}

export interface BatchIC50PreviewItem {
  smiles: string;
  time: number;
  dose: number;
  cellLineId: string;
  lnIc50: number;
  ic50: number;
  unit: string;
  confidence?: number;
  toolVersion?: string;
  drugId?: string;
}

export interface BatchIC50Response {
  success: boolean;
  data: {
    processedRows: number;
    preview: BatchIC50PreviewItem[];
    downloadUrl: string;
  };
}

export interface SweepIC50Request {
  smiles: string;
  cellLineId: string;
  sweepVariable: SweepVariable;
  fixedParamValue: number;
  range: {
    start: number;
    end: number;
    steps: number;
  };
  previewCount?: number;
}

export interface SweepIC50Response {
  success: boolean;
  data: {
    processedRows: number;
    preview: BatchIC50PreviewItem[];
    downloadUrl: string;
  };
}

// Structure for cell line options
export interface CellLineOption {
  id: string;
  name: string;
  tissue: string;
  description: string;
}

export enum AppState {
  IDLE = 'IDLE',
  LOADING = 'LOADING',
  SUCCESS = 'SUCCESS',
  ERROR = 'ERROR'
}

// --- METABOLOMICS (scFEA) TYPES ---

export interface MetabolomicsTaskResponse {
  job_id: string;
  task_id: string;
  message: string;
}

export type MetabolomicsTaskStatus =
  | 'PENDING'
  | 'TRAINING'
  | 'ANALYZING'
  | 'SUCCESS'
  | 'FAILURE'
  | 'UNKNOWN';

export interface MetabolomicsTaskResult {
  files: Record<string, string>;
  images: {
    umap_cluster: string;
    umap_stress: string;
    heatmap_module: string;
    heatmap_pathway: string;
    volcano: string;
    stress_corr: string;
  };
}

export interface MetabolomicsTaskStatusResponse {
  task_id: string;
  status: MetabolomicsTaskStatus;
  progress: number;
  message: string;
  result?: MetabolomicsTaskResult;
  error?: string;
}

export interface MetabolomicsAnalyzeParams {
  file: File;
  epochs: number;
  imputation: boolean;
  nClusters: number;
}

export interface MetabolomicsDiffPlotResponse {
  status: 'success' | 'error';
  image: string;
  csv: string;
  error?: string;
}

// --- Gene Enrichment TYPES ---

export interface EnrichmentRequest {
  geneIds: string[];
  library: string;
  topN?: number;
}

export interface EnrichmentTerm {
  term: string;
  pValue: number;
  adjP?: number;
  combinedScore?: number;
  genes?: string[];
}

export interface EnrichmentResponse {
  success: boolean;
  library: string;
  data: EnrichmentTerm[];
  unmapped?: string[];
}
