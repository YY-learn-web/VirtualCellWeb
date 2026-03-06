import { CellLineOption } from './types';

// Default Cell Lines available in the system
export const PRESET_CELL_LINES: CellLineOption[] = [
  { id: 'PC3', name: 'PC3', tissue: 'Prostate', description: 'Prostate adenocarcinoma (Grade IV)' },
  { id: 'A375', name: 'A375', tissue: 'Skin', description: 'Malignant melanoma' },
  { id: 'MCF7', name: 'MCF7', tissue: 'Breast', description: 'Invasive breast ductal carcinoma' },
  { id: 'A549', name: 'A549', tissue: 'Lung', description: 'Lung carcinoma' },
  { id: 'HEPG2', name: 'HepG2', tissue: 'Liver', description: 'Hepatocellular carcinoma' },
  { id: 'HT29', name: 'HT-29', tissue: 'Colon', description: 'Colorectal adenocarcinoma' },
];

export const MODEL_VERSION = "ASCEND-v2.1-Beta";

// Constraints
export const MAX_TIME_HOURS = 360;
export const MAX_DOSE_UM = 30000;
export const DEFAULT_GENE_COUNT = 978; // L1000 landmark genes