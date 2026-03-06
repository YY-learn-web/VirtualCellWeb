import { PredictionRequest, PredictionResponse, GeneExpressionResult, BatchSweepRequest, BatchResponse, SweepPoint } from '../types';
import { DEFAULT_GENE_COUNT } from '../constants';

/**
 * MOCK API SERVICE
 */

// Helper to generate realistic-looking random data for the demo
const generateMockData = (geneCount: number, scalar: number = 1): GeneExpressionResult[] => {
  return Array.from({ length: geneCount }).map((_, i) => {
    // Generate a bell curve-ish distribution
    const u1 = Math.random();
    const u2 = Math.random();
    const z = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
    
    // Simulate some genes being up-regulated (positive) or down-regulated (negative)
    const expression = z * 2 * scalar; 

    return {
      geneId: `GENE_${i + 1}`,
      expressionLevel: parseFloat(expression.toFixed(3)),
    };
  });
};

export const ASCENDService = {
  predictExpression: async (request: PredictionRequest): Promise<PredictionResponse> => {
    console.log("Sending request to ASCEND Core...", request);

    // Simulate network latency (2 seconds)
    await new Promise(resolve => setTimeout(resolve, 2000));

    // Basic Validation Mock
    if (!request.smiles) {
      throw new Error("Invalid SMILES string provided.");
    }

    // Return Mock Data
    return {
      success: true,
      data: generateMockData(DEFAULT_GENE_COUNT),
      metadata: {
        inferenceTime: "0.42s",
        modelVersion: "ASCEND-v2.1"
      }
    };
  },

  predictBatchSweep: async (request: BatchSweepRequest): Promise<BatchResponse> => {
    console.log("Sending Batch Sweep Request...", request);
    await new Promise(resolve => setTimeout(resolve, 2500));

    const points: SweepPoint[] = [];
    const stepSize = (request.range.end - request.range.start) / (request.range.steps - 1);
    
    // Select 5 "Genes of Interest" to track across the sweep
    const trackedGenes = ['GENE_1', 'GENE_4', 'GENE_12', 'GENE_42', 'GENE_88'];

    for (let i = 0; i < request.range.steps; i++) {
      const xVal = request.range.start + (i * stepSize);
      
      // Create a fake curve effect based on the xVal
      // Normalized 0-1 factor
      const factor = (i + 1) / request.range.steps; 

      const geneData: { [key: string]: number } = {};
      
      trackedGenes.forEach((gene, idx) => {
        // Different curves for different genes (some go up, some go down)
        const direction = idx % 2 === 0 ? 1 : -1;
        const curve = Math.log(xVal + 1) * factor * 2 * direction;
        geneData[gene] = parseFloat(curve.toFixed(3));
      });

      points.push({
        xValue: Math.round(xVal * 100) / 100,
        genes: geneData
      });
    }

    return {
      success: true,
      type: 'sweep',
      data: points
    };
  },

  processBatchUpload: async (file: File): Promise<BatchResponse> => {
    console.log("Processing Batch Upload...", file.name);
    await new Promise(resolve => setTimeout(resolve, 3000));
    
    return {
      success: true,
      type: 'upload',
      data: {
        processedRows: 1240,
        downloadUrl: "#"
      }
    };
  }
};
