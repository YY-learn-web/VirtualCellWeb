import axios from 'axios';

// 后端地址 (如果前后端在不同端口，请修改这里)
export const API_BASE_URL = 'http://localhost:8000';

// 定义接口类型
export interface AnalyzeParams {
  file: File;
  epochs: number;
  imputation: boolean;
}

export interface TaskResponse {
  job_id: string;
  task_id: string;
  message: string;
}

export interface TaskStatus {
  task_id: string;
  status: 'PENDING' | 'TRAINING' | 'ANALYZING' | 'SUCCESS' | 'FAILURE' | 'UNKNOWN';
  progress: number;
  message: string;
  result?: {
    files: Record<string, string>;
    images: {
        umap_cluster: string;
        umap_stress: string;
        heatmap_module: string;
        heatmap_pathway: string;
        volcano: string;
        stress_corr: string;
    };
  };
  error?: string;
}

// 1. 发起分析请求
export const startAnalysis = async (params: AnalyzeParams): Promise<TaskResponse> => {
  const formData = new FormData();
  formData.append('file', params.file);
  formData.append('epochs', params.epochs.toString());
  formData.append('imputation', params.imputation ? 'true' : 'false');

  const response = await axios.post(`${API_BASE_URL}/api/analyze`, formData, {
    headers: {
      'Content-Type': 'multipart/form-data',
    },
  });
  return response.data;
};

// 2. 查询任务状态
export const getTaskStatus = async (taskId: string): Promise<TaskStatus> => {
  const response = await axios.get(`${API_BASE_URL}/api/status/${taskId}`);
  return response.data;
};

// 3. 获取结果图片完整URL
export const getImageUrl = (jobId: string, imageName: string) => {
    // 后端返回的 images 字典里存的是文件名，我们需要拼上静态资源路径
    return `${API_BASE_URL}/results/${jobId}/${imageName}`;
}