<template>
  <div class="scfea-container">
    <el-card class="box-card">
      <template #header>
        <div class="card-header">
          <h2>🧬 scFEA 代谢通量分析平台</h2>
          <el-tag type="info">Deep Learning Flux Prediction</el-tag>
        </div>
      </template>

      <div v-if="taskStatus !== 'SUCCESS'" class="upload-section">
        <el-form :model="form" label-width="140px">

          <el-form-item label="训练轮数 (Epochs)">
            <el-input-number v-model="form.epochs" :min="10" :max="500" :step="10" />
            <span class="help-text">建议 100 轮以获得稳定结果 (本地测试可设为 10)</span>
          </el-form-item>

          <el-form-item label="聚类簇数 (Clusters)">
            <el-input-number v-model="form.nClusters" :min="2" :max="15" :step="1" />
            <span class="help-text">建议 2~8 之间，将用于下游的细胞分群与 UMAP 展示</span>
          </el-form-item>

          <el-form-item label="数据插补 (MAGIC)">
            <el-switch v-model="form.imputation" active-text="启用" inactive-text="关闭" />
            <span class="help-text">针对高度稀疏的单细胞数据推荐开启 (计算较慢)</span>
          </el-form-item>

          <el-upload
            class="upload-demo"
            drag
            action="#"
            :auto-upload="false"
            :limit="1"
            :on-change="handleFileChange"
            :on-remove="() => form.file = null"
          >
            <el-icon class="el-icon--upload"><upload-filled /></el-icon>
            <div class="el-upload__text">
              拖拽 CSV 文件到此处 或 <em>点击上传</em>
            </div>
            <template #tip>
              <div class="el-upload__tip">
                支持 .csv 格式单细胞矩阵 (行=基因, 列=细胞)
              </div>
            </template>
          </el-upload>
        </el-form>

        <div class="action-btn">
          <el-button
            type="primary"
            size="large"
            @click="submitAnalysis"
            :loading="isLoading || isPolling"
            :disabled="!form.file"
          >
            开始分析 (Start Analysis)
          </el-button>
        </div>
      </div>

      <div v-if="isPolling || taskStatus === 'FAILURE'" class="progress-section">
        <el-divider content-position="center">任务状态</el-divider>

        <div class="status-display">
          <el-progress
            type="circle"
            :percentage="progress"
            :status="progressStatus"
          />
          <div class="status-text">
            <h3>{{ currentMessage }}</h3>
            <p v-if="taskStatus === 'TRAINING'">正在训练神经网络，请耐心等待...</p>
            <p v-if="taskStatus === 'ANALYZING'">正在生成降维聚类图及通路热图...</p>
            <p v-if="taskStatus === 'FAILURE'" class="error-text">错误信息: {{ errorMessage }}</p>
          </div>
        </div>
      </div>

      <div v-if="taskStatus === 'SUCCESS' && resultData" class="result-section">
        <el-divider content-position="center">分析结果概览</el-divider>

        <el-tabs v-model="activeTab" class="demo-tabs">

          <el-tab-pane label="降维与聚类 (Dim Reduction)" name="dim">
            <el-row :gutter="20">
              <el-col :span="12">
                <el-card shadow="hover">
                  <template #header>UMAP - By Cluster (K={{ form.nClusters }})</template>
                  <el-image
                    :src="getImageUrl(currentJobId, resultData.images.umap_cluster)"
                    :preview-src-list="[getImageUrl(currentJobId, resultData.images.umap_cluster)]"
                    fit="contain"
                  />
                </el-card>
              </el-col>
              <el-col :span="12">
                <el-card shadow="hover">
                  <template #header>UMAP - By Metabolic Stress</template>
                  <el-image
                    :src="getImageUrl(currentJobId, resultData.images.umap_stress)"
                    :preview-src-list="[getImageUrl(currentJobId, resultData.images.umap_stress)]"
                    fit="contain"
                  />
                </el-card>
              </el-col>
            </el-row>
          </el-tab-pane>

          <el-tab-pane label="代谢通路 (Pathways)" name="pathway">
            <el-card shadow="hover" class="heatmap-card">
              <template #header>Cluster-wise Pathway Activity</template>
              <el-image
                :src="getImageUrl(currentJobId, resultData.images.heatmap_pathway)"
                :preview-src-list="[getImageUrl(currentJobId, resultData.images.heatmap_pathway)]"
                fit="contain"
              />
            </el-card>
            <div style="margin-top: 20px"></div>
            <el-card shadow="hover" class="heatmap-card">
               <template #header>Cluster-wise Module Activity</template>
               <el-image
                 :src="getImageUrl(currentJobId, resultData.images.heatmap_module)"
                 :preview-src-list="[getImageUrl(currentJobId, resultData.images.heatmap_module)]"
                 fit="contain"
               />
            </el-card>
          </el-tab-pane>

          <el-tab-pane label="差异分析 (Differential)" name="diff">
             <el-row :gutter="20">
              <el-col :span="12">
                <el-card shadow="hover">
                  <template #header>
                    <div style="display: flex; justify-content: space-between; align-items: center;">
                      <span>Volcano Plot (Diff Analysis)</span>
                      <div style="display: flex; gap: 5px; align-items: center;">
                        <el-select v-model="diffC1" size="small" style="width: 80px">
                          <el-option v-for="i in form.nClusters" :key="i-1" :label="`C${i-1}`" :value="i-1" />
                        </el-select>
                        <span style="font-size: 12px; font-weight: bold; color: #606266;">VS</span>
                        <el-select v-model="diffC2" size="small" style="width: 80px">
                          <el-option v-for="i in form.nClusters" :key="i-1" :label="`C${i-1}`" :value="i-1" />
                        </el-select>
                        <el-button type="primary" size="small" @click="updateVolcano" :loading="isUpdatingDiff">
                          重绘
                        </el-button>
                      </div>
                    </div>
                  </template>
                  <el-image
                    v-loading="isUpdatingDiff"
                    :src="currentVolcanoUrl"
                    :preview-src-list="[currentVolcanoUrl]"
                    fit="contain"
                  />
                </el-card>
              </el-col>

              <el-col :span="12">
                <el-card shadow="hover">
                  <template #header>Top Flux-Stress Correlation</template>
                  <el-image
                    :src="getImageUrl(currentJobId, resultData.images.stress_corr)"
                    :preview-src-list="[getImageUrl(currentJobId, resultData.images.stress_corr)]"
                    fit="contain"
                  />
                </el-card>
              </el-col>
            </el-row>
          </el-tab-pane>

        </el-tabs>

        <div class="download-area">
           <h3>数据下载</h3>
           <el-button-group>
             <el-button type="success" :icon="Download" @click="downloadFile('flux')">预测通量表 (Flux)</el-button>
             <el-button type="success" :icon="Download" @click="downloadFile('cluster')">细胞聚类表 (Cluster)</el-button>
             <el-button type="success" :icon="Download" @click="downloadFile('pca')">PCA 坐标表</el-button>
           </el-button-group>
           <el-button type="warning" @click="resetForm" style="margin-left: 20px">开始新任务</el-button>
        </div>

      </div>
    </el-card>
  </div>
</template>

<script setup lang="ts">
import { ref, computed, reactive } from 'vue';
import { UploadFilled, Download } from '@element-plus/icons-vue';
import { ElMessage } from 'element-plus';
import { startAnalysis, getTaskStatus, getImageUrl, getCustomDiffPlot, API_BASE_URL } from '@/api/scfea';
import type { UploadFile } from 'element-plus';

// --- State ---
const form = reactive({
  epochs: 100,
  imputation: false,
  nClusters: 4,  // 默认聚类数
  file: null as File | null
});

const isLoading = ref(false);
const isPolling = ref(false);
const currentJobId = ref('');
const progress = ref(0);
const currentMessage = ref('等待开始...');
const errorMessage = ref('');
const taskStatus = ref('PENDING'); // PENDING, TRAINING, ANALYZING, SUCCESS, FAILURE
const resultData = ref<any>(null);
const activeTab = ref('dim');

// 火山图自定义状态
const diffC1 = ref(0);
const diffC2 = ref(1);
const isUpdatingDiff = ref(false);
const currentVolcanoUrl = ref('');

// --- Computed ---
const progressStatus = computed(() => {
  if (taskStatus.value === 'FAILURE') return 'exception';
  if (taskStatus.value === 'SUCCESS') return 'success';
  return '';
});

// --- Methods ---

// 1. 文件选择处理
const handleFileChange = (uploadFile: UploadFile) => {
  if (uploadFile.raw) {
    form.file = uploadFile.raw;
  }
};

// 2. 提交主任务
const submitAnalysis = async () => {
  if (!form.file) {
    ElMessage.warning('请先上传 CSV 文件');
    return;
  }

  try {
    isLoading.value = true;
    taskStatus.value = 'PENDING';
    progress.value = 0;
    resultData.value = null;

    const res = await startAnalysis({
      file: form.file,
      epochs: form.epochs,
      imputation: form.imputation,
      n_clusters: form.nClusters
    });

    currentJobId.value = res.job_id;
    ElMessage.success('任务提交成功，开始分析...');

    startPolling(res.job_id);

  } catch (error: any) {
    console.error(error);
    ElMessage.error('提交失败: ' + (error.response?.data?.detail || error.message));
    isLoading.value = false;
  }
};

// 3. 轮询状态
const startPolling = (jobId: string) => {
  isLoading.value = false;
  isPolling.value = true;

  const timer = setInterval(async () => {
    try {
      const statusRes = await getTaskStatus(jobId);

      taskStatus.value = statusRes.status;
      progress.value = statusRes.progress;
      currentMessage.value = statusRes.message;

      if (statusRes.status === 'SUCCESS') {
        clearInterval(timer);
        isPolling.value = false;
        resultData.value = statusRes.result;

        // 设置初始火山图 URL (默认 C0 vs C1)
        currentVolcanoUrl.value = getImageUrl(jobId, statusRes.result.images.volcano);
        // 初始化下拉框选项
        diffC1.value = 0;
        diffC2.value = 1;

        ElMessage.success('分析完成！');
      } else if (statusRes.status === 'FAILURE') {
        clearInterval(timer);
        isPolling.value = false;
        errorMessage.value = statusRes.error || '未知错误';
        ElMessage.error('分析失败');
      }
    } catch (error) {
      console.error('Polling error', error);
      clearInterval(timer);
      isPolling.value = false;
      ElMessage.error('无法连接服务器查询状态');
    }
  }, 2000);
};

// 4. 重绘火山图
const updateVolcano = async () => {
  if (diffC1.value === diffC2.value) {
    ElMessage.warning('请选择两个不同的 Cluster 进行对比');
    return;
  }
  isUpdatingDiff.value = true;
  try {
    const res = await getCustomDiffPlot(currentJobId.value, diffC1.value, diffC2.value);
    // 强制刷新图片缓存 (加上时间戳后缀，防止浏览器使用旧缓存)
    currentVolcanoUrl.value = getImageUrl(currentJobId.value, res.image) + `?t=${new Date().getTime()}`;
    ElMessage.success(`生成 C${diffC1.value} vs C${diffC2.value} 对比图成功`);
  } catch (err: any) {
    ElMessage.error('重绘火山图失败: ' + (err.response?.data?.detail || err.message));
  } finally {
    isUpdatingDiff.value = false;
  }
};

// 5. 下载文件
const downloadFile = (fileKey: string) => {
  if (!resultData.value || !resultData.value.files[fileKey]) return;
  const originalPath = resultData.value.files[fileKey];
  const filename = originalPath.split(/[\\/]/).pop();
  const downloadUrl = `${API_BASE_URL}/results/${currentJobId.value}/${filename}`;
  window.open(downloadUrl, '_blank');
};

const resetForm = () => {
  taskStatus.value = 'PENDING';
  resultData.value = null;
  form.file = null;
};
</script>

<style scoped>
.scfea-container {
  max-width: 1200px;
  margin: 20px auto;
  padding: 20px;
}

.card-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
}

.upload-section {
  padding: 20px 0;
  max-width: 800px;
  margin: 0 auto;
}

.help-text {
  margin-left: 10px;
  color: #909399;
  font-size: 12px;
}

.action-btn {
  margin-top: 30px;
  text-align: center;
}

.progress-section {
  padding: 40px 0;
  text-align: center;
}

.status-display {
  display: flex;
  flex-direction: column;
  align-items: center;
  margin-top: 20px;
}

.status-text {
  margin-top: 20px;
}

.error-text {
  color: #f56c6c;
  font-weight: bold;
}

.result-section {
  margin-top: 20px;
}

.heatmap-card {
  margin-bottom: 20px;
}

.download-area {
  margin-top: 40px;
  text-align: center;
  padding: 20px;
  background-color: #f5f7fa;
  border-radius: 4px;
}
</style>