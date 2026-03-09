import os
import sys
import random
import numpy as np
import onnxruntime as ort
import pandas as pd
import joblib
from typing import List, Tuple, Optional
from tqdm import tqdm
import re

# 设置随机种子以确保结果一致性
def set_seed(seed: int = 42) -> None:
    """设置所有随机种子以确保结果可重现"""
    np.random.seed(seed)
    random.seed(seed)
    # 设置环境变量
    os.environ["PYTHONHASHSEED"] = str(seed)
    print(f"随机种子已设置为: {seed}")

# 在导入时设置随机种子
set_seed(42)

# 添加otter-knowledge模块路径
OTTER_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "otter-knowledge-main")
if OTTER_PATH not in sys.path:
    sys.path.append(OTTER_PATH)

# 尝试导入otter-knowledge模块
try:
    from inference import get_embedding
    OTTER_AVAILABLE = True
    print("otter-knowledge module loaded successfully.")
except ImportError as e:
    OTTER_AVAILABLE = False
    print(f"Failed to import otter-knowledge module: {str(e)}")

# 全局变量存储模型和数据
_onnx_session = None
_drug_dict = None
_gene_dict = None
_gene_ids = None
_time_scaler = None
_dose_scaler = None

# 文件路径配置
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_DIR = os.path.join(os.path.dirname(os.path.dirname(BASE_DIR)), "webData/model")
DATA_DIR = os.path.dirname(os.path.dirname(BASE_DIR))
ONNX_MODEL_PATH = os.path.join(DATA_DIR, "webData/xCUDO_final.onnx")
GENE_FILE = os.path.join(DATA_DIR, "webData/geneid2vec128_978.csv")
TIME_SCALER_PATH = os.path.join(MODEL_DIR, "time_scaler.pkl")
DOSE_SCALER_PATH = os.path.join(MODEL_DIR, "dose_scaler.pkl")


class ONNXModel:
    """ONNX模型包装类
    
    根据ONNX模型的实际结构：
    输入:
    - drug: [batch, 128] - 药物特征向量
    - pert_time: [batch] - 时间（小时）
    - cell_id: [batch, 978] - 细胞系特征向量
    - pert_idose: [batch] - 剂量（μM）
    
    输出:
    - prediction: [batch, 978] - 978个基因的表达值
    """
    
    def __init__(self, session, gene_dict, gene_ids, time_scaler, dose_scaler):
        self.session = session
        self.gene_dict = gene_dict
        self.gene_ids = gene_ids
        self.time_scaler = time_scaler
        self.dose_scaler = dose_scaler

    def get_drug_feature(self, smile):
        """根据药物smile获取特征向量"""
        
        # 使用otter-knowledge计算特征向量
        embedding = get_embedding(smile)
        if embedding is not None:
            print("Drug embedding computed with otter-knowledge.")
            """for i in range(len(embedding)):
                print(f"向量维度: {len(embedding)}, 数字: {embedding[i]}")"""
            return embedding
        raise ValueError(f"药物smile {smile} 无法计算")  
    
    def get_cell_feature(self, cell_id):
        """根据细胞系id获取特征向量"""
        cell_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "Sample_data",
            "cell.csv",
        )
        cell = pd.read_csv(cell_path, header=None)
        cell_id = re.sub(r'\([^)]*\)', '', cell_id)
        cell_vector = None
        for i in range(len(cell)):
            if cell.iloc[i, 0] == cell_id:
                cell_vector = cell.iloc[i, 1:].values  # 跳过第一列（ID），取后面的978个值
                break
        if cell_vector is None:
            raise ValueError(f"找不到细胞系ID: {cell_id}。可用的细胞系ID: {cell.iloc[:, 0].tolist()}")
        return cell_vector
     
    def predict_single(self, smiles: str, time: float, dose: float, cell_id) -> List[Tuple[str, float]]:
        """
        函数介绍：预测基因表达-单个smile
        入参:
            smiles: 小分子SMILES
            time: 时间
            dose: 剂量
            cell_id:细胞系id

        出参:
            List of (gene_id, expression_value) tuples
        """
        # 1. 获取药物特征（128维）
        drug_feature = self.get_drug_feature(smiles)
        if drug_feature is None:
            raise ValueError(f"无法找到SMILE '{smiles}' 对应的药物特征")
        
        # 2. 标准化时间和剂量
        time_scaled = self.time_scaler.transform([[time]])[0, 0]
        dose_log = np.log(dose + 1e-6)
        dose_scaled = self.dose_scaler.transform([[dose_log]])[0, 0]
        
        # 3. 获取细胞系特征向量
        cell_vector = self.get_cell_feature(cell_id)

        # 4. 准备ONNX模型输入
        # 确保drug_feature是numpy数组
        if isinstance(drug_feature, list):
            drug_feature = np.array(drug_feature)
        drug_input = drug_feature.astype(np.float64).reshape(1, -1)  # [1, 128]
        pert_time_input = np.array([[time_scaled]], dtype=np.float64)  # [1, 1] -> [1]
        cell_vec_input = cell_vector.astype(np.float64).reshape(1, -1)  # [1, 978]
        pert_idose_input = np.array([[dose_scaled]], dtype=np.float64)  # [1, 1] -> [1]
        # 5. 运行ONNX模型推理
        try:
            outputs = self.session.run(
                None,
                {
                    'drug': drug_input,
                    'pert_time': pert_time_input.flatten(),  
                    'cell_id': cell_vec_input,
                    'pert_idose': pert_idose_input.flatten() 
                }
            )
            predictions = outputs[0]  # [978]
            if predictions.ndim == 1:
                predictions = predictions[None, :]  # (978,) -> (1, 978)            
        except Exception as e:
            raise RuntimeError(f"ONNX模型推理失败: {str(e)}")
        
        # 6. 处理输出（转换为基因ID和表达值的列表）
        predictions_flat = predictions.flatten()
        cell_vector = np.asarray(cell_vector, dtype=np.float64)
        if cell_vector.shape[0] != predictions_flat.shape[0]:
            raise ValueError(
                f"Cell vector length {cell_vector.shape[0]} does not match prediction length {predictions_flat.shape[0]}"
            )
        delta_values = predictions_flat - cell_vector
        results = []
        for gene_id, pred_value in zip(self.gene_ids, delta_values):
            gene_id = "GENE_" + str(gene_id)
            results.append((gene_id, float(pred_value)))
        
        return results
    
    def predict_from_file(self, file_path: str, file_id: str, cell_id: str = "PC3") -> bool:
        """
        函数介绍:从文件读取smiles,dose,time,对每行进行预测
        
        输入文件包含3列
        - smiles: 小分子SMILES
        - time: 时间
        - dose: 剂量

        入参:
            file_path: 用户上传的文件的路径
            
        出参:
            results_df: 包含simles,dose,time,gene_id,expression_value的二维表格
        """
        try:
            # 1. 读取输入文件
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"输入文件不存在: {file_path}")

            if file_path.endswith('.csv'):
                df = pd.read_csv(file_path)
            elif file_path.endswith('.xlsx'):
                df = pd.read_excel(file_path)
 
            df_time = df['time']
            df_dose = df['dose']
            cell_id_input = cell_id or "PC3"
            
            all_results = []
            # 2. 对每一行进行预测
            for index, row in tqdm(df.iterrows(), total=len(df), desc="预测药物"):
                smiles = row['smiles']
                time_value = row.get('time', df_time)
                dose_value = row.get('dose', df_dose)
                cell_id = cell_id_input #模型输入需要 后续需调整
                
                try:
                    # 获取药物特征
                    drug_feature = self.get_drug_feature(smiles)
                    
                    # 标准化时间和剂量
                    time_scaled = self.time_scaler.transform([[time_value]])[0, 0]
                    dose_log = np.log(dose_value + 1e-6)
                    dose_scaled = self.dose_scaler.transform([[dose_log]])[0, 0]
                    
                    # 获取细胞系特征向量
                    cell_vector = self.get_cell_feature(cell_id)
             
                    # 准备ONNX模型输入
                    # 确保drug_feature是numpy数组
                    if isinstance(drug_feature, list):
                        drug_feature = np.array(drug_feature)
                    drug_input = drug_feature.astype(np.float64).reshape(1, -1)
                    pert_time_input = np.array([[time_scaled]], dtype=np.float64)
                    cell_vec_input = cell_vector.astype(np.float64).reshape(1, -1)
                    pert_idose_input = np.array([[dose_scaled]], dtype=np.float64)
                    
                    # 运行ONNX模型推理
                    outputs = self.session.run(
                        None,
                        {
                            'drug': drug_input,
                            'pert_time': pert_time_input.flatten(),
                            'cell_id': cell_vec_input,
                            'pert_idose': pert_idose_input.flatten()
                        }
                    )
                    predictions = outputs[0]
                    if predictions.ndim == 1:
                        predictions = predictions[None, :]
                    
                    # 处理输出
                    predictions_flat = predictions.flatten()
                    cell_vector = np.asarray(cell_vector, dtype=np.float64)
                    if cell_vector.shape[0] != predictions_flat.shape[0]:
                        raise ValueError(
                            f"Cell vector length {cell_vector.shape[0]} does not match prediction length {predictions_flat.shape[0]}"
                        )
                    delta_values = predictions_flat - cell_vector
                    result = {
                        'smiles': smiles,
                        'time': time_value,
                        'dose': dose_value,
                        'cellLineId': cell_id,
                    }
                    
                    # 添加基因表达值
                    for gene_id, pred_value in zip(self.gene_ids, delta_values):
                        result[f"GENE_{gene_id}"] = float(pred_value)
                    
                    all_results.append(result)
                    
                except Exception as e:
                    raise ValueError(f"处理SMILES {smiles} 时出错: {str(e)}")
            
            # 5. 生成输出文件名
            output_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "batch/results", f"result_{file_id}.csv")
            os.makedirs(os.path.dirname(output_path), exist_ok=True)

            # 6. 转换为DataFrame并保存
            results_df = pd.DataFrame(all_results)
            gene_cols = [col for col in results_df.columns if col.startswith("GENE_")]
            for col in gene_cols:
                col_mean = float(results_df[col].mean())
                col_std = float(results_df[col].std(ddof=0))
                if col_std == 0:
                    col_std = 1.0
                results_df[col] = (results_df[col] - col_mean) / col_std
            results_df.to_csv(output_path, index=False)
            print(f"成功保存预测结果到: {output_path}")
            
            return results_df
            
        except Exception as e:
            raise ValueError(f"批量预测失败: {str(e)}")
    

def load_onnx_model():
    """加载ONNX模型和相关数据"""
    global _onnx_session, _gene_dict, _gene_ids, _time_scaler, _dose_scaler
    
    print("正在加载ONNX模型...")
    
    # 1. 加载ONNX模型
    if not os.path.exists(ONNX_MODEL_PATH):
        raise FileNotFoundError(f"ONNX模型文件不存在: {ONNX_MODEL_PATH}")
    
    # 尝试使用不同的执行提供者
    available_providers = ort.get_available_providers()
    providers = [provider for provider in ['CoreMLExecutionProvider', 'CPUExecutionProvider'] if provider in available_providers]
    if not providers:
        providers = ['CPUExecutionProvider']
    session_options = ort.SessionOptions()
    
    try:
        _onnx_session = ort.InferenceSession(
            ONNX_MODEL_PATH,
            sess_options=session_options,
            providers=providers
        )
        print(f"ONNX model loaded successfully: {ONNX_MODEL_PATH}")
    except Exception as e:
        print(f"CoreML provider unavailable; retrying with CPU only: {str(e)}")
        try:
            _onnx_session = ort.InferenceSession(
                ONNX_MODEL_PATH,
                sess_options=session_options,
                providers=['CPUExecutionProvider']
            )
            print(f"ONNX model loaded successfully with CPU: {ONNX_MODEL_PATH}")
        except Exception as e2:
            raise RuntimeError(f"无法加载ONNX模型: {str(e2)}")
    
    # 3. 加载基因数据
    _gene_dict = pd.read_csv(GENE_FILE, header=None, index_col=0).to_dict('index')
    _gene_ids = list(_gene_dict.keys())
    print(f"Gene metadata loaded successfully: {len(_gene_dict)} genes")
    # 4. 加载标准化器
    _time_scaler = joblib.load(TIME_SCALER_PATH)
    _dose_scaler = joblib.load(DOSE_SCALER_PATH)

    # 此警告不影响结果 消除警告需要重新导出onnx模型
    # [W:onnxruntime:, execution_frame.cc:874 VerifyOutputSizes] Expected shape from model of {-1,978} does not match actual shape of {978} for output prediction

    return ONNXModel(_onnx_session, _gene_dict, _gene_ids, _time_scaler, _dose_scaler)

