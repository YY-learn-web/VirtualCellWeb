import argparse
import json
import os.path
import time
import pandas as pd
import numpy as np
import random
from typing import Union, Optional, List
import torch
from torch_geometric.data import DataLoader
from tqdm import tqdm
from data_utils.dataset import InputDataset
from embeddings.esm import ESMProtein
from embeddings.morgan_fingerprint import MorganFingerprint
# 该文件用于获取单个smiles的特征向量

# 设置随机种子以确保结果一致性
def set_seed(seed: int = 42) -> None:
    """设置所有随机种子以确保结果可重现"""
    np.random.seed(seed)
    random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    # 确保CUDA操作是确定性的
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    # 设置环境变量
    os.environ["PYTHONHASHSEED"] = str(seed)
    print(f"随机种子已设置为: {seed}")

# 在导入时设置随机种子
set_seed(42)

# 全局变量存储模型
_net = None
_relation_map = None
_device = None
_initial_model = None
_model_loaded = False
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_MODEL_DIR = os.path.join(BASE_DIR, "models")

def init_model(model_path: str = DEFAULT_MODEL_DIR, 
               no_cuda: bool = False) -> bool:
    """
    初始化模型，加载模型文件和关系映射
    
    Args:
        model_path: 本地模型路径
        no_cuda: 是否禁用CUDA
        
    Returns:
        bool: 是否成功初始化
    """
    global _net, _relation_map, _device, _initial_model, _model_loaded
    
    try:
        # 设置设备
        _device = torch.device('cuda' if torch.cuda.is_available() and not no_cuda else 'cpu')
        
        # 初始化Morgan指纹模型
        try:
            from embeddings.morgan_fingerprint import MorganGenerator
            _initial_model = MorganGenerator()
        except ImportError:
            _initial_model = MorganFingerprint()
            print("⚠️ MorganFingerprint已弃用，建议使用MorganGenerator")
        
        # 加载模型
        model_file = os.path.join(model_path, "model.pt")
        relation_map_path = os.path.join(model_path, "relation_map.json")
        
        if not os.path.exists(model_file):
            raise FileNotFoundError(f"模型文件不存在: {model_file}")
        if not os.path.exists(relation_map_path):
            raise FileNotFoundError(f"关系映射文件不存在: {relation_map_path}")
            
        _net = torch.load(model_file, map_location=torch.device(_device))
        
        # 从JSON文件读取关系映射
        with open(relation_map_path) as f:
            _relation_map = json.load(f)
        
        # 设置模型为评估模式
        _net.eval()
        
        _model_loaded = True
        print(f"✅ 模型初始化成功，使用设备: {_device}")
        return True
    except Exception as e:
        print(f"❌ 模型初始化失败: {str(e)}")
        _model_loaded = False
        return False

def get_embedding(smiles: str, model_path: str = DEFAULT_MODEL_DIR, 
                            no_cuda: bool = False) -> Optional[List[float]]:
    """
    函数介绍:获取单个SMILES的特征向量
    
    入参:
        smiles: 小分子SMILES
        model_path: 本地模型路径（如果模型未初始化，将使用此路径初始化）
        no_cuda: 是否禁用CUDA
        
    出参:
        List[float]: 特征向量,如果失败返回None
    """
    global _net, _relation_map, _device, _initial_model, _model_loaded
    
    
    try:
        # 如果模型未加载，尝试初始化
        if not _model_loaded:
            if not init_model(model_path, no_cuda):
                return None
        
        # 获取初始嵌入
        initial_embeddings = _initial_model.get_embedding([smiles])
        
        # 设置模态和实体信息
        modality = 'morgan-fingerprint'
        entity_name = 'Drug'
        rel_id = _relation_map['smiles']
        
        # 处理每个嵌入
        output_embeddings = []
        for embeddings in initial_embeddings:
            # 创建节点用于评分
            nodes = {
                modality: {
                    'embeddings': embeddings.unsqueeze(0).to(_device),
                    'node_indices': torch.tensor([0]).to(_device)
                },
                entity_name: {
                    'embeddings': [None],
                    'node_indices': torch.tensor([1]).to(_device)
                }
            }

            triples = torch.tensor([[0], [rel_id], [1]]).to(_device)
            with torch.no_grad():
                node_output_embeddings = _net.encoder(nodes, triples)
                output_embeddings.append(node_output_embeddings[1])
        
        # 合并初始嵌入和输出嵌入
        # 使用与get_embeddings相同的处理方式
        final_embedding = []
        for output_emb, init_emb in zip(output_embeddings, initial_embeddings):
            # 转换为numpy并合并
            output_np = output_emb.cpu().numpy()
            init_np = init_emb.cpu().numpy()
            
            # 确保都是一维数组
            if output_np.ndim > 1:
                output_np = output_np.flatten()
            if init_np.ndim > 1:
                init_np = init_np.flatten()
                
            # 合并数组
            combined = np.concatenate([output_np, init_np])
            combined = combined[:128]
            final_embedding.extend(combined.tolist())
        
        return final_embedding
        
    except Exception as e:
        print(f"获取药物嵌入失败: {str(e)}")
        return None
