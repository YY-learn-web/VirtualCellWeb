import os
import sys
import shutil
import traceback
import json
from urllib import request as urlrequest, parse as urlparse, error as urlerror
import math
from functools import lru_cache
import base64
# 确保当前目录在 Python 路径中
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
METABO_BACKEND = os.path.join(REPO_ROOT, 'MetabolomicsPred', 'backend')
if METABO_BACKEND not in sys.path:
    sys.path.append(METABO_BACKEND)


METABO_IMPORT_ERROR = None

try:
    from app.algorithms.scfea_core import run_scfea_training
    from app.algorithms.scfea_analysis import run_downstream_analysis, generate_custom_diff_plot
except Exception as exc:
    METABO_IMPORT_ERROR = str(exc)
    run_scfea_training = None
    run_downstream_analysis = None
    generate_custom_diff_plot = None
    print(f'Metabolomics modules not loaded: {exc}')

from ic50_tool import get_ic50_predictor, expression_to_vector
from fastapi import FastAPI, HTTPException, UploadFile, File, Form, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Union, Literal, Any
import uvicorn
import numpy as np
import pandas as pd
import uuid
import zipfile
import io
import time
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from onnx_model import ONNXModel,load_onnx_model



app = FastAPI(title="Virtual Cell Predictor")

class GeneExpressionPayload(BaseModel):
    geneId: str
    expressionLevel: float

class IC50Request(BaseModel):
    cellLineId: str
    smiles: Optional[str] = None
    expression: List[GeneExpressionPayload]
    metadata: Optional[Dict[str, Any]] = None

class BatchIC50Request(BaseModel):
    fileId: str = Field(..., min_length=1)
    cellLineId: Optional[str] = None
    previewCount: int = Field(5, ge=1, le=50)

class SweepRange(BaseModel):
    start: float
    end: float
    steps: int = Field(..., ge=1)

class SweepIC50Request(BaseModel):
    smiles: str
    cellLineId: str
    sweepVariable: Literal["time", "dose"]
    fixedParamValue: float
    range: SweepRange
    previewCount: int = Field(5, ge=1, le=50)


class EnrichmentRequest(BaseModel):
    geneIds: List[str]
    library: str = Field("GO_Biological_Process_2021", min_length=1)
    topN: int = Field(20, ge=1, le=200)


class StringNetworkRequest(BaseModel):
    genes: List[GeneExpressionPayload]
    minExpression: float = Field(1.5)
    species: int = Field(9606, ge=1)
    requiredScore: int = Field(400, ge=0, le=1000)
    networkType: Literal["functional", "physical"] = "functional"
    maxGenes: int = Field(200, ge=5, le=1000)

ic50_predictor = get_ic50_predictor()

@app.post("/api/ic50/predict")
async def predict_ic50(request: IC50Request):
    if not request.expression:
        raise HTTPException(status_code=400, detail="Expression payload is empty")
    try:
        vector = expression_to_vector(
            [gene.model_dump() for gene in request.expression]
        )
        result = ic50_predictor.predict(vector, request.model_dump())
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"IC50 prediction failed: {exc}")
    return {"success": True, "data": result}


@app.post("/api/ic50/batch/predict")
async def predict_ic50_batch(request: BatchIC50Request):
    file_id = request.fileId
    if not file_id:
        raise HTTPException(status_code=400, detail="Missing fileId")

    result_csv = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "batch/results",
        f"result_{file_id}.csv",
    )
    if not os.path.exists(result_csv):
        raise HTTPException(status_code=404, detail="Batch result file not found")

    try:
        df = pd.read_csv(result_csv)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Failed to read batch result file: {exc}")

    if df.empty:
        raise HTTPException(status_code=400, detail="Batch result file is empty")

    gene_ids = getattr(ic50_predictor, "_gene_ids", None)
    if not gene_ids:
        raise HTTPException(status_code=500, detail="IC50 model gene list not loaded")

    gene_cols = [f"GENE_{int(gene_id)}" for gene_id in gene_ids]
    missing_cols = [col for col in gene_cols if col not in df.columns]
    if missing_cols:
        raise HTTPException(status_code=400, detail=f"Missing gene columns: {missing_cols[:5]}...")

    expression_matrix = df[gene_cols].to_numpy(dtype=np.float64)
    results = []
    for idx, row in df.iterrows():
        smiles = row.get("smiles")
        time_val = row.get("time")
        dose_val = row.get("dose")
        cell_line = request.cellLineId or row.get("cellLineId") or "PC3"
        context = {
            "cellLineId": cell_line,
            "smiles": smiles,
            "metadata": {
                "timeHours": float(time_val),
                "doseUm": float(dose_val),
            },
        }
        prediction = ic50_predictor.predict(expression_matrix[idx], context)
        results.append({
            "smiles": str(smiles),
            "time": float(time_val),
            "dose": float(dose_val),
            "cellLineId": str(cell_line),
            "lnIc50": prediction.get("lnIc50"),
            "ic50": prediction.get("ic50"),
            "unit": prediction.get("unit", "uM"),
            "confidence": prediction.get("confidence"),
            "toolVersion": prediction.get("toolVersion"),
            "drugId": prediction.get("drugId"),
        })

    preview_count = min(request.previewCount, len(results))
    preview = results[:preview_count]

    try:
        csv_filename = f"ic50_{file_id}.csv"
        zip_filename = f"ic50_{file_id}.zip"
        output_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "batch/results",
            zip_filename,
        )
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        output_df = pd.DataFrame(results)
        csv_buffer = io.StringIO()
        output_df.to_csv(csv_buffer, index=False)
        csv_content = csv_buffer.getvalue()

        with zipfile.ZipFile(output_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            zipf.writestr(csv_filename, csv_content)

        download_url = f"http://localhost:8080/api/download/{zip_filename}"
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Failed to build IC50 batch output: {exc}")

    return {
        "success": True,
        "data": {
            "processedRows": len(results),
            "preview": preview,
            "downloadUrl": download_url,
        },
    }


@app.post("/api/ic50/sweep/predict")
async def predict_ic50_sweep(request: SweepIC50Request):
    start = float(request.range.start)
    end = float(request.range.end)
    steps = int(request.range.steps)
    if steps < 1:
        raise HTTPException(status_code=400, detail="Sweep steps must be >= 1")

    step_size = (end - start) / (steps - 1) if steps > 1 else 0.0
    sweep_points = []
    for i in range(steps):
        x_value = start + i * step_size
        if request.sweepVariable == "time":
            time_val = x_value
            dose_val = float(request.fixedParamValue)
        else:
            time_val = float(request.fixedParamValue)
            dose_val = x_value
        sweep_points.append((x_value, time_val, dose_val))

    model = load_onnx_model()
    ic50_cache: Dict[tuple, Dict[str, Any]] = {}
    results = []
    for x_value, time_val, dose_val in sweep_points:
        cache_key = (request.smiles, request.cellLineId)
        if cache_key not in ic50_cache:
            predictions = model.predict_single(
                request.smiles, float(time_val), float(dose_val), request.cellLineId
            )
            vector = expression_to_vector(
                [{"geneId": gid, "expressionLevel": val} for gid, val in predictions]
            )
            context = {
                "cellLineId": request.cellLineId,
                "smiles": request.smiles,
                "metadata": {"timeHours": float(time_val), "doseUm": float(dose_val)},
            }
            ic50_cache[cache_key] = ic50_predictor.predict(vector, context)

        prediction = ic50_cache[cache_key]
        results.append(
            {
                "smiles": request.smiles,
                "time": float(time_val),
                "dose": float(dose_val),
                "cellLineId": request.cellLineId,
                "lnIc50": prediction.get("lnIc50"),
                "ic50": prediction.get("ic50"),
                "unit": prediction.get("unit", "uM"),
                "confidence": prediction.get("confidence"),
                "toolVersion": prediction.get("toolVersion"),
                "drugId": prediction.get("drugId"),
                "xValue": round(float(x_value), 2),
            }
        )

    preview_count = min(request.previewCount, len(results))
    preview = results[:preview_count]

    try:
        file_id = uuid.uuid4().hex[:8]
        csv_filename = f"ic50_sweep_{file_id}.csv"
        zip_filename = f"ic50_sweep_{file_id}.zip"
        output_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "batch/results",
            zip_filename,
        )
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        output_df = pd.DataFrame(results)
        csv_buffer = io.StringIO()
        output_df.to_csv(csv_buffer, index=False)
        csv_content = csv_buffer.getvalue()

        with zipfile.ZipFile(output_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            zipf.writestr(csv_filename, csv_content)

        download_url = f"http://localhost:8080/api/download/{zip_filename}"
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Failed to build IC50 sweep output: {exc}")

    return {
        "success": True,
        "data": {
            "processedRows": len(results),
            "preview": preview,
            "downloadUrl": download_url,
        },
    }


# CORS配置
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

'''
    函数介绍：选择分析模式
    入参：
        AnalysisMode:分析模式
        data:分析数据
    出参：
        anser:预测结果
''' 

# 根据不同分析模式路由到不同的处理函数
@app.post("/api/Analysis")
async def Analysis(request: Dict[str, Any]):
    
    AnalysisMode = request.get("AnalysisMode")
    print(f"AnalysisMode: {AnalysisMode}")
    if AnalysisMode == 'single':
        anser = await SingleAnalysis(request)
    else:
        BatchMethod = request.get("BatchMethod")
        if BatchMethod == 'sweep':
            anser = await BatchSweep(request)
        else:
            anser = await BatchUpload(request)

    return anser
        
'''
    函数介绍：单个分析
    入参：
        data:分析数据
    出参：
        anser:预测结果
'''
# 单个分析接口
@app.post("/api/Analysis/Single")
async def SingleAnalysis(request: Dict[str, Any]):
    success = True
    start_time = time.time()
    
    # 1. 读取数据
    smiles = request.get("smiles")
    time1 = request.get("time")
    dose = request.get("dose")
    cellLineMode = request.get("cellLineMode")    
    cellLineId = request.get("cellLineId")

    print(f"Single Analysis - smiles: {smiles}, time: {time1}, dose: {dose}, cellLineMode: {cellLineMode}, cellLineId: {cellLineId}")

    # 2. 使用模型进行预测
    model = load_onnx_model()

    try:
        if(cellLineMode == "preset"):
            anser = model.predict_single(smiles, float(time1), float(dose), cellLineId)
        print(f"✅ ONNX模型预测成功，返回 {len(anser)} 个基因结果")
    except RuntimeError as e:
        if "模型尚未加载" in str(e):
            raise HTTPException(
                status_code=503,
                detail="模型尚未加载，无法进行预测。请检查模型加载状态。"
            )
        raise
    
    # 打印所有预测结果
    print(f"\n=== 基因表达预测结果 ===")
    print(f"SMILES: {smiles}")
    print(f"时间: {time1}h, 剂量: {dose}μM")
    print(f"{'基因ID':<10} {'表达值':<12} {'状态':<8}")
    print("-" * 35)
    for i, (gene_id, expression) in enumerate(anser):
        if expression > 0:
            status = "上调"
        elif expression < 0:
            status = "下调"
        else:
            status = "不变"
        print(f"{gene_id:<10} {expression:>10.4f} {status:<8}")
        
        # 只显示前20个，避免输出过长
        if i >= 19 and len(anser) > 20:
            print(f"... 还有 {len(anser) - 20} 个基因结果")
            break
    print("=" * 35 + "\n")

    # 3. 返回结果
    return {
        'success': success,
        'data': anser,
        'metadata': {            
            'inferenceTime': f"{(time.time() - start_time):.3f}s",
            'modelVersion': 'xCUDO-final'
        }
    }

'''
    函数介绍：批量分析
    入参：
        data:分析数据
    出参：
        anser:预测结果
'''
# 参数扫描接口
@app.post("/api/Analysis/Batch/Sweep")
async def BatchSweep(request: Dict[str, Any]):
    success = True
    print(f"完整请求内容: {request}")

    # 1. 读取数据
    BatchMethod = "sweep"
    SweepVariable = request.get("SweepVariable")
    # 获取范围
    range_config = request.get('range')
    
    # 从范围配置中获取参数
    start = range_config.get('start')
    end = range_config.get('end')
    steps = range_config.get('steps')
    fixedParamValue = request.get("fixedParamValue")
    smiles = request.get("smiles")
    cellLineId = request.get("cellLineId")
        
    print(f"Sweep Analysis - smiles: {smiles}, variable: {SweepVariable}, "
              f"range: {start}-{end}, steps: {steps}"
              f"fixedParamValue: {fixedParamValue}"
              f"cellLineId: {cellLineId}")
    
    res = []
    
    # 生成扫描点
    step_size = (end - start) / (steps - 1) if steps > 1 else 0
    
    for i in range(steps):
        x_value = start + i * step_size
        
        fixed_params = {'time': float(fixedParamValue)}  
        # 根据扫描变量设置参数
        if SweepVariable == 'time':
            current_time = x_value
            current_dose = float(fixed_params.get('dose', 10))
        else:            
            current_time = float(fixed_params.get('time', 24))
            current_dose = x_value
         
        # 进行预测
        model = load_onnx_model()

        try:
            predictions = model.predict_single(smiles, current_time, current_dose, cellLineId)
            
            # 转换为扫描格式（只取前5个基因用于显示）
            genes = {}
            for gene_id, expression in predictions[:5]:
                genes[gene_id] = expression
            
            res.append({
                "xValue": round(x_value, 2),
                "genes": genes
            })
            
            # 打印当前扫描点的详细信息
            print(f"  扫描点 {i+1}/{steps}: x={x_value:.2f}, time={current_time:.1f}h, dose={current_dose:.1f}μM")
            print(f"    前5个基因: {[(gid, f'{exp:.3f}') for gid, exp in predictions[:5]]}")
            
        except Exception as e:
            success = False
            print(f"❌ 扫描点 {x_value} 预测失败: {str(e)}")
    
    print(f"✅ 批量扫描完成，生成 {len(res)} 个数据点")

    # 3. 返回结果
    anser = {
        "success": success,
        "type": 'sweep',
        "data": res
    }
    return anser

# 文件上传接口
@app.post("/api/analysis/batch/upload")
async def BatchUpload(file: UploadFile = File(...), cellLineId: Optional[str] = Form(None)):
    success = True
    # 1. 读取数据
    BatchMethod = "upload"
    # 2.验证文件
    if not file.filename:
        raise HTTPException(status_code=400, detail="未选择文件")
    
    # 3.检查文件大小
    file_content = await file.read()
    max_size = 50 * 1024 * 1024  # MAX 50MB
    
    if len(file_content) == 0:
        raise HTTPException(status_code=400, detail="文件为空")
    
    if len(file_content) > max_size:
        raise HTTPException(status_code=413, detail=f"文件大小超过限制({max_size/1024/1024:.1f}MB)，当前文件大小: {len(file_content)/1024/1024:.1f}MB")
    
    # 4.重置文件指针
    await file.seek(0)
    
    # 5.FastAPI中文件上传需要不同的处理方式
    form = None
    try:
        if file.filename.endswith('.csv') :
            df = pd.read_csv(file.file)
            form = ".csv"
        elif file.filename.endswith('.xlsx'):
            # 对于Excel文件，需要先将内容保存到临时文件
            import tempfile
            with tempfile.NamedTemporaryFile(delete=False, suffix='.xlsx') as tmp:
                tmp.write(file_content)
                tmp_path = tmp.name
            
            try:
                df = pd.read_excel(tmp_path)
                form = ".xlsx"
            finally:
                # 清理临时文件
                os.unlink(tmp_path)
        else:
            raise HTTPException(status_code=400, detail="只支持CSV或Excel文件")
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"文件读取失败: {str(e)}")
    
    # 6.验证文件格式（必须三列）
    if len(df.columns) != 3:
        raise HTTPException(status_code=400, detail="文件必须包含恰好三列")
    
    # 7.验证列名
    expected_columns = {'smiles', 'dose', 'time'}
    if not set(df.columns).issuperset(expected_columns):
        raise HTTPException(status_code=400, detail=f"文件必须包含列名: {expected_columns}")
    
    # 8.预测
    file_id = str(uuid.uuid4())[:8]
    file_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"batch/upload", f"upload_{file_id}{form}")
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    with open(file_path, "wb") as f:
        f.write(file_content)
    
    print(f"文件已保存到: {file_path}")
        
    model = load_onnx_model()
    try:
        pred_df = model.predict_from_file(file_path, file_id, cellLineId)
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"预测失败: {str(e)}")
    
    try:
        csv_filename = f"result_{file_id}.csv"
        zip_filename = f"result_{file_id}.zip"
        output_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "batch/results", zip_filename)
            
        # 确保目录存在
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
            
        # 创建ZIP文件
        with zipfile.ZipFile(output_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            # 将DataFrame转换为CSV字符串
            csv_buffer = io.StringIO()
            pred_df.to_csv(csv_buffer, index=False)
            csv_content = csv_buffer.getvalue()
                
            # 将CSV内容添加到ZIP文件
            zipf.writestr(csv_filename, csv_content)
            
        # 生成下载URL
        download_url = f"http://localhost:8080/api/download/{zip_filename}"
                
        print(f"✅ ZIP文件已生成: {output_path}")
        print(f"🔗 下载链接: {download_url}")
                
    except Exception as e:
        print(f"❌ 生成zip文件失败: {str(e)}")
        raise HTTPException(status_code=500, detail=f"生成结果文件失败: {str(e)}")
    
    # 10.构建返回结果
    result = {
        'success': success,
        'type': 'upload',
        'data': {
            'processedRows': len(pred_df),
            'downloadUrl': download_url,
            'fileId': file_id,
        }
    }
    
    # 打印处理结果
    print(f"\n=== 处理结果 ===")
    print(f"处理状态: {'成功' if result['success'] else '失败'}")
    print(f"处理类型: {result['type']}")
    print(f"处理行数: {result['data']['processedRows']}")
    print(f"完整返回数据: {result}")
    print("================\n")
    
    return result

# 文件下载接口
# ?????????????
# ???????????????
@app.post("/api/analysis/batch/sweep/export")
async def export_sweep_for_metabolomics(request: Dict[str, Any]):
    sweep_variable = request.get("sweepVariable") or request.get("SweepVariable")
    range_config = request.get("range") or {}
    start = range_config.get("start")
    end = range_config.get("end")
    steps = range_config.get("steps")
    fixed_param = request.get("fixedParamValue")
    smiles = request.get("smiles")
    cell_line = request.get("cellLineId")

    if not smiles or start is None or end is None or steps is None or not sweep_variable:
        raise HTTPException(status_code=400, detail="Missing sweep parameters")

    steps = int(steps)
    if steps < 1:
        raise HTTPException(status_code=400, detail="Sweep steps must be >= 1")

    step_size = (float(end) - float(start)) / (steps - 1) if steps > 1 else 0.0
    model = load_onnx_model()

    rows = []
    for i in range(steps):
        x_value = float(start) + i * step_size
        if sweep_variable == "time":
            current_time = x_value
            current_dose = float(fixed_param)
        else:
            current_time = float(fixed_param)
            current_dose = x_value

        predictions = model.predict_single(smiles, current_time, current_dose, cell_line)
        row = {"cell_id": f"sweep_{i + 1}"}
        for gene_id, expression in predictions:
            row[gene_id] = float(expression)
        rows.append(row)

    export_df = pd.DataFrame(rows)
    file_id = uuid.uuid4().hex[:8]
    export_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "batch/results",
        f"metabolomics_sweep_{file_id}.csv",
    )
    os.makedirs(os.path.dirname(export_path), exist_ok=True)
    export_df.to_csv(export_path, index=False)

    return FileResponse(
        path=export_path,
        filename=f"metabolomics_sweep_{file_id}.csv",
        media_type="text/csv",
    )
@app.get("/api/analysis/batch/export/{file_id}")
async def export_batch_for_metabolomics(file_id: str):
    result_csv = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "batch/results",
        f"result_{file_id}.csv",
    )
    if not os.path.exists(result_csv):
        raise HTTPException(status_code=404, detail="Batch result file not found")

    try:
        df = pd.read_csv(result_csv)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Failed to read batch result file: {exc}")

    gene_cols = [col for col in df.columns if col.startswith("GENE_")]
    if not gene_cols:
        raise HTTPException(status_code=400, detail="No gene columns found in batch result")

    export_df = df[gene_cols].copy()
    export_df.insert(0, "cell_id", [f"cell_{idx + 1}" for idx in range(len(export_df))])

    export_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "batch/results",
        f"metabolomics_{file_id}.csv",
    )
    export_df.to_csv(export_path, index=False)

    return FileResponse(
        path=export_path,
        filename=f"metabolomics_{file_id}.csv",
        media_type="text/csv",
    )
@app.get("/api/download/{filename}")
async def download_file(filename: str):
    if not filename.endswith(".zip"):
        raise HTTPException(status_code=400, detail="仅支持 .zip 文件")
    
    file_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "batch/results",
        filename
    )
    
    if not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail="文件不存在")
    
    # 返回文件响应
    return FileResponse(
        path=file_path,
        filename=filename,
        media_type="application/zip",
    )

# --- Metabolomics (scFEA) integration ---
METABO_UPLOAD_DIR = os.path.join(METABO_BACKEND, "upload_files")
METABO_RESULTS_DIR = os.path.join(METABO_BACKEND, "results")
METABO_ASSETS_DIR = os.path.join(METABO_BACKEND, "app", "assets")

os.makedirs(METABO_UPLOAD_DIR, exist_ok=True)
os.makedirs(METABO_RESULTS_DIR, exist_ok=True)

app.mount("/api/metabolomics/results", StaticFiles(directory=METABO_RESULTS_DIR), name="metabolomics-results")

METABO_TASK_STORE = {}


def _ensure_metabolomics_ready():
    if run_scfea_training is None or run_downstream_analysis is None:
        detail = "Metabolomics backend not available"
        if METABO_IMPORT_ERROR:
            detail = f"{detail}: {METABO_IMPORT_ERROR}"
        raise HTTPException(status_code=503, detail=detail)


def process_metabolomics_task(job_id: str, input_file_path: str, params: dict):
    output_dir = os.path.join(METABO_RESULTS_DIR, job_id)

    try:
        METABO_TASK_STORE[job_id] = {
            "status": "TRAINING",
            "progress": 10,
            "message": "Training scFEA metabolomics model..."
        }

        epochs = params.get("epochs", 100)
        use_imputation = params.get("imputation", False)
        n_clusters = params.get("n_clusters", 4)

        flux_file, balance_file = run_scfea_training(
            input_file=input_file_path,
            output_dir=output_dir,
            assets_dir=METABO_ASSETS_DIR,
            sc_imputation=use_imputation,
            epochs=epochs
        )

        METABO_TASK_STORE[job_id] = {
            "status": "ANALYZING",
            "progress": 60,
            "message": "Running downstream metabolomics analysis..."
        }

        result_paths = run_downstream_analysis(
            flux_file=flux_file,
            balance_file=balance_file,
            output_dir=output_dir,
            assets_dir=METABO_ASSETS_DIR,
            n_clusters=n_clusters
        )

        METABO_TASK_STORE[job_id] = {
            "status": "SUCCESS",
            "progress": 100,
            "message": "Analysis complete",
            "result": result_paths
        }

    except Exception as exc:
        traceback.print_exc()
        METABO_TASK_STORE[job_id] = {
            "status": "FAILURE",
            "progress": 0,
            "message": "Analysis failed",
            "error": str(exc)
        }


@app.post("/api/metabolomics/analyze")
async def start_metabolomics_analysis(
    background_tasks: BackgroundTasks,
    file: UploadFile = File(...),
    epochs: int = Form(100),
    imputation: bool = Form(False),
    n_clusters: int = Form(4)
):
    _ensure_metabolomics_ready()
    job_id = str(uuid.uuid4())

    file_ext = os.path.splitext(file.filename)[1]
    saved_filename = f"{job_id}{file_ext}"
    file_path = os.path.join(METABO_UPLOAD_DIR, saved_filename)

    try:
        with open(file_path, "wb+") as buffer:
            shutil.copyfileobj(file.file, buffer)
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"File save failed: {exc}")

    METABO_TASK_STORE[job_id] = {
        "status": "PENDING",
        "progress": 0,
        "message": "Task queued"
    }

    background_tasks.add_task(
        process_metabolomics_task,
        job_id,
        file_path,
        {"epochs": epochs, "imputation": imputation, "n_clusters": n_clusters}
    )

    return {
        "job_id": job_id,
        "task_id": job_id,
        "message": "Metabolomics analysis started"
    }


@app.get("/api/metabolomics/status/{task_id}")
async def get_metabolomics_status(task_id: str):
    task_info = METABO_TASK_STORE.get(task_id)

    if not task_info:
        return {
            "task_id": task_id,
            "status": "UNKNOWN",
            "message": "Task not found"
        }

    return {
        "task_id": task_id,
        "status": task_info["status"],
        "progress": task_info.get("progress", 0),
        "message": task_info.get("message", ""),
        "result": task_info.get("result", None),
        "error": task_info.get("error", None)
    }


@app.get("/api/metabolomics/diff_plot")
async def get_metabolomics_diff_plot(job_id: str, c1: int, c2: int):
    _ensure_metabolomics_ready()
    output_dir = os.path.join(METABO_RESULTS_DIR, job_id)
    try:
        img_name, csv_name = generate_custom_diff_plot(output_dir, c1, c2)
        return {
            "status": "success",
            "image": img_name,
            "csv": csv_name
        }
    except Exception as exc:
        traceback.print_exc()
        raise HTTPException(status_code=400, detail=str(exc))


# --- Gene Enrichment (Enrichr) ---
ENRICHMENT_API_BASE = "https://maayanlab.cloud/Enrichr"
ENRICHMENT_LIB_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets", "enrichment")
ENRICHMENT_GO_LIBS = {
    "GO_Biological_Process_2021": "GO_Biological_Process_2021.gmt",
    "GO_Molecular_Function_2021": "GO_Molecular_Function_2021.gmt",
    "GO_Cellular_Component_2021": "GO_Cellular_Component_2021.gmt",
}
os.makedirs(ENRICHMENT_LIB_DIR, exist_ok=True)

STRING_API_BASE = "https://string-db.org/api"


def _load_entrez_symbol_map(assets_dir: str) -> Dict[str, str]:
    candidates = [
        "entrez_to_symbol.csv",
        "entrez_to_symbol.tsv",
        "geneid_to_symbol.csv",
        "geneid_to_symbol.tsv",
        "geneinfo_beta.txt",
        "gene_info.txt",
    ]
    for name in candidates:
        path = os.path.join(assets_dir, name)
        if not os.path.exists(path):
            continue
        sep = "\t" if name.endswith(".tsv") or name.endswith(".txt") else ","
        try:
            df = pd.read_csv(path, sep=sep)
        except Exception:
            continue
        if df.empty:
            continue
        cols = [c.lower() for c in df.columns]
        if "symbol" not in cols and "gene_symbol" not in cols:
            continue
        entrez_col = None
        for candidate in ["entrez_id", "gene_id", "geneid", "entrez"]:
            if candidate in cols:
                entrez_col = df.columns[cols.index(candidate)]
                break
        if entrez_col is None:
            entrez_col = df.columns[0]
        symbol_col = df.columns[cols.index("symbol")] if "symbol" in cols else df.columns[cols.index("gene_symbol")]
        mapping = (
            df[[entrez_col, symbol_col]]
            .dropna()
            .astype(str)
            .set_index(entrez_col)[symbol_col]
            .to_dict()
        )
        if mapping:
            return mapping
    return {}


def _normalize_enrichment_genes(gene_ids: List[str], assets_dir: str):
    mapping = _load_entrez_symbol_map(assets_dir)
    symbols: List[str] = []
    unmapped: List[str] = []
    for gene in gene_ids:
        if gene is None:
            continue
        raw = str(gene).strip()
        if not raw:
            continue
        key = raw[5:] if raw.startswith("GENE_") else raw
        symbol = None
        if mapping:
            symbol = mapping.get(key)
        if symbol:
            symbols.append(symbol)
        else:
            if any(ch.isalpha() for ch in key):
                symbols.append(key)
            else:
                unmapped.append(raw)
    seen = set()
    unique_symbols = []
    for sym in symbols:
        if sym not in seen:
            seen.add(sym)
            unique_symbols.append(sym)
    return unique_symbols, unmapped


def _map_gene_id_to_symbol(gene_id: str, mapping: Dict[str, str]) -> Optional[str]:
    if gene_id is None:
        return None
    raw = str(gene_id).strip()
    if not raw:
        return None
    key = raw[5:] if raw.startswith("GENE_") else raw
    symbol = mapping.get(key) if mapping else None
    if symbol:
        return symbol
    if any(ch.isalpha() for ch in key):
        return key
    return None


def _enrichr_add_list(genes: List[str]) -> str:
    payload = urlparse.urlencode({
        "list": "\n".join(genes),
        "description": "ASCEND transcriptomics enrichment"
    }).encode("utf-8")
    req = urlrequest.Request(f"{ENRICHMENT_API_BASE}/addList", data=payload, method="POST")
    req.add_header("Content-Type", "application/x-www-form-urlencoded")
    with urlrequest.urlopen(req, timeout=20) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    user_list_id = data.get("userListId")
    if not user_list_id:
        raise ValueError("Enrichr did not return a userListId")
    return str(user_list_id)


def _enrichr_fetch_results(user_list_id: str, library: str):
    query = urlparse.urlencode({"userListId": user_list_id, "backgroundType": library})
    with urlrequest.urlopen(f"{ENRICHMENT_API_BASE}/enrich?{query}", timeout=20) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    return data.get(library, [])


def _download_go_library(library: str, dest_path: str):
    url = f"{ENRICHMENT_API_BASE}/geneSetLibrary?mode=text&libraryName={urlparse.quote(library)}"
    with urlrequest.urlopen(url, timeout=30) as resp:
        content = resp.read()
    with open(dest_path, "wb") as f:
        f.write(content)


def _ensure_go_library(library: str) -> str:
    filename = ENRICHMENT_GO_LIBS.get(library)
    if not filename:
        raise ValueError(f"Offline GO library not supported: {library}")
    path = os.path.join(ENRICHMENT_LIB_DIR, filename)
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        _download_go_library(library, path)
    return path


@lru_cache(maxsize=8)
def _load_gmt_library(path: str) -> Dict[str, List[str]]:
    gene_sets = {}
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            term = parts[0]
            genes = [g for g in parts[2:] if g]
            if genes:
                gene_sets[term] = genes
    if not gene_sets:
        raise ValueError("GMT library is empty")
    return gene_sets


def _log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _logsumexp(log_values: List[float]) -> float:
    if not log_values:
        return float("-inf")
    max_log = max(log_values)
    if max_log == float("-inf"):
        return max_log
    total = sum(math.exp(v - max_log) for v in log_values)
    return max_log + math.log(total)


def _hypergeom_sf(k: int, M: int, n: int, N: int) -> float:
    # Survival function for hypergeometric: P[X >= k]
    if k <= 0:
        return 1.0
    max_k = min(n, N)
    if k > max_k:
        return 0.0
    log_denom = _log_comb(M, N)
    log_terms = []
    for i in range(k, max_k + 1):
        log_num = _log_comb(n, i) + _log_comb(M - n, N - i)
        log_terms.append(log_num - log_denom)
    log_p = _logsumexp(log_terms)
    if log_p == float("-inf"):
        return 0.0
    return min(max(math.exp(log_p), 0.0), 1.0)


def _bh_fdr(p_values: List[float]) -> List[float]:
    n = len(p_values)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: p_values[i])
    adj = [0.0] * n
    prev = 1.0
    for rank, idx in enumerate(order, start=1):
        p = p_values[idx]
        val = min(p * n / rank, 1.0)
        prev = min(prev, val)
        adj[idx] = prev
    # ensure monotonicity
    for i in range(n - 2, -1, -1):
        adj[i] = min(adj[i], adj[i + 1])
    return adj


def _run_go_enrichment(genes: List[str], library: str, top_n: int) -> List[Dict[str, Any]]:
    gmt_path = _ensure_go_library(library)
    gene_sets = _load_gmt_library(gmt_path)
    gene_universe = set()
    for gset in gene_sets.values():
        gene_universe.update(gset)
    if not gene_universe:
        raise ValueError("GO library gene universe is empty")

    query = [g for g in genes if g in gene_universe]
    if len(query) < 5:
        raise ValueError("Not enough mapped genes for GO enrichment after filtering to library")

    M = len(gene_universe)
    N = len(set(query))
    results = []
    for term, gset in gene_sets.items():
        gs = set(gset)
        overlap = list(set(query).intersection(gs))
        k = len(overlap)
        if k == 0:
            continue
        n = len(gs)
        p_val = _hypergeom_sf(k, M, n, N)
        expected = (N * n) / M
        combined = -math.log(max(p_val, 1e-300)) * (k / max(expected, 1e-9))
        results.append({
            "term": term,
            "pValue": p_val,
            "adjP": None,
            "combinedScore": combined,
            "genes": overlap
        })

    pvals = [r["pValue"] for r in results]
    adj = _bh_fdr(pvals)
    for i, val in enumerate(adj):
        results[i]["adjP"] = val

    results.sort(key=lambda r: r["adjP"] if r["adjP"] is not None else r["pValue"])
    return results[:top_n]


@app.post("/api/enrichment/run")
async def run_gene_enrichment(request: EnrichmentRequest):
    if not request.geneIds:
        raise HTTPException(status_code=400, detail="geneIds is required")

    genes, unmapped = _normalize_enrichment_genes(request.geneIds, METABO_ASSETS_DIR)
    if len(genes) < 5:
        raise HTTPException(
            status_code=400,
            detail="Not enough mapped genes for enrichment. Check gene ID mapping."
        )

    results = []
    if request.library in ENRICHMENT_GO_LIBS:
        try:
            results = _run_go_enrichment(genes, request.library, request.topN)
        except urlerror.URLError as exc:
            raise HTTPException(status_code=503, detail=f"GO library download failed: {exc}")
        except Exception as exc:
            raise HTTPException(status_code=500, detail=f"Offline GO enrichment failed: {exc}")
    else:
        try:
            user_list_id = _enrichr_add_list(genes)
            raw_results = _enrichr_fetch_results(user_list_id, request.library)
        except urlerror.URLError as exc:
            raise HTTPException(status_code=503, detail=f"Enrichment service unavailable: {exc}")
        except Exception as exc:
            raise HTTPException(status_code=500, detail=f"Enrichment failed: {exc}")

        for row in raw_results:
            if not row or len(row) < 6:
                continue
            term = row[1]
            p_value = float(row[2]) if row[2] is not None else 1.0
            combined_score = float(row[4]) if len(row) > 4 and row[4] is not None else None
            overlap_genes = row[5] if len(row) > 5 else ""
            adj_p = float(row[6]) if len(row) > 6 and row[6] is not None else None
            gene_list = [g for g in str(overlap_genes).replace(";", ",").split(",") if g]
            results.append({
                "term": term,
                "pValue": p_value,
                "adjP": adj_p,
                "combinedScore": combined_score,
                "genes": gene_list
            })

        results = results[:request.topN]
    return {
        "success": True,
        "library": request.library,
        "data": results,
        "unmapped": unmapped
    }


def _get_gene_columns(df: pd.DataFrame) -> List[str]:
    gene_cols = [col for col in df.columns if str(col).startswith("GENE_")]
    if gene_cols:
        return gene_cols
    meta_cols = {"smiles", "time", "dose", "cellLineId", "cell_id"}
    candidate_cols = [col for col in df.columns if str(col) not in meta_cols]
    numeric_cols = []
    for col in candidate_cols:
        series = pd.to_numeric(df[col], errors="coerce")
        if series.notna().any():
            numeric_cols.append(col)
    return numeric_cols


def _summarize_gene_matrix(df: pd.DataFrame) -> List[List[Union[str, float]]]:
    gene_cols = _get_gene_columns(df)
    if not gene_cols:
        raise ValueError("No gene columns found in batch result")
    gene_df = df[gene_cols].apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0)
    means = gene_df.mean(axis=0)
    return [[str(gene_id), float(means[gene_id])] for gene_id in gene_cols]


@app.get("/api/enrichment/batch/summary/{file_id}")
async def get_batch_enrichment_summary(file_id: str):
    result_csv = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "batch/results",
        f"result_{file_id}.csv",
    )
    if not os.path.exists(result_csv):
        raise HTTPException(status_code=404, detail="Batch result file not found")
    try:
        df = pd.read_csv(result_csv)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Failed to read batch result file: {exc}")
    if df.empty:
        raise HTTPException(status_code=400, detail="Batch result file is empty")
    try:
        data = _summarize_gene_matrix(df)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    return {
        "success": True,
        "data": data,
        "processedRows": len(df)
    }


@app.post("/api/enrichment/sweep/summary")
async def get_sweep_enrichment_summary(request: Dict[str, Any]):
    sweep_variable = request.get("sweepVariable") or request.get("SweepVariable")
    range_config = request.get("range") or {}
    start = range_config.get("start")
    end = range_config.get("end")
    steps = range_config.get("steps")
    fixed_param = request.get("fixedParamValue")
    smiles = request.get("smiles")
    cell_line = request.get("cellLineId")

    if not smiles or start is None or end is None or steps is None or not sweep_variable:
        raise HTTPException(status_code=400, detail="Missing sweep parameters")

    steps = int(steps)
    if steps < 1:
        raise HTTPException(status_code=400, detail="Sweep steps must be >= 1")

    step_size = (float(end) - float(start)) / (steps - 1) if steps > 1 else 0.0
    model = load_onnx_model()

    rows = []
    for i in range(steps):
        x_value = float(start) + i * step_size
        if sweep_variable == "time":
            current_time = x_value
            current_dose = float(fixed_param)
        else:
            current_time = float(fixed_param)
            current_dose = x_value

        predictions = model.predict_single(smiles, current_time, current_dose, cell_line)
        row = {"cell_id": f"sweep_{i + 1}"}
        for gene_id, expression in predictions:
            row[gene_id] = float(expression)
        rows.append(row)

    df = pd.DataFrame(rows)
    try:
        data = _summarize_gene_matrix(df)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    return {
        "success": True,
        "data": data,
        "processedRows": len(df)
    }


@app.post("/api/enrichment/string_network")
async def get_string_network(request: StringNetworkRequest):
    if not request.genes:
        raise HTTPException(status_code=400, detail="genes is required")

    mapping = _load_entrez_symbol_map(METABO_ASSETS_DIR)
    filtered = [g for g in request.genes if g.expressionLevel is not None and g.expressionLevel > request.minExpression]
    if not filtered:
        raise HTTPException(status_code=400, detail="No genes passed the expression filter")

    filtered.sort(key=lambda g: g.expressionLevel, reverse=True)
    filtered = filtered[: request.maxGenes]
    symbols = []
    for item in filtered:
        sym = _map_gene_id_to_symbol(item.geneId, mapping)
        if sym:
            symbols.append(sym)
    if not symbols:
        raise HTTPException(status_code=400, detail="No valid gene symbols found after mapping")

    identifiers = "\r".join(symbols)
    payload = urlparse.urlencode({
        "identifiers": identifiers,
        "species": str(request.species),
        "required_score": str(request.requiredScore),
        "network_type": request.networkType,
        "caller_identity": "virtual-cell-ai",
    }).encode("utf-8")

    try:
        svg_url = f"{STRING_API_BASE}/svg/network"
        req = urlrequest.Request(svg_url, data=payload, method="POST")
        req.add_header("Content-Type", "application/x-www-form-urlencoded")
        with urlrequest.urlopen(req, timeout=30) as resp:
            svg_data = resp.read()

        link_url = f"{STRING_API_BASE}/tsv/get_link"
        link_req = urlrequest.Request(link_url, data=payload, method="POST")
        link_req.add_header("Content-Type", "application/x-www-form-urlencoded")
        link = None
        with urlrequest.urlopen(link_req, timeout=30) as resp:
            text = resp.read().decode("utf-8", errors="ignore").strip()
            if text:
                # Skip header if present
                lines = [line for line in text.splitlines() if line.strip()]
                if len(lines) > 1:
                    parts = lines[1].split("\t")
                    if parts:
                        link = parts[0]

        image_b64 = base64.b64encode(svg_data).decode("utf-8")
        image_uri = f"data:image/svg+xml;base64,{image_b64}"

    except urlerror.HTTPError as exc:
        raise HTTPException(status_code=503, detail=f"STRING API error: {exc}")
    except urlerror.URLError as exc:
        raise HTTPException(status_code=503, detail=f"STRING API unavailable: {exc}")
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"STRING network failed: {exc}")

    return {
        "success": True,
        "filteredCount": len(filtered),
        "mappedCount": len(symbols),
        "image": image_uri,
        "link": link
    }

if __name__ == "__main__":
    print("Loading xCUDO ONNX model and FastAPI starting server ...")
    print("Please wait until server has fully started")
    print("Server will be available at: http://localhost:8080")
    print("Frontend can connect from: http://localhost:3000")
    
    # 加载ONNX模型
    try:
        load_onnx_model()
        print("✅ xCUDO ONNX model loaded successfully!")
    except Exception as e:
        print(f"❌ Failed to load xCUDO ONNX model: {str(e)}")
    
    # 启动服务器
    uvicorn.run(app, host="0.0.0.0", port=8080, reload=False)
