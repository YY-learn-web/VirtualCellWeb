import os
import sys
# 确保当前目录在 Python 路径中
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ic50_tool import get_ic50_predictor, expression_to_vector
from fastapi import FastAPI, HTTPException, UploadFile, File, Form
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
