import os
import shutil
import uuid
import traceback
from fastapi import FastAPI, UploadFile, File, Form, HTTPException, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

# 引入第一阶段写好的核心算法
from app.algorithms.scfea_core import run_scfea_training
from app.algorithms.scfea_analysis import run_downstream_analysis
from app.algorithms.scfea_analysis import run_downstream_analysis, generate_custom_diff_plot
app = FastAPI(title="scFEA Web API (Local No-Redis Version)")

# ================= 配置区 =================
# 1. CORS 设置
origins = [
    "http://localhost",
    "http://localhost:3000",
    "http://localhost:5173",
    "http://127.0.0.1:5173",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# 2. 目录配置
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
UPLOAD_DIR = os.path.join(BASE_DIR, "upload_files")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
ASSETS_DIR = os.path.join(BASE_DIR, "app", "assets")

os.makedirs(UPLOAD_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# 3. 挂载静态文件
app.mount("/results", StaticFiles(directory=RESULTS_DIR), name="results")

# ================= 内存数据库 (替代 Redis) =================
# 在本地测试时，我们用一个全局字典存任务状态
# 格式: { "job_id": { "status": "...", "progress": 0, "message": "..." } }
LOCAL_TASK_STORE = {}


# ================= 后台处理逻辑 (替代 Worker) =================
def process_scfea_task(job_id: str, input_file_path: str, params: dict):
    """
    这是一个将在后台线程运行的函数
    """
    output_dir = os.path.join(RESULTS_DIR, job_id)

    try:
        # 1. 更新状态：开始训练
        LOCAL_TASK_STORE[job_id] = {
            "status": "TRAINING",
            "progress": 10,
            "message": "正在进行深度学习代谢通量预测..."
        }
        print(f"[{job_id}] 开始训练...")

        # 获取参数
        epochs = params.get("epochs", 100)
        use_imputation = params.get("imputation", False)
        n_clusters = params.get("n_clusters", 4)

        # === 调用核心训练代码 ===
        flux_file, balance_file = run_scfea_training(
            input_file=input_file_path,
            output_dir=output_dir,
            assets_dir=ASSETS_DIR,
            sc_imputation=use_imputation,
            epochs=epochs
        )

        # 2. 更新状态：开始分析
        LOCAL_TASK_STORE[job_id] = {
            "status": "ANALYZING",
            "progress": 60,
            "message": "正在进行下游生信分析与绘图..."
        }
        print(f"[{job_id}] 开始分析...")

        # === 调用下游分析代码 ===
        result_paths = run_downstream_analysis(
            flux_file=flux_file,
            balance_file=balance_file,
            output_dir=output_dir,
            assets_dir=ASSETS_DIR,
            n_clusters=n_clusters

        )

        # 3. 任务完成
        LOCAL_TASK_STORE[job_id] = {
            "status": "SUCCESS",
            "progress": 100,
            "message": "分析完成",
            "result": result_paths
        }
        print(f"[{job_id}] 任务完成！")

    except Exception as e:
        error_msg = str(e)
        traceback.print_exc()  # 在控制台打印报错详情
        LOCAL_TASK_STORE[job_id] = {
            "status": "FAILURE",
            "progress": 0,
            "message": "任务失败",
            "error": error_msg
        }


# ================= API 接口 =================

@app.get("/")
def read_root():
    return {"message": "scFEA Local Server Running"}


@app.post("/api/analyze")
async def start_analysis(
        background_tasks: BackgroundTasks,  # <--- FastAPI 核心黑魔法
        file: UploadFile = File(...),
        epochs: int = Form(100),
        imputation: bool = Form(False),
        n_clusters: int = Form(4)
):
    """
    接收文件，开启后台任务
    """
    job_id = str(uuid.uuid4())

    # 保存文件
    file_ext = os.path.splitext(file.filename)[1]
    saved_filename = f"{job_id}{file_ext}"
    file_path = os.path.join(UPLOAD_DIR, saved_filename)

    try:
        with open(file_path, "wb+") as buffer:
            shutil.copyfileobj(file.file, buffer)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"文件保存失败: {str(e)}")

    # 初始化任务状态
    LOCAL_TASK_STORE[job_id] = {
        "status": "PENDING",
        "progress": 0,
        "message": "任务已提交，等待开始..."
    }

    # 添加到后台任务 (这里不需要 Celery，FastAPI 会自己处理)
    background_tasks.add_task(
        process_scfea_task,
        job_id,
        file_path,
        {"epochs": epochs, "imputation": imputation, "n_clusters": n_clusters}
    )

    # 这里的 task_id 和 job_id 设为一样，方便前端逻辑
    return {
        "job_id": job_id,
        "task_id": job_id,
        "message": "Analysis started locally"
    }


@app.get("/api/status/{task_id}")
async def get_status(task_id: str):
    """
    查询任务状态 (从内存字典中读取)
    """
    task_info = LOCAL_TASK_STORE.get(task_id)

    if not task_info:
        # 如果找不到任务 (可能重启了服务器或者ID错)
        return {
            "task_id": task_id,
            "status": "UNKNOWN",
            "message": "任务不存在 (服务器重启后任务会清空)"
        }

    # 构造返回结构，保持和之前 Celery 版本一致，这样前端不用改
    return {
        "task_id": task_id,
        "status": task_info["status"],
        "progress": task_info.get("progress", 0),
        "message": task_info.get("message", ""),
        "result": task_info.get("result", None),
        "error": task_info.get("error", None)
    }

# === 3. 自定义火山图接口  ===
@app.get("/api/diff_plot")
async def get_custom_diff_plot(job_id: str, c1: int, c2: int):
    output_dir = os.path.join(RESULTS_DIR, job_id)
    try:
        img_name, csv_name = generate_custom_diff_plot(output_dir, c1, c2)
        return {
            "status": "success",
            "image": img_name,
            "csv": csv_name
        }
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=400, detail=str(e))