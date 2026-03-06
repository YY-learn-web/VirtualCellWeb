import os
from celery import states
from app.core.celery_app import celery_app
from app.algorithms.scfea_core import run_scfea_training
from app.algorithms.scfea_analysis import run_downstream_analysis

# 定义资源目录 (请确保这些目录存在)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ASSETS_DIR = os.path.join(BASE_DIR, "app", "assets")


@celery_app.task(bind=True)
def task_scfea_pipeline(self, job_id: str, input_file_path: str, params: dict):
    """
    Celery 异步任务：执行 scFEA 全流程
    """
    output_dir = os.path.join(BASE_DIR, "results", job_id)

    try:
        # 1. 状态更新：开始训练
        self.update_state(state='TRAINING', meta={'progress': 10, 'message': '正在进行深度学习代谢通量预测...'})

        epochs = params.get("epochs", 100)
        use_imputation = params.get("imputation", False)

        # 调用 Phase 1 的训练函数
        flux_file, balance_file = run_scfea_training(
            input_file=input_file_path,
            output_dir=output_dir,
            assets_dir=ASSETS_DIR,
            sc_imputation=use_imputation,
            epochs=epochs
        )

        # 2. 状态更新：开始分析
        self.update_state(state='ANALYZING', meta={'progress': 60, 'message': '正在进行下游生信分析与绘图...'})

        # 调用 Phase 1 的分析函数
        result_paths = run_downstream_analysis(
            flux_file=flux_file,
            balance_file=balance_file,
            output_dir=output_dir,
            assets_dir=ASSETS_DIR
        )

        # 3. 完成
        return {
            "status": "completed",
            "job_id": job_id,
            "results": result_paths  # 包含所有生成的图片和CSV路径
        }

    except Exception as e:
        # 错误处理
        import traceback
        error_msg = str(e)
        traceback.print_exc()
        self.update_state(state=states.FAILURE, meta={'error': error_msg})
        # 抛出异常以便 Celery 标记为 Failed
        raise e