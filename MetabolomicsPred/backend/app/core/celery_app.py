import os
from celery import Celery

# 如果使用 Docker，这里通常是 "redis://redis:6379/0"
# 本地开发通常是 "redis://localhost:6379/0"
REDIS_URL = os.getenv("REDIS_URL", "redis://localhost:6379/0")

celery_app = Celery(
    "scfea_worker",
    broker=REDIS_URL,
    backend=REDIS_URL
)

celery_app.conf.update(
    task_serializer="json",
    accept_content=["json"],
    result_serializer="json",
    timezone="Asia/Shanghai",
    enable_utc=True,
)