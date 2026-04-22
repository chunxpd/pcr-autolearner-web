"""
메인 페이지 라우트
- 대시보드: 데이터셋 수, 모델 수, 최근 활동 요약
"""
from flask import Blueprint, render_template
from app.models.database import Dataset, MLModel, AppLog

main_bp = Blueprint("main", __name__)


@main_bp.route("/")
def index():
    """대시보드 페이지"""
    dataset_count = Dataset.query.count()
    model_count = MLModel.query.count()
    recent_logs = AppLog.query.order_by(AppLog.created_at.desc()).limit(10).all()

    return render_template(
        "index.html",
        dataset_count=dataset_count,
        model_count=model_count,
        recent_logs=recent_logs,
    )