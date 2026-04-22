"""
학습 + 모델 관리 + 예측 라우트
"""
from flask import Blueprint, render_template, request, redirect, url_for, flash, jsonify
from app.models.database import Dataset, MLModel, TrainingRun, Prediction
from app.services.training_service import train_rf_model, get_model_list, delete_model
from app.services.prediction_service import predict_with_rf

training_bp = Blueprint("training", __name__)


# ─── 모델 목록 ───
@training_bp.route("/models")
def model_list():
    models = get_model_list()
    return render_template("training/models.html", models=models)


# ─── RF 학습 페이지 ───
@training_bp.route("/training/rf", methods=["GET", "POST"])
def train_rf():
    if request.method == "POST":
        dataset_id = request.form.get("dataset_id", type=int)
        model_name = request.form.get("model_name", "").strip()
        normalization_mode = request.form.get("normalization_mode", "rotation")
        ddcq = request.form.get("ddcq") == "true"

        if not dataset_id or not model_name:
            flash("데이터셋과 모델 이름을 입력해주세요.", "error")
            return redirect(url_for("training.train_rf"))

        result = train_rf_model(dataset_id, model_name, normalization_mode, ddcq)

        if result["success"]:
            flash("학습 완료! 정확도: " + str(result["accuracy"]) + " (" + str(result["sample_count"]) + "개 샘플)", "success")
            return redirect(url_for("training.model_list"))
        else:
            flash("학습 실패: " + result["message"], "error")
            return redirect(url_for("training.train_rf"))

    datasets = Dataset.query.order_by(Dataset.created_at.desc()).all()
    return render_template("training/train_rf.html", datasets=datasets)


# ─── 모델 삭제 ───
@training_bp.route("/models/<int:model_id>/delete", methods=["POST"])
def delete_model_route(model_id):
    success = delete_model(model_id)
    if success:
        flash("모델이 삭제되었습니다.", "success")
    else:
        flash("모델을 찾을 수 없습니다.", "error")
    return redirect(url_for("training.model_list"))


# ─── 예측 페이지 ───
@training_bp.route("/prediction", methods=["GET", "POST"])
def predict_page():
    if request.method == "POST":
        model_id = request.form.get("model_id", type=int)
        dataset_id = request.form.get("dataset_id", type=int)
        normalization_mode = request.form.get("normalization_mode", "rotation")

        if not model_id or not dataset_id:
            flash("모델과 데이터셋을 선택해주세요.", "error")
            return redirect(url_for("training.predict_page"))

        result = predict_with_rf(model_id, dataset_id, normalization_mode)

        if result["success"]:
            return render_template(
                "prediction/results.html",
                results=result["results"],
                summary=result["summary"],
            )
        else:
            flash("예측 실패: " + result["message"], "error")
            return redirect(url_for("training.predict_page"))

    datasets = Dataset.query.order_by(Dataset.created_at.desc()).all()
    models = MLModel.query.filter_by(model_type="rf").order_by(MLModel.created_at.desc()).all()
    return render_template("prediction/form.html", datasets=datasets, models=models)


# ─── 학습 이력 ───
@training_bp.route("/training/history")
def training_history():
    runs = TrainingRun.query.order_by(TrainingRun.started_at.desc()).limit(50).all()
    return render_template("training/history.html", runs=runs)