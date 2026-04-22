"""
분석 라우트
- 타겟별 곡선 그래프 페이지
- Baseline / Normalize 전처리 API
- Anomaly Detection API
"""
from flask import Blueprint, render_template, request, jsonify
from app.models.database import Dataset, Target
from app.services.analysis_service import get_graph_data, get_anomaly_data

analysis_bp = Blueprint("analysis", __name__)


@analysis_bp.route("/datasets/<int:dataset_id>/analysis")
def graph_page(dataset_id):
    """그래프 분석 페이지"""
    dataset = Dataset.query.get_or_404(dataset_id)
    targets = Target.query.filter_by(dataset_id=dataset_id).order_by(Target.target_name).all()
    selected_target = request.args.get("target", "")

    return render_template(
        "analysis/graph.html",
        dataset=dataset,
        targets=targets,
        selected_target=selected_target,
    )


@analysis_bp.route("/api/datasets/<int:dataset_id>/graph")
def graph_data_api(dataset_id):
    """곡선 데이터 JSON API"""
    target_name = request.args.get("target", "")
    preprocess = request.args.get("preprocess", "raw")

    if not target_name:
        return jsonify({"error": "target 파라미터가 필요합니다."}), 400

    data = get_graph_data(dataset_id, target_name, preprocess)
    if "error" in data:
        return jsonify(data), 404

    return jsonify(data)


@analysis_bp.route("/api/datasets/<int:dataset_id>/anomaly")
def anomaly_api(dataset_id):
    """Anomaly Detection JSON API"""
    target_name = request.args.get("target", "")
    threshold = request.args.get("threshold", 0.01, type=float)

    if not target_name:
        return jsonify({"error": "target 파라미터가 필요합니다."}), 400

    data = get_anomaly_data(dataset_id, target_name, threshold)
    if "error" in data:
        return jsonify(data), 400

    return jsonify(data)

@analysis_bp.route("/api/datasets/<int:dataset_id>/train-anomaly", methods=["POST"])
def train_anomaly_api(dataset_id):
    """Anomaly 모델 학습 API"""
    from app.services.analysis_service import train_anomaly_model

    data = request.get_json() or {}
    model_type = data.get("model_type", "gru")
    window_size = data.get("window_size", 3)
    epochs = data.get("epochs", 50)

    result = train_anomaly_model(dataset_id, model_type, window_size, epochs)
    return jsonify(result)

@analysis_bp.route("/datasets/<int:dataset_id>/clustering", methods=["GET", "POST"])
def clustering_page(dataset_id):
    """클러스터링 페이지"""
    dataset = Dataset.query.get_or_404(dataset_id)
    return render_template("analysis/clustering.html", dataset=dataset)


@analysis_bp.route("/api/datasets/<int:dataset_id>/clustering", methods=["POST"])
def clustering_api(dataset_id):
    """클러스터링 API"""
    from app.services.analysis_service import run_clustering

    data = request.get_json() or {}
    n_clusters = data.get("n_clusters", 3)
    method = data.get("method", "kmeans")
    early_weight = data.get("early_weight", 1.0)

    result = run_clustering(dataset_id, n_clusters, method, early_weight)
    return jsonify(result)