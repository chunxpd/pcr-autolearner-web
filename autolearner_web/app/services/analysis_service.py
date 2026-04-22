"""
분석 서비스 레이어
- 기존 AutoLearner를 재사용하여 타겟별 곡선 데이터 생성
- Baseline Subtraction, Normalization 처리
- Anomaly Detection 지원
"""
import os
import copy
import numpy as np
from flask import current_app
from app.models.database import Dataset, DatasetFile, PreprocessingRun
from app.extensions import db


def _build_autolearner(dataset_id):
    """
    데이터셋의 정상 CSV 파일들을 AutoLearner에 로딩하여 반환한다.
    """
    from src.pcr_autolearner import AutoLearner

    learner = AutoLearner()
    normal_files = DatasetFile.query.filter_by(
        dataset_id=dataset_id, is_abnormal=False
    ).all()

    for df in normal_files:
        try:
            learner.read_input(df.file_path)
        except Exception as e:
            print(f"파일 읽기 실패: {df.filename} - {e}")

    return learner


# 다양한 색상 팔레트 (기존 Qt 버전처럼 각 곡선마다 다른 색)
COLOR_PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
    "#393b79", "#5254a3", "#6b6ecf", "#9c9ede", "#637939",
    "#8ca252", "#b5cf6b", "#cedb9c", "#8c6d31", "#bd9e39",
    "#e7ba52", "#e7cb94", "#843c39", "#ad494a", "#d6616b",
    "#e7969c", "#7b4173", "#a55194", "#ce6dbd", "#de9ed6",
]


def get_graph_data(dataset_id, target_name, preprocess="raw"):
    """
    특정 타겟의 곡선 데이터를 Plotly.js용 JSON으로 반환한다.
    """
    learner = _build_autolearner(dataset_id)

    try:
        learner.set_target(target_name)
    except ValueError as e:
        return {"error": str(e)}

    # 전처리 적용
    if preprocess == "baseline":
        learner.baseline_substraction(smoothing_size=3, mode="rotation")
        _save_preprocess_run(dataset_id, target_name, "baseline", "rotation")
    elif preprocess == "normalize":
        learner.normalization(smoothing_size=3, mode="rotation")
        _save_preprocess_run(dataset_id, target_name, "normalize", "rotation")

    # Plotly용 데이터 구성 — 각 곡선마다 고유 색상
    traces = []
    keys = sorted(learner.Cq_values.keys())
    for idx, key in enumerate(keys):
        pos_neg = learner.meta_info.get(key, {}).get("Pos/Neg", "")
        timestamp = learner.gen_timestamps.get(key, "")
        label = str(key) + "_" + timestamp if timestamp else str(key)
        color = COLOR_PALETTE[idx % len(COLOR_PALETTE)]

        traces.append({
            "key": key,
            "label": label,
            "values": learner.Cq_values[key],
            "pos_neg": pos_neg,
            "color": color,
        })

    return {
        "cycles": learner.cycles,
        "traces": traces,
        "target_name": target_name,
        "preprocess": preprocess,
        "sample_count": len(traces),
    }


def get_anomaly_data(dataset_id, target_name, threshold=0.01):
    """
    Anomaly Detection 실행.
    Raw 데이터에서 primitive normalization 후 GRU/TFLite 모델로 이상치 탐지.
    
    Returns:
        dict: {
            "anomalies": { "label": [cycle_indices], ... },
            "target_name": "...",
            "threshold": 0.01,
            "graph_data": { ... }  # raw 그래프 데이터 포함
        }
    """
    from src.pcr_autolearner import AutoLearner
    from src.pcr_analyzer import PCRmachine

    # 모델 파일 찾기
    model_path = _find_anomaly_model()
    if not model_path:
        return {"error": "Anomaly 모델 파일을 찾을 수 없습니다. (.tflite 또는 .keras)"}

    try:
        from src.pcr_anomaly_detector import PCRAnomalyDetector
        detector = PCRAnomalyDetector(seq_len=3)
        detector.load_model(model_path)
    except Exception as e:
        return {"error": "Anomaly 모델 로드 실패: " + str(e)}

    # 데이터 로딩
    learner = _build_autolearner(dataset_id)
    try:
        learner.set_target(target_name)
    except ValueError as e:
        return {"error": str(e)}

    # Raw 그래프 데이터도 함께 반환 (anomaly는 raw에서만 동작)
    graph_data = get_graph_data(dataset_id, target_name, preprocess="raw")

    # Anomaly detection 실행
    pcr_analyzer = PCRmachine()
    anomalies_dict = {}

    keys = sorted(learner.Cq_values.keys())
    for key in keys:
        timestamp = learner.gen_timestamps.get(key, "")
        label = str(key) + "_" + timestamp if timestamp else str(key)

        y = learner.Cq_values.get(key)
        x = learner.cycles
        if y is None:
            continue

        # pseudo PCR runner: primitive normalization
        pcr_analyzer.reset()
        for cq_value in y:
            pcr_analyzer.add_one_cq_value(cq_value, set_winsize=3, start=x[0])
        pcr_analyzer.pcr_machine_cal_dCq()
        adj_y = pcr_analyzer.get_cq_values()
        baseline_val = adj_y["Cq_smoothing"][0]
        adj_y_list = [value - baseline_val for value in adj_y["Cq_smoothing"]]
        if y[0] > 0:
            adj_y_list = [value / y[0] for value in adj_y_list]

        anomalies = detector.detect_anomaly(adj_y_list, threshold=threshold)
        if anomalies:
            anomalies_dict[label] = anomalies

    return {
        "anomalies": anomalies_dict,
        "target_name": target_name,
        "threshold": threshold,
        "graph_data": graph_data,
    }


def _find_anomaly_model():
    """Anomaly 모델 파일을 찾는다 (saved_models/anomaly/ 또는 기존 repo models/)
    .keras 파일을 우선 사용 (tflite는 Select TF Ops 문제가 있을 수 있음)
    """
    # 1) 웹 프로젝트의 saved_models/anomaly/
    anomaly_dir = current_app.config.get("ANOMALY_MODEL_FOLDER", "")
    if anomaly_dir and os.path.isdir(anomaly_dir):
        # keras 우선
        for fn in sorted(os.listdir(anomaly_dir), reverse=True):
            if fn.endswith(".keras"):
                return os.path.join(anomaly_dir, fn)
        for fn in sorted(os.listdir(anomaly_dir), reverse=True):
            if fn.endswith(".tflite"):
                return os.path.join(anomaly_dir, fn)

    # 2) 기존 repo의 models/ 디렉토리
    repo_models = os.path.join(current_app.config.get("REPO_ROOT", ""), "models")
    if os.path.isdir(repo_models):
        # keras 우선 (anomaly 관련)
        for fn in sorted(os.listdir(repo_models), reverse=True):
            if fn.endswith(".keras") and "anomaly" in fn.lower():
                return os.path.join(repo_models, fn)
        # 그 다음 tflite
        for fn in sorted(os.listdir(repo_models), reverse=True):
            if fn.endswith(".tflite"):
                return os.path.join(repo_models, fn)

    return None


def _save_preprocess_run(dataset_id, target_name, method, mode):
    """전처리 실행 이력을 DB에 저장"""
    run = PreprocessingRun(
        dataset_id=dataset_id,
        target_name=target_name,
        method=method,
        mode=mode,
        smoothing_size=3,
    )
    db.session.add(run)
    db.session.commit()
    
def train_anomaly_model(dataset_id, model_type="gru", window_size=3, epochs=50):
    """
    Anomaly Detection 모델을 학습하고 저장한다.
    
    Returns:
        dict: {"success": True/False, "message": "...", "model_path": "..."}
    """
    from src.pcr_autolearner import AutoLearner
    from src.pcr_autolearner_ai import load_data
    from src.pcr_anomaly_detector import PCRAnomalyDetector
    import numpy as np

    try:
        dataset = Dataset.query.get(dataset_id)
        if not dataset:
            return {"success": False, "message": "데이터셋을 찾을 수 없습니다."}

        data_dir = dataset.source_path
        if not data_dir or not os.path.isdir(data_dir):
            return {"success": False, "message": "데이터셋 경로가 유효하지 않습니다."}

        # 데이터 로딩
        runner = AutoLearner()
        input_rfus, input_labels, targets, source_files = load_data(
            runner, data_dir, normalization_mode="primitive", ddcq=False
        )

        if len(input_rfus) == 0:
            return {"success": False, "message": "학습할 데이터가 없습니다."}

        # 시퀀스 생성 (Positive 샘플만)
        X = []
        y = []
        for cq_values, label in zip(input_rfus, input_labels):
            if label == "Negative":
                continue
            for i in range(len(cq_values) - window_size):
                X.append(cq_values[i:i + window_size])
                y.append(cq_values[i + window_size])

        if not X:
            return {"success": False, "message": "학습 가능한 Positive 샘플이 없습니다."}

        X = np.array(X).reshape(-1, window_size, 1)
        y = np.array(y).reshape(-1, 1)

        # 모델 학습
        save_dir = current_app.config["ANOMALY_MODEL_FOLDER"]
        os.makedirs(save_dir, exist_ok=True)
        model_path = os.path.join(save_dir, "anomaly_detector.keras")

        detector = PCRAnomalyDetector(seq_len=window_size)
        detector.fit(X, y, epochs=epochs, model_path=model_path, model_type=model_type)

        return {
            "success": True,
            "message": "Anomaly 모델 학습 완료: " + model_path,
            "model_path": model_path,
        }

    except Exception as e:
        return {"success": False, "message": str(e)}    
    
def run_clustering(dataset_id, n_clusters=3, method="kmeans", early_weight=1.0):
    """
    전체 타겟 데이터를 클러스터링한다.
    기존 PCRCluster 클래스를 재사용.

    Returns:
        dict: {
            "success": True/False,
            "pca": {"x": [...], "y": [...], "labels": [...], "sample_labels": [...]},
            "profiles": [{"cluster": 0, "n": 312, "mean": [...], "curves": [[...], ...]}, ...],
            "n_clusters": 3,
            "method": "kmeans",
        }
    """
    from src.pcr_autolearner import AutoLearner
    from src.pcr_autolearner_ai import PCRCluster
    from sklearn.decomposition import PCA

    try:
        learner = _build_autolearner(dataset_id)
        targets = learner.get_target_names()

        if not targets:
            return {"success": False, "message": "타겟이 없습니다."}

        input_data = []
        input_labels = []

        for target in targets:
            learner.set_target(target)
            learner.normalization(smoothing_size=3, mode="primitive")
            keys = list(learner.meta_info.keys())
            input_data.extend([learner.Cq_values[key] for key in keys])
            for key in keys:
                timestamp = learner.gen_timestamps.get(key, "")
                label = target + "_" + timestamp if timestamp else target
                input_labels.append(label)

        input_data = np.array(input_data)
        if len(input_data) == 0:
            return {"success": False, "message": "클러스터링할 데이터가 없습니다."}

        # 클러스터링 실행
        clusterer = PCRCluster(
            n_clusters=n_clusters,
            method=method,
            early_weight=early_weight
        )

        normalized_data, valid_indices = clusterer.preprocess(input_data, skip_normalization=True)
        if len(normalized_data) == 0:
            return {"success": False, "message": "전처리 후 유효한 데이터가 없습니다."}

        filtered_labels = [input_labels[i] for i in valid_indices]
        cluster_labels = clusterer.fit_predict(normalized_data)

        # PCA 2D
        pca = PCA(n_components=2)
        reduced = pca.fit_transform(normalized_data)

        # Cluster Profiles
        profiles = []
        unique_labels = sorted(set(cluster_labels.tolist()))
        for cl in unique_labels:
            mask = cluster_labels == cl
            cluster_data = normalized_data[mask]
            mean_curve = np.mean(cluster_data, axis=0).tolist()

            # 개별 곡선은 최대 50개까지만 (브라우저 성능)
            sample_curves = cluster_data[:50].tolist()

            profiles.append({
                "cluster": int(cl),
                "n": int(mask.sum()),
                "mean": mean_curve,
                "curves": sample_curves,
            })

        return {
            "success": True,
            "pca": {
                "x": reduced[:, 0].tolist(),
                "y": reduced[:, 1].tolist(),
                "labels": cluster_labels.tolist(),
                "sample_labels": filtered_labels,
            },
            "profiles": profiles,
            "n_clusters": n_clusters,
            "method": method,
        }

    except Exception as e:
        return {"success": False, "message": str(e)}    