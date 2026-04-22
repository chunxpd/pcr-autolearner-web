"""
학습 서비스 레이어
- RF 모델 학습 (기존 PCRPredicter 재사용)
- 모델 저장/목록/삭제
"""
import os
import json
from datetime import datetime
from flask import current_app

from app.extensions import db
from app.models.database import Dataset, DatasetFile, MLModel, TrainingRun, AppLog


def _log(level, category, message):
    log_entry = AppLog(level=level, category=category, message=message)
    db.session.add(log_entry)
    db.session.commit()


def train_rf_model(dataset_id, model_name, normalization_mode="rotation", ddcq=False):
    """
    RF Pos/Neg 모델을 학습하고 저장한다.

    Args:
        dataset_id: 데이터셋 ID
        model_name: 모델 표시 이름
        normalization_mode: 'rotation' 또는 'primitive'
        ddcq: ddCq 적용 여부

    Returns:
        dict: {"success": True/False, "message": "...", "model_id": N, "accuracy": 0.xx}
    """
    from src.pcr_autolearner import AutoLearner
    from src.pcr_autolearner_ai import PCRPredicter, load_data

    # 1) 학습 이력 시작
    training_run = TrainingRun(
        dataset_id=dataset_id,
        status="running",
        params_json=json.dumps({
            "model_type": "rf",
            "normalization_mode": normalization_mode,
            "ddcq": ddcq,
        }),
        started_at=datetime.utcnow(),
    )
    db.session.add(training_run)
    db.session.commit()

    try:
        # 2) 데이터 로딩
        dataset = Dataset.query.get(dataset_id)
        if not dataset:
            raise ValueError("데이터셋을 찾을 수 없습니다.")

        # 데이터셋 경로 결정
        data_dir = dataset.source_path
        if not data_dir or not os.path.isdir(data_dir):
            raise ValueError("데이터셋 경로가 유효하지 않습니다: " + str(data_dir))

        runner = AutoLearner()
        input_data, input_labels, targets, source_files = load_data(
            runner, data_dir,
            normalization_mode=normalization_mode,
            ddcq=ddcq,
        )

        if len(input_data) == 0:
            raise ValueError("학습할 데이터가 없습니다.")

        # 3) RF 학습
        save_dir = current_app.config["RF_MODEL_FOLDER"]
        os.makedirs(save_dir, exist_ok=True)

        pcr_learner = PCRPredicter()
        model = pcr_learner.train_model(
            input_data,
            input_labels,
            save_dir=save_dir,
            ai_model="rf",
        )

        # 4) 정확도 계산 (train_model 내부에서 출력되지만 여기서도 계산)
        import numpy as np
        from sklearn.metrics import accuracy_score
        processed_data, processed_labels = pcr_learner.preprocess_data(input_data, input_labels)
        predictions = model.predict(processed_data)
        accuracy = accuracy_score(processed_labels, predictions)

        # 5) 모델 파일 경로
        model_file = os.path.join(save_dir, "best_model_rf.pkl")

        # 타임스탬프 붙인 파일명으로 복사 (덮어쓰기 방지)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        final_filename = "rf_" + timestamp + ".pkl"
        final_path = os.path.join(save_dir, final_filename)

        import shutil
        shutil.copy2(model_file, final_path)

        # 6) DB에 모델 등록
        ml_model = MLModel(
            name=model_name,
            model_type="rf",
            file_path=final_path,
            accuracy=accuracy,
            dataset_id=dataset_id,
            params_json=json.dumps({
                "normalization_mode": normalization_mode,
                "ddcq": ddcq,
                "sample_count": len(input_data),
                "targets": targets,
            }),
        )
        db.session.add(ml_model)
        db.session.flush()

        # 7) 학습 이력 업데이트
        training_run.model_id = ml_model.id
        training_run.accuracy = accuracy
        training_run.sample_count = len(input_data)
        training_run.status = "success"
        training_run.finished_at = datetime.utcnow()
        db.session.commit()

        _log("INFO", "train", "RF 모델 학습 완료: " + model_name + " (정확도: " + str(round(accuracy, 4)) + ")")

        return {
            "success": True,
            "message": "학습 완료",
            "model_id": ml_model.id,
            "accuracy": round(accuracy, 4),
            "sample_count": len(input_data),
        }

    except Exception as e:
        training_run.status = "failed"
        training_run.error_message = str(e)
        training_run.finished_at = datetime.utcnow()
        db.session.commit()

        _log("ERROR", "train", "RF 학습 실패: " + str(e))

        return {
            "success": False,
            "message": str(e),
        }


def get_model_list():
    """저장된 모델 목록 반환"""
    return MLModel.query.order_by(MLModel.created_at.desc()).all()


def delete_model(model_id):
    """모델 삭제"""
    model = MLModel.query.get(model_id)
    if not model:
        return False

    # 파일 삭제
    if os.path.exists(model.file_path):
        os.remove(model.file_path)

    db.session.delete(model)
    db.session.commit()
    _log("INFO", "train", "모델 삭제: " + model.name)
    return True