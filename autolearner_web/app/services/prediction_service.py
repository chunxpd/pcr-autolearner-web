"""
예측 서비스 레이어
- RF 모델로 Pos/Neg 예측
- 예측 결과 DB 저장
"""
import os
import json
import numpy as np
import joblib
from datetime import datetime
from flask import current_app

from app.extensions import db
from app.models.database import Dataset, DatasetFile, MLModel, Prediction, AppLog


def _log(level, category, message):
    log_entry = AppLog(level=level, category=category, message=message)
    db.session.add(log_entry)
    db.session.commit()


def predict_with_rf(model_id, dataset_id, normalization_mode="rotation", ddcq=False):
    """
    저장된 RF 모델로 데이터셋을 예측한다.

    Returns:
        dict: {"success": True/False, "results": [...], "summary": {...}}
    """
    from src.pcr_autolearner import AutoLearner
    from src.pcr_autolearner_ai import PCRPredicter, load_data_target

    try:
        # 1) 모델 로드
        ml_model = MLModel.query.get(model_id)
        if not ml_model:
            raise ValueError("모델을 찾을 수 없습니다.")
        if not os.path.exists(ml_model.file_path):
            raise ValueError("모델 파일이 존재하지 않습니다: " + ml_model.file_path)

        model = joblib.load(ml_model.file_path)

        # 2) 데이터 로드
        dataset = Dataset.query.get(dataset_id)
        if not dataset:
            raise ValueError("데이터셋을 찾을 수 없습니다.")

        data_dir = dataset.source_path
        if not data_dir or not os.path.isdir(data_dir):
            raise ValueError("데이터셋 경로가 유효하지 않습니다.")

        runner = AutoLearner()
        runner.read_dir(data_dir)
        targets = runner.get_target_names()

        pcr_learner = PCRPredicter()
        results = []
        total_correct = 0
        total_count = 0

        for target in targets:
            input_data, input_labels, target_name, keys = load_data_target(
                runner, target,
                normalization_mode=normalization_mode,
                ddcq=ddcq,
            )

            # 전처리
            processed_data, processed_labels = pcr_learner.preprocess_data(
                input_data, input_labels
            )

            # 예측
            predictions = model.predict(processed_data)

            for i, key in enumerate(keys):
                actual = "Positive" if processed_labels[i] == 1 else "Negative"
                predicted = "Positive" if predictions[i] == 1 else "Negative"
                file_name = os.path.basename(runner.file_names.get(key, "unknown"))

                is_correct = (actual == predicted)
                if actual != "Negative" or predicted != "Negative":
                    total_count += 1
                    if is_correct:
                        total_correct += 1

                # DB에 예측 결과 저장
                pred_record = Prediction(
                    model_id=model_id,
                    dataset_id=dataset_id,
                    target_name=target_name,
                    sample_key=key,
                    actual_label=actual,
                    predicted_label=predicted,
                    file_name=file_name,
                )
                db.session.add(pred_record)

                results.append({
                    "target": target_name,
                    "key": key,
                    "actual": actual,
                    "predicted": predicted,
                    "correct": is_correct,
                    "file_name": file_name,
                })

        db.session.commit()

        accuracy = round(total_correct / total_count, 4) if total_count > 0 else 0

        _log("INFO", "predict", "예측 완료: 모델=" + ml_model.name + ", 데이터셋=" + dataset.name + ", 정확도=" + str(accuracy))

        return {
            "success": True,
            "results": results,
            "summary": {
                "total": len(results),
                "correct": total_correct,
                "accuracy": accuracy,
                "model_name": ml_model.name,
                "dataset_name": dataset.name,
            },
        }

    except Exception as e:
        _log("ERROR", "predict", "예측 실패: " + str(e))
        return {
            "success": False,
            "message": str(e),
        }