"""
SQLAlchemy 데이터베이스 모델
- 모든 테이블 정의를 이 파일에 모아둠
- 테이블: datasets, dataset_files, targets, preprocessing_runs,
          ml_models, training_runs, predictions, app_logs
"""
from datetime import datetime
from app.extensions import db


class Dataset(db.Model):
    """업로드/등록된 데이터셋 단위"""
    __tablename__ = "datasets"

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(200), nullable=False)
    source_type = db.Column(db.String(20), nullable=False)     # 'upload' 또는 'local_dir'
    source_path = db.Column(db.String(500), nullable=True)
    description = db.Column(db.Text, nullable=True)
    created_at = db.Column(db.DateTime, default=datetime.utcnow)

    # 관계
    files = db.relationship("DatasetFile", backref="dataset", lazy=True, cascade="all, delete-orphan")
    targets = db.relationship("Target", backref="dataset", lazy=True, cascade="all, delete-orphan")

    def __repr__(self):
        return f"<Dataset {self.id}: {self.name}>"


class DatasetFile(db.Model):
    """데이터셋 내 개별 CSV 파일"""
    __tablename__ = "dataset_files"

    id = db.Column(db.Integer, primary_key=True)
    dataset_id = db.Column(db.Integer, db.ForeignKey("datasets.id"), nullable=False)
    filename = db.Column(db.String(300), nullable=False)
    file_path = db.Column(db.String(500), nullable=False)
    file_version = db.Column(db.String(20), nullable=True)
    is_abnormal = db.Column(db.Boolean, default=False)
    row_count = db.Column(db.Integer, nullable=True)
    created_at = db.Column(db.DateTime, default=datetime.utcnow)

    def __repr__(self):
        return f"<DatasetFile {self.filename}>"


class Target(db.Model):
    """데이터셋에서 추출된 Target Name"""
    __tablename__ = "targets"

    id = db.Column(db.Integer, primary_key=True)
    dataset_id = db.Column(db.Integer, db.ForeignKey("datasets.id"), nullable=False)
    target_name = db.Column(db.String(200), nullable=False)
    sample_count = db.Column(db.Integer, default=0)
    positive_count = db.Column(db.Integer, default=0)
    negative_count = db.Column(db.Integer, default=0)
    created_at = db.Column(db.DateTime, default=datetime.utcnow)

    def __repr__(self):
        return f"<Target {self.target_name}>"


class PreprocessingRun(db.Model):
    """전처리(Baseline/Normalize) 실행 이력"""
    __tablename__ = "preprocessing_runs"

    id = db.Column(db.Integer, primary_key=True)
    dataset_id = db.Column(db.Integer, db.ForeignKey("datasets.id"), nullable=False)
    target_name = db.Column(db.String(200), nullable=False)
    method = db.Column(db.String(30), nullable=False)
    mode = db.Column(db.String(30), default="rotation")
    smoothing_size = db.Column(db.Integer, default=3)
    created_at = db.Column(db.DateTime, default=datetime.utcnow)


class MLModel(db.Model):
    """저장된 ML 모델 (RF, Anomaly 등)"""
    __tablename__ = "ml_models"

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(200), nullable=False)
    model_type = db.Column(db.String(20), nullable=False)
    file_path = db.Column(db.String(500), nullable=False)
    accuracy = db.Column(db.Float, nullable=True)
    params_json = db.Column(db.Text, nullable=True)
    dataset_id = db.Column(db.Integer, db.ForeignKey("datasets.id"), nullable=True)
    created_at = db.Column(db.DateTime, default=datetime.utcnow)

    training_runs = db.relationship("TrainingRun", backref="model", lazy=True)
    predictions = db.relationship("Prediction", backref="model", lazy=True)

    def __repr__(self):
        return f"<MLModel {self.name} ({self.model_type})>"


class TrainingRun(db.Model):
    """학습 실행 이력"""
    __tablename__ = "training_runs"

    id = db.Column(db.Integer, primary_key=True)
    model_id = db.Column(db.Integer, db.ForeignKey("ml_models.id"), nullable=True)
    dataset_id = db.Column(db.Integer, db.ForeignKey("datasets.id"), nullable=True)
    accuracy = db.Column(db.Float, nullable=True)
    sample_count = db.Column(db.Integer, nullable=True)
    params_json = db.Column(db.Text, nullable=True)
    status = db.Column(db.String(20), default="pending")
    error_message = db.Column(db.Text, nullable=True)
    started_at = db.Column(db.DateTime, default=datetime.utcnow)
    finished_at = db.Column(db.DateTime, nullable=True)


class Prediction(db.Model):
    """예측 결과 개별 레코드"""
    __tablename__ = "predictions"

    id = db.Column(db.Integer, primary_key=True)
    model_id = db.Column(db.Integer, db.ForeignKey("ml_models.id"), nullable=False)
    dataset_id = db.Column(db.Integer, db.ForeignKey("datasets.id"), nullable=False)
    target_name = db.Column(db.String(200), nullable=True)
    sample_key = db.Column(db.Integer, nullable=True)
    actual_label = db.Column(db.String(20), nullable=True)
    predicted_label = db.Column(db.String(20), nullable=True)
    file_name = db.Column(db.String(300), nullable=True)
    created_at = db.Column(db.DateTime, default=datetime.utcnow)


class AppLog(db.Model):
    """애플리케이션 로그"""
    __tablename__ = "app_logs"

    id = db.Column(db.Integer, primary_key=True)
    level = db.Column(db.String(10), default="INFO")
    category = db.Column(db.String(30), nullable=True)
    message = db.Column(db.Text, nullable=False)
    created_at = db.Column(db.DateTime, default=datetime.utcnow)