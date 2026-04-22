"""
Flask Application Factory
- 앱 생성, 설정 로드, 확장 초기화, Blueprint 등록을 한 곳에서 처리
- 기존 repo의 src/ 모듈을 import할 수 있도록 sys.path 설정
"""
import os
import sys
from flask import Flask
from app.config import Config
from app.extensions import db


def create_app(config_class=Config):
    """Flask 앱을 생성하고 반환하는 팩토리 함수"""

    app = Flask(__name__)
    app.config.from_object(config_class)

    # ──────────────────────────────────────────────
    # 기존 repo 루트를 sys.path에 추가
    # 이렇게 하면 from src.pcr_autolearner import AutoLearner 가 가능
    # ──────────────────────────────────────────────
    repo_root = app.config["REPO_ROOT"]
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)

    # ──────────────────────────────────────────────
    # 확장 초기화
    # ──────────────────────────────────────────────
    db.init_app(app)

    # ──────────────────────────────────────────────
    # 필요한 디렉토리 자동 생성
    # ──────────────────────────────────────────────
    for folder_key in [
        "UPLOAD_FOLDER",
        "PROCESSED_FOLDER",
        "TEMP_FOLDER",
        "RF_MODEL_FOLDER",
        "ANOMALY_MODEL_FOLDER",
    ]:
        folder_path = app.config[folder_key]
        os.makedirs(folder_path, exist_ok=True)

    # instance 폴더 (DB 파일 위치)
    instance_dir = os.path.join(
        os.path.abspath(os.path.dirname(os.path.dirname(__file__))),
        "instance"
    )
    os.makedirs(instance_dir, exist_ok=True)

    # ──────────────────────────────────────────────
    # Blueprint 등록
    # ──────────────────────────────────────────────
    from app.routes.main import main_bp
    app.register_blueprint(main_bp)

     # ★ 추가: 데이터셋 라우트
    from app.routes.datasets import datasets_bp
    app.register_blueprint(datasets_bp)

     # ★ 추가: 분석 라우트
    from app.routes.analysis import analysis_bp
    app.register_blueprint(analysis_bp)

    # ★ 추가: 학습/모델/예측 라우트
    from app.routes.training import training_bp
    app.register_blueprint(training_bp)

    # ──────────────────────────────────────────────
    # DB 테이블 자동 생성 (개발용)
    # ──────────────────────────────────────────────
    with app.app_context():
        from app.models import database  # noqa: F401
        db.create_all()

    return app