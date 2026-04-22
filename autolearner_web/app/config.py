"""
Flask 애플리케이션 설정
- DB 경로, 업로드 경로, 모델 저장 경로 등을 중앙 관리
"""
import os

# 프로젝트 루트: autolearner_web/ 디렉토리
BASE_DIR = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))


class Config:
    """Flask 기본 설정"""

    # Flask 시크릿 키 (세션 등에 사용)
    SECRET_KEY = os.environ.get("SECRET_KEY", "autolearner-dev-secret-key")

    # SQLite DB 경로
    SQLALCHEMY_DATABASE_URI = (
        f"sqlite:///{os.path.join(BASE_DIR, 'instance', 'autolearner.db')}"
    )
    SQLALCHEMY_TRACK_MODIFICATIONS = False

    # 파일 업로드 설정
    UPLOAD_FOLDER = os.path.join(BASE_DIR, "data", "uploads")
    PROCESSED_FOLDER = os.path.join(BASE_DIR, "data", "processed")
    TEMP_FOLDER = os.path.join(BASE_DIR, "data", "temp")
    MAX_CONTENT_LENGTH = 500 * 1024 * 1024  # 최대 500MB

    # 모델 저장 경로
    RF_MODEL_FOLDER = os.path.join(BASE_DIR, "saved_models", "rf")
    ANOMALY_MODEL_FOLDER = os.path.join(BASE_DIR, "saved_models", "anomaly")

    # 허용 확장자
    ALLOWED_EXTENSIONS = {"csv", "zip"}

    # ★ 기존 repo 루트 경로 (src/ 모듈 import용)
    # autolearner_web 의 상위 디렉토리 = D:\program\igenetech\ctlabs-primerwrapper-main
    REPO_ROOT = os.path.abspath(os.path.join(BASE_DIR, ".."))