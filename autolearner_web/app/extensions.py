"""
Flask 확장 모듈 초기화
- SQLAlchemy 등 확장을 여기서 한 번만 생성하고, app factory에서 init_app() 호출
"""
from flask_sqlalchemy import SQLAlchemy

# SQLAlchemy 인스턴스 (아직 app에 바인딩되지 않은 상태)
db = SQLAlchemy()