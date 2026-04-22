"""
Flask 앱 실행 진입점
사용법: python run.py
브라우저에서 http://localhost:5000 접속
"""
from app import create_app

app = create_app()

if __name__ == "__main__":
    print("=" * 50)
    print("  PCR Autolearner Web MVP")
    print("  http://localhost:5000")
    print("=" * 50)
    app.run(debug=True, host="0.0.0.0", port=5000)