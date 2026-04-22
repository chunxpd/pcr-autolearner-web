"""
데이터셋 서비스 레이어
- CSV/ZIP 업로드 처리
- 서버 로컬 디렉토리 등록
- 기존 src.pcr_autolearner.AutoLearner를 재사용하여 CSV 파싱
"""
import os
import zipfile
import shutil
from datetime import datetime
from werkzeug.utils import secure_filename
from flask import current_app

from app.extensions import db
from app.models.database import Dataset, DatasetFile, Target, AppLog


def allowed_file(filename):
    """허용된 확장자인지 확인"""
    allowed = current_app.config["ALLOWED_EXTENSIONS"]
    return "." in filename and filename.rsplit(".", 1)[1].lower() in allowed


def _log(level, category, message):
    """앱 로그를 DB에 저장"""
    log_entry = AppLog(level=level, category=category, message=message)
    db.session.add(log_entry)
    db.session.commit()


def upload_files(files, dataset_name, description=""):
    """
    업로드된 파일(CSV 또는 ZIP)을 처리하여 데이터셋을 생성한다.
    
    Args:
        files: FileStorage 객체 리스트 (request.files.getlist)
        dataset_name: 데이터셋 이름
        description: 설명 (선택)
    
    Returns:
        Dataset 객체
    """
    # 1) 데이터셋 레코드 생성
    dataset = Dataset(
        name=dataset_name,
        source_type="upload",
        description=description,
    )
    db.session.add(dataset)
    db.session.flush()  # id 확보

    # 2) 데이터셋 전용 디렉토리 생성
    dataset_dir = os.path.join(
        current_app.config["UPLOAD_FOLDER"],
        f"dataset_{dataset.id}"
    )
    os.makedirs(dataset_dir, exist_ok=True)
    dataset.source_path = dataset_dir

    # 3) 파일 저장
    saved_csv_paths = []

    for file in files:
        if not file or not file.filename:
            continue
        
        filename = secure_filename(file.filename)

        if filename.lower().endswith(".zip"):
            # ZIP 파일 처리: 임시 저장 후 압축 해제
            zip_path = os.path.join(current_app.config["TEMP_FOLDER"], filename)
            file.save(zip_path)
            
            try:
                with zipfile.ZipFile(zip_path, "r") as zf:
                    for member in zf.namelist():
                        # CSV 파일만 추출
                        if member.lower().endswith(".csv") and not member.startswith("__MACOSX"):
                            # 폴더 구조 무시하고 파일명만 사용
                            csv_filename = secure_filename(os.path.basename(member))
                            if not csv_filename:
                                continue
                            target_path = os.path.join(dataset_dir, csv_filename)
                            
                            # ZIP에서 추출하여 저장
                            with zf.open(member) as src, open(target_path, "wb") as dst:
                                dst.write(src.read())
                            
                            saved_csv_paths.append((csv_filename, target_path))
            finally:
                # 임시 ZIP 파일 삭제
                if os.path.exists(zip_path):
                    os.remove(zip_path)

        elif filename.lower().endswith(".csv"):
            # 개별 CSV 파일 저장
            target_path = os.path.join(dataset_dir, filename)
            file.save(target_path)
            saved_csv_paths.append((filename, target_path))

    # 4) 각 CSV 파일을 AutoLearner로 파싱하여 버전 확인 + DB 등록
    from src.pcr_autolearner import AutoLearner
    learner = AutoLearner()

    for filename, file_path in saved_csv_paths:
        try:
            version = learner.check_file_version(file_path)
            is_abnormal = (version == "abnormal")
            
            # 행 수 (cycle 수) 대략 계산
            row_count = 0
            if not is_abnormal:
                with open(file_path, "r") as f:
                    row_count = sum(1 for line in f) - 1  # 헤더 제외

            df = DatasetFile(
                dataset_id=dataset.id,
                filename=filename,
                file_path=file_path,
                file_version=version,
                is_abnormal=is_abnormal,
                row_count=row_count,
            )
            db.session.add(df)

        except Exception as e:
            # 파싱 실패해도 파일 레코드는 생성 (abnormal 표시)
            df = DatasetFile(
                dataset_id=dataset.id,
                filename=filename,
                file_path=file_path,
                file_version="abnormal",
                is_abnormal=True,
                row_count=0,
            )
            db.session.add(df)

    db.session.commit()

    # 5) 타겟 이름 추출
    _extract_targets(dataset)

    # 6) 로그
    _log("INFO", "upload", f"데이터셋 '{dataset_name}' 생성 완료 (파일 {len(saved_csv_paths)}개)")

    return dataset


def register_local_directory(dir_path, dataset_name, description=""):
    """
    서버 내 로컬 디렉토리를 데이터셋으로 등록한다.
    파일을 복사하지 않고, 원본 경로를 그대로 참조한다.
    
    Args:
        dir_path: 서버 내 디렉토리 경로
        dataset_name: 데이터셋 이름
        description: 설명 (선택)
    
    Returns:
        Dataset 객체 또는 에러 메시지 문자열
    """
    if not os.path.isdir(dir_path):
        return f"경로가 존재하지 않습니다: {dir_path}"

    # CSV 파일 목록 확인
    csv_files = sorted([
        f for f in os.listdir(dir_path)
        if f.lower().endswith(".csv")
    ])

    if not csv_files:
        return f"해당 디렉토리에 CSV 파일이 없습니다: {dir_path}"

    # 데이터셋 생성
    dataset = Dataset(
        name=dataset_name,
        source_type="local_dir",
        source_path=dir_path,
        description=description,
    )
    db.session.add(dataset)
    db.session.flush()

    # 파일 등록
    from src.pcr_autolearner import AutoLearner
    learner = AutoLearner()

    for filename in csv_files:
        file_path = os.path.join(dir_path, filename)
        try:
            version = learner.check_file_version(file_path)
            is_abnormal = (version == "abnormal")
            row_count = 0
            if not is_abnormal:
                with open(file_path, "r") as f:
                    row_count = sum(1 for line in f) - 1

            df = DatasetFile(
                dataset_id=dataset.id,
                filename=filename,
                file_path=file_path,
                file_version=version,
                is_abnormal=is_abnormal,
                row_count=row_count,
            )
            db.session.add(df)
        except Exception:
            df = DatasetFile(
                dataset_id=dataset.id,
                filename=filename,
                file_path=file_path,
                file_version="abnormal",
                is_abnormal=True,
                row_count=0,
            )
            db.session.add(df)

    db.session.commit()

    # 타겟 추출
    _extract_targets(dataset)

    _log("INFO", "upload", f"로컬 디렉토리 '{dir_path}' 등록 완료 (파일 {len(csv_files)}개)")

    return dataset


def _extract_targets(dataset):
    """
    데이터셋의 CSV 파일들을 읽어서 Target Name을 추출하고 DB에 저장한다.
    """
    from src.pcr_autolearner import AutoLearner
    learner = AutoLearner()

    # 정상 파일만 읽기
    normal_files = DatasetFile.query.filter_by(
        dataset_id=dataset.id, is_abnormal=False
    ).all()

    for df in normal_files:
        try:
            learner.read_input(df.file_path)
        except Exception as e:
            print(f"파일 읽기 실패: {df.filename} - {e}")

    # 타겟 이름 추출
    try:
        target_names = learner.get_target_names()
    except Exception:
        target_names = []

    for tname in target_names:
        try:
            learner.set_target(tname)
            sample_count = learner.get_sample_size()
            pos_count = learner.get_pos_num()
            neg_count = sample_count - pos_count
        except Exception:
            sample_count = 0
            pos_count = 0
            neg_count = 0

        target = Target(
            dataset_id=dataset.id,
            target_name=tname,
            sample_count=sample_count,
            positive_count=pos_count,
            negative_count=neg_count,
        )
        db.session.add(target)

    db.session.commit()


def delete_dataset(dataset_id):
    """데이터셋과 관련 파일을 삭제한다."""
    dataset = Dataset.query.get(dataset_id)
    if not dataset:
        return False

    # 업로드된 파일이면 디렉토리 삭제
    if dataset.source_type == "upload" and dataset.source_path:
        if os.path.isdir(dataset.source_path):
            shutil.rmtree(dataset.source_path, ignore_errors=True)

    # DB 레코드 삭제 (cascade로 files, targets도 삭제됨)
    db.session.delete(dataset)
    db.session.commit()

    _log("INFO", "upload", f"데이터셋 '{dataset.name}' (ID:{dataset_id}) 삭제")
    return True