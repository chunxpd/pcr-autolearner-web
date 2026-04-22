"""
데이터셋 관리 라우트
- 업로드 페이지
- 목록 페이지
- 상세 페이지
- 삭제
- 로컬 디렉토리 등록
"""
from flask import Blueprint, render_template, request, redirect, url_for, flash, jsonify
from app.models.database import Dataset, DatasetFile, Target
from app.services.dataset_service import (
    upload_files,
    register_local_directory,
    delete_dataset,
)

datasets_bp = Blueprint("datasets", __name__, url_prefix="/datasets")


@datasets_bp.route("/")
def list_datasets():
    """데이터셋 목록 페이지"""
    datasets = Dataset.query.order_by(Dataset.created_at.desc()).all()
    return render_template("datasets/list.html", datasets=datasets)


@datasets_bp.route("/upload", methods=["GET", "POST"])
def upload():
    """CSV/ZIP 업로드 페이지"""
    if request.method == "POST":
        # 업로드 방식 확인
        upload_type = request.form.get("upload_type", "file")
        dataset_name = request.form.get("dataset_name", "").strip()
        description = request.form.get("description", "").strip()

        if not dataset_name:
            flash("데이터셋 이름을 입력해주세요.", "error")
            return redirect(url_for("datasets.upload"))

        if upload_type == "file":
            # 파일 업로드 방식
            files = request.files.getlist("files")
            if not files or all(f.filename == "" for f in files):
                flash("파일을 선택해주세요.", "error")
                return redirect(url_for("datasets.upload"))

            dataset = upload_files(files, dataset_name, description)
            flash(f"데이터셋 '{dataset.name}'이 생성되었습니다.", "success")
            return redirect(url_for("datasets.detail", dataset_id=dataset.id))

        elif upload_type == "local_dir":
            # 로컬 디렉토리 등록 방식
            dir_path = request.form.get("dir_path", "").strip()
            if not dir_path:
                flash("디렉토리 경로를 입력해주세요.", "error")
                return redirect(url_for("datasets.upload"))

            result = register_local_directory(dir_path, dataset_name, description)
            if isinstance(result, str):
                # 에러 메시지
                flash(result, "error")
                return redirect(url_for("datasets.upload"))
            else:
                flash(f"데이터셋 '{result.name}'이 등록되었습니다.", "success")
                return redirect(url_for("datasets.detail", dataset_id=result.id))

    return render_template("datasets/upload.html")


@datasets_bp.route("/<int:dataset_id>")
def detail(dataset_id):
    """데이터셋 상세 페이지"""
    dataset = Dataset.query.get_or_404(dataset_id)
    files = DatasetFile.query.filter_by(dataset_id=dataset_id).order_by(DatasetFile.filename).all()
    targets = Target.query.filter_by(dataset_id=dataset_id).order_by(Target.target_name).all()

    return render_template(
        "datasets/detail.html",
        dataset=dataset,
        files=files,
        targets=targets,
    )


@datasets_bp.route("/<int:dataset_id>/delete", methods=["POST"])
def delete(dataset_id):
    """데이터셋 삭제"""
    success = delete_dataset(dataset_id)
    if success:
        flash("데이터셋이 삭제되었습니다.", "success")
    else:
        flash("데이터셋을 찾을 수 없습니다.", "error")
    return redirect(url_for("datasets.list_datasets"))

@datasets_bp.route("/browse-dir")
def browse_directory():
    """서버 내 디렉토리 탐색 API"""
    import os
    path = request.args.get("path", "")

    # 기본 경로 (빈 값이면 드라이브 목록)
    if not path:
        # Windows 드라이브 목록
        import string
        drives = []
        for letter in string.ascii_uppercase:
            drive = letter + ":\\"
            if os.path.exists(drive):
                drives.append({"name": drive, "path": drive, "type": "drive"})
        return jsonify({"items": drives, "current": ""})

    # 해당 경로의 하위 폴더 목록
    if not os.path.isdir(path):
        return jsonify({"items": [], "current": path, "error": "경로가 존재하지 않습니다."})

    items = []
    try:
        for name in sorted(os.listdir(path)):
            full_path = os.path.join(path, name)
            if os.path.isdir(full_path):
                items.append({"name": name, "path": full_path, "type": "folder"})
    except PermissionError:
        return jsonify({"items": [], "current": path, "error": "접근 권한이 없습니다."})

    # 상위 폴더
    parent = os.path.dirname(path.rstrip("\\/"))
    if parent == path.rstrip("\\/"):
        parent = ""  # 루트면 드라이브 목록으로

    # CSV 파일 수 표시
    csv_count = 0
    try:
        csv_count = len([f for f in os.listdir(path) if f.lower().endswith(".csv")])
    except Exception:
        pass

    return jsonify({
        "items": items,
        "current": path,
        "parent": parent,
        "csv_count": csv_count,
    })