# save this as app.py
import json
import logging
import os
import os.path
import time
import uuid

# from importlib.metadata import requires
from time import sleep

from flask import Flask, jsonify, render_template, request, send_file
from werkzeug.utils import secure_filename

from src.call_blastn import CallBlast
from src.cd_hit import CDHIT
from src.check_status import ChkStatus
from src.common import Common
from src.config import Config
from src.download_genebank import DwGeneBank
from src.multiple_sequence_alignment import MSA
from src.parse_blastxml import Parser
from src.taxid_generator import TxGenerator
from src.hsp_pair_score import Score

app = Flask(__name__)
app.secret_key = "ctlabs seceret key"
config = Config()
TAXINFO = TxGenerator(config.config)
logging.basicConfig(
    level=logging.INFO,
    datefmt="%m/%d/%Y %I:%M:%S %p",
    format="%(asctime)s:%(levelname)s:%(message)s",
)


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/blastn", methods=["GET", "POST"])
def run_blastn():
    if request.method == "POST":
        data = request.form
    else:
        data = request.args
    logging.info("=============run_blastn start================")
    tt_num_primers = 0
    for kkk, vvv in data.items():
        if (
            "in_primer_sequence" in kkk
            and "reverseq" not in kkk
            and "megablastq" not in kkk
            and vvv != ""
        ):
            tt_num_primers += 1
    logging.info(f"input: {data}")
    jobid = uuid.uuid1()  # for running system
    if config.config["TESTQ"]:
        jobid = "715781dc-2fd8-11ed-861d-00155d9efa9a"  # for test
    config.config["jobid"] = jobid
    logging.info("==================run_blastn end======================")
    return render_template(
        "blastn.html", data=data, jobid=jobid, total_num_primers=tt_num_primers
    )


@app.route("/blastn_all", methods=["POST"])
def run_sub_blastn_all():
    res = dict()
    data = request.get_json()
    ## step 1 : blastn all
    print(data)
    config.config["jobid"] = data["jobid"]
    query_blast = CallBlast(config.config)
    print(f"starting: blastn_all. semaphore: {query_blast.semaphore}")
    strand = -1 if data["in_primer_sequence_reverseq"] == "on" else 1
    megablastq = True if data["in_primer_sequence_megablastq"] == "on" else False
    xmlfn = query_blast.query_blast_standalone(
        data["in_primer_sequence"], strand=strand, megablastq=megablastq
    )
    print(f"ending: blastn_all. semaphore: {query_blast.semaphore}")
    ## step 2: parse xml file
    parser = Parser(xmlfn)
    hsps_list = parser.get_hsps()
    res["res"] = parser.gen_text(hsps_list[:100])  # 결과를 보여줄 때 너무 느려서 앞의 100개만 가지고 온다
    res["cnt"] = len(hsps_list)
    return jsonify(res)


@app.route("/blastn_exclude", methods=["POST"])
def run_sub_blastn_exclude():
    sleep(
        int(config.config["SLEEPTIME"])
    )  # blastn all 과 동시 시행 될 경우 작업 거부 될 수 있어 10초 쉰다.
    res = dict()
    data = request.get_json()
    print(data)
    config.config["jobid"] = data["jobid"]
    query_blast = CallBlast(config.config)
    print(f"starting: blastn_exclude. semaphore: {query_blast.semaphore}")
    strand = -1 if data["in_primer_sequence_reverseq"] == "on" else 1
    megablastq = True if data["in_primer_sequence_megablastq"] == "on" else False
    xmlexcfn = query_blast.query_blast_exclude_taxa_standalone(
        data["in_primer_sequence"],
        data["in_exclude_txids"],
        strand=strand,
        megablastq=megablastq,
    )
    print(f"ending: blastn_exclude. semaphore: {query_blast.semaphore}")
    # res["res"] = xmlexcfn
    ## step 2: parse xml file
    parser = Parser(xmlexcfn)
    hsps_list = parser.get_hsps()
    res["res"] = parser.gen_text(hsps_list[:100])  # 결과를 보여줄 때 너무 느려서 앞의 100개만 가지고 온다
    res["cnt"] = len(hsps_list)
    return jsonify(res)


@app.route("/genebank", methods=["GET", "POST"])
def run_genebank():
    data = request.args
    config.config["jobid"] = data["jobid"]
    chk_primers = []
    for ppp in [1, 2, 3, 4]:
        if f"dw_primer{ppp}" not in data:
            continue
        chk_primers.append(data[f"dw_primer{ppp}"])
    return render_template("genebank.html", data=chk_primers, config=config.config)


@app.route("/download_genebank", methods=["GET", "POST"])
def download_genebank():
    res = dict()
    data = request.get_json()
    print(f"from ajax:{data}")
    if data["testq"] == True:
        time.sleep(10)
    config.config["jobid"] = data["jobid"]
    print("genebank 내려받기 시작")
    print(f"config:{config.config}")
    for primer_seq in data["primer_seqs"]:
        genebank = DwGeneBank(config.config, primer_seq)
        outfn = genebank.download()
        res[primer_seq] = outfn
    print("genebank 내려받기 끝")
    return jsonify(res)


@app.route("/download_genebank_cnt", methods=["GET", "POST"])
def sub_download_genebank_cnt():
    """현재 내려받고 있는 레코드의 갯수를 리턴한다.
    genebank 파일 내려받을 때.
    """
    res = {}
    data = request.get_json()
    config.config["jobid"] = data["jobid"]
    checker = ChkStatus(config.config, data["jobid"])
    res.update(checker.check_progress_log())
    logging.info(f"Check status: {res}")
    # for primer in data["primer_seqs"]:
    #     checker = ChkStatus(config.config, primer, data["jobid"])
    #     res[primer] = checker.check_num_fasta()
    return jsonify(res)


@app.route("/local_download/<jobid>/<primerseq>", methods=["GET", "POST"])
def local_download(jobid, primerseq):
    """서버에 저장된 파일을 로컬 컴퓨터로 가지고 온다"""
    assert config.config["jobid"] == jobid
    file_name = Common.get_genbank_fn(config.config, primerseq)
    return send_file(
        file_name,
        mimetype="text/csv",
        as_attachment=True,
        download_name=os.path.basename(file_name),
    )


@app.route("/multiplesequencealignment", methods=["GET", "POST"])
def multiple_sequence_alignment():
    """
    1. fasta 파일들을 upload 할 수 있도록 form 을 제공한다.
    2. clustalo 를 실행하고
    3. high homology region (HHR)을 찾고
    4. aln 파일을 내려받을 수 있게 하고
    5. 어떤 부위가 HHR 인지를 표 형식으로 화면에 출력한다.
    """
    # data = request.args
    # assert config.config["jobid"] == data["jobid"]
    # config.config["genebankfn"] = Common.get_genbank_fn(config.config)
    # msa = MSA(config.config)
    return render_template("multiplesequencealignment.html", config=config.config)


@app.route("/upload", methods=["POST"])
def upload_fasta_files():
    tmp_jobid = uuid.uuid1()
    res = {"jobid": tmp_jobid}
    up_list = []
    if request.method == "POST":
        # check if the post request has the file part
        if len(request.files) == 0:
            return "업로드할 파일이 없습니다."
        files = request.files
        # If the user does not select a file, the browser submits an
        # empty file without a filename.
    for fa_file in files:
        file_storage = files[fa_file]
        if file_storage.filename == "":
            continue
        if not (file_storage.filename and allowed_file(file_storage.filename)):
            continue
        ofilename = secure_filename(file_storage.filename)
        os.makedirs(Common.get_upload_path(config.config, tmp_jobid), exist_ok=True)
        file_storage.save(
            Common.get_upload_file_path(config.config, tmp_jobid, ofilename)
        )
        up_list.append(file_storage.filename)
    res["size"] = len(up_list)
    res["data"] = up_list
    return jsonify(res)


@app.route("/cd_hit", methods=["GET", "POST"])
def cd_hit():
    """
    input ::
        {
            'files': '{
                "data":[fa1, fa2, ...],
                "jobid":uuid,
                "size":2
            }',
            'lower_bound': '0.4'
        }

    Returns:
        json: {"res": final_fasta, "count": fasta내 record 수}
    """
    data = request.get_json()
    logging.info(data)
    cd_hit_obj = CDHIT()
    logging.info("START:: CONCAT input fastas to ")
    cat_fa = cd_hit_obj.concat_files(
        config.config, data["files"]["jobid"], data["files"]["data"]
    )
    logging.info(f"END:: CONCAT input fastas to {cat_fa}")
    logging.info("START:: CD-HIT")
    logging.info(config)
    res = cd_hit_obj.main(cat_fa, float(data["lower_bound"]))
    logging.info(f"END:: CD-HIT {res}")
    return jsonify(res)


def allowed_file(filename, ftype="fa"):
    if ftype == "fa":
        return "." in filename and filename.rsplit(".", 1)[1].lower() in ["fa", "fasta"]
    elif ftype == "xml":
        return "." in filename and filename.rsplit(".", 1)[1].lower() in ["xml"]


@app.route("/run_msa", methods=["GET", "POST"])
def run_msa():
    res = {}
    data = request.get_json()
    logging.info("START:: run_msa")
    logging.info(f"    run_msa config: {data['config']}")
    logging.info(f"    run_msa in_fa: {data['in_fa']}")
    msa = MSA(data["config"])
    aln_fn = msa.run_msa(data["in_fa"])
    hhr_regs = msa.find_hrr(aln_fn)
    hhr_table = msa.hrr_list_to_html(hhr_regs)
    res["msa_res"] = aln_fn
    res["hhr_res"] = hhr_table
    logging.info("END:: run_msa")
    return jsonify(res)


@app.route("/local_download_msa/<filename>", methods=["GET", "POST"])
def local_download_msa(filename):
    """서버에 저장된 파일을 로컬 컴퓨터로 가지고 온다"""
    # assert config.config["jobid"] == jobid
    # file_name = Common.get_msa_fn(config.config)
    file_name = filename
    return send_file(
        file_name,
        mimetype="text/csv",
        as_attachment=True,
        download_name="clustalo.aln",
    )


@app.route("/read_clustalo_log", methods=["GET", "POST"])
def read_clustalo_log():
    data = request.get_json()
    try:
        with open(data["log_fn"], "r") as handle:
            cont = handle.readlines()
            cont = "".join(cont[-10:])
    except:
        cont = ""
    return jsonify({"res": cont})


@app.route("/taxid_test", methods=["GET", "POST"])
def taxid_test():
    return render_template("taxid_test.html")


@app.route("/taxid_search", methods=["GET"])
def taxid_search():
    qstr = request.args.get("qstr")
    logging.info(f"autocomplete str: {qstr}")
    res = TAXINFO.get_taxinfo_query(qstr)
    return jsonify(json_list=res)


@app.route("/blastnpairscores.html", methods=["GET", "POST"])
def blastn_pair_scores():
    if request.method == "POST":
        data = request.form
    else:
        data = request.args
    return render_template("blastnpairscores.html")


@app.route("/upload_blastn", methods=["POST"])
def upload_blastn_files():
    """blastn xml 파일들을 업로드 한다."""
    tmp_jobid = uuid.uuid1()
    res = {"jobid": tmp_jobid}
    up_list = []
    if request.method == "POST":
        # check if the post request has the file part
        if len(request.files) == 0:
            return "업로드할 파일이 없습니다."
        files = request.files
        # If the user does not select a file, the browser submits an
        # empty file without a filename.
    for fa_file in files:
        file_storage = files[fa_file]
        if file_storage.filename == "":
            continue
        if not (
            file_storage.filename and allowed_file(file_storage.filename, ftype="xml")
        ):
            continue
        ofilename = secure_filename(file_storage.filename)
        os.makedirs(Common.get_upload_path(config.config, tmp_jobid), exist_ok=True)
        file_storage.save(
            Common.get_upload_file_path(config.config, tmp_jobid, ofilename)
        )
        up_list.append(file_storage.filename)
    res["size"] = len(up_list)
    res["data"] = up_list
    return jsonify(res)


@app.route("/cal_blastn_xml_pair", methods=["POST"])
def cal_blastn_xml_pair():
    res = {"jobid": "test", "filename": "test"}
    data = request.get_json()
    xml_files = data["files"]
    jobid = data["jobid"]
    runner = Score()
    hspses = []
    for fn in xml_files:
        filepath = Common.get_upload_file_path(config.config, jobid, fn)
        parser = Parser(filepath)
        prc_hsps = runner.preprocessing(parser.get_hsps())
        hspses.append(prc_hsps)
    cmb_hspes = runner.combine_hsp_sets(*hspses)
    outstr = runner.mk_csv_file(cmb_hspes, delim="\t")
    outfn = Common.get_blastn_pair_score_fn(config.config, jobid, "hsp_score_table.tsv")
    with open(outfn, "w") as handle:
        handle.write(outstr)
    res["jobid"] = jobid
    res["filename"] = "hsp_score_table.tsv"
    print(res)
    return jsonify(res)


@app.route("/local_download_general/<jobid>/<filename>", methods=["GET", "POST"])
def local_download_general(jobid, filename):
    """서버에 저장된 파일을 로컬 컴퓨터로 가지고 온다"""
    file_name = Common.get_blastn_pair_score_fn(config.config, jobid, filename)
    return send_file(
        file_name,
        mimetype="text/csv",
        as_attachment=True,
        download_name=os.path.basename(file_name),
    )


if __name__ == "__main__":

    app.run(host="0.0.0.0")
