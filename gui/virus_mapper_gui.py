import os
import shutil
import subprocess
import sys
from PyQt5.QtWidgets import *
from PyQt5.QtGui import QFont
from PyQt5.Qt import Qt
from PyQt5 import uic
from openpyxl import Workbook
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import pysam
from src.config import Config
from Bio import SeqIO
from src.virus_mapper import Mapper
from Bio.SeqRecord import SeqRecord
from PyQt5.QtGui import QStandardItemModel, QStandardItem
from PyQt5.QtCore import pyqtSlot


# UI파일 연결 코드
UI_class = uic.loadUiType("gui/virus_mapper.ui")[0]


class MyWindow(QMainWindow, UI_class):
    def __init__(self):
        super().__init__()
        self.setupUi(self)

        self.config = Config().config
        self.mapper = Mapper()
        self.REFGAPDB = {}
        self.sam_mapped_pos = None
        self.open_ref_fasta.clicked.connect(self.open_ref_fasta_clicked)
        self.open_tag_fasta.clicked.connect(self.open_query_fasta_clicked)
        self.run_bowtie2.clicked.connect(self.run_bowtie2_clicked)
        self.download_excel.clicked.connect(self.download_excel_clicked)

    def open_ref_fasta_clicked(self):
        """이 이벤트가 발생하면
        1. event 가 발생하면
        2. 파일 입력 창이 뜨고 fasta 파일 선택
        3. tempdir/bowtie2 아래에 fasta 파일을 넣고
        4. bowtie2 index 실행하고
        5. list_ref_seqs 창에 내용을 표시한다.
        그리고
        self.config 에 'BT2_IDX_BASE'를 저장하고 다른 method 에서 참고할 수 있게 한다.
        self.REFGAPDB 에 ref 서열들의 gap(-) 의 index 정보들을 저장해서 차후 사용 할 수 있게 한다.
        """
        in_filen = self.open_file_dialog()
        if not in_filen:
            self.show_alert(
                "No File Selected",
                "No reference FASTA file was selected. Please select a file to proceed.",
            )
            return
        ref_inx_basename = self.run_bowtie2_index(in_filen)
        self.config["BT2_IDX_BASE"] = ref_inx_basename

    def open_query_fasta_clicked(self):
        """이 이벤트가 발생하면
        1. event가 발생
        2. 파일 입력 창이 뜨고 fasta 파일 선택
        3. 열린 파일을 파싱하고 degenerancy(regular expression)를 처리하여 서열 추가 생성
        4. tempdir 아래에 추가된 서열드를 fastq 파일로 저장
        5. 화면 결과창에 입력 된 서열 갯수, 생성된 서열 갯수, 생성 서열 들을 표시한다.
        저장된 파일 경로 및 이름을 self.config['TMP_QRY_SEQ']에 저장한다.
        """
        in_filen = self.open_file_dialog()
        if not in_filen:
            self.show_alert(
                "No File Selected",
                "No query(tag) FASTA file was selected. Please select a file to proceed.",
            )
            return
        in_seqrecs = list(SeqIO.parse(open(in_filen), "fasta"))
        query_seqs = self.mapper.gen_fasta_from_degenerancy(in_seqrecs)
        # tar_fn = f"{self.config['TEMPDIR']}/{os.path.basename(in_filen)}"

        # Define the output FASTQ file path
        tar_fn = f"{self.config['TEMPDIR']}/{os.path.splitext(os.path.basename(in_filen))[0]}.fastq"

        # Convert sequences to FASTQ format with dummy quality scores
        fastq_seqs = []
        for seq in query_seqs:
            qual = "I" * len(
                seq.seq
            )  # Using a dummy quality score of 'I' (high quality)
            fastq_rec = SeqRecord(
                seq.seq,
                id=seq.id,
                description=seq.description,
                letter_annotations={"phred_quality": [40] * len(seq.seq)},
            )
            fastq_seqs.append(fastq_rec)
        # SeqIO.write(query_seqs, tar_fn, "fasta")
        # Write the FASTQ file
        SeqIO.write(fastq_seqs, tar_fn, "fastq")

        self.query_seq_res.append(f"{len(in_seqrecs)} -> {len(query_seqs)}")
        self.query_seq_res.append(
            "\n".join(["\n".join([f">{v.id}", f"{v.seq}"]) for v in query_seqs])
        )
        self.config["TMP_QRY_SEQ"] = tar_fn

    def open_file_dialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file, _ = QFileDialog.getOpenFileName(
            self,
            "Open File",
            "",
            "Fasta Files (*.fasta *.fa)",
            options=options,
        )
        return file

    def run_bowtie2_index(self, in_fasta, debugq=False):
        """1. Bowtie 2 version 2.4.4 by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)
            Usage: bowtie2-build [options]* <reference_in> <bt2_index_base>
                reference_in            comma-separated list of files with ref sequences
                bt2_index_base          write bt2 data to files with this dir/basename
        2. Save the gap positions into self.REFGAPDB (dict)
        Options:
            -f                      reference files are Fasta (default)

            Return: The path of bowtie2 index basename
        """
        tar_dir = f"{self.config['TEMPDIR']}/bowtie2"
        # Copy the file to the target directory
        if not os.path.exists(tar_dir):
            os.makedirs(tar_dir)
        destination_path = os.path.join(tar_dir, os.path.basename(in_fasta))
        shutil.copy(in_fasta, destination_path)

        # gap info db and remove gap
        gapped_fasta = list(SeqIO.parse(open(destination_path), "fasta"))
        self.REFGAPDB = self.mapper.mk_gap_dic(gapped_fasta)
        ungapped_fasta = self.mapper.rm_gaps(gapped_fasta)

        # Save the ungapped fasta to destination_path
        SeqIO.write(ungapped_fasta, destination_path, "fasta")

        # Run Bowtie2 indexing
        bowtie2_index_command = ["bowtie2-build", destination_path, destination_path]
        result = subprocess.run(bowtie2_index_command, capture_output=True, text=True)

        if result.returncode == 0:
            self.ref_indexing_res.append("Bowtie2 indexing completed successfully.")
            if debugq:
                self.ref_indexing_res.append(result.stdout)
            self.ref_indexing_res.append(f"indexing file locaton: {destination_path}")

            # Display the contents of REFGAPDB
            self.ref_indexing_res.append("REFGAPDB contents:")
            self.ref_indexing_res.append(str(self.REFGAPDB))
        else:
            self.ref_indexing_res.append("Bowtie2 indexing failed.")
            self.ref_indexing_res.append(result.stderr)

        return f"{destination_path}"

    def run_bowtie2_clicked(self):
        """1. reference path 와 query path 가 잘 저장되어있는지 확인
        2. bowtied2 running
        3. sam file 열기
        4. mapping 위치 gap 포함해서 변경
        5. 다음의 형식으로 결과 출력
                query1  query2  ... queryN
        ref1    map_pos map_pos ... map_pos
        ref2
        ref3
        ...
        refM
        """
        ## 1. check the prerequisites: reference and queries
        checkq = all(k in self.config for k in ["BT2_IDX_BASE", "TMP_QRY_SEQ"])
        if checkq:
            pass
        else:
            self.show_alert(
                "Error", "A Reference sequence or a query sequence are missing."
            )
            return
        ## 2. bowtie2 running
        sam_out = self.local_run_bowtie2()

        ## 3. SAM file parsing
        sam_mapped_pos: dict = self.mapper.parse_sam_file(
            pysam.AlignmentFile(sam_out, "r")
        )

        ## 4. adjustment of mapped postion based on the GAP data
        for ref_id, query_seq in sam_mapped_pos.items():
            for query_id, mapped_pos in query_seq.items():
                sam_mapped_pos[ref_id][query_id] = self.mapper.adj_map_pos(
                    mapped_pos, ref_id, self.REFGAPDB
                )

        ## 5. formatting the mapped dictionary
        self.populate_table_view(sam_mapped_pos, self.table_map_res)

        ## 6. Set the sam_mapped_pos as a class variable
        self.sam_mapped_pos = sam_mapped_pos

        # ## 6. save excel file
        # excel_file = os.path.join(
        #     self.config["TEMPDIR"],
        #     os.path.basename(self.config["TMP_QRY_SEQ"]).split(".")[0] + ".xlsx",
        # )
        # self.save_to_excel(sam_mapped_pos, excel_file)

    def show_alert(self, title, message):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setWindowTitle(title)
        msg.setText(message)
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()

    def local_run_bowtie2(self) -> str:
        """run bowtie2

        Return:
            SAM file output
        """
        # Retrieve the file paths from the config
        reference_file = self.config["BT2_IDX_BASE"]
        query_file = self.config["TMP_QRY_SEQ"]

        # Create output path using TEMP directory and basename of query file
        output_file = os.path.join(
            self.config["TEMPDIR"], os.path.basename(query_file).split(".")[0] + ".sam"
        )

        # Construct the Bowtie2 command
        bowtie2_command = [
            "bowtie2",
            "-a",
            "-p",
            "2",
            "-x",
            reference_file,  # Reference index base
            "-U",
            query_file,  # Query file (unpaired reads)
            "-S",
            output_file,  # Output SAM file
        ]
        print(" ".join(bowtie2_command))

        try:
            # Run the Bowtie2 command
            result = subprocess.run(
                bowtie2_command, check=True, capture_output=True, text=True
            )
            print("Bowtie2 output:")
            print(result.stdout)

        except subprocess.CalledProcessError as e:
            print("An error occurred while running Bowtie2:")
            print(e.stderr)

        return output_file

    def populate_table_view(self, data, table_view):
        # Transform data into matrix format
        query_ids, ref_ids, matrix = self.mapper.transform_data_to_matrix(data)

        # Create QStandardItemModel
        model = QStandardItemModel()

        # Set the row headers (Reference IDs)
        model.setVerticalHeaderLabels(ref_ids)

        # Set the column headers (Query IDs)
        model.setHorizontalHeaderLabels(query_ids)

        # Populate the model with data
        for row_index, ref_id in enumerate(ref_ids):
            for col_index, query_id in enumerate(query_ids):
                mapped_pos = matrix[ref_id].get(query_id, "")
                item = QStandardItem(str(mapped_pos) if mapped_pos is not None else "")
                model.setItem(row_index, col_index, item)
            # Set the model to the table view
        table_view.setModel(model)

    def save_to_excel(self, data, filepath):
        # Transform data into matrix format
        query_ids, ref_ids, matrix = self.mapper.transform_data_to_matrix(data)

        # Create a new workbook and select the active sheet
        wb = Workbook()
        ws = wb.active

        # Set the column headers (Query IDs)
        for col_index, query_id in enumerate(
            query_ids, start=2
        ):  # Start at 2 to leave space for row headers
            ws.cell(row=1, column=col_index, value=query_id)

        # Set the row headers (Reference IDs)
        for row_index, ref_id in enumerate(
            ref_ids, start=2
        ):  # Start at 2 to leave space for column headers
            ws.cell(row=row_index, column=1, value=ref_id)

        # Populate the sheet with data
        for row_index, ref_id in enumerate(ref_ids, start=2):
            for col_index, query_id in enumerate(query_ids, start=2):
                mapped_pos = matrix[ref_id].get(query_id, "")
                ws.cell(row=row_index, column=col_index, value=mapped_pos)

        # Save the workbook to the specified filepath
        wb.save(filepath)

    @pyqtSlot()
    def download_excel_clicked(self):
        if self.sam_mapped_pos is None:
            self.show_alert("Error", "Please run mapping first")
            return
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Excel File",
            "/mnt/c/Users",
            "Excel Files (*.xlsx);;All Files (*)",
            options=options,
        )
        if file_path and self.sam_mapped_pos is not None:
            self.save_to_excel(self.sam_mapped_pos, file_path)


if __name__ == "__main__":
    app = QApplication(sys.argv)

    Window = MyWindow()

    Window.show()

    app.exec_()
