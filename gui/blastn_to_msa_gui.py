import glob
import os
import shlex
import subprocess
import sys
import threading
from PyQt5.QtWidgets import *
from PyQt5.QtGui import QFont
from PyQt5.Qt import Qt
from PyQt5 import uic
from openpyxl import Workbook
from src.common import Common  # ui 파일을 사용하기 위한 모듈 import
from PyQt5.QtGui import QMovie

# UI파일 연결 코드
UI_class = uic.loadUiType("gui/blastn_to_msa.ui")[0]


class MyWindow(QMainWindow, UI_class):
    def __init__(self):
        super().__init__()
        self.setupUi(self)

        self.in_xml_files = []

        self.open_files_button.clicked.connect(self.open_files_button_clicked)
        self.blastn_to_MSA_btn.clicked.connect(self.blastn_to_MSA_button_clicked)
        self.add_group_button.setEnabled(False)  # not yet implemented
        self.download_msa_files.setEnabled(False)  # not yet implemented

    def open_files_button_clicked(self):
        """이 이벤트가 발생하면
        1. file 열기 창을 띄우고
        2. xml 파일들을 복수개 선택해서
        3. file_list_window 에 표시한다.
        """
        self.file_list_widget.clear()
        self.MSA_file_list_widget.clear()
        self.in_xml_files = []
        in_files = self.open_file_dialog()
        self.populate_file_list_window(in_files)
        self.in_xml_files = in_files

    def open_file_dialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file, _ = QFileDialog.getOpenFileNames(
            self, "Open File", "", "XML Files (*.xml);;All Files (*)", options=options
        )
        return file

    def populate_file_list_window(self, in_files):
        for file in in_files:
            item = QListWidgetItem(file)
            self.file_list_widget.addItem(item)

    def blastn_to_MSA_button_clicked(self):
        """blastn_to_MSA button 이 눌리면
        background_xml_to_msa 를 호출하고
        background_xml_to_msa 는 결과가 나오면 populate_msa_file_list_window 를 호출해서 결과를 표시한다.
        """
        self.MSA_file_list_widget.clear()
        self.blastn_to_MSA_btn.setEnabled(False)
        self.loading_label.setText(f"total {len(self.in_xml_files)} files")
        QApplication.processEvents()  # flush out all event result immediately after clicking
        self.background_xml_to_msa()
        self.loading_label.setText("Job done")
        self.blastn_to_MSA_btn.setEnabled(True)

    def background_xml_to_msa(self):
        custom_env = {
            "JOBID": f"job{Common.get_date_string()}",
        }

        if "VIRTUAL_ENV" in os.environ:
            custom_env["VIRTUAL_ENV"] = os.environ["VIRTUAL_ENV"]
        cnt = 0
        ttcnt = len(self.in_xml_files)
        for one_xml in self.in_xml_files:
            cnt += 1
            command_string = f"scripts/xml_to_msa.sh {one_xml}"
            self.loading_label.setText(f"{cnt}/{ttcnt}: (running: {one_xml})")
            QApplication.processEvents()  # flush out all event result immediately after clicking
            subprocess.run(
                command_string, env=custom_env, capture_output=False, text=True, shell=True
            )
            self.populate_msa_file_list_window(custom_env["JOBID"], one_xml)

    def populate_msa_file_list_window(self, jobid: str, one_xml: str):
        xml_file_name = os.path.splitext(os.path.basename(one_xml))[0]

        out_files = glob.glob(f"tempdir/{jobid}*{xml_file_name}*.aln")
        for file in out_files:
            item = QListWidgetItem(file)
            self.MSA_file_list_widget.addItem(item)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    Window = MyWindow()
    Window.show()
    app.exec_()
