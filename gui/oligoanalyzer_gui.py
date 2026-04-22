import sys
from PyQt5.QtWidgets import *
from PyQt5.QtGui import QFont
from PyQt5.Qt import Qt
from PyQt5 import uic
from openpyxl import Workbook  # ui 파일을 사용하기 위한 모듈 import
from src.oligoanalyzer import OligoAnalyzer
from Bio.Seq import Seq

try:
    from Bio.SeqUtils import GC
except ImportError:
    from Bio.SeqUtils import gc_fraction

    GC = gc_fraction
# from PyQt5.Qt import gc

# UI파일 연결 코드
UI_class = uic.loadUiType("gui/oligoanalyzer.ui")[0]


class MyWindow(QMainWindow, UI_class):
    def __init__(self):
        super().__init__()
        self.setupUi(self)

        self.initialize_class_variables_and_ui("all")
        self.uploadButton.clicked.connect(self.load_btn_clicked)
        self.step1_run_button.clicked.connect(self.step1_tm_filter)
        self.step1_excel_export.clicked.connect(self.step1_excel_export_run)
        self.step2_run_button.clicked.connect(self.step2_homodimer_hairpin_filter)
        self.step2_select_done.clicked.connect(self.step2_select_done_filter)
        self.step2_set_original_button.clicked.connect(self.step2_set_original)
        self.step2_excel_export.clicked.connect(self.step2_excel_export_run)
        self.step3_run_button.clicked.connect(self.step3_run)
        self.step3_excel_export.clicked.connect(self.step3_excel_export_run)

    def initialize_class_variables_and_ui(self, step_name: str):
        """클래스 변수 초기화 및 UI 내 화면 초기화. 각 클래스 변수는 앞 단계의 결과에 의존적이기 때문에
        all 은 전체
        step0 는 step0~step3
        step1 은 step1~step3
        step2 는 step2~step3
        step3 는 step3
        를 초기화 한다.

        Args:
            step_name (str): 다음 중 하나 ["all", "step0", "step1", "step2", "step3"]
        """

        def step0():
            self.oligo_analyzer = None  # initialize at step0
            self.in_sequences = None  # initialize at step0

        def step1():
            self.step1_tm_filter_res = {}
            self.step1_res_tt_seqn = 0
            self.step1_result_tabWidget.clear()  # UI clear

        def step2():
            self.step2_homodimer_hairpin_filter_res = {}
            self.step2_homodimer_hairpin_filter_res_ori = {}  # for recovery
            self.step2_result_tabWidget.clear()  # UI clear
            self.step2_select_done.hide()
            self.step2_set_original_button.hide()
            self.step2_status_label.setText("status: 0 sequence(s) done")
            self.step2_select_label.setText("selection results will be dispalyed here")

        def step3():
            self.step3_heterodimer_res = {}
            self.step3_label.setText("please run step1 and 2")
            self.step3_result_tabWidget.clear()

        init_func = [step0, step1, step2, step3]
        idx_dic = {"all": 0, "step0": 0, "step1": 1, "step2": 2, "step3": 3}
        for iii in range(idx_dic[step_name], len(init_func)):
            init_func[iii]()

    def load_btn_clicked(self):
        """이 이벤트가 발생하면
        1. file 열기 창을 띄우고
        2. oligoanalyzer 입력용 excel 파일을 입력받는다
        3. 그 다음 그 파일을 열여서 필요한 입력 값들을 가지고 온 다음
        4. 표를 만들어서
        5. uploadReaultTable 에 넣는다.
        """
        in_filen = self.open_file_dialog()
        self.step0_populate_table(in_filen)
        self.uploadResultTable.horizontalHeader().setSectionResizeMode(
            2, QHeaderView.Stretch
        )

    def open_file_dialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file, _ = QFileDialog.getOpenFileName(
            self, "Open File", "", "All Files (*);;Text Files (*.txt)", options=options
        )
        return file

    def step0_populate_table(self, in_filen: str):
        """step0:
        setup the first input
        ex::
        {
            "pan candida": {
                "forward": {"seq": "CGCGCTACACTGACGGA", "length": 17, "GC%": 64.7, "Tm": 57.7},
                "probe": {"seq": "TGCGCTTACTAGGAATTCCTCGTTGAAG", "length": 28, "GC%": 46.4, "Tm": 60.2},
                "reverse": {"seq": "CAGGGACGTAATCAACGCA", "length": 19, "GC%": 52.6, "Tm": 55.1},
                "template": "CGCGCTACACTGACGGAGCCAGCGAGTATTAACCTTGGCCGAGAGGTCTGGGAAATCTTGTGAAACTCCGTCGTGCTGGGGATAGAGCATTGTAATTATTGCTCTTCAACGAGGAATTCCTAGTAAGCGCAAGTCATCAGCTTGCGTTGATTACGTCCCTG",
            },
            <strain2>:{},...
        }
        Args:
            in_filen (str): _description_
        """
        self.oligo_analyzer = OligoAnalyzer()
        in_wb = self.oligo_analyzer.set_workbook(in_filen)
        ws_names = in_wb.sheetnames
        selected_ws, _ = QInputDialog.getItem(
            self, "Select String", "Choose a string:", ws_names
        )
        self.in_sequences = self.oligo_analyzer.read_sheet(selected_ws)
        in_seqlist = self.oligo_analyzer.read_sheet_to_list(self.in_sequences)
        n_rows = len(in_seqlist)
        n_cols = len(in_seqlist[0])
        self.uploadResultTable.setRowCount(n_rows)
        self.uploadResultTable.setColumnCount(n_cols)
        for iii in range(n_rows):
            for jjj in range(n_cols):
                try:
                    item = QTableWidgetItem(str(in_seqlist[iii][jjj]))
                    self.uploadResultTable.setItem(iii, jjj, item)
                except:
                    continue

    def step1_tm_filter(self):
        if self.oligo_analyzer is None:
            QMessageBox.information(self, "Alert", "Please open an input excel file")
            return None
        try:
            tm_from = float(self.tm_from_inbox.text())
            tm_to = float(self.tm_to_inbox.text())
            tm_from_probe = float(self.tm_from_inbox_probe.text())
            tm_to_probe = float(self.tm_to_inbox_probe.text())
        except (
            ValueError
        ) as e:  # Tm 값 범위가 지정 되지 않았을 때는 0~100도로 설정하고 진행
            tm_from = 0
            tm_to = 100
            tm_from_probe = 0
            tm_to_probe = 100
            QMessageBox.information(
                self, "Alert", f"Default Tm range: {tm_from} ~ {tm_to}"
            )
            self.tm_from_inbox.setText(str(tm_from))
            self.tm_to_inbox.setText(str(tm_to))
            self.tm_from_inbox_probe.setText(str(tm_from_probe))
            self.tm_to_inbox_probe.setText(str(tm_to_probe))
        # | try
        na_conc = float(self.step1_na_conc.text())
        self.initialize_class_variables_and_ui("step1")
        strains = self.in_sequences.keys()
        tt_seqn = 0
        for each_strain in strains:
            _tmp_step1_res = self.oligo_analyzer.step1_tm_filter_main(
                self.in_sequences[each_strain],
                [tm_from, tm_to, tm_from_probe, tm_to_probe],
                salt_concentration=na_conc,
            )
            self.step1_tm_filter_res[each_strain] = _tmp_step1_res
            _seqn = self.step1_tm_filter_populate_table(_tmp_step1_res, each_strain)
            tt_seqn += _seqn
        self.step2_label_2.setText(f"total {len(strains)} strains\n{tt_seqn} sequences")
        self.step1_res_tt_seqn = tt_seqn

    def step1_tm_filter_populate_table(
        self, step1_tm_filter_res: dict, strain_name: str
    ) -> int:
        """population of the input dictionary in table format
        Args:
            step1_tm_filter_res (dict):
                        {
                        "forward": {
                            "TGCACCACATGTGTTTTTTATTG": 57.51,
                            "TGCACCACATGTGTTTTTTATTGA": 58.4,
                            "TGCACCACATGTGTTTTTTATTGAA": 59.22,
                        },
                        "probe": {
                            "AATTTCTTTGGTGGCGGGAGCAATCCTACC": 72.07,
                            "ATTTCTTTGGTGGCGGGAGCAATCCTACC": 71.83,
                            "ATTTCTTTGGTGGCGGGAGCAATCCTACCG": 73.43,
                            "CTTTGGTGGCGGGAGCAATCCTACC": 70.7,
                            "CTTTGGTGGCGGGAGCAATCCTACCG": 72.59,
                            "CTTTGGTGGCGGGAGCAATCCTACCGC": 74.34,
                            "CTTTGGTGGCGGGAGCAATCCTACCGCC": 75.97,
                            "CTTTGGTGGCGGGAGCAATCCTACCGCCA": 76.07,
                            "CTTTGGTGGCGGGAGCAATCCTACCGCCAG": 77.53,
                            "TCTTTGGTGGCGGGAGCAATCCTACC": 71.02,
                            "TCTTTGGTGGCGGGAGCAATCCTACCG": 72.83,
                            "TCTTTGGTGGCGGGAGCAATCCTACCGC": 74.51,
                            "TCTTTGGTGGCGGGAGCAATCCTACCGCC": 76.07,
                            "TCTTTGGTGGCGGGAGCAATCCTACCGCCA": 76.17,
                            "TTCTTTGGTGGCGGGAGCAATCCTACC": 71.31,
                            "TTCTTTGGTGGCGGGAGCAATCCTACCG": 73.04,
                            "TTCTTTGGTGGCGGGAGCAATCCTACCGC": 74.66,
                            "TTCTTTGGTGGCGGGAGCAATCCTACCGCC": 76.17,
                            "TTTCTTTGGTGGCGGGAGCAATCCTACC": 71.58,
                            "TTTCTTTGGTGGCGGGAGCAATCCTACCG": 73.24,
                            "TTTCTTTGGTGGCGGGAGCAATCCTACCGC": 74.8,
                        },
                        "reverse": {"TTTGACTGTAAATAAAAAGTTTGG": 54.98, "TTTGACTGTAAATAAAAAGTTTGGT": 55.94},
                    }
        Return:
            int : 하나의 종에 대한 전체 서열 갯수
        """
        tab_header = ["primer type", "ori.seq.", "adj.seq.", "Tm", "GC%", "ID"]
        _tmp_widget, _tmp_layout, _tmp_table = self._step1_populate_table_set_layout(
            tab_header, step1_tm_filter_res
        )
        iii = 0
        for p_type, seqdict in step1_tm_filter_res.items():
            for seq, tm in seqdict.items():
                jjj = 0
                _tmp_table.setItem(iii, jjj, QTableWidgetItem(p_type))
                jjj += 1
                _tmp_seq = self.in_sequences[strain_name][p_type]["seq"]
                _tmp_table.setItem(iii, jjj, QTableWidgetItem(_tmp_seq))
                jjj += 1
                _tmp_table.setItem(iii, jjj, QTableWidgetItem(seq))
                jjj += 1
                _tmp_table.setItem(iii, jjj, QTableWidgetItem(str(tm)))
                jjj += 1
                _item = QTableWidgetItem()
                _tmp_gc = GC(_tmp_seq)
                _item.setData(Qt.DisplayRole, _tmp_gc)
                _tmp_table.setItem(iii, jjj, _item)
                jjj += 1
                _tmp_table.setItem(
                    iii, jjj, QTableWidgetItem(self.gen_id(strain_name, p_type, seq))
                )
                iii += 1
        _tmp_layout.addWidget(_tmp_table)
        _tmp_table.resizeColumnsToContents()
        _tmp_table.resizeRowsToContents()
        self.step1_result_tabWidget.addTab(_tmp_widget, strain_name)
        return iii

    def _step1_populate_table_set_layout(self, tab_header, step1_tm_filter_res: dict):
        ncols = len(tab_header)  # 3
        nrows = 0
        for primertype in step1_tm_filter_res:
            nrows += len(step1_tm_filter_res[primertype])
        _tmp_widget = QWidget()
        _tmp_layout = QVBoxLayout(_tmp_widget)
        _tmp_table = QTableWidget()
        _tmp_table.setSortingEnabled(True)
        _tmp_table.setColumnCount(ncols)
        _tmp_table.setRowCount(nrows)
        _tmp_table.setHorizontalHeaderLabels(tab_header)
        # _tmp_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        return _tmp_widget, _tmp_layout, _tmp_table

    def gen_id(self, strain_name: str, p_type: str, tar_seq: str) -> str:
        """주어진 strain_name 과 tar_seq (primer or probe 서열)을 바탕으로
        <strain>_{F,R,P}_<position>_<mer> 형식의 id 를 만드는 기능

        Args:
            strain_name (str): strain name of the sequence
            p_type (str): a sequence primer type like one of primer, probe, reverse
            tar_seq (str): primer/probe sequence itself

        Returns:
            str: idenfier string
        """
        template_seq = self.in_sequences[strain_name]["template"]
        if p_type != "reverse":
            position = self.oligo_analyzer.step1_2_find_align_index(
                tar_seq, template_seq
            )
        else:
            position = self.oligo_analyzer.step1_2_find_align_index(
                str(Seq(tar_seq).reverse_complement()), template_seq
            )
        seq_type = "F" if p_type == "forward" else ("R" if p_type == "reverse" else "P")
        mer = len(tar_seq)
        res = f"{strain_name}_{seq_type}_{position}_{mer}"
        return res

    def step2_homodimer_hairpin_filter(self):
        """step1_tm_filter 의 결과물인 self.step1_tm_filter_res 에 저장된 것을 input 으로 하여
        oligoanalyzer::step2_run_mfold 와
        oligoanalyzer::step2_homodimer_main 을 순차적으로 호출한다.
        """
        if len(self.step1_tm_filter_res) == 0:
            QMessageBox.information(self, "Alert", "Please run step1 Tm filter first")
            return None
        self.step2_run_button.hide()
        self.step2_status_label.setText("Processing... Wait a moment")
        QApplication.processEvents()  # flush out all event result immediately after clicking
        self.initialize_class_variables_and_ui("step2")
        step2_seq_ttnum = 0
        for each_strain, tmp_res in self.step1_tm_filter_res.items():
            _step2_mfold_res = self.oligo_analyzer.step2_run_mfold(
                tmp_res, qt_label=[each_strain, self.step2_status_label, QApplication]
            )
            _step2_final = self.oligo_analyzer.step2_homodimer_main(_step2_mfold_res)
            step2_seq_ttnum += self.step2_populate_table(_step2_final, each_strain)
            self.step2_homodimer_hairpin_filter_res[each_strain] = _step2_final
        if len(self.step2_homodimer_hairpin_filter_res_ori) == 0:
            # 처음으로 step1 에서 결과가 넘어온 상태일 경우는 selection 에서 원복하기 위해서 첫 결과를 저장해 놓는다
            self.step2_homodimer_hairpin_filter_res_ori = (
                self.step2_homodimer_hairpin_filter_res
            )
        self._step2_homodimer_hairpin_filter_post_ui_layout(step2_seq_ttnum)

    def _step2_homodimer_hairpin_filter_post_ui_layout(self, step2_seq_ttnum: int):
        self.step2_select_label.setText(f"{step2_seq_ttnum} sequences are remaining")
        self.step3_label.setText(f"total {step2_seq_ttnum} sequences")
        self.step2_run_button.show()
        self.step2_select_done.show()
        self.step2_set_original_button.show()

    def step2_populate_table(
        self, step2_final: dict, strain_name: str, default_check_status: int = 0
    ):
        """step2_final 의 형태는 아래와 같음
        {
            "forward": {
                "TGCACCACATGTGTTTTTTATTG": {
                    "Tm": 49.46,
                    "homodimer": {"Tm": 24.0, "deltaG": -5.7}
                    "mfold": {
                        "1": {
                            "Tm(Celsius)": 85.0,
                            "dG": -2.92,
                            "dH": -21.8,
                            "dS": -60.87,
                            "fold": (
                                "\n"
                                " \n"
                                "T--------|   CA\n"
                                "         GCAC  \\\n"
                                "         TGTG  C\n"
                                "GTTATTTTT^   TA\n"
                                "  20        10\n"
                                " \n"
                                " \n"
                                " "
                            ),
                        },
                        "Foldingat(Celsius)": 37.0,
                        "Mg++": 0.0,
                        "Na+": 50.0,
                    },
                },
                "TGCACCACATGTGTTTTTTATTGA": {"Tm": 50.72, "mfold": {}},
                "TGCACCACATGTGTTTTTTATTGAA": {"Tm": 51.44, "mfold": {}},
            }
        }

        Args:
            step2_final (dict): _description_
            strain_name (str): _description_
            default_check_sttus (int) : 0 is for unchecked, 1 is for checked
        """
        tab_header = [
            "select",
            "type",
            "original seq",
            "adjusted seq",
            "homodimer property",
            "hairpin structure",
            "hairpin property",
            "position to template",
            "ori. seq. == adj. seq",
            "ID",
        ]
        _tmp_widget, _tmp_layout, _tmp_table = self._step2_populate_table_set_layout(
            tab_header, step2_final
        )
        iii = 0
        for p_type, seqdict in step2_final.items():
            for seq, step2_detail in seqdict.items():
                # 0 column : selection box
                self._step2_populate_table_0_column(
                    _tmp_table, default_check_status, iii
                )
                # 1 column: primer seq type
                _item = QTableWidgetItem(p_type)
                _tmp_table.setItem(iii, 1, _item)
                # 2 column: original input sequence. Monospace font
                _item = QTableWidgetItem(self.in_sequences[strain_name][p_type]["seq"])
                _item.setFont(QFont("Monospace"))
                _tmp_table.setItem(iii, 2, _item)
                # 3. column: adjusted seq. Monospace font
                _item = QTableWidgetItem(seq)
                _item.setFont(QFont("Monospace"))
                _tmp_table.setItem(iii, 3, _item)
                # 4 column: homodimer properties
                _tmp_table.setItem(
                    iii,
                    4,
                    QTableWidgetItem(
                        self._concat_dict(step2_detail["homodimer"], sep="\n")
                    ),
                )
                # 5 and 6 columns
                _mfold_res = step2_detail["mfold"]
                self._step2_populate_table_5and6_columns(_mfold_res, _tmp_table, iii)
                # 7 column : position to template
                if p_type == "reverse":
                    # index 를 찾기위해 reverse primer 일 경우 임시로 seq 를 reverse
                    _tmp_p_seq = str(Seq(seq).reverse_complement())
                else:
                    _tmp_p_seq = seq
                _tmp_position = self.oligo_analyzer.step1_2_find_align_index(
                    _tmp_p_seq, self.in_sequences[strain_name]["template"]
                )
                _item = QTableWidgetItem()
                _item.setData(Qt.DisplayRole, _tmp_position)
                _tmp_table.setItem(iii, 7, _item)
                # 8 column : ori. seq == adj. seq
                _item = self.in_sequences[strain_name][p_type]["seq"] == seq
                _item = QTableWidgetItem(str(_item))
                _tmp_table.setItem(iii, 8, _item)
                # 9 column: identifier or 일련번호
                _item = self.gen_id(strain_name, p_type, seq)
                _item = QTableWidgetItem(_item)
                _tmp_table.setItem(iii, 9, _item)
                iii += 1
        _tmp_table.resizeColumnsToContents()
        _tmp_table.resizeRowsToContents()
        _tmp_layout.addWidget(_tmp_table)
        self.step2_result_tabWidget.addTab(_tmp_widget, strain_name)
        return iii

    def _step2_populate_table_set_layout(self, tab_header, step2_final: dict):
        ncols = len(tab_header)
        nrows = 0
        for primertype in step2_final:
            nrows += len(step2_final[primertype])
        _tmp_widget = QWidget()
        _tmp_layout = QVBoxLayout(_tmp_widget)
        _tmp_table = QTableWidget()
        _tmp_table.setSortingEnabled(True)
        _tmp_table.setColumnCount(ncols)
        _tmp_table.setRowCount(nrows)
        _tmp_table.setHorizontalHeaderLabels(tab_header)
        _tmp_table.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        _tmp_table.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        return _tmp_widget, _tmp_layout, _tmp_table

    def _step2_populate_table_0_column(self, _tmp_table, default_check_status, iii):
        _item = QTableWidgetItem()
        _item.setFlags(_item.flags() | 0x00001000)
        _item.setCheckState(Qt.Checked if default_check_status else 0)
        _tmp_table.setItem(iii, 0, _item)

    def _step2_populate_table_5and6_columns(self, _mfold_res, _tmp_table, iii):
        # _mfold_res = step2_detail["mfold"]
        if "1" in _mfold_res:
            # 5 column : hairpin structure. monospace font
            _item = QTableWidgetItem(_mfold_res["1"]["fold"])
            _item.setFont(QFont("Monospace"))
            _tmp_table.setItem(iii, 5, _item)
            # 6 column: hairpin properties
            _tmp_table.setItem(
                iii,
                6,
                QTableWidgetItem(
                    self._concat_dict(_mfold_res["1"], exc_key=["fold"], sep="\n")
                ),
            )
        else:
            _tmp_table.setItem(iii, 5, QTableWidgetItem("no result"))
            _tmp_table.setItem(iii, 6, QTableWidgetItem("no result"))

    def _concat_dict(self, tmp_dict: dict, exc_key: list = [], sep: str = ",") -> str:
        """_summary_

        Args:
            tmp_dict (dict): _description_
            exc_key (list): _description_

        Returns:
            str: _description_
        """
        res = []
        for kkk, ttt in tmp_dict.items():
            if kkk in exc_key:
                continue
            res.append(f"{kkk}:{ttt}")
        return sep.join(res)

    def step2_select_done_filter(self):
        """_summary_"""
        selection_res = {}
        final_res = {}
        num_tabs = self.step2_result_tabWidget.count()
        # 1. get the list of checked sequences
        for index in range(num_tabs):
            tab = self.step2_result_tabWidget.widget(index)
            strain_name = self.step2_result_tabWidget.tabText(index)
            each_step2_table = tab.findChild(QTableWidget)
            self._step2_select_done_filter_let_remain(
                strain_name, each_step2_table, selection_res
            )
        for strain_name in selection_res:
            final_res[strain_name] = self.oligo_analyzer.step2_selection_done(
                self.step2_homodimer_hairpin_filter_res[strain_name],
                selection_res[strain_name],
            )
        self.step2_homodimer_hairpin_filter_res = final_res
        # 2. populate again the table widget
        self.step2_result_tabWidget.clear()  # step2 UI clear
        remaining_seq_ttnum = 0
        for (
            strain_name,
            selection_done_dict,
        ) in self.step2_homodimer_hairpin_filter_res.items():
            remaining_seq_ttnum += self.step2_populate_table(
                selection_done_dict, strain_name, default_check_status=1
            )
        self.step2_select_label.setText(
            f"{remaining_seq_ttnum} sequences are remaining"
        )
        self.step3_label.setText(f"total {remaining_seq_ttnum} sequences")

    def _step2_select_done_filter_let_remain(
        self, strain_name, each_step2_table, selection_res
    ):
        for row in range(each_step2_table.rowCount()):
            item = each_step2_table.item(row, 0)
            if item.checkState() == Qt.Checked:
                seq_type = each_step2_table.item(row, 1).text()
                adj_seq = each_step2_table.item(row, 3).text()
                if strain_name not in selection_res:
                    selection_res[strain_name] = {}
                if seq_type not in selection_res[strain_name]:
                    selection_res[strain_name][seq_type] = []
                selection_res[strain_name][seq_type].append(adj_seq)

    def step2_set_original(self):
        if len(self.step2_homodimer_hairpin_filter_res_ori) == 0:
            QMessageBox.information(self, "Alert", "Please run step1 Tm filter first")
            return None
        if (
            self.step2_homodimer_hairpin_filter_res
            == self.step2_homodimer_hairpin_filter_res_ori
        ):
            QMessageBox.information(
                self, "Alert", "Already original step1 results are given"
            )
            return None
        self.step2_result_tabWidget.clear()  # UI clear
        seq_ttnum = 0
        for (
            strain_name,
            seq_dict,
        ) in self.step2_homodimer_hairpin_filter_res_ori.items():
            seq_ttnum += self.step2_populate_table(
                seq_dict, strain_name, default_check_status=0
            )
        self.step2_homodimer_hairpin_filter_res = (
            self.step2_homodimer_hairpin_filter_res_ori
        )
        self.step2_select_label.setText(f"{seq_ttnum} sequences get to step3")
        self.step3_label.setText(f"total {seq_ttnum} sequences")

    def step3_run_each_strain(self):
        """_summary_"""
        if len(self.step2_homodimer_hairpin_filter_res) == 0:
            QMessageBox.information(
                self, "Alert", "Please run step2 homodimer and hairpin filter first"
            )
            return None
        res = {}
        for strain_name, seq_dict in self.step2_homodimer_hairpin_filter_res.items():
            res[strain_name] = self.oligo_analyzer.step3_heterodimer_main(seq_dict)
        self.step3_heterodimer_res = res
        self.step3_populate_table(res)

    def step3_populate_table_each_strain(self, step3_res: dict):
        ## populate table
        self.step3_result_tabWidget.clear()
        tab_header = [
            "sequence type",
            "sequence",
            "length",
            "Tm",
            "homodimer\n properties",
            "heterodimer\n properties",
            "partner seq. type",
            "partner sequence",
            "partner seq's\n length",
            "partner seq's\n Tm",
            "partner seq's\n homodimer properties",
        ]
        for strain_name in step3_res:
            _tmp_widget, _tmp_layout, _tmp_table = (
                self._step2_populate_table_set_layout(tab_header, step3_res)
            )
            row_num = 0
            for step3_key in step3_res[strain_name]:
                seqtype1, seq1, seqtype2, seq2 = step3_key
                seq_len1 = len(seq1)
                seq_len2 = len(seq2)
                melting_temp1 = self.step2_homodimer_hairpin_filter_res[strain_name][
                    seqtype1
                ][seq1]["Tm"]
                homodimer_prop1 = self.step2_homodimer_hairpin_filter_res[strain_name][
                    seqtype1
                ][seq1]["homodimer"]
                melting_temp2 = self.step2_homodimer_hairpin_filter_res[strain_name][
                    seqtype2
                ][seq2]["Tm"]
                homodimer_prop2 = self.step2_homodimer_hairpin_filter_res[strain_name][
                    seqtype2
                ][seq2]["homodimer"]
                heterodimer_prop = step3_res[strain_name][step3_key]
                ## 0 column: sequence type 1
                _item = QTableWidgetItem(seqtype1)
                _tmp_table.setItem(row_num, 0, _item)
                ## 1 column: sequence 1
                _item = QTableWidgetItem(seq1)
                _item.setFont(QFont("Monospace"))
                _tmp_table.setItem(row_num, 1, _item)
                ## 2 column: len(seq1)
                _item = QTableWidgetItem(str(seq_len1))
                _tmp_table.setItem(row_num, 2, _item)
                ## 3 column: melting temperature
                _item = QTableWidgetItem(f"{melting_temp1: .2f}")
                _tmp_table.setItem(row_num, 3, _item)
                ## 4 column: homodimer properties
                _item = QTableWidgetItem(self._concat_dict(homodimer_prop1, sep="\n"))
                _tmp_table.setItem(row_num, 4, _item)
                ## 5 column: heterodimer properties
                _item = QTableWidgetItem(self._concat_dict(heterodimer_prop, sep="\n"))
                _tmp_table.setItem(row_num, 5, _item)
                ## 6 column: hetero seq. type
                _item = QTableWidgetItem(seqtype2)
                _tmp_table.setItem(row_num, 6, _item)
                ## 7 column: hetero sequence itself
                _item = QTableWidgetItem(seq2)
                _item.setFont(QFont("Monospace"))
                _tmp_table.setItem(row_num, 7, _item)
                ## 8 column : hetero seq.len
                _item = QTableWidgetItem(f"{seq_len2}")
                _tmp_table.setItem(row_num, 8, _item)
                ## 9 column: hetero seq. Tm
                _item = QTableWidgetItem(f"{melting_temp2: .2f}")
                _tmp_table.setItem(row_num, 9, _item)
                ## 10 column: het seq' homodimer prop.
                _item = QTableWidgetItem(self._concat_dict(homodimer_prop2, sep="\n"))
                _tmp_table.setItem(row_num, 10, _item)
                row_num += 1
            _tmp_table.resizeColumnsToContents()
            _tmp_table.resizeRowsToContents()
            _tmp_table.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
            _tmp_table.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
            _tmp_layout.addWidget(_tmp_table)
            self.step3_result_tabWidget.addTab(_tmp_widget, strain_name)

    def step3_run(self):
        """_summary_"""
        if len(self.step2_homodimer_hairpin_filter_res) == 0:
            QMessageBox.information(
                self, "Alert", "Please run step2 homodimer and hairpin filter first"
            )
            return None
        res = {}
        step3_input = self.oligo_analyzer.step3_combine_all_strains(
            self.step2_homodimer_hairpin_filter_res
        )
        res = self.oligo_analyzer.step3_heterodimer_all_strains_main(step3_input)
        self.step3_heterodimer_res = res
        self.step3_populate_table(res)

    def step3_populate_table(self, step3_res: dict):
        ## populate table
        DELIM = ":"
        self.step3_result_tabWidget.clear()
        tab_header = [
            "sequence type",  # 0
            "length",  # 1
            "Tm",  # 2
            "homodimer\n properties",  # 3
            "alignment",  # 4
            "heterodimer\n properties",  # 5
            "total\n binding",  # 6
            "partner seq. type",  # 7
            "partner seq's\n length",  # 8
            "partner seq's\n Tm",  # 9
            "partner seq's\n homodimer properties",  # 10
        ]
        _tmp_widget, _tmp_layout, _tmp_table = self._step3_populate_table_set_layout(
            tab_header, step3_res
        )
        row_num = 0
        for step3_key in step3_res:
            seqtype1, seq1, seqtype2, seq2 = step3_key
            seq_len1 = len(seq1)
            seq_len2 = len(seq2)
            strain_name1, p_type1 = seqtype1.split(DELIM)
            strain_name2, p_type2 = seqtype2.split(DELIM)
            melting_temp1 = self.step2_homodimer_hairpin_filter_res[strain_name1][
                p_type1
            ][seq1]["Tm"]
            homodimer_prop1 = self.step2_homodimer_hairpin_filter_res[strain_name1][
                p_type1
            ][seq1]["homodimer"]
            melting_temp2 = self.step2_homodimer_hairpin_filter_res[strain_name2][
                p_type2
            ][seq2]["Tm"]
            homodimer_prop2 = self.step2_homodimer_hairpin_filter_res[strain_name2][
                p_type2
            ][seq2]["homodimer"]
            heterodimer_prop = step3_res[step3_key]
            ## 0 column: sequence type 1
            _item = seqtype1 + DELIM + self.gen_id(strain_name1, p_type1, seq1)
            _item = QTableWidgetItem(_item.replace(DELIM, "\n"))
            # _item = QTableWidgetItem(seqtype1.replace(DELIM, "\n"))
            _tmp_table.setItem(row_num, 0, _item)
            ## 1 column: len(seq1)
            _item = QTableWidgetItem(str(seq_len1))
            _tmp_table.setItem(row_num, 1, _item)
            ## 2 column: melting temperature
            _item = QTableWidgetItem(f"{melting_temp1: .2f}")
            _tmp_table.setItem(row_num, 2, _item)
            ## 3 column: homodimer properties
            _item = QTableWidgetItem(self._concat_dict(homodimer_prop1, sep="\n"))
            _tmp_table.setItem(row_num, 3, _item)
            ## 4 column: pairwise alignment
            _item = QTableWidgetItem(heterodimer_prop["pairwise"])
            _item.setFont(QFont("Monospace"))
            _tmp_table.setItem(row_num, 4, _item)
            ## 5 column: heterodimer properties
            _item = QTableWidgetItem(
                self._concat_dict(
                    heterodimer_prop, sep="\n", exc_key=["pairwise", "tt_binding"]
                )
            )
            _tmp_table.setItem(row_num, 5, _item)
            ## 6 column: total bidining (integer type)
            _item = QTableWidgetItem()
            _item.setData(Qt.DisplayRole, heterodimer_prop["tt_binding"])
            _tmp_table.setItem(row_num, 6, _item)
            ## 7 column: hetero seq. type
            _item = seqtype2 + DELIM + self.gen_id(strain_name2, p_type2, seq2)
            _item = QTableWidgetItem(_item.replace(DELIM, "\n"))
            _tmp_table.setItem(row_num, 7, _item)
            ## 8 column : hetero seq.len
            _item = QTableWidgetItem(f"{seq_len2}")
            _tmp_table.setItem(row_num, 8, _item)
            ## 9 column: hetero seq. Tm
            _item = QTableWidgetItem(f"{melting_temp2: .2f}")
            _tmp_table.setItem(row_num, 9, _item)
            ## 10 column: het seq' homodimer prop.
            _item = QTableWidgetItem(self._concat_dict(homodimer_prop2, sep="\n"))
            _tmp_table.setItem(row_num, 10, _item)
            row_num += 1
        _tmp_table.resizeColumnsToContents()
        _tmp_table.resizeRowsToContents()
        _tmp_table.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        _tmp_table.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        _tmp_layout.addWidget(_tmp_table)
        scroll_area = QScrollArea()
        scroll_area.setWidget(_tmp_widget)
        self.step3_result_tabWidget.addTab(_tmp_widget, "results")

    def _step3_populate_table_set_layout(self, tab_header, step3_final: dict):
        ncols = len(tab_header)
        nrows = len(step3_final)
        # for primertype in step2_final:
        #     nrows += len(step2_final[primertype])
        _tmp_widget = QWidget()
        _tmp_layout = QVBoxLayout(_tmp_widget)
        _tmp_table = QTableWidget()
        _tmp_table.setSortingEnabled(True)
        _tmp_table.setColumnCount(ncols)
        _tmp_table.setRowCount(nrows)
        _tmp_table.setHorizontalHeaderLabels(tab_header)
        return _tmp_widget, _tmp_layout, _tmp_table

    def step1_excel_export_run(self):
        """step1 의 tabWidget 의 결과값을 가지고 와서 excel 로 아웃풋한다.
        self.step1_result_tabWidget 에 결과가 저장되어있다."""
        outfn = "step1.xlsx"
        self.to_excel(self.step1_result_tabWidget, outfn)
        QMessageBox.information(self, "Alert", f"{outfn} is generated")

    def step2_excel_export_run(self):
        outfn = "step2.xlsx"
        self.to_excel(self.step2_result_tabWidget, outfn)
        QMessageBox.information(self, "Alert", f"{outfn} is generated")

    def step3_excel_export_run(self):
        outfn = "step3.xlsx"
        self.to_excel(self.step3_result_tabWidget, outfn)
        QMessageBox.information(self, "Alert", f"{outfn} is generated")

    def to_excel(self, ori_data, wb_fn: str):
        workbook = Workbook()

        # Remove the default sheet created by Workbook
        default_sheet = workbook.active
        workbook.remove(default_sheet)

        for i in range(ori_data.count()):
            tab_name = ori_data.tabText(i)
            sheet = workbook.create_sheet(tab_name)
            tab = ori_data.widget(i)
            table = tab.layout().itemAt(0).widget()
            for row in range(table.rowCount()):
                for col in range(table.columnCount()):
                    item = table.item(row, col)
                    value = item.text() if item else ""
                    sheet.cell(row=row + 1, column=col + 1, value=value)

        workbook.save(wb_fn)


if __name__ == "__main__":
    app = QApplication(sys.argv)

    Window = MyWindow()

    Window.show()

    app.exec_()
