import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QVBoxLayout,
    QWidget,
    QComboBox,
    QSlider,
    QLabel,
    QLineEdit,
)
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from src.pcr_analyzer import PCRmachine
from colorama import Fore, Back, Style


def read_excel(file_name, sheet_name):
    input_data = pd.read_excel(file_name, sheet_name=sheet_name)

    # 빈 행 처리: 모든 값이 NaN인 첫 번째 행을 찾고 데이터 슬라이싱
    if input_data.isnull().all(axis=1).any():
        input_data = input_data.loc[: input_data.isnull().all(axis=1).idxmax() - 1]

    df = input_data.set_index("X")
    df.columns = df.columns.astype(str)
    df_cq = df.apply(pd.to_numeric, errors="coerce")
    return df_cq


def read_csv(file_name):
    """
    Reads a CSV file and processes it into two outputs:
    1. A pandas DataFrame for the first block (before the first empty line).
    2. A list of headers based on the 'Target Name' column in the third block.

    Args:
        file_path (str): Path to the CSV file.

    Returns:
        df_cq: DataFrame for the first block
    """
    first_block = []
    third_block_headers = []
    third_block_started = False

    with open(file_name, "r") as file:
        for line in file:
            line = line.strip()

            # Stop reading the first block at the first empty line
            if not line and not third_block_started:
                break

            first_block.append(line)

        for line in file:
            line = line.strip()
            # Detect the start of the third block
            if "Num, Target Name" in line:
                third_block_started = True
                continue

            # Process the first block
            if not third_block_started:
                continue

            # Process the third block to extract headers
            if third_block_started and line:
                columns = line.split(",")
                if len(columns) > 1:  # Ensure there are at least two columns
                    third_block_headers.append(
                        columns[1].strip()
                    )  # Extract 'Target Name'

    # Convert the first block into a pandas DataFrame
    if first_block:
        first_block_data = [row.split(",") for row in first_block]
        df = pd.DataFrame(first_block_data[1:], columns=first_block_data[0])
        df.columns = df.columns.astype(str)
        df = df.apply(pd.to_numeric, errors="coerce")
    else:
        df = pd.DataFrame()  # Return an empty DataFrame if no data is found
    if "X" in df.columns:
        # Set the first column ("X") as the index
        df_cq = df.set_index("X")
    print(third_block_headers)
    df_cq.columns = third_block_headers

    print(df_cq)
    return df_cq


class PCRAnalyzerApp(QMainWindow):
    def __init__(
        self, df_cq, cycle_tt_num, ddcq_threshold_scale=0.01, update_func_type=1
    ):
        super().__init__()
        self.df_cq = df_cq
        self.sample_index = df_cq.columns[0]
        self.arrow_annotation = None
        self.cycle_tt_num = cycle_tt_num
        self.max_baseline_end = cycle_tt_num
        self.ignore_init_cycles = 3  # Default value for ignore_init_cycles
        self.ddcq_threshold_scal = ddcq_threshold_scale
        self.smoothing_winsize = 3  # Default smoothing window size
        self.update_plot = [self.update_plot2, self.update_plot1][update_func_type]

        self.initUI()

    def initUI(self):
        self.setWindowTitle("PCR Analyzer")

        # Main widget
        central_widget = QWidget(self)
        self.setCentralWidget(central_widget)

        # Layout
        layout = QVBoxLayout(central_widget)

        # Sample selection dropdown
        self.sample_dropdown = QComboBox(self)
        self.sample_dropdown.addItems(self.df_cq.columns)
        self.sample_dropdown.currentTextChanged.connect(self.on_dropdown_change)
        layout.addWidget(self.sample_dropdown)

        # Update function dropdown
        self.plot_dropdown = QComboBox(self)
        self.plot_dropdown.addItems(["update_plot2", "update_plot1"])
        self.plot_dropdown.currentIndexChanged.connect(self.on_plot_dropdown_change)
        layout.addWidget(self.plot_dropdown)

        # Smoothing window size dropdown
        self.smoothing_dropdown = QComboBox(self)
        self.smoothing_dropdown.addItems(["3", "5"])
        self.smoothing_dropdown.currentTextChanged.connect(self.on_smoothing_change)
        layout.addWidget(QLabel("Smoothing Window Size"))
        layout.addWidget(self.smoothing_dropdown)

        # Cycle threshold input
        self.cycle_input = QLineEdit(self)
        self.cycle_input.setText(str(self.cycle_tt_num))
        self.cycle_input.setPlaceholderText("Total cycle count")
        self.cycle_input.editingFinished.connect(self.on_cycle_input_change)
        layout.addWidget(QLabel("Max baseline"))
        layout.addWidget(self.cycle_input)

        # Ignore initial cycles input
        self.ignore_cycles_input = QLineEdit(self)
        self.ignore_cycles_input.setText(str(self.ignore_init_cycles))
        self.ignore_cycles_input.setPlaceholderText("Ignore initial cycles")
        self.ignore_cycles_input.editingFinished.connect(
            self.on_ignore_cycles_input_change
        )
        layout.addWidget(QLabel("Ignore Initial Cycles"))
        layout.addWidget(self.ignore_cycles_input)

        # Matplotlib figure and canvas
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvas(self.fig)
        layout.addWidget(self.canvas)

        # Matplotlib toolbar
        self.toolbar = NavigationToolbar(self.canvas, self)
        layout.addWidget(self.toolbar)

        # Cycle slider
        self.cycle_slider = QSlider(Qt.Horizontal, self)
        self.cycle_slider.setMinimum(1)
        self.cycle_slider.setMaximum(len(self.df_cq))
        self.cycle_slider.valueChanged.connect(self.on_slider_change)
        layout.addWidget(self.cycle_slider)

        # Cycle label
        self.label = QLabel(f"Cycle: {self.cycle_slider.value()}", self)
        layout.addWidget(self.label)

        # ddCq Threshold slider
        self.ddcq_slider = QSlider(Qt.Horizontal, self)
        self.ddcq_slider.setMinimum(0)
        self.ddcq_slider.setMaximum(1000)
        self.ddcq_slider.setValue(10)
        self.ddcq_slider.valueChanged.connect(self.on_ddcq_slider_change)
        layout.addWidget(self.ddcq_slider)

        # ddCq Threshold label
        self.ddcq_label = QLabel(
            f"ddCq Threshold Scale: {self.ddcq_threshold_scal * 100: .2f}%", self
        )
        layout.addWidget(self.ddcq_label)

        # Initial plot setup
        (self.l,) = self.ax.plot([], [], ".-", label="original")
        (self.l2,) = self.ax.plot([], [], ".-", label="rotation")
        self.ax.legend()

        # Initial plot update
        self.update_plot(self.sample_index, self.cycle_slider.value())

    def on_dropdown_change(self):
        self.sample_index = self.sample_dropdown.currentText()
        self.update_plot(self.sample_index, self.cycle_slider.value())

    def on_plot_dropdown_change(self, index):
        self.update_plot = [self.update_plot2, self.update_plot1][index]
        self.update_plot(self.sample_index, self.cycle_slider.value())

    def on_smoothing_change(self):
        self.smoothing_winsize = int(self.smoothing_dropdown.currentText())
        self.update_plot(self.sample_index, self.cycle_slider.value())

    def on_slider_change(self):
        self.label.setText(f"Cycle: {self.cycle_slider.value()}")
        self.update_plot(self.sample_index, self.cycle_slider.value())

    def on_ddcq_slider_change(self):
        self.ddcq_threshold_scal = self.ddcq_slider.value() / 10000.0
        self.ddcq_label.setText(
            f"ddCq Threshold Scale: {self.ddcq_threshold_scal * 100: .2f}%"
        )
        self.update_plot(self.sample_index, self.cycle_slider.value())

    def on_cycle_input_change(self):
        try:
            self.max_baseline_end = int(self.cycle_input.text())
            self.update_plot(self.sample_index, self.cycle_slider.value())
        except ValueError:
            pass  # Ignore invalid input

    def on_ignore_cycles_input_change(self):
        try:
            self.ignore_init_cycles = int(self.ignore_cycles_input.text())
            self.update_plot(self.sample_index, self.cycle_slider.value())
        except ValueError:
            pass  # Ignore invalid input

    def _draw_arrow(self, sample_index, baseline_end):
        """baseline_end 위치에 화살표
        Args:
            sample_index: 현재 선택된 sample index
            baseline_end: 1-based index of the baseline end. 첫번째라는 뜻
        """
        baseline_end_offset = self.df_cq[sample_index].index[0] - 1
        if baseline_end:  # baseline_end 위치에 화살표
            if self.arrow_annotation is not None:
                try:
                    self.arrow_annotation.remove()
                except ValueError:
                    self.arrow_annotation = None

            self.arrow_annotation = self.ax.annotate(
                f"Baseline End({baseline_end+baseline_end_offset})",
                xy=(
                    self.df_cq[sample_index].index[baseline_end + baseline_end_offset],
                    self.df_cq[sample_index][baseline_end],
                ),
                xytext=(
                    self.df_cq[sample_index].index[
                        baseline_end + baseline_end_offset - 1
                    ],
                    self.df_cq[sample_index][baseline_end + baseline_end_offset] * 1.1,
                ),
                arrowprops=dict(facecolor="red", shrink=0.05),
            )
        else:
            if self.arrow_annotation is not None:
                try:
                    self.arrow_annotation.remove()
                except ValueError:
                    self.arrow_annotation = None

    def _set_plot_axis(self, sample_index):
        self.ax.set_xlim([self.df_cq.index.min(), self.df_cq.index.max()])
        self.ax.set_ylim(
            [self.df_cq[sample_index].min() * 0.9, self.df_cq[sample_index].max() * 1.1]
        )

        self.ax.set_title(f"Sample Index: {sample_index}")

        self.ax.relim()
        self.ax.autoscale_view()
        self.canvas.draw()

    def update_plot2(self, sample_index, cycle_slider_val):
        """slide bar 의 선택된 cycle까지 모두 입력하고 T 계산하고 baseline_end 찾기"""
        self.ax.clear()
        (self.l,) = self.ax.plot([], [], ".-", label="original")
        (self.l2,) = self.ax.plot([], [], ".-", label="rotation")
        self.ax.legend()

        cur_cycle = int(cycle_slider_val)
        ignore_init_cycles = self.ignore_init_cycles
        baseline_end = None

        pcr_analyzer = PCRmachine(max_cycle=self.cycle_tt_num)

        # 일단 현재 cycle 까지 모두 입력
        for cycle, cq_value in enumerate(self.df_cq[sample_index][:cur_cycle]):
            pcr_analyzer.add_one_cq_value(
                cq_value, set_winsize=self.smoothing_winsize, start=self.df_cq.index[0]
            )
            pcr_analyzer.pcr_machine_cal_dCq()
        # 입력 후 플랏팅
        self.l.set_data(
            self.df_cq[sample_index].index[:cur_cycle],
            np.array(self.df_cq[sample_index][:cur_cycle]),
        )
        if cur_cycle > ignore_init_cycles:
            adj_cq_values = pcr_analyzer.rotate_cqs(
                cur_cycle, ignore_init_cycles=ignore_init_cycles, debugq=True
            )
            self.l2.set_data(
                adj_cq_values["cycle"][:cur_cycle],
                np.array(adj_cq_values["Cq_smoothing"][:cur_cycle]),
            )
        else:
            self.l2.set_data([], [])
        ## debug
        self.debug_print(
            sample_index, baseline_end, pcr_analyzer, self.max_baseline_end
        )
        # baseline end 찾기
        baseline_end = pcr_analyzer.find_baseline_end(
            ignore_init_cycles=ignore_init_cycles,
            ddcq_threshold_scale=self.ddcq_threshold_scal,
            max_cycle=self.max_baseline_end,
        )
        if baseline_end:  # baseline_end 를 중심으로 그래프 rotation
            baseline_end = min(baseline_end, self.max_baseline_end)
            adj_cq_values = pcr_analyzer.rotate_cqs(
                baseline_end, ignore_init_cycles=ignore_init_cycles, debugq=True
            )
            self.l2.set_data(
                adj_cq_values["cycle"][:cur_cycle],
                np.array(adj_cq_values["Cq_smoothing"][:cur_cycle]),
            )
            # debug
            self.debug_print(
                sample_index, baseline_end, pcr_analyzer, self.max_baseline_end
            )
        self._draw_arrow(sample_index, baseline_end)
        self._set_plot_axis(sample_index)

    def update_plot1(self, sample_index, cycle_slider_val):
        """Cq 값 입력될때마다 T 계산하고 baseline_end 찾기"""
        self.ax.clear()
        (self.l,) = self.ax.plot([], [], ".-", label="original")
        (self.l2,) = self.ax.plot([], [], ".-", label="rotation")
        self.ax.legend()

        cur_cycle = int(cycle_slider_val)
        ignore_init_cycles = self.ignore_init_cycles
        baseline_end = None

        pcr_analyzer = PCRmachine(max_cycle=self.cycle_tt_num)

        for cycle, cq_value in enumerate(self.df_cq[sample_index][:cur_cycle]):
            pcr_analyzer.add_one_cq_value(
                cq_value, set_winsize=self.smoothing_winsize, start=self.df_cq.index[0]
            )
            pcr_analyzer.pcr_machine_cal_dCq()
            baseline_end = pcr_analyzer.find_baseline_end(
                ddcq_threshold_scale=self.ddcq_threshold_scal,
                ignore_init_cycles=ignore_init_cycles,
            )
            self.l.set_data(
                self.df_cq[sample_index].index[:cur_cycle],
                np.array(self.df_cq[sample_index][:cur_cycle]),
            )
            if cycle > ignore_init_cycles:
                adj_cq_values = pcr_analyzer.rotate_cqs(
                    cycle + 1, ignore_init_cycles=ignore_init_cycles
                )
                self.l2.set_data(
                    adj_cq_values["cycle"][:cur_cycle],
                    np.array(adj_cq_values["Cq_smoothing"][:cur_cycle]),
                )
            else:
                self.l2.set_data([], [])
            if not baseline_end:
                continue
            baseline_end = min(baseline_end, self.max_baseline_end)
            break

        self._draw_arrow(sample_index, baseline_end)

        for cq_value in self.df_cq[sample_index][cycle + 1 : cur_cycle]:
            pcr_analyzer.add_one_cq_value(
                cq_value, set_winsize=self.smoothing_winsize, start=self.df_cq.index[0]
            )
            pcr_analyzer.pcr_machine_cal_dCq()
            adj_cq_values = pcr_analyzer.rotate_cqs(
                baseline_end, ignore_init_cycles=ignore_init_cycles
            )
            self.l2.set_data(
                adj_cq_values["cycle"][:cur_cycle],
                np.array(adj_cq_values["Cq_smoothing"][:cur_cycle]),
            )
        self._set_plot_axis(sample_index)

    def debug_print(self, sample_index, baseline_end, pcr_analyzer, max_cycle):
        ddcq_threshold = pcr_analyzer.cal_ddcq_threshold(
            scale_factor=self.ddcq_threshold_scal, max_cycle=max_cycle
        )
        print(
            f"sample: {sample_index}, baseline end(1-based): {baseline_end}, ddCq thr: {ddcq_threshold}"
        )
        df = pd.DataFrame(pcr_analyzer.get_cq_values())
        selected_cols = ["cycle", "Cq_ori", "Cq_smoothing", "dCq", "ddCq"]

        # Define column widths for alignment
        column_widths = {
            "cycle": max(5, len("cycle") + 1),
            "Cq_ori": max(13, len("Cq_ori") + 1),
            "Cq_smoothing": max(13, len("Cq_smoothing") + 1),
            "dCq": max(13, len("dCq") + 1),
            "ddCq": max(13, len("ddCq") + 1),
        }

        # Generate header with aligned columns
        header = "".join(f"{col:<{column_widths[col]}}" for col in selected_cols)
        print(header)

        # Iterate over rows
        for _, row in df.iterrows():
            ddcq = row["ddCq"]
            formatted_row = []
            for col in selected_cols:
                if col == "cycle":
                    # Cycle column is integer
                    formatted_row.append(f"{int(row[col]):<{column_widths[col]}}")
                elif col == "ddCq":
                    # Apply highlighting only to ddCq column
                    value = (
                        f"{row[col]:<{column_widths[col]}.2f}"
                        if isinstance(row[col], (int, float))
                        else f"{str(row[col]):<{column_widths[col]}}"
                    )
                    if ddcq is not None and ddcq >= ddcq_threshold:
                        formatted_row.append(
                            f"{Back.WHITE}{Fore.BLACK}{value}{Style.RESET_ALL}"
                        )
                    else:
                        formatted_row.append(value)
                else:
                    # Default formatting for other columns
                    formatted_row.append(
                        f"{row[col]:<{column_widths[col]}.2f}"
                        if isinstance(row[col], (int, float))
                        else f"{str(row[col]):<{column_widths[col]}}"
                    )

            # Print the formatted row
            print("".join(formatted_row))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PCR Analyzer")
    parser.add_argument("input_excelfile", type=str, help="Path to input Excel file")
    parser.add_argument("sheet_name", type=str, help="Sheet name in the Excel file")
    parser.add_argument(
        "--update_func",
        type=int,
        choices=[1, 2],
        default=1,
        help="Update function type (1 or 2)",
    )
    args = parser.parse_args()

    try:
        df_cq = read_excel(args.input_excelfile, args.sheet_name)
    except ValueError:
        df_cq = read_csv(args.input_excelfile)

    cycle_tt_num = df_cq.shape[0]

    app = QApplication(sys.argv)
    window = PCRAnalyzerApp(df_cq, cycle_tt_num, update_func_type=args.update_func - 1)
    window.show()
    sys.exit(app.exec_())
