import sys
import os
import json
import argparse
import tempfile
import time
from PyQt5.QtWidgets import (
    QApplication,
    QWidget,
    QLabel,
    QListView,
    QPushButton,
    QVBoxLayout,
    QHBoxLayout,
    QFrame,
    QFileDialog,
    QMenuBar,
    QAction,
    QSplitter,
    QLineEdit,
    QListWidget,
    QTextEdit,
    QTabWidget,
)
from PyQt5.QtCore import Qt, QStringListModel, QThread, pyqtSignal
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import time
import joblib
import numpy as np
from src.pcr_autolearner import AutoLearner
from src.pcr_anomaly_detector import PCRAnomalyDetector
from src.pcr_analyzer import PCRmachine
from src.pcr_autolearner_ai import PCRPredicter, PCRCluster, load_data, load_data_target
from PyQt5.QtWidgets import (
    QDialog,
    QFormLayout,
    QSpinBox,
    QDoubleSpinBox,
    QComboBox,
    QMessageBox,
    QDialogButtonBox,
)


class ClusteringThread(QThread):
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)

    def __init__(self, autolearner, params, parent=None):
        super().__init__(parent)
        self.autolearner = autolearner
        self.params = params

    def run(self):
        try:
            # Extract data
            targets = self.autolearner.get_target_names()
            input_data = []
            input_labels = []
            
            for target in targets:
                self.autolearner.set_target(target)
                self.autolearner.normalization(smoothing_size=3, mode="primitive")
                keys = list(self.autolearner.meta_info.keys())
                input_data.extend([self.autolearner.Cq_values[key] for key in keys])
                for key in keys:
                    timestamp = self.autolearner.gen_timestamps.get(key, "")
                    label = f"{target}_{timestamp}"
                    input_labels.append(label)
            
            input_data = np.array(input_data)
            if len(input_data) == 0:
                self.error.emit("No data found.")
                return

            clusterer = PCRCluster(
                n_clusters=self.params["n_clusters"],
                method=self.params["method"],
                early_weight=self.params["early_weight"]
            )
            
            normalized_data, valid_indices = clusterer.preprocess(input_data, skip_normalization=True)
            if len(normalized_data) == 0:
                self.error.emit("No valid data after preprocessing.")
                return

            labels_to_plot = [input_labels[i] for i in valid_indices]
            labels = clusterer.fit_predict(normalized_data)
            
            # Just return data and labels, visualization will be handled in the UI thread
            # because Matplotlib figures are not thread-safe for GUI backend interaction.
            self.finished.emit({
                "data": normalized_data,
                "labels": labels,
                "input_labels": labels_to_plot,
                "method": self.params["method"]
            })
        except Exception as e:
            self.error.emit(str(e))


class ClusteringResultWindow(QWidget):
    def __init__(self, cluster_results, parent=None):
        super().__init__(parent)
        self.cluster_results = cluster_results
        self.setWindowTitle("Clustering Results")
        self.setMinimumSize(900, 700)
        
        self.layout = QVBoxLayout(self)
        
        self.tabs = QTabWidget()
        self.layout.addWidget(self.tabs)
        
        # Initialize PCRCluster to use its visualization logic
        self.clusterer = PCRCluster(method=cluster_results["method"])
        self.figures = self.clusterer.visualize_clusters(
            cluster_results["data"], 
            cluster_results["labels"], 
            sample_labels=cluster_results["input_labels"]
        )
        
        # Tab 1: PCA Plot
        self.pca_tab = QWidget()
        pca_layout = QVBoxLayout(self.pca_tab)
        self.pca_canvas = FigureCanvas(self.figures['pca'])
        pca_layout.addWidget(self.pca_canvas)
        self.tabs.addTab(self.pca_tab, "PCA Scatter Plot")
        
        # Tab 2: Cluster Profiles
        self.profiles_tab = QWidget()
        profiles_main_layout = QVBoxLayout(self.profiles_tab)
        
        # Create a scroll area or just a simple splitter if few clusters
        # For simplicity, let's use a sub-tab widget for each cluster profile
        self.profile_tabs = QTabWidget()
        for label, fig in self.figures['profiles']:
            canvas = FigureCanvas(fig)
            tab = QWidget()
            tab_layout = QVBoxLayout(tab)
            tab_layout.addWidget(canvas)
            self.profile_tabs.addTab(tab, f"Cluster {label}")
            
        profiles_main_layout.addWidget(self.profile_tabs)
        self.tabs.addTab(self.profiles_tab, "Cluster Profiles")
        
        # Export Button
        self.export_btn = QPushButton("Export All Plots to Directory")
        self.export_btn.clicked.connect(self.export_plots)
        self.layout.addWidget(self.export_btn)

    def export_plots(self):
        save_dir = QFileDialog.getExistingDirectory(self, "Select Save Directory")
        if save_dir:
            try:
                self.clusterer.visualize_clusters(
                    self.cluster_results["data"],
                    self.cluster_results["labels"],
                    sample_labels=self.cluster_results["input_labels"],
                    save_dir=save_dir
                )
                QMessageBox.information(self, "Success", f"All plots saved to:\n{save_dir}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to save plots: {e}")


class AnomalyTrainingThread(QThread):
    finished = pyqtSignal(str)
    error = pyqtSignal(str)

    def __init__(self, autolearner, params):
        super().__init__()
        self.autolearner = autolearner
        self.params = params

    def run(self):
        try:
            data_dir = self.params["data_dir"]
            initial_model = self.params.get("initial_model", None)
            
            # Filtering logic for incremental training
            history_file = os.path.join(os.path.dirname(os.path.abspath(self.params["save_dir"])), "training_history.json")
            history = {}
            if os.path.exists(history_file):
                with open(history_file, 'r') as f:
                    history = json.load(f)

            norm_initial_model = os.path.abspath(initial_model) if initial_model else None
            files_to_exclude = []
            if norm_initial_model and norm_initial_model in history:
                files_to_exclude = history[norm_initial_model]
            
            # Use a fresh AutoLearner instance
            training_autolearner = AutoLearner()

            # Match the standalone CLI pattern (load_data handles filtering):
            input_rfus, input_labels, _, source_files = load_data(
                training_autolearner, 
                data_dir, 
                normalization_mode="primitive", 
                ddcq=False,
                exclude_files=files_to_exclude
            )

            seq_len = self.params["window_size"]
            X = []
            y = []
            processed_files = []
            for cq_values, label, src_file in zip(input_rfus, input_labels, source_files):
                if label == "Negative":
                    continue
                for i in range(len(cq_values) - seq_len):
                    X.append(cq_values[i : i + seq_len])
                    y.append(cq_values[i + seq_len])
                processed_files.append(src_file)

            if not X:
                if files_to_exclude:
                    self.error.emit("No new positive samples found (all may have been excluded).")
                else:
                    self.error.emit("No positive samples found in the data directory.")
                return

            print(f"Anomaly Incremental Training: Using {len(processed_files)} new files.")

            X = np.array(X).reshape(-1, seq_len, 1)
            y = np.array(y).reshape(-1, 1)

            # Build and train detector
            detector = PCRAnomalyDetector(seq_len=seq_len)
            model_save_path = os.path.abspath(os.path.join(self.params["save_dir"], "anomaly_detector.keras"))
            
            detector.fit(
                X, y, 
                epochs=self.params["epochs"], 
                model_path=model_save_path, 
                model_type=self.params["model_type"],
                initial_model_path=initial_model if initial_model and os.path.exists(initial_model) else None
            )

            # Export to TFLite float16
            tflite_path = os.path.abspath(os.path.join(self.params["save_dir"], "anomaly_detector_float16.tflite"))
            detector.export_tflite_float16(tflite_path)

            # Update history
            combined_history = list(set(files_to_exclude + processed_files))
            history[model_save_path] = combined_history
            history[tflite_path] = combined_history
            
            os.makedirs(os.path.dirname(history_file), exist_ok=True)
            with open(history_file, 'w') as f:
                json.dump(history, f, indent=4)

            self.finished.emit(tflite_path)
        except Exception as e:
            self.error.emit(str(e))


class AnomalyTrainingDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Train Anomaly Detection Model")
        self.setFixedWidth(500)
        
        layout = QFormLayout(self)
        
        self.data_dir_input = QLineEdit()
        data_dir_layout = QHBoxLayout()
        data_dir_layout.addWidget(self.data_dir_input)
        browse_btn = QPushButton("Browse")
        browse_btn.clicked.connect(self.browse_data_dir)
        data_dir_layout.addWidget(browse_btn)
        layout.addRow("Data Directory:", data_dir_layout)
        
        self.model_type_combo = QComboBox()
        self.model_type_combo.addItems(["lstm", "gru", "rnn"])
        self.model_type_combo.setCurrentText("gru")
        layout.addRow("Model Type:", self.model_type_combo)
        
        self.window_size_input = QSpinBox()
        self.window_size_input.setRange(2, 10)
        self.window_size_input.setValue(3)
        layout.addRow("Window Size:", self.window_size_input)
        
        self.epochs_input = QSpinBox()
        self.epochs_input.setRange(1, 1000)
        self.epochs_input.setValue(100)
        layout.addRow("Epochs:", self.epochs_input)
        
        self.save_dir_input = QLineEdit()
        save_dir_layout = QHBoxLayout()
        save_dir_layout.addWidget(self.save_dir_input)
        save_browse_btn = QPushButton("Browse")
        save_browse_btn.clicked.connect(self.browse_save_dir)
        save_dir_layout.addWidget(save_browse_btn)
        self.save_dir_input.setText(os.path.join(os.getcwd(), "models"))
        layout.addRow("Save Directory:", save_dir_layout)
        
        self.initial_model_input = QLineEdit()
        self.initial_model_input.setPlaceholderText("Optional: Existing Keras model to continue training")
        initial_model_layout = QHBoxLayout()
        initial_model_layout.addWidget(self.initial_model_input)
        initial_browse_btn = QPushButton("Browse")
        initial_browse_btn.clicked.connect(self.browse_initial_model)
        initial_model_layout.addWidget(initial_browse_btn)
        layout.addRow("Initial Model:", initial_model_layout)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addRow(buttons)

    def browse_data_dir(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Data Directory")
        if dir_path:
            self.data_dir_input.setText(dir_path)

    def browse_save_dir(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Save Directory")
        if dir_path:
            self.save_dir_input.setText(dir_path)

    def browse_initial_model(self):
        model_path, _ = QFileDialog.getOpenFileName(
            self, "Select Initial Anomaly Model", "", "Keras Models (*.keras)"
        )
        if model_path:
            self.initial_model_input.setText(model_path)

    def get_params(self):
        return {
            "data_dir": self.data_dir_input.text(),
            "model_type": self.model_type_combo.currentText(),
            "window_size": self.window_size_input.value(),
            "epochs": self.epochs_input.value(),
            "save_dir": self.save_dir_input.text(),
            "initial_model": self.initial_model_input.text()
        }


class ClusteringDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Clustering Parameters")
        layout = QFormLayout(self)

        self.n_clusters = QSpinBox()
        self.n_clusters.setRange(2, 10)
        self.n_clusters.setValue(3)
        layout.addRow("Number of Clusters:", self.n_clusters)

        self.method = QComboBox()
        self.method.addItems(["kmeans", "gmm", "bgmm"])
        layout.addRow("Clustering Method:", self.method)

        self.early_weight = QDoubleSpinBox()
        self.early_weight.setRange(0.1, 10.0)
        self.early_weight.setSingleStep(0.1)
        self.early_weight.setValue(1.0)
        layout.addRow("Early Cycle Weight:", self.early_weight)

        self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, self)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        layout.addRow(self.buttons)

    def get_params(self):
        return {
            "n_clusters": self.n_clusters.value(),
            "method": self.method.currentText(),
            "early_weight": self.early_weight.value()
        }


class TrainingThread(QThread):
    finished = pyqtSignal(str)
    error = pyqtSignal(str)

    def __init__(self, autolearner, params, parent=None):
        super().__init__(parent)
        self.autolearner = autolearner
        self.params = params

    def run(self):
        try:
            data_dir = self.params["data_dir"]
            initial_model = self.params.get("initial_model", None)
            
            # Filtering logic for incremental training
            history_file = os.path.join(os.path.dirname(os.path.abspath(self.params["save_dir"])), "training_history.json")
            history = {}
            if os.path.exists(history_file):
                with open(history_file, 'r') as f:
                    history = json.load(f)

            norm_initial_model = os.path.abspath(initial_model) if initial_model else None
            files_to_exclude = []
            if norm_initial_model and norm_initial_model in history:
                files_to_exclude = history[norm_initial_model]
            
            # Use a fresh AutoLearner instance for training to avoid shared state issues
            training_autolearner = AutoLearner()
            
            # Use load_data from pcr_autolearner_ai with pre-filtering
            input_data, input_labels, _, source_files = load_data(
                training_autolearner,
                data_dir,
                normalization_mode=self.params["normalization_mode"],
                ddcq=self.params["ddcq"],
                exclude_files=files_to_exclude
            )
            
            if len(input_data) == 0:
                self.error.emit("No new training data found in the selected directory (all files may have been excluded).")
                return

            print(f"Incremental Training: Using {len(input_data)} new samples.")
            
            train_data = input_data
            train_labels = input_labels
            used_new_files = list(set(source_files))

            pcr_learner = PCRPredicter()
            pcr_learner.train_model(
                train_data,
                train_labels,
                learning_rate=self.params["learning_rate"],
                epochs=self.params["epochs"],
                save_dir=self.params["save_dir"],
                ai_model=self.params["ai_model"],
                initial_model_path=initial_model if initial_model and os.path.exists(initial_model) else None
            )
            
            # Update history with absolute paths
            model_ext = ".keras" if self.params['ai_model'] == 'cnn' else ".pkl"
            new_model_path = os.path.abspath(os.path.join(self.params["save_dir"], f"best_model_{self.params['ai_model']}{model_ext}"))
            
            # Combine history from initial model if it exists
            combined_history = list(set(files_to_exclude + used_new_files))
            history[new_model_path] = combined_history
            
            os.makedirs(os.path.dirname(history_file), exist_ok=True)
            with open(history_file, 'w') as f:
                json.dump(history, f, indent=4)
            
            self.finished.emit(self.params["save_dir"])
        except Exception as e:
            self.error.emit(str(e))


class TrainingDialog(QDialog):
    def __init__(self, current_dir="", parent=None):
        super().__init__(parent)
        self.setWindowTitle("AI Model Training Parameters")
        self.setMinimumWidth(400)
        layout = QFormLayout(self)

        self.data_dir_input = QLineEdit(current_dir)
        self.browse_btn = QPushButton("Browse")
        self.browse_btn.clicked.connect(self.browse_directory)
        dir_layout = QHBoxLayout()
        dir_layout.addWidget(self.data_dir_input)
        dir_layout.addWidget(self.browse_btn)
        layout.addRow("Data Directory:", dir_layout)

        self.ai_model = QComboBox()
        self.ai_model.addItems(["cnn", "rf"])
        layout.addRow("AI Model Type:", self.ai_model)

        self.learning_rate = QDoubleSpinBox()
        self.learning_rate.setRange(0.0001, 0.1)
        self.learning_rate.setDecimals(4)
        self.learning_rate.setSingleStep(0.001)
        self.learning_rate.setValue(0.001)
        layout.addRow("Learning Rate:", self.learning_rate)

        self.epochs = QSpinBox()
        self.epochs.setRange(1, 1000)
        self.epochs.setValue(50)
        layout.addRow("Epochs:", self.epochs)

        self.normalization_mode = QComboBox()
        self.normalization_mode.addItems(["rotation", "primitive"])
        layout.addRow("Normalization Mode:", self.normalization_mode)

        self.ddcq = QComboBox()
        self.ddcq.addItems(["False", "True"])
        layout.addRow("Apply ddcq:", self.ddcq)

        self.save_dir_input = QLineEdit("models/")
        layout.addRow("Save Directory:", self.save_dir_input)

        self.initial_model_input = QLineEdit()
        self.initial_model_input.setPlaceholderText("Optional: Select an existing model to continue training")
        initial_model_layout = QHBoxLayout()
        initial_model_layout.addWidget(self.initial_model_input)
        initial_browse_btn = QPushButton("Browse")
        initial_browse_btn.clicked.connect(self.browse_initial_model)
        initial_model_layout.addWidget(initial_browse_btn)
        layout.addRow("Initial Model:", initial_model_layout)

        self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, self)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        layout.addRow(self.buttons)

    def browse_directory(self):
        directory = QFileDialog.getExistingDirectory(self, "Select Data Directory")
        if directory:
            self.data_dir_input.setText(directory)

    def browse_initial_model(self):
        model_path, _ = QFileDialog.getOpenFileName(
            self, "Select Initial Model File", "", "Model Files (*.keras *.pkl)"
        )
        if model_path:
            self.initial_model_input.setText(model_path)

    def get_params(self):
        return {
            "data_dir": self.data_dir_input.text(),
            "ai_model": self.ai_model.currentText(),
            "learning_rate": self.learning_rate.value(),
            "epochs": self.epochs.value(),
            "normalization_mode": self.normalization_mode.currentText(),
            "ddcq": self.ddcq.currentText() == "True",
            "save_dir": self.save_dir_input.text(),
            "initial_model": self.initial_model_input.text()
        }


class ParameterSearchThread(QThread):
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)

    def __init__(self, autolearner, parent=None):
        super().__init__(parent)
        self.autolearner = autolearner

    def run(self):
        try:
            start_time = time.time()
            scale_res = self.autolearner.iteration_training_scaling_factor(
                cycle_range=(3, 27), epsilon=1e-4, max_iteration=100, step_size=1e-3
            )
            scaling_factor_min, scaling_factor_max = scale_res["scaling_factor"]
            scaling_loss_min, scaling_loss_max = scale_res["loss_value"]
            cycle_res = self.autolearner.iteration_training_cycle_range(
                (scaling_factor_min, scaling_factor_max)
            )
            cycle_min, cycle_max = cycle_res["cycle_range"]
            cycle_loss_min, cycle_loss_max = cycle_res["loss_value"]
            end_time = time.time()
            print(f"Parameter search time: {end_time - start_time:.2f} seconds")
            self.finished.emit(
                {
                    "scaling_factor_min": scaling_factor_min,
                    "scaling_factor_max": scaling_factor_max,
                    "scaling_loss_min": scaling_loss_min,
                    "scaling_loss_max": scaling_loss_max,
                    "cycle_min": cycle_min,
                    "cycle_max": cycle_max,
                    "cycle_loss_min": cycle_loss_min,
                    "cycle_loss_max": cycle_loss_max,
                }
            )
        except Exception as e:
            self.error.emit(str(e))


class UIWindow(QWidget):
    def __init__(self, anomaly_args: dict = None):
        super().__init__()
        self.initUI()
        self.autolearner = AutoLearner()
        # 인스턴스명과 메서드명이 충돌하지 않도록 이름 변경
        self.anomaly_detector = PCRAnomalyDetector(seq_len=3)
        self.target_name = ""
        self.ai_model = None  # AI 모델을 저장할 변수
        self.ai_model_type = None  # AI 모델 유형 (cnn 또는 rf)
        self.baseline_substracted = (
            False  # baseline substraction이 되었는지 여부를 추적
        )
        self.normalized = False  # 데이터가 normalization되었는지 여부를 추적
        self.anomaly_args = anomaly_args

    def initUI(self):
        main_layout = QVBoxLayout()
        # 메뉴바
        self.menu_bar = self.create_menu_bar()
        main_layout.setMenuBar(self.menu_bar)
        # AI 모델 및 라벨
        main_layout.addLayout(self.create_ai_model_layout())
        # Y축 범위
        main_layout.addLayout(self.create_yscale_layout())
        # 이하 기존 splitter/프레임 조립 코드 유지

        # Vertical QSplitter to manage height
        vertical_splitter = QSplitter(Qt.Vertical)

        # Horizontal QSplitter to manage width
        horizontal_splitter = QSplitter(Qt.Horizontal)

        ## 타겟 이름 영역
        target_frame, target_layout = self.add_target_frame()
        self.initial_target_list(target_frame, target_layout)
        horizontal_splitter.addWidget(target_frame)

        ## 그래프 영역 (Matplotlib 그래프)
        graph_frame = self.add_graph_frame()
        horizontal_splitter.addWidget(graph_frame)

        # Set stretch factors for horizontal splitter
        # 라벨 리스트 박스 (오른쪽)
        self.label_list_widget = QListWidget()
        self.label_list_widget.itemClicked.connect(self.on_label_box_clicked)  # 추가
        horizontal_splitter.addWidget(self.label_list_widget)
        horizontal_splitter.setStretchFactor(0, 1)
        horizontal_splitter.setStretchFactor(1, 3)
        horizontal_splitter.setStretchFactor(2, 1)

        # Add horizontal splitter to vertical splitter
        vertical_splitter.addWidget(horizontal_splitter)

        # Add parameter search frame to the bottom of the vertical splitter
        param_layout_frame = QFrame()
        param_layout = self.add_paramsearch_frame()
        param_layout_frame.setLayout(param_layout)
        vertical_splitter.addWidget(param_layout_frame)

        # Set stretch factors for vertical splitter
        vertical_splitter.setStretchFactor(0, 4)  # Target name + Plot 영역
        vertical_splitter.setStretchFactor(1, 1)  # Parameter search 영역

        # Add vertical splitter to the main layout
        main_layout.addWidget(vertical_splitter)

        # 계산 중 어두운 화면 추가
        self.add_blocklayer()

        # UI 최종 마무리
        self.setLayout(main_layout)
        self.setWindowTitle("PCR parameter autolearner v0.0")
        self.setGeometry(300, 300, 800, 600)
        self.show()

    def create_menu_bar(self):
        menu_bar = QMenuBar(self)
        # File 메뉴
        file_menu = menu_bar.addMenu("File")
        open_dir_action = QAction("Open Directory", self)
        open_dir_action.triggered.connect(self.open_directory)
        file_menu.addAction(open_dir_action)

        # Calculation 메뉴
        calc_menu = menu_bar.addMenu("Calculation")
        logistic_action = QAction("Logistic Fitting", self)
        logistic_action.triggered.connect(self.logistic_fitting)
        calc_menu.addAction(logistic_action)

        param_search_action = QAction("Parameter Search", self)
        param_search_action.triggered.connect(self.parameter_search)
        calc_menu.addAction(param_search_action)

        # AI 메뉴
        ai_menu = menu_bar.addMenu("AI")

        anomaly_detect_action = QAction("Anomaly Detection", self)
        anomaly_detect_action.triggered.connect(self.anomaly_detection)
        ai_menu.addAction(anomaly_detect_action)
        
        train_anomaly_action = QAction("Train Anomaly Detection Model", self)
        train_anomaly_action.triggered.connect(self.train_anomaly_detection)
        ai_menu.addAction(train_anomaly_action)

        cluster_action = QAction("Cluster RFU Patterns", self)
        cluster_action.triggered.connect(self.cluster_data)
        ai_menu.addAction(cluster_action)
        
        train_action = QAction("Train Pos/Neg AI Model", self)
        train_action.triggered.connect(self.train_ai_model)
        ai_menu.addAction(train_action)

        import_pred_action = QAction("Import Pos/Neg AI Model", self)
        import_pred_action.triggered.connect(self.import_ai_model)
        ai_menu.addAction(import_pred_action)

        return menu_bar

    def create_button_layout(self):
        button_layout = QHBoxLayout()
        basesub_button = self.add_basesub_button()
        button_layout.addWidget(basesub_button)
        normalize_button = self.add_normalize_button()
        button_layout.addWidget(normalize_button)
        return button_layout

    def create_ai_model_layout(self):
        ai_model_layout = QHBoxLayout()
        # Baseline Sub / Normalize buttons added back here
        basesub_button = self.add_basesub_button()
        ai_model_layout.addWidget(basesub_button)
        normalize_button = self.add_normalize_button()
        ai_model_layout.addWidget(normalize_button)

        self.ai_model_label = QLabel("No Pos/Neg AI Model Loaded")
        self.ai_model_label.setStyleSheet("color: red; font-weight: bold;")
        ai_model_layout.addWidget(self.ai_model_label)
        self.normalization_warning = QLabel("Pos/Neg AI Model works only with normalized data.")
        self.normalization_warning.setStyleSheet("color: orange; font-style: italic;")
        ai_model_layout.addWidget(self.normalization_warning)
        return ai_model_layout

    def create_yscale_layout(self):
        yscale_layout = QHBoxLayout()
        yscale_label = QLabel("Y-Axis Range:")
        yscale_layout.addWidget(yscale_label)
        self.yscale_min_input = QLineEdit()
        self.yscale_min_input.setPlaceholderText("Min")
        yscale_layout.addWidget(self.yscale_min_input)
        self.yscale_max_input = QLineEdit()
        self.yscale_max_input.setPlaceholderText("Max")
        yscale_layout.addWidget(self.yscale_max_input)
        yscale_button = QPushButton("Set Y-Axis Range")
        yscale_button.clicked.connect(self.set_yscale_range)
        yscale_layout.addWidget(yscale_button)
        return yscale_layout

    def add_directory_button(self):
        # 디렉토리 라벨
        directory_button = QPushButton("Open directory: read PCR output files")
        directory_button.clicked.connect(self.open_directory)  # 클릭 이벤트 연결
        return directory_button

    def add_basesub_button(self):
        """baseline substraction"""
        basesub_button = QPushButton("Baseline Sub.")
        basesub_button.clicked.connect(self.baseline_substraction)
        return basesub_button

    def add_normalize_button(self):
        """normalize"""
        normalize_button = QPushButton("Normalize")
        normalize_button.clicked.connect(self.normalize)
        return normalize_button

    def add_target_frame(self):
        # 타겟 이름 영역
        target_frame = QFrame()
        target_frame.setFrameShape(QFrame.Box)
        target_layout = QVBoxLayout()

        target_name_label = QLabel("Target Name")
        target_layout.addWidget(target_name_label)

        self.target_list = QListView()
        self.target_model = QStringListModel()
        self.target_list.setModel(self.target_model)
        target_layout.addWidget(self.target_list)

        return target_frame, target_layout

    def add_graph_frame(self):
        graph_frame = QFrame()
        graph_frame.setFrameShape(QFrame.Box)
        graph_layout = QVBoxLayout()
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        graph_layout.addWidget(self.canvas)
        graph_frame.setLayout(graph_layout)
        return graph_frame

    def add_paramsearch_frame(self):
        param_layout = QVBoxLayout()

        param_header = QLabel("Parameter Search Results:")
        param_header.setStyleSheet("font-weight: bold;")
        param_layout.addWidget(param_header)

        self.param_label = QTextEdit()
        self.param_label.setReadOnly(True)
        self.param_label.setPlaceholderText("Scaling factor and cycle range results will appear here...")
        param_layout.addWidget(self.param_label)
        return param_layout

    def add_blocklayer(self):
        self.overlay = QWidget(self)
        self.overlay.setParent(self)  # 부모 위젯 명시적으로 설정
        self.overlay.setLayout(QVBoxLayout())  # 레이아웃 설정
        self.overlay.setStyleSheet("background-color: rgba(0, 0, 0, 150);")
        self.overlay.hide()
        self.overlay.lower()

    def initial_target_list(self, target_frame, target_layout):
        # 초기 문자열 리스트 설정
        self.set_initial_target_list(
            [
                "Here will be",
                " displayed target names",
                ", after opening",
                " the csv directory",
            ]
        )
        # 리스트 아이템 클릭 이벤트 연결
        self.target_list.clicked.connect(self.target_item_clicked)
        target_frame.setLayout(target_layout)

    def set_initial_target_list(self, items):
        # 초기 문자열 리스트 설정
        self.target_model.setStringList(items)

    def set_yscale_range(self):
        """Set the y-axis range based on user input."""
        try:
            y_min = float(self.yscale_min_input.text())
            y_max = float(self.yscale_max_input.text())

            if y_min >= y_max:
                raise ValueError("Minimum value must be less than maximum value.")

            if self.figure.axes:  # Ensure there is at least one subplot
                ax = self.figure.axes[0]
                ax.set_ylim(y_min, y_max)  # Update the y-axis range
                self.canvas.draw()  # Redraw the canvas to apply changes
        except ValueError as e:
            print(f"Invalid input for y-axis range: {e}")

    def target_item_clicked(self, index):
        # 리스트 아이템 클릭 이벤트 처리
        self.target_name = self.target_model.stringList()[index.row()]
        try:
            self.autolearner.set_target(self.target_name)
        except ValueError:
            print("No target name selected.")
            return
        print(f"Clicked target_name: {self.target_name}")
        self.baseline_substracted = False  # baseline substraction이 되지 않았음을 표시
        self.normalized = False  # 데이터가 normalization되지 않았음을 표시
        self.update_graph(self.target_name)  # 그래프 업데이트

    def open_directory(self):
        """Open a directory and load PCR output files."""
        # 디렉토리 선택 다이얼로그 열기
        directory = QFileDialog.getExistingDirectory(self, "Select Directory")
        if directory:
            # 기존 데이터 초기화
            self.autolearner.reset()  # AutoLearner의 데이터를 초기화하는 메서드가 있다고 가정
            self.target_model.setStringList([])  # 타겟 리스트 초기화
            self.lines = {}  # 그래프 라인 초기화
            self.figure.clear()  # 그래프 초기화
            self.canvas.draw()  # 캔버스 업데이트

            # 새 데이터를 불러오기
            print(f"Selected directory: {directory}")
            self.update_target_list(directory)

    def update_target_list(self, directory):
        # 선택된 디렉토리의 파일 목록을 읽어 target_list에 아이템 추가
        self.autolearner.read_dir(directory)
        target_names = self.autolearner.get_target_names()
        target_names.sort()
        self.set_initial_target_list(target_names)
        print(f"abnormal files: {self.autolearner.abnormal_files}")

    def update_graph(self, target_name, anomalies_dict=None):
        """Update the graph with the current target data."""
        # Clear the existing figure
        self.figure.clear()
        ax = self.figure.add_subplot(111)

        # Store line objects for toggling visibility
        self.lines = {}
        self.anomaly_markers = {}
        self.label_list_widget.clear()  # 라벨 리스트 박스 초기화

        # Plot each dataset
        for kkk in self.autolearner.Cq_values:
            line = self._plot_dataset(ax, kkk)

            # Add AI model predictions if applicable
            if self.ai_model and self.normalized:
                self._add_ai_predictions(line, kkk)
            else:
                line_key = f"{kkk}_{self.autolearner.gen_timestamps[kkk]}"
                self.lines[line_key] = line
                self.label_list_widget.addItem(line_key)  # 라벨을 리스트에 추가

            # anomalies_dict가 있으면 해당 위치에 X표시 추가
            label = f"{kkk}_{self.autolearner.gen_timestamps[kkk]}"
            y = self.autolearner.Cq_values[kkk]
            x = self.autolearner.cycles
            marker_objs = []
            if anomalies_dict and label in anomalies_dict:
                for idx in anomalies_dict[label]:
                    (marker_obj,) = ax.plot(
                        x[idx],
                        y[idx],
                        marker="x",
                        color="black",
                        markersize=10,
                        linestyle="None",
                    )
                    marker_objs.append(marker_obj)
            self.anomaly_markers[label] = marker_objs

        # Add legend and configure interactivity
        self._configure_legend(ax)

        # Set the title and redraw the canvas
        ax.set_title(f"{target_name}")
        self.canvas.draw()

    def _plot_dataset(self, ax, kkk):
        """Plot a single dataset and return the line object."""
        (line,) = ax.plot(
            self.autolearner.cycles,
            self.autolearner.Cq_values[kkk],
            ".-",
            # label=str(kkk),  # Ensure the label is a string
            label=f"{kkk}_{self.autolearner.gen_timestamps[kkk]}",
        )
        return line

    def _add_ai_predictions(self, line, kkk):
        """Add AI model predictions and ground truth to the legend."""
        sample_data = np.array(self.autolearner.Cq_values[kkk]).reshape(1, -1)
        prediction = self._predict_with_ai_model(sample_data)

        # Get the ground truth (actual value) from meta_info
        ground_truth = self.autolearner.meta_info.get(kkk, {}).get("Pos/Neg", "Unknown")
        ground_truth = "P" if ground_truth == "Positive" else "N"
        # Determine if the prediction is positive or negative
        is_positive = (
            prediction == 1 if np.isscalar(prediction) else prediction.item() == 1
        )
        predicted_label = "P" if is_positive else "N"

        # Format the legend label with ground truth and prediction
        line_key = f"{kkk}_{self.autolearner.gen_timestamps[kkk]}"
        label = f"{line_key} (tr_{ground_truth}:pr_{predicted_label})"
        line.set_label(label)
        self.lines[label] = line  # Use the legend label as the key
        self.label_list_widget.addItem(label)  # 라벨을 리스트에 추가

        # ground_truth와 predicted_label이 다르면 붉은색으로 표시
        if ground_truth != predicted_label:
            # QListWidget.addItem은 None을 반환하므로, 마지막 아이템을 가져와서 색상 지정
            last_row = self.label_list_widget.count() - 1
            item_widget = self.label_list_widget.item(last_row)
            item_widget.setForeground(Qt.red)

    def _predict_with_ai_model(self, sample_data):
        """Make predictions using the loaded AI model."""
        if self.ai_model_type == "rf":  # Random Forest 모델
            return self.ai_model.predict(sample_data)[0]
        elif self.ai_model_type == "cnn":  # CNN 모델
            sample_data = sample_data.reshape(
                (sample_data.shape[0], sample_data.shape[1], 1)
            )
            return np.argmax(self.ai_model.predict(sample_data), axis=1)[0]
        else:
            raise ValueError("Unknown AI model type.")

    def _configure_legend(self, ax):
        # """Configure the legend and interactivity."""
        # # Calculate the number of columns for the legend
        # num_items = len(self.lines)
        # max_items_per_column = 30
        # num_columns = (num_items // max_items_per_column) + (
        #     1 if num_items % max_items_per_column > 0 else 0
        # )

        # # Create the legend
        # legend = ax.legend(
        #     loc="center left",  # Position the legend to the left of the bounding box
        #     bbox_to_anchor=(1.0, 0.5),  # Move the legend outside the plot on the right
        #     ncol=num_columns,  # Dynamically set the number of columns
        # )
        # for legend_line in legend.get_lines():
        #     legend_line.set_picker(10)  # Set picker sensitivity to 10 pixels

        # # Adjust the layout to make space for the legend
        # self.figure.subplots_adjust(left=0.1, right=0.75, top=0.9, bottom=0.1)

        # Connect the pick event to toggle visibility
        self.canvas.mpl_connect("pick_event", self.on_legend_click)

        # Connect the button press event to toggle visibility of all lines
        self.canvas.mpl_connect("button_press_event", self.on_plot_click)

        # Connect the motion notify event to display tooltips
        self.canvas.mpl_connect("motion_notify_event", self.on_plot_hover)

    def on_plot_click(self, event):
        """Toggle the visibility of all lines when the plot is clicked."""
        # Check if the event occurred inside the plot area
        if event.inaxes is not None:
            # Determine the current visibility state of the first line
            if self.lines:
                first_line = next(iter(self.lines.values()))
                new_visibility = not first_line.get_visible()

                # Toggle visibility for all lines
                for line in self.lines.values():
                    line.set_visible(new_visibility)

                # anomaly_marker도 같이 토글
                if hasattr(self, "anomaly_markers"):
                    for marker_list in self.anomaly_markers.values():
                        for marker in marker_list:
                            marker.set_visible(new_visibility)

                # Update the legend to reflect the visibility state
                legend = event.inaxes.get_legend()
                if legend:
                    for legend_line in legend.get_lines():
                        legend_line.set_alpha(1.0 if new_visibility else 0.2)
                # Redraw the canvas to apply changes
                self.canvas.draw()

    def on_legend_click(self, event):
        """Toggle the visibility of a plot when its legend label is clicked.
        deprecated, use on_label_box_clicked instead"""
        legend_line = event.artist
        label = legend_line.get_label()

        # Ensure the label matches the keys in self.lines
        if label in self.lines:
            # Toggle the visibility of the corresponding plot line
            line = self.lines[label]
            visible = not line.get_visible()
            line.set_visible(visible)

            # Update the legend to reflect the visibility state
            legend_line.set_alpha(1.0 if visible else 0.2)

            # Redraw the canvas to apply changes
            self.canvas.draw()
        else:
            print(f"Warning: Label '{label}' not found in self.lines")

    def on_plot_hover(self, event):
        """Display the label of a plot line when the mouse hovers over it."""
        if event.inaxes is not None:  # Ensure the event occurs inside the plot area
            for label, line in self.lines.items():
                if line.contains(event)[0]:  # Check if the mouse is near the line
                    # Display the label as a tooltip
                    self.canvas.setToolTip(label)
                    return
            # If no line is hovered, clear the tooltip
            self.canvas.setToolTip("")

    def on_label_box_clicked(self, item):
        """라벨 박스에서 라벨 클릭 시 해당 라인의 표시/숨김 토글"""
        label = item.text()
        if label in self.lines:
            line = self.lines[label]
            visible = not line.get_visible()
            line.set_visible(visible)
            if hasattr(self, "anomaly_markers") and label in self.anomaly_markers:
                for marker in self.anomaly_markers[label]:
                    marker.set_visible(visible)
            # 선택된 항목의 배경색으로 상태 표시(선택 시 회색, 숨김 시 연한 회색)
            if visible:
                item.setBackground(Qt.white)
            else:
                item.setBackground(Qt.lightGray)
                item.setSelected(False)  # 선택 해제
            self.canvas.draw()
        else:
            print(f"Label '{label}' not found in self.lines")

    def parameter_search(self):
        self.param_label.setText("Auto learning...")
        self.overlay.setGeometry(0, 0, self.width(), self.height())
        self.overlay.show()
        QApplication.processEvents()  # 즉시 블록 레이어 작동
        print("Overlay shown")
        self.overlay.raise_()

        self.search_thread = ParameterSearchThread(self.autolearner)
        self.search_thread.finished.connect(self.on_search_finished)
        self.search_thread.error.connect(self.on_search_error)
        self.search_thread.start()

    def on_search_finished(self, result):
        scaling_factor_min = result["scaling_factor_min"]
        scaling_factor_max = result["scaling_factor_max"]
        scaling_loss_min = result["scaling_loss_min"]
        scaling_loss_max = result["scaling_loss_max"]
        cycle_min = result["cycle_min"]
        cycle_max = result["cycle_max"]
        cycle_loss_min = result["cycle_loss_min"]
        cycle_loss_max = result["cycle_loss_max"]
        output_str = (
            f"Scaling Factor Range: {scaling_factor_min:.4f} ~ {scaling_factor_max:.4f}\n"
            f"Scaling Loss: {scaling_loss_min:.4f} ~ {scaling_loss_max:.4f}\n"
            f"Cycle Range: {cycle_min} ~ {cycle_max}\n"
            f"Cycle Loss: {cycle_loss_min:.4f} ~ {cycle_loss_max:.4f}"
        )
        self.param_label.setText(output_str)
        print(output_str)
        self.overlay.hide()
        print("Overlay hidden")

    def on_search_error(self, message):
        print(f"Error during parameter search: {message}")
        self.param_label.setText(f"Error: {message}")
        self.overlay.hide()
        print("Overlay hidden")

    def resizeEvent(self, event):
        print(f"Window size: {self.geometry()}")  # 창 크기 확인
        self.overlay.setGeometry(
            0, 0, self.width(), self.height()
        )  # overlay 크기 명시적으로 설정
        super().resizeEvent(event)

    def baseline_substraction(self):
        # 현재 선택된 target name을 가지고 진행
        try:
            self.autolearner.set_target(self.target_name)
        except ValueError:
            print("No target name selected.")
            return
        self.autolearner.baseline_substraction(mode="rotation")
        self.baseline_substracted = True  # baseline substraction만 되었음을 표시
        self.normalized = False  # 데이터가 normalization되지 않았음을 표시
        self.update_graph(self.target_name)  # 그래프 업데이트

    def normalize(self):
        # 현재 선택된 target name을 가지고 진행
        self.autolearner.set_target(self.target_name)
        self.autolearner.normalization(mode="rotation")
        self.baseline_substracted = False
        self.normalized = True  # 데이터가 normalization되었음을 표시
        self.update_graph(self.target_name)

    def cluster_data(self):
        """Perform clustering on current data."""
        dialog = ClusteringDialog(self)
        if dialog.exec_():
            params = dialog.get_params()
            self.param_label.setText("Clustering data...")
            self.overlay.setGeometry(0, 0, self.width(), self.height())
            self.overlay.show()
            self.overlay.raise_()
            
            self.cluster_thread = ClusteringThread(self.autolearner, params)
            self.cluster_thread.finished.connect(self.on_clustering_finished)
            self.cluster_thread.error.connect(self.on_search_error)
            self.cluster_thread.start()

    def on_clustering_finished(self, result):
        self.overlay.hide()
        self.param_label.setText("Clustering finished.")
        
        # open the result window
        self.result_window = ClusteringResultWindow(result)
        self.result_window.show()

    def train_anomaly_detection(self):
        """이상치 탐지 모델 학습"""
        dialog = AnomalyTrainingDialog(self)
        if dialog.exec_() == QDialog.Accepted:
            params = dialog.get_params()
            if not params["data_dir"]:
                QMessageBox.warning(self, "Warning", "Please select a data directory.")
                return

            self.overlay.show()
            self.param_label.setText("Training Anomaly Detection Model...")
            
            self.anomaly_train_thread = AnomalyTrainingThread(self.autolearner, params)
            self.anomaly_train_thread.finished.connect(self.on_anomaly_training_finished)
            self.anomaly_train_thread.error.connect(self.on_search_error)
            self.anomaly_train_thread.start()

    def on_anomaly_training_finished(self, tflite_path):
        self.overlay.hide()
        self.param_label.setText(f"Anomaly detection training finished. Exported to: {tflite_path}")
        QMessageBox.information(
            self, 
            "Training Complete", 
            f"Anomaly detection model training and TFLite export complete.\n\nSaved to: {tflite_path}"
        )

    def train_ai_model(self):
        """AI 모델 학습"""
        # 현재 로드된 디렉토리가 있으면 기본값으로 사용
        current_dir = ""
        if hasattr(self.autolearner, "input_dir"):
             current_dir = self.autolearner.input_dir
        
        dialog = TrainingDialog(current_dir, self)
        if dialog.exec_():
            params = dialog.get_params()
            if not params["data_dir"]:
                QMessageBox.warning(self, "Warning", "Data directory is required.")
                return

            self.param_label.setText("Training AI model...")
            self.overlay.setGeometry(0, 0, self.width(), self.height())
            self.overlay.show()
            self.overlay.raise_()
            
            self.train_thread = TrainingThread(self.autolearner, params)
            self.train_thread.finished.connect(self.on_training_finished)
            self.train_thread.error.connect(self.on_search_error)
            self.train_thread.start()

    def on_training_finished(self, save_dir):
        self.overlay.hide()
        self.param_label.setText(f"Training finished. Model saved in {save_dir}")
        QMessageBox.information(self, "Training Finished", f"AI Model training completed successfully.\nSaved to: {save_dir}")

    def import_ai_model(self):
        """AI 모델을 import하는 메서드"""
        model_path, _ = QFileDialog.getOpenFileName(
            self, "Select AI Model File", "", "Model Files (*.pkl *.keras)"
        )
        if model_path:
            try:
                if model_path.lower().endswith(".pkl"):
                    self.ai_model = joblib.load(model_path)  # Random Forest 모델 로드
                    self.ai_model_type = "rf"  # 모델 유형 저장
                    self.ai_model_label.setText("Pos/Neg AI Model: RF")  # 모델 유형 표시
                    self.ai_model_label.setStyleSheet(
                        "color: green; font-weight: bold;"
                    )
                elif model_path.lower().endswith(".keras"):
                    from tensorflow.keras.models import load_model

                    self.ai_model = load_model(model_path)  # CNN 모델 로드
                    self.ai_model_type = "cnn"  # 모델 유형 저장
                    self.ai_model_label.setText("Pos/Neg AI Model: CNN")  # 모델 유형 표시
                    self.ai_model_label.setStyleSheet(
                        "color: green; font-weight: bold;"
                    )
                else:
                    raise ValueError(
                        "Unsupported model file format. Use .pkl or .keras."
                    )
                print(f"Pos/Neg AI Model loaded from {model_path} (Type: {self.ai_model_type})")
            except Exception as e:
                print(f"Error loading Pos/Neg AI model: {e}")
                self.ai_model = None
                self.ai_model_type = None
                self.ai_model_label.setText("No Pos/Neg AI Model Loaded")
                self.ai_model_label.setStyleSheet("color: red; font-weight: bold;")

    def logistic_fitting(self):
        """현재 plot된 그래프들을 label별로 autolearner.sigmoid_fitting() 호출"""
        if not hasattr(self, "lines") or not self.lines:
            print("No lines to fit.")
            return
        for label, line in self.lines.items():
            # label에서 kkk 추출 (kkk는 label의 앞부분, 예: "3_20240529_120000" 또는 "3 (tr_P:pr_N)")
            kkk = None
            if "(" in label:
                kkk = label.split(" ")[0]
            else:
                kkk = label.split("_")[0]
            try:
                kkk = int(kkk)
            except Exception:
                print(f"Cannot extract kkk from label: {label}")
                continue
            y = self.autolearner.Cq_values.get(kkk)
            x = self.autolearner.cycles
            if y is not None:
                popt = self.autolearner.sigmoid_fitting(x, y)
                print(f"{label}\t{popt}")
            else:
                print(f"No Cq_values for label: {label}")

    def _check_current_mode_anomaly_detection(self):
        # baseline_substraction 또는 normalization 상태에서는 anomaly_detection을 실행하지 않음
        if self.normalized or self.baseline_substracted:
            print("Anomaly detection은 원본 데이터에서만 작동합니다.")
            if self.baseline_substracted:
                print("현재 Baseline substracted 상태입니다.")
            elif self.normalized:
                print("현재 Normalized 상태입니다.")
            return True
        if not hasattr(self, "lines") or not self.lines:
            print("No lines to detect anomalies")
            return True
        return False

    def _get_threshold_anomaly_detection(self):
        from PyQt5.QtWidgets import QInputDialog

        default_thresh = (
            getattr(self.anomaly_args, "anomaly_threshold", 0.01)
            if self.anomaly_args
            else 0.01
        )
        thresh, ok = QInputDialog.getDouble(
            self,
            "Anomaly threshold",
            "Enter threshold (absolute error):",
            value=float(default_thresh),
            min=0.0,
            max=1.0,
            decimals=2,
        )
        if not ok:
            print("Anomaly detection cancelled by user.")
            return False
        # 선택된 값을 현재 인스턴스/args에 저장
        if self.anomaly_args is not None:
            try:
                self.anomaly_args.anomaly_threshold = float(thresh)
            except Exception:
                pass
        threshold = float(thresh)
        return threshold

    def _pseudo_pcr_runner(self, pcr_analyzer, y, x):
        for cq_value in y:
            pcr_analyzer.add_one_cq_value(cq_value, set_winsize=3, start=x[0])
        pcr_analyzer.pcr_machine_cal_dCq()
        adj_y = pcr_analyzer.get_cq_values()
        baseline = adj_y["Cq_smoothing"][0]
        adj_y = [value - baseline for value in adj_y["Cq_smoothing"]]
        adj_y = [value / y[0] for value in adj_y] if y[0] > 0 else adj_y
        return adj_y

    def anomaly_detection(self):
        """현재 plot된 그래프들을 label별로 GRU 를 이용하여 anomal RFU들을 표시. 표시는 X표를 각 값들 위에 표시한다.
        현재 GRU 모델은 pcr_autolearner::AutoLearner.normalization(smoothing_size=3, mode='primitive') 인 것을 대상으로 만들어져있다.
        normalization이나 baseline substraction에서는 작동하지 않도록 한다.
        이 기능이 trigger되면 popup 창이 뜨고 거기에 threshold를 입력 할 수 있게 한다.
        """
        # baseline_substraction 또는 normalization 상태에서는 anomaly_detection을 실행하지 않음
        if self._check_current_mode_anomaly_detection():
            return

        # 사용자로부터 threshold 입력 받기 (popup)
        threshold = self._get_threshold_anomaly_detection()
        if threshold is False:
            return

        pcr_analyzer = PCRmachine()
        self.anomaly_detector.load_model(self.anomaly_args.anomaly_model)

        anomalies_dict = {}
        total = len(self.lines)
        for idx, (label, line) in enumerate(self.lines.items()):
            percent = int((idx + 1) / total * 100)
            print(
                f"Detecting anomalies for {label}...[{percent}%]", end="\r", flush=True
            )
            pcr_analyzer.reset()
            kkk = None
            kkk = self.autolearner.get_key_from_label(label)
            try:
                kkk = int(kkk)
            except Exception:
                print(f"Cannot extract kkk from label: {label}")
                continue
            y = self.autolearner.Cq_values.get(kkk)
            x = self.autolearner.cycles
            if y is not None:
                adj_y = self._pseudo_pcr_runner(pcr_analyzer, y, x)
                anomalies = self.anomaly_detector.detect_anomaly(
                    adj_y, threshold=threshold
                )
                anomalies_dict[label] = anomalies
            else:
                print(f"No Cq_values for label: {label}")
        self.update_graph(self.target_name, anomalies_dict=anomalies_dict)
        try:
            self.print_anomaly_csv_files(anomalies_dict)
        except Exception as e:
            print(f"Failed to show anomaly file diaglos:{e}")

    def print_anomaly_csv_files(self, anomalies_dict):
        """콘솔에 이상치가 3개 이상인 샘플의 파일명(디렉토리 제외)을 'label:파일이름' 형식으로 출력."""
        threshold = 3
        if not anomalies_dict:
            return
        print(f"====== files whose number of anomal cycles >= {threshold} =======")
        for label, anomalies in anomalies_dict.items():
            try:
                if not anomalies or len(anomalies) < threshold:
                    continue

                # kkk 추출
                kkk = self.autolearner.get_key_from_label(label)
                if kkk is None:
                    if "(" in label:
                        kkk = label.split(" ")[0]
                    else:
                        kkk = label.split("_")[0]
                try:
                    kkk_int = int(kkk)
                except Exception:
                    kkk_int = kkk

                files = self.autolearner.file_names.get(
                    kkk_int, self.autolearner.file_names.get(kkk, None)
                )

                if not files:
                    print(f"{label}:No file info available")
                    continue

                if isinstance(files, (list, tuple)):
                    for fn in files:
                        print(f"{label}:{os.path.basename(str(fn))}")
                else:
                    print(f"{label}:{os.path.basename(str(files))}")
            except Exception as e:
                print(f"Failed to print files for {label}: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PCR Autolearner GUI")
    parser.add_argument(
        "--anomaly_model",
        dest="anomaly_model",
        default="models/pcr_anomaly_detector_gru_ws3_float16.tflite",
        help="Path to anomaly detector model (.tflite or Keras). Default: models/pcr_anomaly_detector_gru_ws3_float16.tflite",
    )
    parser.add_argument(
        "--anomaly_threshold",
        type=float,
        default=0.01,
        help="Anomaly detection threshold. Default: 0.01",
    )
    args = parser.parse_args()

    app = QApplication(sys.argv)
    # UIWindow accepts optional model_path (will try to load; warning printed on failure)
    window = UIWindow(anomaly_args=args)
    sys.exit(app.exec_())
