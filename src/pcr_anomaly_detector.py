import os
import numpy as np
import src.pcr_autolearner_ai as auto_ai


class PCRAnomalyDetector:
    def __init__(self, seq_len=3, units=(64, 32), use_mixed_precision=False):
        """
        units: tuple(int, int) - 유닛 수 (첫 RNN 레이어, 두번째 RNN 레이어)
        use_mixed_precision: bool - True면 mixed_float16 policy 사용 (훈련 메모리 감소 가능)
        """
        self.seq_len = seq_len
        self.model = None
        self.units = tuple(units)
        self.use_mixed_precision = use_mixed_precision
        if self.use_mixed_precision:
            try:
                from tensorflow.keras import mixed_precision

                mixed_precision.set_global_policy("mixed_float16")
            except Exception:
                pass

    def build_model(self, model_type="lstm"):
        if model_type == "lstm":
            self.build_model_lstm()
        elif model_type == "gru":
            self.build_model_gru()
        elif model_type == "rnn":
            self.build_model_rnn()
        else:
            raise ValueError(f"Unsupported model type: {model_type}")

    def build_model_lstm(self):
        from tensorflow.keras.models import Sequential
        from tensorflow.keras.layers import LSTM, Dense, Dropout

        u1, u2 = self.units
        model = Sequential()
        model.add(LSTM(u1, return_sequences=True, input_shape=(self.seq_len, 1)))
        model.add(Dropout(0.2))
        model.add(LSTM(u2))
        model.add(Dense(1))
        model.compile(optimizer="adam", loss="mse")
        self.model = model

    def build_model_gru(self):
        from tensorflow.keras.models import Sequential
        from tensorflow.keras.layers import GRU, Dense, Dropout

        u1, u2 = self.units
        model = Sequential()
        model.add(GRU(u1, return_sequences=True, input_shape=(self.seq_len, 1)))
        model.add(Dropout(0.2))
        model.add(GRU(u2))
        model.add(Dense(1))
        model.compile(optimizer="adam", loss="mse")
        self.model = model

    def build_model_rnn(self):
        from tensorflow.keras.models import Sequential
        from tensorflow.keras.layers import SimpleRNN, Dense, Dropout

        u1, u2 = self.units
        model = Sequential()
        model.add(SimpleRNN(u1, return_sequences=True, input_shape=(self.seq_len, 1)))
        model.add(Dropout(0.2))
        model.add(SimpleRNN(u2))
        model.add(Dense(1))
        model.compile(optimizer="adam", loss="mse")
        self.model = model

    def fit(self, cq_sequences, y, epochs=100, model_path=None, model_type="lstm", initial_model_path=None):
        if initial_model_path and os.path.exists(initial_model_path):
            print(f"Loading initial model for continued training: {initial_model_path}")
            self.load_model(initial_model_path)
        else:
            self.build_model(model_type)
        callbacks = []
        if model_path:
            from tensorflow.keras.callbacks import ModelCheckpoint

            checkpoint = ModelCheckpoint(
                model_path, monitor="val_loss", save_best_only=True, verbose=0
            )
            callbacks.append(checkpoint)
        self.model.fit(
            cq_sequences,
            y,
            epochs=epochs,
            verbose=1,
            validation_split=0.2,
            callbacks=callbacks,
        )

    def predict_next(self, cq_sequence):
        # cq_sequence: (seq_len, 1)
        inp = cq_sequence.reshape(1, self.seq_len, 1)
        # TFLite interpreter가 있으면 그것을 사용
        if hasattr(self, "tflite_interpreter") and self.tflite_interpreter is not None:
            # 입력 dtype/shape에 맞게 변환
            input_details = self._tflite_input_details
            idx = input_details[0]["index"]
            expected_dtype = input_details[0]["dtype"]
            data = inp.astype(expected_dtype)
            self.tflite_interpreter.set_tensor(idx, data)
            self.tflite_interpreter.invoke()
            output_details = self._tflite_output_details
            out = self.tflite_interpreter.get_tensor(output_details[0]["index"])
            return float(out.reshape(-1)[0])
        # fallback: Keras 모델 사용
        return self.model.predict(inp, verbose=0)[0, 0]

    def detect_anomaly(
        self, cq_values, threshold=0.01, save_plot_path=None, label=None
    ):
        """_summary_

        Args:
            cq_values (_type_): _description_
            threshold (float, optional): _description_. Defaults to 0.01.
            save_plot_path (_type_, optional): _description_. Defaults to None.
            label (_type_, optional): _description_. Defaults to None.

        Returns:
            list: indices of detected anomalies. 0-based index.
        """
        anomalies = []
        preds = []
        indices = []
        for i in range(self.seq_len, len(cq_values)):
            seq = np.array(cq_values[i - self.seq_len : i]).reshape(self.seq_len, 1)
            pred = self.predict_next(seq)
            preds.append(pred)
            indices.append(i)
            error = abs(pred - cq_values[i])
            # 이전 값보다 감소하는 경우에만 이상치로 간주
            if error > threshold and cq_values[i - 1] > cq_values[i]:
                anomalies.append(i)

        # 시각화 및 저장
        if save_plot_path is not None:
            import matplotlib
            import matplotlib.pyplot as plt

            matplotlib.use("Agg")  # GUI 없이 파일로만 저장하는 백엔드

            plt.figure(figsize=(10, 5))
            plt.plot(
                range(len(cq_values)), cq_values, label="True Cq values", marker="o"
            )
            plt.plot(indices, preds, label="Predicted", marker="x")
            if anomalies:
                plt.plot(
                    anomalies,
                    np.array(cq_values)[anomalies],
                    color="red",
                    marker="o",
                    linestyle="None",
                    label="Anomaly",
                )
            if label:
                plt.title(f"Anomaly Detection - {label}")
            else:
                plt.title("Anomaly Detection")
            plt.xlabel("Cycle")
            plt.ylabel("Cq Value")
            plt.ylim(min(-0.1, min(cq_values)), max(1.0, max(cq_values)))
            plt.legend()
            plt.tight_layout()
            plt.savefig(save_plot_path)
            plt.close()

        return anomalies

    def get_error(self, cq_values):
        errors = []
        for i in range(self.seq_len, len(cq_values)):
            seq = np.array(cq_values[i - self.seq_len : i]).reshape(self.seq_len, 1)
            pred = self.predict_next(seq)
            error = abs(pred - cq_values[i])
            errors.append(error)
        return errors

    def save_model(self, filepath="models/pcr_anomaly_detector_lstm.keras"):
        self.model.save(filepath)

    def export_tflite_float16(
        self, out_path="models/pcr_anomaly_detector_float16.tflite"
    ):
        """
        Keras 모델을 임시 SavedModel로 저장한 뒤 float16 양자화된 TFLite로 변환하여 저장합니다.
        변환 중 TensorList/TensorArray 관련 오류를 우회하기 위해 Select TF ops 허용,
        _experimental_lower_tensor_list_ops 비활성화 및 resource variables 사용 플래그를 설정합니다.
        """
        if self.model is None:
            raise RuntimeError("No model to export. Build or load a model first.")
        try:
            import tensorflow as tf
            import tempfile
            import os
            from pathlib import Path
        except Exception as e:
            raise RuntimeError(f"TensorFlow import failed: {e}")

        try:
            with tempfile.TemporaryDirectory() as td:
                # SavedModel 형식으로 저장: Keras 3에서는 model.export()가 있으면 사용
                if hasattr(self.model, "export"):
                    self.model.export(td)
                else:
                    # fallback: model.save는 일부 환경에서 Keras3와 충돌할 수 있음.
                    try:
                        self.model.save(td)
                    except Exception as e:
                        raise RuntimeError(
                            "Failed to save model as SavedModel. "
                            "Use a TF/Keras version that supports model.export() or model.save() for SavedModel."
                            f" Underlying error: {e}"
                        )

                # SavedModel로부터 TFLite 변환기 생성
                converter = tf.lite.TFLiteConverter.from_saved_model(td)
                converter.optimizations = [tf.lite.Optimize.DEFAULT]
                # float16 타겟 타입 설정
                converter.target_spec.supported_types = [tf.float16]

                # TensorList/TensorArray 관련 오류가 발생하는 모델(예: GRU)이 있을 수 있으므로
                # Select TF ops 허용 및 tensor list lowering 비활성화 권장
                try:
                    converter.target_spec.supported_ops = [
                        tf.lite.OpsSet.TFLITE_BUILTINS,
                        tf.lite.OpsSet.SELECT_TF_OPS,
                    ]
                except Exception:
                    # 일부 TF 버전에서는 target_spec.supported_ops 가 없을 수 있음 -> 무시
                    pass

                # 권장 실험적 플래그 설정 (버전에 따라 속성명이 다를 수 있음)
                # resource 변수 상수 폴딩 문제 회피
                if hasattr(converter, "experimental_enable_resource_variables"):
                    converter.experimental_enable_resource_variables = True
                # tensor list lowering 비활성화
                try:
                    # 일부 TF 버전에서 private 속성명 사용
                    converter._experimental_lower_tensor_list_ops = False
                except Exception:
                    # 다른 이름이 있으면 무시
                    pass

                tflite_model = converter.convert()

                # 출력 디렉터리 생성
                out_dir = os.path.dirname(out_path)
                if out_dir:
                    Path(out_dir).mkdir(parents=True, exist_ok=True)
                with open(out_path, "wb") as f:
                    f.write(tflite_model)
        except Exception as e:
            raise RuntimeError(f"TFLite float16 export failed: {e}")
        return out_path

    def load_model(self, filepath):
        # allow loading either Keras (.keras/.h5) or TFLite (.tflite)
        # reset current models/interpreters
        self.model = None
        if hasattr(self, "tflite_interpreter"):
            try:
                del self.tflite_interpreter
            except Exception:
                pass
        if filepath.lower().endswith(".tflite"):
            try:
                import tensorflow as tf

                self.tflite_interpreter = tf.lite.Interpreter(model_path=filepath)
                self.tflite_interpreter.allocate_tensors()
                self._tflite_input_details = self.tflite_interpreter.get_input_details()
                self._tflite_output_details = (
                    self.tflite_interpreter.get_output_details()
                )
                return
            except Exception as e:
                raise RuntimeError(f"Failed to load TFLite model: {e}")
        else:
            from tensorflow.keras.models import load_model

            self.model = load_model(filepath)
            # ensure no leftover tflite interpreter
            if hasattr(self, "tflite_interpreter"):
                try:
                    del self.tflite_interpreter
                except Exception:
                    pass
            return


def train_mode(args, runner):
    input_rfus, input_labels, targets = auto_ai.load_data(
        runner, args.data_dir, normalization_mode="primitive", ddcq=False
    )

    # LSTM 입력 시퀀스 생성 (예: sliding window)
    seq_len = args.window_size
    X = []
    y = []
    for cq_values, label in zip(input_rfus, input_labels):
        if label == "Negative":
            continue  # 음성 샘플은 제외
        for i in range(len(cq_values) - seq_len):
            X.append(cq_values[i : i + seq_len])
            y.append(cq_values[i + seq_len])

    X = np.array(X).reshape(-1, seq_len, 1)  # (num_samples, seq_len, 1)
    y = np.array(y).reshape(-1, 1)  # (num_samples,)

    # LSTM 학습
    detector = PCRAnomalyDetector(seq_len=seq_len)
    detector.fit(
        X, y, epochs=100, model_path=args.model_path, model_type=args.model_type
    )

    # ModelCheckpoint로 best 모델을 저장했다면 해당 파일을 불러와 float16 TFLite로 변환하여 저장
    if args.model_path:
        try:
            # 강제로 파일에 저장된 best 모델을 로드하도록 현재 메모리 모델을 제거
            detector.model = None
            detector.load_model(args.model_path)
        except Exception as e:
            print(f"Warning: failed to load model from {args.model_path}: {e}")

        try:
            base = os.path.splitext(args.model_path)[0]
            tflite_path = base + "_float16.tflite"
            detector.export_tflite_float16(tflite_path)
            print(f"Exported float16 TFLite: {tflite_path}")
        except Exception as e:
            print(f"Failed to export float16 TFLite: {e}")


def train_mode_saturation(args, runner):
    """
    다양한 수의 정상 샘플로 학습시켜 validation loss 기준으로 성능 포화(saturation) 지점을 탐색.
    - args: 기존 parse_arguments()에서 가져오는 args (window_size 등 사용)
    - runner: auto_ai.AutoLearner 인스턴스

    Returns: dict {sample_size: (mean_val_loss, std_val_loss)}
    """
    import csv
    import random
    import numpy as np

    """ 아래 변수 설명:
    - sample_sizes: 리스트로 테스트할 샘플 수들. None이면 자동 생성.
    - epochs: 각 실험의 학습 epoch 수
    - model_type: "lstm" | "gru" | "rnn"
    - repeats: 각 샘플 수에 대해 반복 실험 횟수 (무작위 샘플링)
    - out_csv: 결과를 저장할 CSV 경로"""
    sample_sizes = args.sample_sizes
    epochs = args.epochs
    model_type = args.model_type
    repeats = args.repeats
    out_csv = args.out_csv

    # 데이터 로드 (정규화 모드 등 train_mode와 동일)
    input_rfus, input_labels, targets = auto_ai.load_data(
        runner, args.data_dir, normalization_mode="primitive", ddcq=False
    )

    seq_len = args.window_size
    # 정상(Positive 또는 비빈 레이블을 정상으로 간주)의 시퀀스 리스트 수집
    positive_sequences = []
    for cq_values, label in zip(input_rfus, input_labels):
        if label == "Negative":
            continue
        # sliding-window로부터 한 시퀀스 전체(샘플)로 취급: 각 sample은 cq_values 전체 시퀀스
        positive_sequences.append(np.array(cq_values))

    n_pos = len(positive_sequences)
    if n_pos == 0:
        raise RuntimeError("사용 가능한 정상 샘플이 없습니다.")

    # sample_sizes 기본 생성: 5,10,20,... up to n_pos
    if sample_sizes is None:
        sizes = []
        s = 5
        while s < n_pos:
            sizes.append(s)
            s = int(s * 2) if s * 2 < n_pos else s + max(1, int((n_pos - s) / 5))
            if s > n_pos - 1:
                break
        sizes.append(n_pos)
        sample_sizes = sorted(list(set(sizes)))
    results = {}

    print(
        f"Start saturation test: total positive samples={n_pos}, seq_len={seq_len}, model_type={model_type}"
    )

    for sample_count in sample_sizes:
        if sample_count > n_pos:
            print(f"skip {sample_count} (only {n_pos} available)")
            continue
        val_losses = []
        for r in range(repeats):
            # 무작위로 sample_count개의 샘플 선택
            chosen = random.sample(positive_sequences, sample_count)
            # sliding window으로 X,y 생성 (모든 선택된 샘플을 합침)
            X = []
            y = []
            for cq_values in chosen:
                for i in range(len(cq_values) - seq_len):
                    X.append(cq_values[i : i + seq_len])
                    y.append(cq_values[i + seq_len])
            X = np.array(X).reshape(-1, seq_len, 1)
            y = np.array(y).reshape(-1, 1)
            if len(X) == 0:
                # 윈도우가 너무 커서 샘플 생성 불가
                print(
                    f" sample_count={sample_count}, repeat={r}: no windows (seq_len too large). skip trial."
                )
                continue

            # 새 모델 생성 및 학습 (validation_split 사용)
            detector = PCRAnomalyDetector(seq_len=seq_len)
            detector.build_model(model_type)
            history = detector.model.fit(
                X, y, epochs=epochs, verbose=0, validation_split=0.2
            )
            val_loss_min = float(np.min(history.history.get("val_loss", [np.nan])))
            val_losses.append(val_loss_min)
            print(
                f"\rsample_count={sample_count} repeat={r+1}/{repeats} val_loss_min={val_loss_min:.6f}",
                end="",
                flush=True,
            )
        print("")  # 한 줄 마무리

        if len(val_losses) == 0:
            mean_loss = float("nan")
            std_loss = float("nan")
        else:
            mean_loss = float(np.mean(val_losses))
            std_loss = float(np.std(val_losses))
        results[sample_count] = (mean_loss, std_loss)

        # 중간 저장(append) — CSV 헤더 처리
        header = ["sample_count", "mean_val_loss", "std_val_loss"]
        write_header = False
        try:
            if not os.path.exists(out_csv):
                write_header = True
        except Exception:
            write_header = True
        with open(out_csv, "a", newline="") as f:
            writer = csv.writer(f)
            if write_header:
                writer.writerow(header)
            writer.writerow([sample_count, mean_loss, std_loss])

    print("Saturation test finished. Results saved to", out_csv)
    return results


def simulation_mode(args, runner):
    # 모델 로드
    detector = PCRAnomalyDetector(seq_len=args.window_size)
    detector.load_model(args.model_path)
    print(f"Model loaded from {args.model_path}")

    # 데이터 로드
    runner.read_dir(args.data_dir)
    targets = runner.get_target_names()
    for target in targets:
        input_rfus, input_labels, target, keys = auto_ai.load_data_target(
            runner, target, normalization_mode="primitive", ddcq=False
        )
        if args.debug:
            input_labels = list(input_labels)
            _n_pos = input_labels.count("Positive")
            _n_neg = input_labels.count("Negative")
            _n_void = input_labels.count("")
            print(
                f"Total samples {target}: {len(input_rfus)}={len(input_labels)}={len(keys)}"
            )
            print(f"Positive: {_n_pos}, Negative: {_n_neg}, Void: {_n_void}")
        for cq_values, label, key in zip(input_rfus, input_labels, keys):
            # if label != "Positive":  # Positive 또는 Negative 또는 빈칸
            #     continue
            file_name = runner.file_names[key]
            detector.detect_anomaly(
                cq_values,
                threshold=args.threshold,
                save_plot_path=f"{args.imgpath}/{target}_{key}_{label}.png",
                label=f"{target}_{key}:{label}",
            )


def export_tflite_mode(args):
    """Keras(.keras/.h5) 모델을 float16 TFLite로 변환해서 저장."""
    detector = PCRAnomalyDetector(seq_len=getattr(args, "window_size", 3))
    # 로드(keras 파일이어야 함)
    try:
        detector.load_model(args.model_path)
    except Exception as e:
        print(f"Failed to load model '{args.model_path}': {e}")
        return

    # 입력이 이미 tflite면 작업 없음
    if args.model_path.lower().endswith(".tflite"):
        print("Input model is already a TFLite file. Nothing to export.")
        return

    # 출력 경로 결정
    out_path = args.out_path
    if not out_path:
        import os

        base = os.path.splitext(args.model_path)[0]
        out_path = base + "_float16.tflite"

    try:
        tpath = detector.export_tflite_float16(out_path)
        print(f"Exported float16 TFLite: {tpath}")
    except Exception as e:
        print(f"Export failed: {e}")


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description="PCR Anomaly Detector Training")
    subparsers = parser.add_subparsers(
        dest="mode", required=True, help="Subcommands: train"
    )
    ## Subcommand: train
    train_parser = subparsers.add_parser("train", help="Train the anomaly detector")
    train_parser.add_argument(
        "data_dir",
        type=str,
        help="Directory containing PCR csv files",
    )
    train_parser.add_argument(
        "model_path",
        type=str,
        help="Path to save the trained model (e.g. models/pcr_anomaly_detector_lstm.keras)",
    )
    train_parser.add_argument(
        "--window_size",
        type=int,
        default=3,
        help="Sliding window size for RFU cycles (default: 3)",
    )
    train_parser.add_argument(
        "--model_type",
        type=str,
        default="lstm",
        choices=["lstm", "gru", "rnn"],
        help="Type of RNN model to use (default: lstm)",
    )

    ## Subcommand: train_saturation
    sat_parser = subparsers.add_parser(
        "train_saturation",
        help="Train the anomaly detector with saturation test",
    )
    sat_parser.add_argument(
        "data_dir",
        type=str,
        help="Directory containing PCR csv files",
    )
    sat_parser.add_argument(
        "--window_size",
        type=int,
        default=3,
        help="Sliding window size for RFU cycles (default: 3)",
    )
    sat_parser.add_argument(
        "--sample_sizes",
        type=int,
        nargs="+",
        default=None,
        help="List of sample sizes to test (default: auto-generated)",
    )
    sat_parser.add_argument(
        "--epochs",
        type=int,
        default=50,
        help="Number of epochs for each training (default: 50)",
    )
    sat_parser.add_argument(
        "--model_type",
        type=str,
        default="gru",
        choices=["lstm", "gru", "rnn"],
        help="Type of RNN model to use (default: gru)",
    )
    sat_parser.add_argument(
        "--repeats",
        type=int,
        default=3,
        help="Number of repeats for each sample size (default: 3)",
    )
    sat_parser.add_argument(
        "--out_csv",
        type=str,
        default="saturation_results.csv",
        help="Output CSV file to save saturation test results (default: saturation_results.csv)",
    )

    ## Subcommand: simulation
    sim_parser = subparsers.add_parser(
        "simulation", help="Run simulation with the trained model"
    )
    sim_parser.add_argument(
        "data_dir",
        type=str,
        help="Directory containing PCR csv files",
    )
    sim_parser.add_argument(
        "model_path",
        type=str,
        help="Path to the trained model (e.g. models/pcr_anomaly_detector_lstm_wsN.keras)",
    )
    sim_parser.add_argument(
        "--window_size",
        type=int,
        default=3,
        help="Sliding window size for RFU cycles (default: 3)",
    )
    sim_parser.add_argument(
        "--threshold",
        type=float,
        default=0.01,
        help="Anomaly detection threshold (default: 0.01)",
    )
    sim_parser.add_argument(
        "--imgpath",
        type=str,
        default="results",
        help="Directory to save anomaly detection plots (default: results)",
    )
    sim_parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug mode to print additional information",
    )

    ## Subcommand: export_tflite
    exp_parser = subparsers.add_parser(
        "export_tflite", help="Convert a Keras model (.keras/.h5) to float16 TFLite"
    )
    exp_parser.add_argument(
        "model_path",
        type=str,
        help="Path to the Keras model to convert (e.g. models/model.keras)",
    )
    exp_parser.add_argument(
        "out_path",
        nargs="?",
        default=None,
        help="Output tflite path (optional). If omitted, '<model>_float16.tflite' is used.",
    )
    exp_parser.add_argument(
        "--window_size",
        type=int,
        default=3,
        help="(optional) seq_len for detector instance (default: 3)",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    runner = auto_ai.AutoLearner()
    if args.mode == "train":
        train_mode(args, runner)
    elif args.mode == "simulation":
        simulation_mode(args, runner)
    elif args.mode == "train_saturation":
        train_mode_saturation(args, runner)
    elif args.mode == "export_tflite":
        export_tflite_mode(args)
