# import sys
import argparse
import os
import numpy as np

# import joblib  # Ensure joblib is imported at the top of the file
from src.pcr_autolearner import AutoLearner  # Import the Autolearner class

# from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split


class PCRPredicter:
    def __init__(self):
        self.autolearner = AutoLearner()  # Use the Autolearner class
        self.target_id = None

    def preprocess_data(self, data, labels):
        """Preprocess the data to normalize and set Positive/Negative labels."""

        # Set Positive (1) and Negative (0) labels
        labels = np.array([1 if label == "Positive" else 0 for label in labels])

        # Ensure all samples have 45 cycles
        processed_data = []
        for sample in data:
            if len(sample) < 45:
                # Span the last value to make it 45 cycles
                sample = np.pad(sample, (0, 45 - len(sample)), mode="edge")
            processed_data.append(sample)
        processed_data = np.array(processed_data)

        return processed_data, labels

    def build_cnn_model(self, input_shape, learning_rate=0.001):
        """Build and compile a CNN model for classification."""
        from tensorflow.keras.models import Sequential
        from tensorflow.keras.layers import (
            Conv1D,
            Dense,
            Flatten,
            MaxPooling1D,
            Dropout,
            BatchNormalization,
        )
        from tensorflow.keras.optimizers import Adam

        # Define the CNN architecture
        model = Sequential(
            [
                Conv1D(32, kernel_size=3, activation="relu", input_shape=input_shape),
                MaxPooling1D(pool_size=2),
                Conv1D(64, kernel_size=3, activation="relu"),
                MaxPooling1D(pool_size=2),
                Flatten(),
                Dense(128, activation="relu"),
                Dropout(0.5),  # Dropout with 50% rate
                BatchNormalization(),
                Dense(
                    2, activation="softmax"
                ),  # Output layer for binary classification
            ]
        )

        # Compile the model
        model.compile(
            optimizer=Adam(learning_rate=learning_rate),
            loss="categorical_crossentropy",
            metrics=["accuracy"],
        )

        return model

    def train_model(
        self,
        data,
        labels,
        learning_rate=0.001,
        epochs=50,
        save_dir="models/",
        ai_model="cnn",
        initial_model_path=None,
    ):
        """Train a model (CNN or Random Forest) using the provided data and labels."""
        # Preprocess the data
        data, labels = self.preprocess_data(data, labels)

        if ai_model == "cnn":
            model = self._train_cnn(data, labels, learning_rate, epochs, save_dir, initial_model_path)
        elif ai_model == "rf":
            if initial_model_path:
                print("Warning: Incremental training for Random Forest is not fully supported via weight continuation. Re-training from scratch if needed.")
            model = self._train_random_forest(data, labels, save_dir)
        else:
            raise ValueError("Invalid ai_model. Choose 'cnn' or 'rf'.")

        return model

    def _train_cnn(self, data, labels, learning_rate, epochs, save_dir, initial_model_path=None):
        """Train a CNN model."""
        from tensorflow.keras.utils import to_categorical
        from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping
        from tensorflow.keras.models import load_model

        # Convert labels to categorical format
        labels = to_categorical(labels, num_classes=2)

        # Build or Load the CNN model
        input_shape = (data.shape[1], 1)  # (45 cycles, 1 feature per cycle)
        if initial_model_path and os.path.exists(initial_model_path):
            print(f"Loading initial model weights from {initial_model_path}")
            model = load_model(initial_model_path)
            from tensorflow.keras.optimizers import Adam
            model.compile(
                optimizer=Adam(learning_rate=learning_rate),
                loss="categorical_crossentropy",
                metrics=["accuracy"],
            )
        else:
            model = self.build_cnn_model(input_shape, learning_rate=learning_rate)

        # Reshape data for CNN input
        data = data.reshape((data.shape[0], data.shape[1], 1))

        # Define callbacks
        checkpoint_callback = ModelCheckpoint(
            filepath=os.path.join(save_dir, "best_model_cnn.keras"),
            monitor="val_loss",
            save_best_only=True,
            save_weights_only=False,
            verbose=1,
        )
        early_stopping = EarlyStopping(
            monitor="val_loss",
            patience=10,
            verbose=1,
            restore_best_weights=True,
        )

        # Train the CNN model
        model.fit(
            data,
            labels,
            epochs=epochs,
            batch_size=32,
            validation_split=0.2,
            callbacks=[checkpoint_callback, early_stopping],
        )

        return model

    def _train_random_forest(self, data, labels, save_dir):
        """Train a Random Forest model."""
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.metrics import accuracy_score
        from sklearn.model_selection import train_test_split
        import joblib

        # Split data into training and validation sets
        X_train, X_val, y_train, y_val = train_test_split(
            data, labels, test_size=0.2, random_state=42
        )

        # Train a Random Forest model
        model = RandomForestClassifier(n_estimators=100, random_state=42)
        model.fit(X_train, y_train)

        # Evaluate the model on the validation set
        val_predictions = model.predict(X_val)
        accuracy = accuracy_score(y_val, val_predictions)
        print(f"Validation Accuracy (Random Forest): {accuracy:.4f}")

        # Save the Random Forest model
        model_path = os.path.join(save_dir, "best_model_rf.pkl")
        joblib.dump(model, model_path)
        print(f"Random Forest model saved to {model_path}")

        return model

    def predict(self, data, model_path):
        """Make predictions using the trained model."""
        import numpy as np

        import tensorflow as tf
        from tensorflow.keras.models import load_model

        # Set random seeds for reproducibility
        np.random.seed(42)
        tf.random.set_seed(42)

        # Preprocess the data
        data, _ = self.preprocess_data(data, labels=[])

        # Load the model
        model = load_model(model_path)
        model.trainable = False  # Ensure the model is not trainable

        # Reshape data for CNN input
        data = data.reshape((data.shape[0], data.shape[1], 1))

        # Make predictions
        predictions = model.predict(data, batch_size=32, verbose=0)
        predicted_labels = np.argmax(predictions, axis=1)

        return predicted_labels

    def save_rfu_plot_one(self, data, label, target, target_num, save_dir="tempdir/"):
        """Save RFU plot for the given data and labels without axes."""
        import matplotlib.pyplot as plt
        import os

        # Set the QT_QPA_PLATFORM environment variable. 이게 없으면 X 연결해야 함
        os.environ["QT_QPA_PLATFORM"] = "offscreen"

        # Ensure the save directory exists
        os.makedirs(save_dir, exist_ok=True)

        # Calculate the y-axis scale
        y_min = data[0]
        y_max = max(max(data), 2 * data[0])

        # Create a plot
        plt.figure(figsize=(10, 6))
        plt.plot(data.T, alpha=0.5)
        plt.ylim(y_min - 0.1 * y_min, y_max + 0.1 * y_min)  # Set y-axis limits
        plt.axis("off")  # Turn off both x-axis and y-axis
        plt.savefig(
            f"{save_dir}/{target}_{target_num}_{label}.png",
            bbox_inches="tight",
            pad_inches=0,
        )
        plt.close()

    def save_rfu_plot_all(self, data, labels, target, keys, save_dir="tempdir/"):
        """Save RFU plots for all data and labels."""
        import os

        # Ensure the save directory exists
        os.makedirs(save_dir, exist_ok=True)

        # Create a plot for each sample
        for i in range(len(data)):
            self.save_rfu_plot_one(
                data[i], labels[i], target, keys[i], save_dir=save_dir
            )

    def find_horizontal_lines_hough(
        self,
        image_path,
        rho=1,
        theta=np.pi / 180,
        threshold=100,
        min_line_length=50,
        max_line_gap=10,
    ):
        """
        허프 변환을 사용하여 이미지에서 수평선을 검출합니다.

        Args:
            image_path (str): 이미지 파일 경로.
            rho (int): 허프 공간에서 직선의 거리에 대한 해상도 (픽셀 단위).
            theta (float): 허프 공간에서 직선의 각도에 대한 해상도 (라디안 단위).
            threshold (int): 허프 공간에서 직선으로 간주될 최소 교차점 수.
            min_line_length (int): 검출될 선분의 최소 길이 (픽셀 단위).
            max_line_gap (int): 동일한 선분으로 간주될 수 있는 최대 간격 (픽셀 단위).

        Returns:
            list: 검출된 수평선의 좌표 (x1, y1, x2, y2) 리스트.
                이미지 로드 실패 시 None 반환.
        """
        import cv2

        img = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
        if img is None:
            return None

        # 에지 검출
        edges = cv2.Canny(img, 50, 150, apertureSize=3)

        # 허프 선 변환 적용
        lines = cv2.HoughLinesP(
            edges,
            rho,
            theta,
            threshold,
            minLineLength=min_line_length,
            maxLineGap=max_line_gap,
        )

        horizontal_lines = []
        if lines is not None:
            for line in lines:
                x1, y1, x2, y2 = line[0]
                # 수평에 가까운 선분 필터링 (각도 기반)
                angle = np.arctan2(abs(y2 - y1), abs(x2 - x1))
                if angle < np.pi / 18:  # 약 10도 이내를 수평으로 간주 (조정 가능)
                    horizontal_lines.append((x1, y1, x2, y2))

        if horizontal_lines:
            print("검출된 수평선:")
            for x1, y1, x2, y2 in horizontal_lines:
                print(f"({x1}, {y1}) - ({x2}, {y2})")
            # img = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)  # Load as grayscale
            if img is None:
                raise FileNotFoundError(f"Image not found at {image_path}")
            img_color = cv2.cvtColor(
                img, cv2.COLOR_GRAY2BGR
            )  # Convert grayscale to BGR
            for x1, y1, x2, y2 in horizontal_lines:
                cv2.line(img_color, (x1, y1), (x2, y2), (0, 0, 255), 2)
            cv2.imshow("Horizontal Lines Detected", img_color)
            cv2.waitKey(0)
            cv2.destroyAllWindows()

        else:
            print("수평선을 찾지 못했습니다.")
        return horizontal_lines



class PCRCluster:
    def __init__(self, n_clusters=3, method="kmeans", early_weight=1.0):
        self.n_clusters = n_clusters
        self.method = method
        self.early_weight = early_weight
        self.model = None

    def preprocess(self, data, skip_normalization=False):
        """Normalize data by dividing by the first cycle value.
        Rows with 0 in the first cycle are removed.
        
        Args:
            data (np.ndarray): Shape (n_samples, n_cycles)
            skip_normalization (bool): If True, skip normalization and filtering.
            
        Returns:
            tuple: (normalized_data, valid_indices)
        """
        if skip_normalization:
            # If skipping normalization, we assume data is already preprocessed (e.g. primitive normalization)
            # We return all indices as valid.
            return data, np.arange(len(data))

        # Filter out samples with 0 in the first cycle
        first_cycle_vals = data[:, 0]
        valid_indices = np.where(first_cycle_vals != 0)[0]
        
        if len(valid_indices) < len(data):
            print(f"Removed {len(data) - len(valid_indices)} samples with 0 value at the first cycle.")
            
        filtered_data = data[valid_indices]
        
        if len(filtered_data) == 0:
            return np.array([]), valid_indices

        # Division by first cycle value
        # reshaping for broadcasting: (N, 1)
        divisors = filtered_data[:, 0:1] 
        normalized = filtered_data / divisors
        
        # Apply weighting to first 10 cycles if early_weight != 1.0
        if self.early_weight != 1.0 and normalized.shape[1] >= 10:
            normalized[:, :10] *= self.early_weight
            
        return normalized, valid_indices

    def fit_predict(self, data):
        """Fit clustering model and predict cluster labels."""
        if self.method == "kmeans":
            from sklearn.cluster import KMeans
            self.model = KMeans(n_clusters=self.n_clusters, random_state=42, n_init=10)
            labels = self.model.fit_predict(data)
        elif self.method == "gmm":
            from sklearn.mixture import GaussianMixture
            self.model = GaussianMixture(n_components=self.n_clusters, random_state=42, n_init=10)
            labels = self.model.fit_predict(data)
        elif self.method == "bgmm":
            from sklearn.mixture import BayesianGaussianMixture
            # weight_concentration_prior_type='dirichlet_process' allows inferring effective number of components
            # from a larger n_components upper bound.
            self.model = BayesianGaussianMixture(
                n_components=self.n_clusters,
                random_state=42,
                n_init=10,
                weight_concentration_prior_type='dirichlet_process'
            )
            labels = self.model.fit_predict(data)
        else:
            raise ValueError(f"Unknown method: {self.method}")
            
        return labels

    def visualize_clusters(self, data, labels, sample_labels=None, save_dir=None):
        """Visualize clustering results.
        
        Args:
            data (np.ndarray): The normalized data.
            labels (np.ndarray): Cluster labels.
            sample_labels (list, optional): Labels for the PCA scatter plot.
            save_dir (str, optional): If provided, saves plots to this directory.
            
        Returns:
            dict: A dictionary containing the PCA figure and a list of cluster profile figures.
        """
        import matplotlib.pyplot as plt
        import os
        from sklearn.decomposition import PCA

        results = {}

        # 1. PCA Scatter Plot
        pca = PCA(n_components=2)
        reduced_data = pca.fit_transform(data)

        fig_pca = plt.figure(figsize=(12, 8))
        ax_pca = fig_pca.add_subplot(111)
        scatter = ax_pca.scatter(
            reduced_data[:, 0], reduced_data[:, 1], c=labels, cmap="viridis", alpha=0.6
        )
        
        if sample_labels and len(sample_labels) == len(data):
            for i, label in enumerate(sample_labels):
                ax_pca.text(
                    reduced_data[i, 0], 
                    reduced_data[i, 1], 
                    label, 
                    fontsize=8, 
                    alpha=0.7
                )

        ax_pca.set_title(f"PCR Clustering Results (PCA) - {self.method}")
        ax_pca.set_xlabel("PC1")
        ax_pca.set_ylabel("PC2")
        fig_pca.colorbar(scatter, label="Cluster Label")
        
        results['pca'] = fig_pca

        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            fig_pca.savefig(os.path.join(save_dir, "pca_scatter.png"))

        # 2. Cluster Profiles
        unique_labels = np.unique(labels)
        profile_figs = []
        for label in unique_labels:
            cluster_data = data[labels == label]
            fig_profile = plt.figure(figsize=(10, 6))
            ax_profile = fig_profile.add_subplot(111)
            
            # Plot individual curves faintly
            ax_profile.plot(cluster_data.T, color="gray", alpha=0.1)
            # Plot mean curve
            mean_curve = np.mean(cluster_data, axis=0)
            ax_profile.plot(mean_curve, color="red", linewidth=2, label="Mean Pattern")
            
            ax_profile.set_title(f"Cluster {label} Profile (n={len(cluster_data)})")
            ax_profile.set_xlabel("Cycle")
            ax_profile.set_ylabel("Normalized RFU")
            ax_profile.legend()
            
            profile_figs.append((label, fig_profile))
            
            if save_dir:
                fig_profile.savefig(os.path.join(save_dir, f"cluster_{label}_profile.png"))

        results['profiles'] = profile_figs

        if save_dir:
            print(f"Clustering visualizations saved to {save_dir}")
        
        return results


def parse_arguments():
    """Parse and validate command-line arguments."""
    parser = argparse.ArgumentParser(
        description="PCR Autolearner AI for PCR prediction"
    )
    subparsers = parser.add_subparsers(
        dest="mode",
        required=True,
        help="Subcommands: train, predict, save_rfu_plot, or cluster",
    )

    # Subcommand: cluster
    cluster_parser = subparsers.add_parser("cluster", help="Cluster RFU patterns")
    cluster_parser.add_argument(
        "data_dir", nargs="+", help="Directory(s) containing the PCR data files"
    )
    cluster_parser.add_argument(
        "--n_clusters",
        type=int,
        default=3,
        help="Number of clusters (or components) for K-Means, GMM, BGMM (default: 3)",
    )
    cluster_parser.add_argument(
        "--output_dir",
        default="cluster_results",
        help="Directory to save clustering results (default: cluster_results)",
    )
    cluster_parser.add_argument(
        "--method",
        choices=["kmeans", "gmm", "bgmm"],
        default="kmeans",
        help="Clustering method: kmeans (default), gmm, or bgmm",
    )
    cluster_parser.add_argument(
        "--early_weight",
        type=float,
        default=1.0,
        help="Weight factor for the first 10 cycles (default: 1.0)",
    )
    cluster_parser.add_argument(
        "--no_labels",
        action="store_true",
        help="Do not show labels on the plot",
    )


    # Subcommand: train
    train_parser = subparsers.add_parser("train", help="Train a model")
    train_parser.add_argument(
        "data_dir", help="Directory containing the PCR data files"
    )
    # train_parser.add_argument("model_path", help="Directory to save the trained model")
    train_parser.add_argument(
        "--learning_rate",
        type=float,
        default=0.001,
        help="Learning rate for training (default: 0.001)",
    )
    train_parser.add_argument(
        "--epochs",
        type=int,
        default=50,
        help="Number of epochs for training (default: 50)",
    )
    train_parser.add_argument(
        "--ai_model",
        choices=["cnn", "rf"],
        default="cnn",
        help="AI model to use for training (cnn or rf, default: cnn)",
    )
    train_parser.add_argument(
        "--normalization_mode",
        choices=["primitive", "rotation"],
        default="rotation",
        help="Normalization mode (primitive or rotation, default: rotation)",
    )
    train_parser.add_argument(
        "--ddcq",
        action="store_true",
        help="Apply ddcq (default: False)",
    )
    train_parser.add_argument(
        "--save_dir",
        default="models/",
        help="Directory to save the trained model (default: models/)",
    )

    # Subcommand: predict
    predict_parser = subparsers.add_parser(
        "predict", help="Make predictions using a trained model"
    )
    predict_parser.add_argument(
        "data_dir", help="Directory containing the PCR data files"
    )
    predict_parser.add_argument("model_path", help="Path to the trained model file")
    predict_parser.add_argument(
        "--ai_model",
        choices=["cnn", "rf"],
        default="cnn",
        help="AI model to use for prediction (cnn or rf, default: cnn)",
    )
    predict_parser.add_argument(
        "--normalization_mode",
        choices=["primitive", "rotation"],
        default="rotation",
        help="Normalization mode (primitive or rotation, default: rotation)",
    )
    predict_parser.add_argument(
        "--ddcq",
        action="store_true",
        help="Apply ddcq (default: False)",
    )

    # Subcommand: save_rfu_plot
    save_rfu_plot = subparsers.add_parser(
        "save_rfu_plot", help="Save RFU plots for PCR data"
    )
    save_rfu_plot.add_argument(
        "data_dir", help="Directory containing the PCR data files"
    )
    save_rfu_plot.add_argument(
        "--output_dir", help="Output directory. default=tempdir/"
    )

    args = parser.parse_args()

    # Validate model_path based on the mode
    import os

    # if args.mode == "train" and not os.path.isdir(args.model_path):
    #     parser.error("In train mode, model_path must be a directory.")
    if args.mode == "predict" and not os.path.isfile(args.model_path):
        parser.error("In predict mode, model_path must be a file.")

    return args


def load_data(runner, data_dir, normalization_mode="primitive", ddcq=False, exclude_files=None):
    """Load and preprocess data from the specified directory.
    Returns: (input_data, input_labels, targets, source_files)
    """
    runner.reset()
    if exclude_files is None:
        exclude_files = []
    
    # Normalize exclude_files to absolute paths for robust comparison
    exclude_files = [os.path.abspath(f) for f in exclude_files]
    
    # runner.read_dir replaced with manual iteration to support exclusion before loading/normalization
    file_list = sorted(os.listdir(data_dir))
    for fn in file_list:
        if not fn.lower().endswith(".csv"):
            continue
        full_path = os.path.abspath(os.path.join(data_dir, fn))
        if full_path in exclude_files:
            # print(f"Skipping already trained file: {fn}")
            continue
        
        try:
            runner.read_input(full_path)
        except ValueError as e:
            print(f"file name: {fn}, Error: {e}")

    targets = runner.get_target_names()
    input_data, input_labels, source_files = [], [], []
    for target in targets:
        runner.set_target(target)
        runner.normalization(smoothing_size=3, mode=normalization_mode)
        keys = list(runner.meta_info.keys())
        input_data.extend([runner.Cq_values[key] for key in keys])
        input_labels.extend([runner.meta_info[key]["Pos/Neg"] for key in keys])
        # Track which file each sample came from
        for key in keys:
            source_files.append(os.path.abspath(runner.file_names.get(key, "unknown")))

    if ddcq:
        input_data = np.diff(np.diff(input_data, axis=1), axis=1)
    
    # Handle case where no new data is loaded
    if len(input_data) == 0:
        return np.array([]), np.array([]), targets, []
        
    return np.array(input_data), np.array(input_labels), targets, source_files


def load_data_target(runner, target, normalization_mode="rotation", ddcq=False):
    """Load and preprocess data for a specific target."""
    runner.set_target(target)
    runner.normalization(smoothing_size=3, mode=normalization_mode)
    keys = list(runner.meta_info.keys())  # key는 샘플 번호임
    input_data = np.array([runner.Cq_values[key] for key in keys])
    input_labels = np.array([runner.meta_info[key]["Pos/Neg"] for key in keys])
    if ddcq:
        input_data = np.diff(np.diff(input_data, axis=1), axis=1)
    return (input_data, input_labels, target, keys)


def train_mode(args, runner):
    """Handle the training mode."""
    input_data, input_labels, _, _ = load_data(
        runner,
        args.data_dir,
        normalization_mode=args.normalization_mode,
        ddcq=args.ddcq,
    )
    pcr_learner = PCRPredicter()
    trained_model = pcr_learner.train_model(
        input_data,
        input_labels,
        learning_rate=args.learning_rate,
        epochs=args.epochs,
        save_dir=args.save_dir,
        ai_model=args.ai_model,
    )


def predict_mode(args, runner):
    """Handle the prediction mode."""
    pcr_learner = PCRPredicter()

    if args.ai_model == "rf":
        import joblib

        model = joblib.load(args.model_path)  # Load the Random Forest model
    elif args.ai_model == "cnn":
        pass
    else:
        raise ValueError("Invalid ai_model. Choose 'cnn' or 'rf'.")

    runner.read_dir(args.data_dir)
    targets = runner.get_target_names()
    for target in targets:
        input_data, input_labels, target, keys = load_data_target(
            runner, target, normalization_mode=args.normalization_mode, ddcq=args.ddcq
        )
        if args.ai_model == "cnn":
            predictions = pcr_learner.predict(input_data, args.model_path)
        elif args.ai_model == "rf":
            predictions = model.predict(input_data)
        for iii, key in enumerate(keys):
            actual = (
                "Positive"
                if runner.meta_info[key]["Pos/Neg"] == "Positive"
                else "Negative"
            )
            predicted = "Positive" if predictions[iii] == 1 else "Negative"
            file_name = runner.file_names[key]
            print(
                f"{target}-{key}\tActual: {actual}\tPredicted: {predicted}\tfile: {file_name}"
            )


def save_rfu_plot(args, runner):
    """Handle the find_baseline_start mode."""
    runner.read_dir(args.data_dir)
    targets = runner.get_target_names()

    # Generation of RFU plots for each target
    for target in targets:
        runner.set_target(target)
        keys = list(runner.meta_info.keys())

        input_data = np.array([runner.Cq_values[key] for key in keys])
        input_labels = np.array(
            [
                (
                    "Positive"
                    if "Positive" == runner.meta_info[key]["Pos/Neg"]
                    else "Negative"
                )
                for key in keys
            ]
        )

        pcr_learner = PCRPredicter()
        pcr_learner.save_rfu_plot_all(
            input_data, input_labels, target, keys, save_dir=args.output_dir
        )


def cluster_mode(args, runner):
    """Handle the clustering mode."""
    # Data Loading for multiple directories
    # Iterate over all provided directories and read them into the runner
    for d in args.data_dir:
        runner.read_dir(d)

    # Extract data (similar logic to load_data but post-read)
    targets = runner.get_target_names()
    input_data = []
    input_labels = [] # Stores strings like "No_Target_Time"
    
    for target in targets:
        runner.set_target(target)
        runner.normalization(smoothing_size=3, mode="primitive")
        
        keys = list(runner.meta_info.keys())
        input_data.extend([runner.Cq_values[key] for key in keys])
        
        # Generate labels: No. + Target + Date/Time
        for key in keys:
            # original_no = runner.meta_info[key].get("Original_No", key)
            timestamp = runner.gen_timestamps.get(key, "")
            # Clean timestamp if needed (e.g., regex already did some parsing, but filename might be better)
            # Actually timestamp comes from filename parsing in AutoLearner
            
            # label = f"{original_no}_{target}_{timestamp}"
            label = f"{target}_{timestamp}"
            input_labels.append(label)
    
    input_data = np.array(input_data)
    
    if len(input_data) == 0:
        print("No data found in the provided directories.")
        return
    
    # Perform Clustering
    clusterer = PCRCluster(
        n_clusters=args.n_clusters,
        method=args.method,
        early_weight=args.early_weight
    )
    
    # Preprocess
    # Since we used "primitive" normalization above in runner.normalization, we skip it here
    normalized_data, valid_indices = clusterer.preprocess(input_data, skip_normalization=True)
    
    if len(normalized_data) == 0:
        print("No valid data after preprocessing.")
        return

    # Filter labels if they exist and we have labels for all input
    if labels_to_plot := (input_labels if not args.no_labels else None):
        # input_labels is a list
        labels_to_plot = [input_labels[i] for i in valid_indices]
    
    # Fit & Predict
    labels = clusterer.fit_predict(normalized_data)
    
    # Visualize
    clusterer.visualize_clusters(normalized_data, labels, sample_labels=labels_to_plot, save_dir=args.output_dir)


def main():
    args = parse_arguments()
    runner = AutoLearner()

    if args.mode == "train":
        train_mode(args, runner)
    elif args.mode == "predict":
        predict_mode(args, runner)
    elif args.mode == "save_rfu_plot":
        save_rfu_plot(args, runner)
    elif args.mode == "cluster":
        cluster_mode(args, runner)


if __name__ == "__main__":
    main()
