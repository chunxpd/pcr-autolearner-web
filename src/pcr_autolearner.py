import os
import numpy as np
import pandas as pd
from src.pcr_analyzer import PCRmachine
from typing import Dict, Optional, Tuple
import copy


class CqValue:
    def __init__(self, fn: str):
        self.fn = fn
        self.gen_timestamp = ""
        self.cycles = []
        self.sample_names = {}
        self.Cq_values = {}
        self.meta_info = {}


class AutoLearner:
    def __init__(self):
        self.reset()  # Initialize all variables using the reset method

    def reset(self):
        """Reset all global variables in the AutoLearner class."""
        self.cq_container = []
        self.cycles = []
        self.sample_names = {}
        self.Cq_values = {}  # RFU values
        self.meta_info = {}
        self.Cq_init_value = {}
        self.gen_timestamps = []
        self.abnormal_files = []  # CSV 포맷이 이상한 파일들

    def parse_timestamp(self, input_fn: str) -> Optional[str]:
        """Parse the timestamp from the input file name.

        Args:
            input_fn (str): Path to the input CSV file.

        Returns:
            Optional[str]: Parsed timestamp or None if not found.
        """
        import re

        filename = os.path.basename(input_fn)
        # YYYYMMDD_HHMMSS 패턴 우선
        match = re.search(r"(\d{8}_\d{6})", filename)
        if match:
            return match.group(1)
        # 6자리 숫자 패턴 (예: 240601)
        match = re.search(r"(\d{6})", filename)
        if match:
            return match.group(1)
        return ""

    def check_file_version(self, input_fn: str) -> str:
        """Check the version of the input file based on its header.

        Args:
            input_fn (str): Path to the input CSV file.

        Returns:
            str: Version of the input file. (version1 or version2)
        """
        with open(input_fn, "r") as handle:
            # 첫 줄이 "X,"로 시작하지 않으면 이상한 파일로 간주
            first_line = handle.readline()
            if not first_line.lstrip().startswith("X,"):
                return "abnormal"

            # 이후 전체 내용 검사 (기존 동작 유지)
            content = first_line + handle.read()
            if "Threshold" in content and "RESULT" in content:
                return "version3"
            if any(
                keyword in content
                for keyword in ["-----", "=====", "SUBSTRACTION", "NORMALIZATION"]
            ):
                return "version2"
            else:
                return "version1"

    def read_input(self, input_fn: str):
        """Read the input file and determine its version.
        Then call the appropriate method to read the input data.

        Args:
            input_fn (str): Path to the input CSV file.
        """
        version = self.check_file_version(input_fn)
        if version == "version1":
            self.read_input_version1(input_fn)
        elif version == "version2":
            self.read_input_version2(input_fn)
        elif version == "version3":
            self.read_input_version3(input_fn)
        else:
            self.abnormal_files.append(input_fn)
            # raise ValueError(f"Unknown file version: {version}")

    def read_input_version1(self, input_fn: str):
        """After parsing the input_fn (csv file format),
        storing the Cq values and meta info
        into self.Cq_values={1:[], 2:[], ...}, and
        self.meta_info={1:{"Target Name":.., "Pos/Neg": ...}}

        Args:
            input_fn (str): Path to the input CSV file.
        """
        cq_value = CqValue(input_fn)
        processed_memo_entries = set()  # Track processed entries in the "Memo" section

        with open(input_fn, "r") as handle:
            header = handle.readline().strip()
            headers = header.split(",")
            cq_value.gen_timestamp = self.parse_timestamp(input_fn)
            cq_value.sample_names = {
                k: v for k, v in zip(range(1, len(headers)), headers[1:])
            }

            # First block: Read cycles and Cq values
            self._read_first_block(handle, cq_value)

            # Read meta info and ignore duplicates
            for line in handle:
                line = line.strip()
                if self._should_skip_line(line):
                    continue

                if "Num, Target Name" in line:
                    headers = line.split(",")
                    continue

                tmp_values = line.split(",")
                # Assuming the first column is the "Memo" entry
                memo_entry = tmp_values[0]
                if memo_entry in processed_memo_entries:
                    continue  # Skip duplicate entries
                processed_memo_entries.add(memo_entry)

                cq_value.meta_info[int(tmp_values[0])] = {
                    k.lstrip(): (
                        f"{tmp_values[0].zfill(2)}_{v.lstrip()}"
                        if k.lstrip() == "Target Name"
                        else v.lstrip()
                    )
                    for k, v in zip(headers[1:], tmp_values[1:])
                }

        self.cq_container.append(cq_value)

    def _read_first_block(self, handle, cq_value):
        """Read the first block of cycles and Cq values from the file handle."""
        for line in handle:
            line = line.strip()
            if line == "":
                break
            tmp_cq_values = line.split(",")
            if tmp_cq_values[0] == "":
                break  # Stop reading if the first value is empty (excel dummy line)
            cq_value.cycles.append(int(tmp_cq_values[0]))
            for iii, cur_cq in enumerate(tmp_cq_values[1:], start=1):
                try:
                    cq_value.Cq_values[iii].append(float(cur_cq))
                except KeyError:
                    cq_value.Cq_values[iii] = [float(cur_cq)]
        # Ensure the data has exactly 45 lines: Moim은 45 cycle 을 기본으로 사용. 하지만 많은 경우 40cycle에서 끝남
        while len(cq_value.cycles) < 45:
            last_cycle = cq_value.cycles[-1]
            cq_value.cycles.append(last_cycle + 1)  # Increment the cycle number
            for iii in cq_value.Cq_values:
                last_value = cq_value.Cq_values[iii][-1]
                cq_value.Cq_values[iii].append(last_value)  # Repeat the last Cq value

    def read_input_version2(self, input_fn: str):
        """After parsing the input_fn (csv file format),
        storing the Cq values and meta info
        into self.Cq_values={1:[], 2:[], ...}, and
        self.meta_info={1:{"Target Name":.., "Pos/Neg": ...}}

        Args:
            input_fn (str): Path to the input CSV file.
        """
        cq_value = CqValue(input_fn)
        processed_memo_entries = set()  # Track processed entries in the "Memo" section

        with open(input_fn, "r") as handle:
            header = handle.readline().strip()
            headers = header.split(",")
            cq_value.gen_timestamp = self.parse_timestamp(input_fn)
            cq_value.sample_names = {
                k: v for k, v in zip(range(1, len(headers)), headers[1:])
            }

            # First block: Read cycles and Cq values
            self._read_first_block(handle, cq_value)

            # skipping lines until meeting a keyword "RESULT"
            for line in handle:
                line = line.strip()
                if self._should_skip_line_v2(line):
                    continue
                if "No., Target, Result" not in line:
                    continue  # Skip until we find the header line
                if "No., Target, Result" in line:
                    headers = line.split(",")
                    headers = [
                        h.replace("No.", "Num")
                        .replace("Target", "Target Name")
                        .replace("Ct", "Ct Value")
                        .replace("Result", "Pos/Neg")
                        .lstrip()
                        for h in headers
                    ]
                    break  # Exit the loop after reading the header
            # Read meta info and ignore duplicates
            for line in handle:
                line = line.strip()
                if self._should_skip_line_v3(line):
                    continue
                if line == "":
                    break

                tmp_values = line.split(",")
                if "".join(tmp_values) == "":
                    break  # Stop reading if the line is empty (excel dummy line)
                memo_entry = tmp_values[
                    0
                ]  # Assuming the first column is the "Memo" entry

                if memo_entry in processed_memo_entries:
                    continue  # Skip duplicate entries
                processed_memo_entries.add(memo_entry)

                cq_value.meta_info[int(tmp_values[0])] = {
                    k.lstrip(): (
                        f"{tmp_values[0].zfill(2)}_{v.lstrip()}"
                        if k.lstrip() == "Target Name"
                        else v.lstrip()
                    )
                    for k, v in zip(headers[1:], tmp_values[1:])
                }

        self.cq_container.append(cq_value)

    def read_input_version3(self, input_fn: str):
        """After parsing the input_fn (csv file format),
        storing the Cq values and meta info
        into self.Cq_values={1:[], 2:[], ...}, and
        self.meta_info={1:{"Target Name":.., "Pos/Neg": ...}}

        Args:
            input_fn (str): Path to the input CSV file.
        """
        cq_value = CqValue(input_fn)
        processed_memo_entries = set()  # Track processed entries in the "Memo" section

        with open(input_fn, "r") as handle:
            header = handle.readline().strip()
            headers = header.split(",")
            cq_value.gen_timestamp = self.parse_timestamp(input_fn)
            cq_value.sample_names = {
                k: v for k, v in zip(range(1, len(headers)), headers[1:])
            }

            # First block: Read cycles and Cq values
            self._read_first_block(handle, cq_value)

            # skipping lines until meeting a keyword "RESULT"
            for line in handle:
                line = line.strip()
                if self._should_skip_line_v2(line):
                    continue
                if "No., Target, Result" not in line:
                    continue  # Skip until we find the header line
                if "No., Target, Result" in line:
                    headers = line.split(",")
                    headers = [
                        h.replace("No.", "Num")
                        .replace("Target", "Target Name")
                        .replace("Ct", "Ct Value")
                        .replace("Result", "Pos/Neg")
                        .lstrip()
                        for h in headers
                    ]
                    break  # Exit the loop after reading the header
            # Read meta info and ignore duplicates
            for line in handle:
                line = line.strip()
                if self._should_skip_line_v3(line):
                    continue
                if line == "":
                    break

                tmp_values = line.split(",")
                if "".join(tmp_values) == "":
                    break  # Stop reading if the line is empty (excel dummy line)
                memo_entry = tmp_values[
                    0
                ]  # Assuming the first column is the "Memo" entry

                if memo_entry in processed_memo_entries:
                    continue  # Skip duplicate entries
                processed_memo_entries.add(memo_entry)

                try:
                    idx = int(tmp_values[0])
                except ValueError:
                    continue

                cq_value.meta_info[idx] = {
                    k.lstrip(): (
                        f"{tmp_values[0].zfill(2)}_{v.lstrip()}"
                        if k.lstrip() == "Target Name"
                        else v.lstrip()
                    )
                    for k, v in zip(headers[1:], tmp_values[1:])
                }

        self.cq_container.append(cq_value)

    def _should_skip_line(self, line: str) -> bool:
        """Check if a line should be skipped based on specific keywords or conditions."""
        skip_keywords = [
            "Criteria Normalization",
            "Criteria Baseline Fitting",
            "Memo",
            "Scenario",
            "Elapsed",
            "Start Time",
            "End   Time",
            "------",
            "======",
            "SUBSTRACTION",
            "NORMALIZATION",
        ]
        return any(keyword in line for keyword in skip_keywords) or line == ""

    def _should_skip_line_v2(self, line: str) -> bool:
        """파일내에서 RESULT 라는 키워드가 나오면 멈춘다"""
        return "RESULT" in line

    def _should_skip_line_v3(self, line: str) -> bool:
        """Check if a line should be skipped based on specific keywords or conditions."""
        skip_keywords = ["Baseline Fitting", "Normalization"]
        return any(keyword in line for keyword in skip_keywords)

    def read_dir(self, input_dir: str):
        """input_dir 에 있는 모든 파일을 읽어서
        self.read_input() 함수를 호출하는 함수

        Args:
            input_dir (str): input directory
        """
        file_list = sorted(os.listdir(input_dir))
        for fn in file_list:
            if not fn.lower().endswith(".csv"):
                continue
            try:
                self.read_input(f"{input_dir}/{fn}")
            except ValueError as e:
                print(f"file name: {fn}, Error: {e}")
                raise ValueError

    def get_target_names(self) -> list:
        """return a list of target names

        Returns:
            list: list of target names
        """
        res = []
        for rec in self.cq_container:
            res.extend([rec.meta_info[i]["Target Name"] for i in rec.meta_info])
        return sorted(list(set(res)))

    def set_target(self, target_name: str):
        """
        Sets the target name and updates the corresponding attributes from original values of self.cq_container

        This method filters the records in `self.cq_container` based on the provided
        `target_name` and updates the following attributes:
        - `self.cycles`: List of cycles for the matching records.
        - `self.sample_names`: Dictionary mapping an index to sample names for the matching records.
        - `self.Cq_values`: Dictionary mapping an index to Cq values for the matching records.
        - `self.meta_info`: Dictionary mapping an index to meta information for the matching records.

        Args:
            target_name (str): The target name to filter records by.
        """
        self.cycles = []
        self.sample_names = {}
        self.Cq_values = {}
        self.meta_info = {}
        self.Cq_init_value = {}
        self.gen_timestamps = {}
        self.file_names = {}
        iii = 1
        for rec in self.cq_container:
            for key in list(rec.meta_info.keys()):
                if rec.meta_info[key]["Target Name"] != target_name:
                    continue
                self.cycles = copy.deepcopy(rec.cycles)
                self.sample_names[iii] = copy.deepcopy(rec.sample_names[key])
                self.Cq_values[iii] = copy.deepcopy(rec.Cq_values[key])
                self.meta_info[iii] = copy.deepcopy(rec.meta_info[key])
                # self.meta_info[iii]["Original_No"] = key  # Store original sample number
                self.Cq_init_value[iii] = self.Cq_values[iii][0]
                self.gen_timestamps[iii] = rec.gen_timestamp
                self.file_names[iii] = rec.fn
                iii += 1
        if len(self.Cq_values) == 0:
            raise ValueError(
                f"Target name '{target_name}' not found in the input data."
            )

    def get_expected_value(self, col_num: int, key: str = None) -> str:
        """from self.meta_info, return a value of a given key of the corresponding colum number

        Args:
            col_num (int): input sample's number. usually start from 1
            key (str): key is one of  "Target Name", "Th_nor", "Th_sub", "Cutoff", "Ct Value", and "Pos/Neg"

        Returns:
            str:
        """
        return self.meta_info[col_num][key]

    def get_sample_size(self):
        return len(self.sample_names)

    def get_pos_num(self):
        res = 0
        for item_key in self.meta_info:
            if "Positive" == self.get_expected_value(item_key, key="Pos/Neg"):
                res += 1
        return res

    def baseline_endq(
        self,
        sample_number: int,
        cycle_range: tuple,
        smoothing_size: int = 3,
        scaling_factor: float = 0.01,
    ) -> int:
        """src.pcr_analyzer 에 있는 find_baseline_end() 함수를 이용하는 wrapping 함수

        Args:
            sample_number (int): 샘플 열번호 (1부터 시작)
            cycle_range (tutple): (ignore_cycle, max_cycle)

        Returns:
            int: baseline_end
        """
        res = 0
        pcr_analyzer = PCRmachine(max_cycle=len(self.Cq_values[1]))
        for cycle, cq_value in enumerate(
            self.Cq_values[sample_number], start=self.cycles[0]
        ):
            pcr_analyzer.add_one_cq_value(
                cq_value, set_winsize=smoothing_size, start=self.cycles[0]
            )
            pcr_analyzer.pcr_machine_cal_dCq()
        adj_cq_values = pcr_analyzer.rotate_cqs(
            cycle_range[1], ignore_init_cycles=cycle_range[0]
        )
        baseline_end = pcr_analyzer.find_baseline_end(
            ignore_init_cycles=cycle_range[0],
            ddcq_threshold_scale=scaling_factor,
            max_cycle=cycle_range[1],
        )
        if cycle_range[0] <= baseline_end <= cycle_range[1]:
            res = 1
        return res

    def cal_loss_function(
        self, cycle_range: tuple, scaling_factor: float = 0.01
    ) -> float:
        """https://github.com/io2oi/ctlabs-primerwrapper/issues/264#issue-2804358209
        calculating loss function
        L = \frac{1}{K} \sum_{i}^{N}{y_i - e_i}
        테스트입력에서 Positive 12개, Negative 1개, no label 3개임.
        그런데 Negative 가 실제로는 Positive
        여기서 Loss function이 0이 되도록 cycle range, scaling factor 등을 조정해 나가야 함

        Args:
            cycle_range (tuple): hyperparameter
            scaling_factor (float, optional): hyperparameter. Defaults to 0.01.

        Returns:
            float: loss function 값
        """
        K = self.get_pos_num()
        N = self.get_sample_size()
        sum = 0
        for iii in self.meta_info:
            y_i = 1 if self.get_expected_value(iii, key="Pos/Neg") == "Positive" else 0
            e_i = self.baseline_endq(
                iii, cycle_range=cycle_range, scaling_factor=scaling_factor
            )
            sum += y_i - e_i
        res = sum / K
        return res

    def iteration_training_scaling_factor(
        self,
        cycle_range: tuple,
        epsilon: float = 0.0001,
        max_iteration: int = 100,
        step_size: float = 1e-3,
    ) -> dict:
        """ "https://github.com/io2oi/ctlabs-primerwrapper/issues/267
        loss function 이 0이 아니면 답이 틀렸다는 이야기.
        - 일 경우: negative 를 positive 로 인식 했다는 뜻. scaling factor 를 늘려본다.
        + 일 경우: Positive 를 negative 로 인식 했다는 뜻. scaling factor 를 줄여본다.
        scaling factor 의 범위는 0.0 (0%) ~ 0.1 (10%)로 고정한다.
        단, scaling factor 가 0.0 이하로 내려가면 안된다. 그리고 scaling factor 가 충분히 작을때에는 증가하는 방향으로만
        학습이 진행되어야 한다. 왜냐하면 baseline end 가 3보다 더 작은 곳에서 잡히기 때문이다.
        Args:
            cycle_range (tuple): hyper parameter. default = (3,35)
            epsilon (float): loss function 이 saturation 되었다고 인정하는 값. default 0.0001
            max_iteration (int): iteration 의 최대 값. default = 100
            step_size (float): scaling_factor 의 변화량

        Returns:
            dict: range of optimal scaling_factor, range of loss fuction for each scaling_facor_min and max
        """
        res = {"scaling_factor": None, "loss_value": None}
        scaling_factor_min = 0.0 + step_size * 2
        scaling_factor_max = 0.11
        iteration_count = 0
        while iteration_count < max_iteration:
            loss_min = self.cal_loss_function(
                cycle_range=cycle_range, scaling_factor=scaling_factor_min
            )
            print(
                f"Iteration {iteration_count}: scaling_factor_min = {scaling_factor_min}, loss_min = {loss_min}"
            )
            if abs(loss_min) < epsilon:
                break
            scaling_factor_min = scaling_factor_min + step_size
            iteration_count += 1  # 추가: 반복 횟수 추적

        iteration_count = 0
        while iteration_count < max_iteration:
            loss_max = self.cal_loss_function(
                cycle_range=cycle_range, scaling_factor=scaling_factor_max
            )
            print(
                f"Iteration {iteration_count}: scaling_factor_max = {scaling_factor_max}, loss_max = {loss_max}"
            )
            if abs(loss_max) < epsilon:
                break
            if loss_max > 0:  # Positive loss: scaling factor 감소
                scaling_factor_max = scaling_factor_max - step_size
            else:  # Negative loss: scaling factor 증가
                scaling_factor_max = scaling_factor_max + step_size
            iteration_count += 1  # 추가: 반복 횟수 추적

        res["scaling_factor"] = (
            max(scaling_factor_min, 0.0),
            min(scaling_factor_max, 0.1),
        )
        loss_min = self.cal_loss_function(
            cycle_range=cycle_range, scaling_factor=scaling_factor_min
        )
        loss_max = self.cal_loss_function(
            cycle_range=cycle_range, scaling_factor=scaling_factor_max
        )
        res["loss_value"] = (loss_min, loss_max)  # 결과 저장

        return res

    def find_max_range(self, cycle_loss_values) -> tuple:
        """[(cycle_range, loss_value), ...]
        가 있을 때 loss_value 의 절대값이 작은 것부터 큰 것까지
        그리고 최대 범위의 cycle_range 를 반환

        Args:
            cycle_loss_values (_type_): _description_

        Returns:
            tuple: _description_
        """
        in_list = sorted(cycle_loss_values, key=lambda x: abs(x[1]))
        min_loss_value = in_list[0][1]
        tmp_cont = []
        for cycle_range, loss_value in in_list:
            if abs(min_loss_value) == abs(loss_value):
                tmp_cont.append(cycle_range)
        tmp_cont = sorted(tmp_cont, key=lambda x: x[1] - x[0])
        return (tmp_cont[-1], min_loss_value)

    def iteration_training_cycle_range(self, given_scaling_factors: tuple):
        """_summary_

        Args:
            given_scaling_factors (tuple): scaling factor min and max
        """
        cont_min = []
        cont_max = []
        for iii in range(3, 5):  ### init cycle: 3~5
            for jjj in range(33, self.cycles[-1]):  ### end cycle: 33~input. usually 40
                tmp_loss_value = self.cal_loss_function(
                    cycle_range=(iii, jjj), scaling_factor=given_scaling_factors[0]
                )
                cont_min.append(((iii, jjj), tmp_loss_value))
                tmp_loss_value = self.cal_loss_function(
                    cycle_range=(iii, jjj), scaling_factor=given_scaling_factors[1]
                )
                cont_max.append(((iii, jjj), tmp_loss_value))
        cycle_range_min, loss_value_min = self.find_max_range(cont_min)
        cycle_range_max, loss_value_max = self.find_max_range(cont_max)

        res = {
            "cycle_range": (cycle_range_min, cycle_range_max),
            "loss_value": (loss_value_min, loss_value_max),
        }
        return res

    def baseline_substraction(self, smoothing_size: int = 3, mode: str = "rotation"):
        """
        Subtracts the baseline Cq value from each sample's Cq values.
        선택된 item target name 에 대한 Cq 값들에 대해서
        1. 임시 baseline 을 3~8 번째 cycle 을 기준으로 잡는다.
        2. pcr_analyzer.add_one_cq_value()를 이용해서 값을 넣은 뒤
        'Cq_smooting' 값을 이용해서 baseline 의 평균 값을 구하고
        3. rotation 을  5번째 cycle 을 중심으로 해준다
        4. 그런다음 self.Cq_values[item]에 substraction 값을 넣는다.

        mode = "rotation" 이면 3~8 CYCLE 의 값을 이용해서 baseline을 잡고 5번째 cycle로 rotation을 하고
        mode = "primitive" 이면 첫 번째 cycle의 값을 이용해서 baseline을 잡고 rotation을 하지 않는다.
        Args:
            smoothing_size (int, optional): Smoothing window size. Defaults to 3.
            mode (str, optional): Mode of processing. Defaults to "rotation".
        """
        # Create a single PCRmachine instance
        pcr_analyzer = PCRmachine(max_cycle=len(self.Cq_values[1]))

        for key in self.Cq_values:
            if mode == "rotation":
                tmp_horizontal_regions = self.find_horizontal_regions(
                    self.Cq_values[key], self.cycles
                )
                if tmp_horizontal_regions == []:
                    tmp_base_start = 3
                    tmp_base_end = 8
                    tmp_rotation_center = tmp_base_end
                    print(
                        f"Sample {key}: No horizontal regions found. horizontal region=({tmp_base_start}~{tmp_base_end}"
                    )
                else:
                    # Use the first horizontal region for baseline calculation
                    tmp_base_start, tmp_base_end = tmp_horizontal_regions[0]
                    tmp_rotation_center = tmp_base_end
                    print(
                        f"Sample {key}: horizontal region=({tmp_base_start}~{tmp_base_end})"
                        f" second region=({tmp_horizontal_regions[1] if len(tmp_horizontal_regions) > 1 else 'N/A'})"
                    )
                self._process_sample(
                    key,
                    pcr_analyzer,
                    tmp_base_start,
                    tmp_base_end,
                    tmp_rotation_center,
                    smoothing_size,
                )
            else:
                self._process_sample_primitive(
                    key,
                    pcr_analyzer,
                    smoothing_size,
                )

    def _process_sample(
        self, key, pcr_analyzer, base_start, base_end, rotation_center, smoothing_size
    ):
        """
        Processes a single sample by resetting the PCRmachine, calculating dCq,
        rotating the Cq values, and subtracting the baseline.

        Args:
            key (int): The sample key.
            pcr_analyzer (PCRmachine): The PCRmachine instance.
            base_start (int): Start cycle for baseline calculation.
            base_end (int): End cycle for baseline calculation.
            rotation_center (int): Cycle for rotation.
            smoothing_size (int): Smoothing window size.
        """
        pcr_analyzer.reset()
        self._add_cq_values_to_analyzer(key, pcr_analyzer, smoothing_size)
        adj_cq_values = pcr_analyzer.rotate_cqs(
            rotation_center, ignore_init_cycles=base_start
        )
        baseline = self._calculate_baseline(adj_cq_values, base_start, base_end)
        self.Cq_values[key] = [
            value - baseline for value in adj_cq_values["Cq_smoothing"]
        ]

    def _process_sample_primitive(self, key, pcr_analyzer, smoothing_size):
        """
        Processes a single sample by resetting the PCRmachine and adding Cq values.
        and subsstracting the baseline. In this primitive version, we do not use rotation and baseline.
        Just substract the first cycle value from all the Cq values.

        Args:
            key (int): The sample key.
            pcr_analyzer (PCRmachine): The PCRmachine instance.
            smoothing_size (int): Smoothing window size.
        """
        pcr_analyzer.reset()
        self._add_cq_values_to_analyzer(key, pcr_analyzer, smoothing_size)
        adj_cq_values = pcr_analyzer.get_cq_values()
        baseline = adj_cq_values["Cq_smoothing"][0]
        self.Cq_values[key] = [
            value - baseline for value in adj_cq_values["Cq_smoothing"]
        ]

    def _add_cq_values_to_analyzer(self, key, pcr_analyzer, smoothing_size):
        """
        Adds Cq values to the PCRmachine instance and calculates dCq.

        Args:
            key (int): The sample key.
            pcr_analyzer (PCRmachine): The PCRmachine instance.
            smoothing_size (int): Smoothing window size.
        """
        for cycle, cq_value in enumerate(self.Cq_values[key], start=self.cycles[0]):
            pcr_analyzer.add_one_cq_value(
                cq_value, set_winsize=smoothing_size, start=self.cycles[0]
            )
        try:
            pcr_analyzer.pcr_machine_cal_dCq()
        except ValueError as e:
            print(self.Cq_values[key])
            print(f"Error file: {self.file_names[key]}")
            raise ValueError(f"Error in sample {key}: {e}")

    def _calculate_baseline(self, adj_cq_values, base_start, base_end):
        """
        Calculates the baseline value using the adjusted Cq values.

        Args:
            adj_cq_values (dict): Adjusted Cq values from the PCRmachine.
            base_start (int): Start cycle for baseline calculation.
            base_end (int): End cycle for baseline calculation.

        Returns:
            float: The calculated baseline value.
        """
        return sum(adj_cq_values["Cq_smoothing"][base_start : base_end + 1]) / (
            base_end - base_start + 1
        )

    def normalization(self, smoothing_size: int = 3, mode="rotation"):
        """
        baseline substraction 이후에 normalization 을 해준다.
        normalization 은 첫번째 cycle의 값을 이용해서 모든 baseline substraction 값들을 나누어 주면 된다.
        첫번째 cycle 의 값이 0이 되는 경우는 normalization 을 하지 않는다.
        """
        self.baseline_substraction(smoothing_size=smoothing_size, mode=mode)
        for key in self.Cq_values:
            if self.Cq_init_value[key] == 0:
                print(
                    f"Sample {key}: Initial Cq value is zero. Skipping normalization."
                )
                continue
            self.Cq_values[key] = [
                value / self.Cq_init_value[key] for value in self.Cq_values[key]
            ]

    def find_horizontal_regions(self, cq_values, cycles, min_length=5):
        """
        Detect horizontal regions in the cq_values.
        cq_values 는 RFU cq_values 를 의미하며, cycles 는 해당 cq_values 의 cycle number 를 의미한다.
        horizontal region 은 연속된 cycle 들의 cq_values 의 변화가 일정 threshold 이하로 유지되는 구간을 의미한다.
        horizontal region 의 시작과 끝은 cycle number 로 표현되며, 연속된 cycle 들의 cq_values 의 변화가 일정 threshold 이하로 유지되는 구간을 의미한다.
        단, 처음 0~3 cycle 은 무시하도록 한다.

        Args:
            cq_values (list or np.ndarray): The RFU cq_values for a single sample.
            cycles (list or np.ndarray): The cycle numbers corresponding to cq_values.
            min_length (int): The minimum number of consecutive cycles to consider a region horizontal.

        Returns:
            list: A list of tuples, where each tuple represents a horizontal region
                as (start_cycle, end_cycle).
        """
        ignore_cycles = 3  # Ignore the first 3 cycles
        horizontal_regions = []
        start_cycle = None
        threshold = 0.025 * (
            max(cq_values[ignore_cycles:]) - min(cq_values[ignore_cycles:])
        )
        for i in range(ignore_cycles, len(cq_values) - 1):
            # Check if the difference between consecutive values is below the threshold
            if abs(cq_values[i + 1] - cq_values[i]) <= threshold:
                if start_cycle is None:
                    start_cycle = i  # Mark the start of a horizontal region
            else:
                if start_cycle is not None:
                    # End the horizontal region
                    if i - start_cycle + 1 >= min_length:
                        horizontal_regions.append((start_cycle, i))
                    start_cycle = None

        # Handle the case where the last region extends to the end of the cq_values
        if start_cycle is not None and len(cq_values) - start_cycle >= min_length:
            horizontal_regions.append((start_cycle, len(cq_values) - 1))

        horizontal_regions = [
            (cycles[start], cycles[end]) for start, end in horizontal_regions
        ]
        return horizontal_regions

    def sigmoid_fitting(
        self, x: list, y: list, function="logistic"
    ) -> Dict[str, float]:
        """Fit a sigmoid function to the data points (x, y).
        function logistic = y0 + a / (1 + (x/x0)**b)
        Args:
            x (list): List of x values (e.g., cycle numbers).
            y (list): List of y values (e.g., Cq values).

        Returns:
            Dict: Parameters of the fitted sigmoid function (y0, a, x0, b).
        """
        from scipy.optimize import curve_fit

        def sigmoid_logistic(x, y0, a, x0, b):
            return y0 + a / (1 + (x / x0) ** b)

        def sigmoid_naive(x, y0, a, x0, b):
            return y0 + a / (1 + np.exp(-b * (x - x0)))

        def sigmoid_gompertz(x, y0, a, x0, b):
            return y0 + a * np.exp(-np.exp(-b * (x - x0)))

        if function == "logistic":
            sigmoid_like = sigmoid_logistic
        elif function == "naive":
            sigmoid_like = sigmoid_naive
        elif function == "gompertz":
            sigmoid_like = sigmoid_gompertz
        else:
            raise ValueError(
                "Invalid function type. Choose 'logistic', 'naive', or 'gompertz'."
            )

        # 초기값 설정: y0=최소값, a=최대-최소, x0=중간 x, b=1
        y0_init = min(y)
        a_init = max(y) - min(y)
        x0_init = np.median(x)
        b_init = 1.0
        p0 = [y0_init, a_init, x0_init, b_init]

        try:
            popt, _ = curve_fit(sigmoid_like, x, y, p0=p0, maxfev=5000)
            y0, a, x0, b = popt
            res = {
                "y0": round(y0, 3),
                "a": round(a, 3),
                "x0": round(x0, 3),
                "b": round(b, 3),
            }
        except RuntimeError as e:
            res = {
                "y0": None,
                "a": None,
                "x0": None,
                "b": None,
            }

        return res

    def get_key_from_label(self, label: str):
        """레이블에서 앞부분의 정수 키(kkk)를 추출하여 int로 반환.
        예: "3_20240529_120000" -> 3, "3 (tr_P:pr_N)" -> 3.
        찾지 못하면 None 반환."""
        if not label:
            return None
        import re

        # 선행 공백 후 바로 오는 숫자(가장 일반적인 케이스)
        m = re.match(r"\s*(\d+)", label)
        if m:
            return int(m.group(1))

        # '_' 또는 공백으로 분할했을 때 첫 토큰이 숫자이면 반환
        first = re.split(r"[_\s]+", label)[0]
        if first.isdigit():
            return int(first)

        # 어떠한 위치의 숫자도 허용하려면 마지막 수를 반환(옵션)
        m2 = re.search(r"(\d+)", label)
        if m2:
            return int(m2.group(1))

        return None


if __name__ == "__main__":
    pass
