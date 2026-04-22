import math
from typing import Tuple, Union
import numpy as np
import pandas as pd
import statsmodels.api as sm


class PCRmachine:
    def __init__(self, min_ddcq_pos_num: int = 4, max_cycle: int = 40):
        """Cq_ori: PCR 기기에서 생산되는 fluorescence 값
        Cq: 회전변환할 때 각도를 미세하게 조절하기 위해서 Cq_ori 의 첫번째 cycle 값으로 나누어준 값
        Args:
            min_ddcq_pos_num (int, optional): baseline 의 끝 cycle 을 잡기 위한 ddCq 값의 연속 횟수. Defaults to 4.
        """
        self.min_ddcq_pos_num = min_ddcq_pos_num
        self.channel = None  # dye 종류
        self.cq_values = {
            "cycle": [],
            "Cq": [],
            "Cq_ori": [],
            "Cq_smoothing": [],
            "dCq": [],
            "ddCq": [],
            "dddCq": [],
        }
        self.cq_smoothing = {"winsize3": [], "winsize5": []}
        self.max_cycle = max_cycle

    def reset(self):
        """Reset the internal state of the PCRmachine."""
        self.cq_values = {
            "cycle": [],
            "Cq": [],
            "Cq_ori": [],
            "Cq_smoothing": [],
            "dCq": [],
            "ddCq": [],
            "dddCq": [],
        }
        self.cq_smoothing = {"winsize3": [], "winsize5": []}

    def get_cq_values(self):
        """각 cycle 별 Cq 값을 반환하는 메서드 (실제 구현)"""
        return self.cq_values

    def moving_average(self, sequence: list, winsize: int = 3) -> float:
        """python 리스트에 값을 하나씩 넣을 때 이전에 들어온 값 2개를 새로 들어오늘 값에 더하고 평균을 한 것을 지금 값으로 하는 코드

        Args:
            sequence (list): "Cq_ori"
            winsize (list): 3 or 5

        Returns:
            float: smoothing 값
        """
        res = None
        # 리스트가 3개 이상일 때만 마지막 3개의 평균을 구함
        if len(sequence) < winsize:
            res = sum(sequence) / len(sequence)  # 값이 3개 미만이면 가능한 값들로 평균
        else:
            # 마지막 3개 값의 평균 계산
            res = sum(sequence[-winsize:]) / winsize
        return res

    def add_one_cq_value(self, value: float, set_winsize: int = 3, start=1) -> bool:
        """PCR 기계에서 Cq 값이 하나씩 들어옴에 따라 값을 저장하고 dCq 값들을 계산한다.
        Cq 값은 항상 흔들림이 있기 때문에 window size 3으로 smoothing 한다.
        Cq_i = [Cq_(i-2) + Cq_(i-1) + Cq_(i)]/3
        첫 cycle의 Cq값이 0이면 tmp_smoothing_cq3 값을 그대로 사용한다.

        Args:
            value (value): Cq 값 한개

        Returns:
            bool: 에러가 발생하면 False, 그 외에는 True
        """
        cycle_offset = start - 1
        try:
            cycle_num = len(self.cq_values["cycle"]) + cycle_offset
            self.cq_values["cycle"].append(cycle_num + 1)
            self.cq_values["Cq_ori"].append(value)
            tmp_smoothing_cq3 = self.moving_average(self.cq_values["Cq_ori"], winsize=3)
            self.cq_smoothing["winsize3"].append(tmp_smoothing_cq3)
            tmp_smoothing_cq5 = self.moving_average(self.cq_values["Cq_ori"], winsize=5)
            self.cq_smoothing["winsize5"].append(tmp_smoothing_cq5)
            if set_winsize == 3:
                self.cq_values["Cq_smoothing"] = self.cq_smoothing["winsize3"]
                tmp_init_value = self.cq_values["Cq_smoothing"][0]
                if tmp_init_value != 0:
                    self.cq_values["Cq"].append(tmp_smoothing_cq3 / tmp_init_value)
                else:
                    self.cq_values["Cq"].append(tmp_smoothing_cq3)
            elif set_winsize == 5:
                self.cq_values["Cq_smoothing"] = self.cq_smoothing["winsize5"]
                tmp_init_value = self.cq_values["Cq_smoothing"][0]
                if tmp_init_value != 0:
                    self.cq_values["Cq"].append(tmp_smoothing_cq5 / tmp_init_value)
                else:
                    self.cq_values["Cq"].append(tmp_smoothing_cq5)
            self.cq_values["dCq"].append(None)
            self.cq_values["ddCq"].append(None)
            self.cq_values["dddCq"].append(None)
        except Exception:
            return False
        return True

    def pcr_machine_cal_dCq(self):
        """Cq 값으로 부터 dCq, ddCq, dddCq 값을 구하는 메소드
        self.cq_values 를 기반으로 계산한다.
        Return:
            res (list): self.cq_values 에 값을 업데이트 하고 그 값을 그대로 반환
        """
        df = pd.DataFrame(self.cq_values)
        df["dCq"] = df["Cq_smoothing"].diff()
        df["ddCq"] = df["dCq"].diff()
        df["dddCq"] = df["ddCq"].diff()
        self.cq_values.update(df.to_dict(orient="list"))
        return self.get_cq_values()

    def cal_ddcq_threshold(
        self, cq_list: list = None, scale_factor: float = 0.001, max_cycle: int = None
    ) -> float:
        """(min-max)*scale_factor

        Args:
            cq_list (list): _description_
            scale_factor (float): 0.1% default

        Returns:
            float: _description_
        """
        if cq_list is None:
            cq_list = self.get_cq_values()["Cq_smoothing"]
        if max_cycle is None:
            diff = max(cq_list) - min(cq_list)
        else:
            diff = max(cq_list[:max_cycle]) - min(
                cq_list[:max_cycle]
            )  # cycle 은 1에서 부터 시작하므로 +1 을 해 줄 필요 없음
        res = diff * scale_factor
        return res

    def find_consecutive_regions(
        self, in_list: list, consecutive_number: int = 4, T: float = 0
    ) -> list:
        """입력 in_list 에서 T 보다 큰 값이 연속으로 consecutive_number 보다 크거나 같은 영역을 찾아주는 함수.
        (start,end) pair 에 대한 list 가 반환값임. python은 기본이 0-based index 이기 때문에
        값들은 in_list[start:end+1] 로 해야 모두 splicing이 가능함

        Args:
            in_list (list): ddCq list
            consecutive_number (int, optional): n(ddCq > T) >= consecutive_number. Defaults to 4.
            T (float, optional): 바로 위의 threshold. Defaults to 0.

        Returns:
            list: [(start, end), ...] 0 based index
        """
        start = None  # 연속 구간 시작 인덱스
        result = []  # 결과 저장 리스트

        for i, val in enumerate(in_list):
            if not math.isnan(val) and val > T:
                if start is None:  # 새로운 연속 구간 시작
                    start = i
            else:
                if (
                    start is not None and i - start >= consecutive_number
                ):  # 조건 만족 구간 저장
                    result.append((start, i - 1))
                start = None  # 연속 구간 종료

        # 리스트 끝까지 탐색 후 처리
        if start is not None and len(in_list) - start >= consecutive_number:
            result.append((start, len(in_list) - 1))

        return result

    def find_all_negative_dcqs(self, dcq_list: list, start: int, end: int) -> bool:
        """from input dcq_list, confirm that all dcq_list[start:end+1] < 0

        Args:
            dcq_list (list): _description_
            start (int): 0-based start index
            end (int): 0-based end index

        Returns:
            bool: all dcq < 0 => return True, else return False
        """
        res = False
        obs_dcqs = dcq_list[start : end + 1]
        res = all(x < 0 for x in obs_dcqs)
        return res

    def find_all_positive_dcqs(self, dcq_list: list, start: int, end: int) -> bool:
        """from input dcq_list, confirm that all dcq_list[start:end+1] >= 0

        Args:
            dcq_list (list): _description_
            start (int): 0-based start index
            end (int): 0-based end index

        Returns:
            bool: all dcq >= 0 --> return True, else return False
        """
        res = False
        obs_dcqs = dcq_list[start : end + 1]
        res = all(x >= 0 for x in obs_dcqs)
        return res

    def find_baseline_end(
        self,
        consecutive_number: int = 4,
        ignore_init_cycles: int = 3,
        ddcq_threshold_scale: float = 0.01,
        max_cycle: int = None,
    ) -> Union[bool, int]:
        """현재까지 업데이트된 self.cq_values 값을 기반으로 계산한다.
        ddCq 가 4연속 나올 때 최초 cycle 로부터 2단계 앞을 baseline end 로 계산해서 반환한다.
        아직 그런 값이 없을 때는 False 반환
        ddCq 가 4연속 나오거나 나오지 않으면
        self.max_cycle - consecutive_number
        즉 max_cycle이 40이라면 baseline_end 의 최대값은 36

        Args:
            consecutive_number (int): ddCq > 0 일 연속 갯수의 값. 사용자가 조정 가능하다. default = 4
            max_cycle (int): PCR cycle 의 최대 값. 보통 40으로 설정하는데 그보다 큰 수도 있고 사용자가 40에 도달하기 전에 인위적으로 멈출 수 있음
            ignore_init_cycles (int): PCR cycle 에서 흔들리는 첫 cyle 수. baseline_end 가 이것보다 작거나 같으면 건너뛴다.
            ddcq_threshold_scle (float): dynamic T 를 정하기 위함. scale* (min(Cq_smoothing)-max(Cq_smoothing), default 0.01 (1%)
            max_cycle (int): 1부터 시작하는 cycle 번호. 여기까지 cycle 만을 고려해서 max(ddCq) 와 min(ddCq) 를 정한다.
        Return:
            res (bool): baseline end (int, 1부터 시작하는 값), 기본은 False
        """
        res = False
        cq_values = self.get_cq_values()
        ddcq_list = cq_values["ddCq"]
        baseline_end_upperbound = self.max_cycle - consecutive_number
        T = self.cal_ddcq_threshold(
            cq_values["Cq_smoothing"],
            scale_factor=ddcq_threshold_scale,
            max_cycle=max_cycle,
        )
        print(f"Threshold: {T}")
        if len(ddcq_list) < consecutive_number:
            return res
        try:
            reg_start, reg_end = self.find_consecutive_regions(
                ddcq_list, consecutive_number=consecutive_number, T=T
            )[-1]
            obs_cons_num = reg_end - reg_start + 1
        except:  # 영역이 없을 때
            obs_cons_num = 0
        if res == False and len(ddcq_list) >= baseline_end_upperbound:
            # baseline_end 가 없고 이미 너무 cycle 이 많이 지났음
            res = baseline_end_upperbound
        if obs_cons_num >= consecutive_number and self.find_all_positive_dcqs(
            cq_values["dCq"], reg_start, reg_end
        ):
            # baseline_end 계산
            # region 시작 -2 가 baseline_end 이고 list 인덱스와 달리 cycle은 1에서 시작하므로 -2 + 1
            res = reg_start - 2 + 1
            if res > baseline_end_upperbound:
                # baseline_end가 너무 크면
                res = baseline_end_upperbound
            if res <= ignore_init_cycles:
                # baseline_end 가 너무 작음
                res = False
        return res

    def fining_baseline_end(
        self, crude_baseline_end: int, ignore_init_cycles: int = 3
    ) -> int:
        """일차적으로 찾은 baseline end 를 가지고 좀 더 정확하게 찾아본다. 현재 제안된 알고리즘
        1. linear regression 을 한다.
        2. baseline +-2 에서 linear regression 을 한다
        3. 가장 작은 p-value 를 가지는 baseline 을 최종 결정한다.

        Args:
            crude_baseline_end (int): self.find_baseline_end 의 return value
            ignore_init_cycles (int, optional): qPCR 의 처음 cycle 들은 버린다. Defaults to 3.

        Returns:
            int: 좀 더 정밀하게 찾은 baseline end (fine end)
        """
        res = []
        cq_values = self.get_cq_values()
        cycles = cq_values["cycle"]
        cq_list = cq_values["Cq"]
        offset_list = [-2, -1, 0, 1, 2]
        for offset in offset_list:
            Y = cq_list[ignore_init_cycles : crude_baseline_end + offset + 1]
            X = cycles[ignore_init_cycles : crude_baseline_end + offset + 1]
            slope, intercept, pvalue = self.lm(X, Y)
            res.append([crude_baseline_end + offset, pvalue])
        sorted_res = sorted(res, key=lambda x: x[1])
        return sorted_res[0][0]

    def lm(self, X: list, Y: list) -> float:
        X = sm.add_constant(X)
        model = sm.OLS(Y, X)
        results = model.fit()
        pvalue = results.pvalues[1]
        intercept, slope = results.params
        return slope, intercept, pvalue

    def baseline_lm(
        self, baseline_end: int, ignore_init_cycles: int = 3, input_cq: dict = None
    ) -> float:
        """linear regression of Cq values from ignore_init_cycles to baseline_end

        Args:
            baseline_end (int): 1-based cycle number. 첫번째 cycle 이 0에서 시작해도 첫번째라는 뜻으로 1부터 시작
            ignore_init_cycles (int): 처음 몇 cycle 을 무시할 것인지. Defaults to 3.
            input_cq (dict): if None, use the self.get_cq_values() output

        Returns:
            float: _description_
        """
        if input_cq is None:
            cq_values = self.get_cq_values()
        else:
            cq_values = input_cq
        cycles = cq_values["cycle"]
        cq_list = cq_values["Cq"]
        baseline_offset = cycles[0] - 1
        try:
            contain_idx = cycles.index(baseline_end + baseline_offset)
        except ValueError as e:
            print(f"cycle list: {cycles}, and looking for {baseline_end}")
            raise ValueError(e)
        Y = cq_list[ignore_init_cycles : contain_idx + 1]
        X = cycles[ignore_init_cycles : contain_idx + 1]
        try:
            slope, intercept, pvalue = self.lm(X, Y)
        except IndexError as e:
            print(f"X:{X}, Y:{Y}")
            raise IndexError(e)
        return slope

    def slope_to_angle(self, slope: float, deg_type: str = "radian") -> float:
        """_summary_

        Args:
            slope (float): _description_

        Returns:
            float: _description_
        """
        res = math.atan(slope)
        if deg_type == "degree":
            res = math.degrees(res)
        return res

    def rotate_cqs(
        self,
        baseline_end: int,
        input_cq: dict = None,
        ignore_init_cycles: int = 3,
        debugq: bool = False,
    ) -> dict:
        """from input cq, calculate the rotated coordinates. baseline point is the rotating center.
        cycle 은 rotate 된것을 사용하지 않고
        Cq 값만 rotate 한 것을 사용한다. 아래 Examples 확인=> 이렇게 하는 것에 사람들이 이견이 없음

        Args:
            baseline_end (int): 1-based cycle number. 첫번째 cycle 이 0에서 시작해도 첫번째라는 뜻으로 1부터 시작
            input_cq (dict): if None, use the self.get_cq_values() output
            ignore_init_cycles (int): 처음 몇 cycle 을 무시할 것인지. Defaults to 3.
            debugq (bool): print out the debug message

        Returns:
            dict:

        Examples:
            예를 들어 -45도 회전을 한다. 그러면 baseline end 이후의 점들은 y축에 평행하게 계산이 되는데
            그 값들을 그대로 이용하면서 cycle은 기존 값을 사용하여 표시한다.
            즉 기울기가 -1 인 그래프를 -45도 회전이동하면 baseline end 이후 값들은 baseline 위로 그대로 올라가게 된다.
            :+     +                    +                                 +
            : +   +    = 회전이동 =>     +  = cycle 을 그대로 이용 =>      +
            :  + +                      +                               +
            :   +                    ++++                           ++++
            :1234567                 1234                           1234567
        """
        res = {"cycle": [], "Cq": [], "Cq_smoothing": []}
        if input_cq is None:
            input_cq = self.get_cq_values()
        baseline_end_offset = input_cq["cycle"][0] - 1
        if debugq:
            print(f"rotate at {baseline_end+baseline_end_offset}")

        if baseline_end - ignore_init_cycles < 2:
            # linear regression 이 가능하려면 baseline_end가 ignore_init_cycle 보다 2이상 커야 함
            slope = self.baseline_lm(
                baseline_end,
                ignore_init_cycles=ignore_init_cycles - 1,
                input_cq=input_cq,
            )
        else:
            slope = self.baseline_lm(
                baseline_end, ignore_init_cycles=ignore_init_cycles, input_cq=input_cq
            )
        # | if
        rad = self.slope_to_angle(-slope, deg_type="radian")
        x_pivot = baseline_end + baseline_end_offset
        y_pivot = self.get_baseline_end_cq(baseline_end, input_cq=input_cq)
        X = np.array(input_cq["cycle"])
        Y = np.array(input_cq["Cq"])
        X, Y = self.rotate_point(*(self.shift_vector(X, Y, x_pivot, y_pivot)), rad)
        res["cycle"] = input_cq["cycle"]
        res["Cq"] = np.round((Y + y_pivot), 4).tolist()
        scale_factor = input_cq["Cq_smoothing"][0]
        res["Cq_smoothing"] = (np.array(res["Cq"]) * scale_factor).tolist()
        return res

    def shift_vector(
        self, X: np.array, Y: np.array, pivot_x: int, pivot_y: float
    ) -> Tuple[np.array, np.array]:
        return X - pivot_x, Y - pivot_y

    def rotate_point(self, x, y, radian):
        angle_radians = radian

        # 회전 행렬 공식 적용
        x_new = x * math.cos(angle_radians) - y * math.sin(angle_radians)
        y_new = x * math.sin(angle_radians) + y * math.cos(angle_radians)

        return x_new, y_new

    def get_baseline_end_cq(self, baseline_end: int, input_cq: dict = None):
        """_summary_

        Args:
            baseline_end (int): 1-based cycle number. 첫번째 cycle 이 0에서 시작해도 첫번째라는 뜻으로 1부터 시작
            input_cq (dict, optional): if None, use the self.get_cq_values() output
        """
        if input_cq is None:
            input_cq = self.get_cq_values()
        baseline_end_offset = input_cq["cycle"][0] - 1
        baseline_cq_idx = input_cq["cycle"].index(baseline_end + baseline_end_offset)
        res_cq = input_cq["Cq"][baseline_cq_idx]
        return res_cq
