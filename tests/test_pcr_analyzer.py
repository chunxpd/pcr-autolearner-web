import math
import numpy as np
import pytest
from src.pcr_analyzer import PCRmachine


# mock PCRmachine fixture
@pytest.fixture(scope="function")
def mock_pcr_machine(mocker):
    """
    cycle	Cq	dCq	ddCq	dddCq
    1	10
    2	11	1
    3	9	-2	-3
    4	13	4	6	9
    5	10	-3	-7	-13
    6	14	4	7	14
    7	11	-3	-7	-14
    8	15	4	7	14
    9	12	-3	-7	-14
    10	16	4	7	14
    11	13	-3	-7	-14
    12	20	7	10	17  <-- 여기가 baseline 임
    13	30	10	3	-7
    14	45	15	5	2
    15	61	16	1	-4
    16	100	39	23	22
    17	200	100	61	38
    18	400	200	100	39
    19	480	80	-120	-220
    20	500	20	-60	60
    21	501	1	-19	41
    """
    # PCRmachine 클래스의 인스턴스 mock 생성

    mock_machine = mocker.Mock(spec=PCRmachine)

    # cycle별로 반환할 Cq 값을 리스트로 정의
    data = {
        "cycle": [
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            20,
            21,
        ],
        "Cq": [
            10,
            11,
            9,
            13,
            10,
            14,
            11,
            15,
            12,
            16,
            13,
            20,
            30,
            45,
            61,
            100,
            200,
            400,
            480,
            500,
            501,
        ],
        "dCq": [
            None,
            1,
            -2,
            4,
            -3,
            4,
            -3,
            4,
            -3,
            4,
            -3,
            7,
            10,
            15,
            16,
            39,
            100,
            200,
            80,
            20,
            1,
        ],
        "ddCq": [
            None,
            None,
            -3,
            6,
            -7,
            7,
            -7,
            7,
            -7,
            7,
            -7,
            10,
            3,
            5,
            1,
            23,
            61,
            100,
            -120,
            -60,
            -19,
        ],
        "dddCq": [
            None,
            None,
            None,
            9,
            -13,
            14,
            -14,
            14,
            -14,
            14,
            -14,
            17,
            -7,
            2,
            -4,
            22,
            38,
            39,
            -220,
            60,
            41,
        ],
    }

    # side_effect로 값을 차례로 반환하고, 다 반환하면 StopIteration 발생
    mock_machine.get_cq_values.return_value = data

    return mock_machine


@pytest.fixture(scope="function")
def empty_obj():
    res = PCRmachine()
    yield res
    del res


@pytest.fixture(scope="module")
def obj2():
    res = PCRmachine()
    for cqval in [
        43704,
        43837,
        43476,
        43775,
        43758,
        43604,
        43241,
        43387,
        43123,
        43286,
    ]:
        res.add_one_cq_value(cqval)
    yield res
    del res


# dictionary 에 nan 이 있을때를 대비해서 만든 assert 용 함수
def compare_dicts(dict1, dict2):
    """nan 이 있을 때 두개의 dict 는 무조건 다르다고 나온다. NaN==NaN => False 임"""
    # 두 dict의 키가 동일한지 확인
    assert dict1.keys() == dict2.keys(), "Keys in dictionaries do not match"

    for key in dict1:
        val1, val2 = dict1[key], dict2[key]

        # 값이 리스트인 경우, 리스트 내부 값 비교
        if isinstance(val1, list) and isinstance(val2, list):
            assert len(val1) == len(val2), f"List lengths for key '{key}' do not match"
            for v1, v2 in zip(val1, val2):
                compare_values(v1, v2)  # 리스트 내부 값 비교
        else:
            compare_values(val1, val2)  # 단일 값 비교


def compare_values(val1, val2):
    # 두 값이 모두 NaN인 경우는 같다고 처리
    if isinstance(val1, float) and isinstance(val2, float):
        if np.isnan(val1) and np.isnan(val2):
            return

        # float 값이 1/1000 수준에서 같으면 assert 통과
        assert (
            abs(val1 - val2) < 1e-3
        ), f"Values do not match within tolerance: {val1} != {val2}"
    else:
        # 값이 같지 않을 경우 assertion error 발생
        assert val1 == val2, f"Values do not match: {val1} != {val2}"


# test function
def test_smoothing(empty_obj):
    obs = [1, 2, 3]
    pred = 2
    assert empty_obj.moving_average(obs, winsize=3) == pred
    obs = [1, 2, 3, 4]
    pred = 3
    assert empty_obj.moving_average(obs, winsize=3) == pred
    obs = [1, 2]
    pred = 1.5
    assert empty_obj.moving_average(obs, winsize=3) == pred
    obs = [1]
    pred = 1
    assert empty_obj.moving_average(obs, winsize=3) == pred
    obs = [1, 2, 3, 4, 5]
    pred = 3
    assert empty_obj.moving_average(obs, winsize=5) == pred


def test_pcr_machine_add_cq(empty_obj):
    """Cq는 Cq_smoothing 과 같은 값으로 Cq_smoothing 의 첫번째 값으로 나누어준 값이다.
    rotation 계산 등에 쓰이는 것이다."""
    # mock_pcr_machine에서 cycle별로 Cq 값을 가져와 테스트
    pred = {
        "cycle": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "Cq_ori": [
            43704,
            43837,
            43476,
            43775,
            43758,
            43604,
            43241,
            43387,
            43123,
            43286,
        ],
        "Cq_smoothing": [
            43704.0,
            (43837 + 43704) / 2,
            (43476 + 43837 + 43704) / 3,
            (43775 + 43476 + 43837) / 3,
            (43758 + 43775 + 43476) / 3,
            (43604 + 43758 + 43775) / 3,
            (43241 + 43604 + 43758) / 3,
            (43387 + 43241 + 43604) / 3,
            (43123 + 43387 + 43241) / 3,
            (43286 + 43123 + 43387) / 3,
        ],
        "Cq": [
            1.0,
            1.0015215998535603,
            0.9992754286411618,
            0.9998169503935567,
            0.9992144121056806,
            1.0001906766733786,
            0.996117822930014,
            0.9932881810970773,
            0.9896195619012753,
            0.9899627799133566,
        ],
        "dCq": [None] * 10,
        "ddCq": [None] * 10,
        "dddCq": [None] * 10,
    }
    obs_values = [
        43704,
        43837,
        43476,
        43775,
        43758,
        43604,
        43241,
        43387,
        43123,
        43286,
    ]
    for obs in obs_values:
        empty_obj.add_one_cq_value(obs)
    assert pred == empty_obj.get_cq_values()


def test_pcr_machine_realtime_dCq(obj2):
    """
    cycle	Cq	dCq	ddCq	dddCq
    1	43704
    2	43837	133
    3	43476	-361	-494
    4	43775	299	660	1154
    5	43758	-17	-316	-976
    6	43604	-154	-137	179
    7	43241	-363	-209	-72
    8	43387	146	509	718
    9	43123	-264	-410	-919
    10	43286	163	427	837
    위와 같이 하였으나 이제는 smoothing 을 기반으로 계산한다.
    """
    pred = {
        "cycle": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "Cq_ori": [
            43704,
            43837,
            43476,
            43775,
            43758,
            43604,
            43241,
            43387,
            43123,
            43286,
        ],
        "Cq_smoothing": [
            43704.0,
            (43837 + 43704) / 2,
            (43476 + 43837 + 43704) / 3,
            (43775 + 43476 + 43837) / 3,
            (43758 + 43775 + 43476) / 3,
            (43604 + 43758 + 43775) / 3,
            (43241 + 43604 + 43758) / 3,
            (43387 + 43241 + 43604) / 3,
            (43123 + 43387 + 43241) / 3,
            (43286 + 43123 + 43387) / 3,
        ],
        "Cq": [
            1.0,
            1.0015215998535603,
            0.9992754286411618,
            0.9998169503935567,
            0.9992144121056806,
            1.0001906766733786,
            0.996117822930014,
            0.9932881810970773,
            0.9896195619012753,
            0.9899627799133566,
        ],
        "dCq": [
            np.nan,
            66.500000,
            -98.166667,
            23.666667,
            -26.333333,
            42.666667,
            -178.000000,
            -123.666667,
            -160.333333,
            15.000000,
        ],
        "ddCq": [
            np.nan,
            np.nan,
            -164.666667,
            121.833333,
            -50.000000,
            69.000000,
            -220.666667,
            54.333333,
            -36.666667,
            175.333333,
        ],
        "dddCq": [
            np.nan,
            np.nan,
            np.nan,
            286.500000,
            -171.833333,
            119.000000,
            -289.666667,
            275.000000,
            -91.000000,
            212.000000,
        ],
    }
    compare_dicts(pred, obj2.pcr_machine_cal_dCq())


def test_pcr_machine_realtime_dCq_exception(empty_obj):
    obs = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    for cq in obs:
        empty_obj.add_one_cq_value(cq)
    pred = {
        "cycle": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "Cq_ori": [0] * 10,
        "Cq_smoothing": [0.0] * 10,
        "Cq": [0.0] * 10,
        "dCq": [np.nan] * 1 + [0.0] * 9,
        "ddCq": [np.nan] * 2 + [0.0] * 8,
        "dddCq": [np.nan] * 3 + [0.0] * 7,
    }
    compare_dicts(pred, empty_obj.pcr_machine_cal_dCq())
    # with pytest.raises(ValueError):
    #     empty_obj.pcr_machine_cal_dCq()


def test_ddcq_threshold(empty_obj):
    """지금까지 들어온 ddCq 값의 min-max 차이의 0.1% 를 넘는 값 T의 ddCq값을 찾기 위해서"""
    obs = [8510]
    pred = 0
    assert empty_obj.cal_ddcq_threshold(obs) == pred
    obs = [8510, 8557]
    pred = (8557 - 8510) * 0.001
    assert empty_obj.cal_ddcq_threshold(obs) == pred
    obs = [8510, 8557, 8355]
    pred = 0
    assert empty_obj.cal_ddcq_threshold(obs, scale_factor=0.0) == pred


def test_ddcq_threshold2(empty_obj):
    """max_cycle개념 도입. 처음을 1cycle 이라고 하고 max_cycle 까지 값 중 최대, 초소 값에 scale factor 곱하기"""
    obs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    pred = (4 - 1) * 1.0
    assert empty_obj.cal_ddcq_threshold(obs, scale_factor=1.0, max_cycle=4) == pred
    pred = (10 - 1) * 1.0
    assert empty_obj.cal_ddcq_threshold(obs, scale_factor=1.0) == pred


def test_add_one_cq_value(empty_obj):
    """Cq 값이 들어가면 cycle 은 1씩 증가하고, Cq 값은 Cq_ori 에, smoothing 값은 Cq_smoothing 에 들어간다.
    start 값이 주어지면 cycle 은 start 값으로 시작한다."""
    obs = [1, 2, 3, 4, 5]
    obs_cycles = [0, 1, 2, 3, 4]
    for cq in obs:
        empty_obj.add_one_cq_value(cq, start=0)
    assert empty_obj.get_cq_values()["Cq_ori"] == obs
    assert empty_obj.get_cq_values()["Cq_smoothing"] == [
        1.0,
        (2 + 1) / 2,
        (3 + 2 + 1) / 3,
        (4 + 3 + 2) / 3,
        (5 + 4 + 3) / 3,
    ]
    assert empty_obj.get_cq_values()["cycle"] == obs_cycles


def test_find_baseline_end(empty_obj):
    """전체 데이터가 확보된 상태에서 찾는 것이 아니라 현재 데이터가 들어간 상태에서 찾는 것이다"""
    obs_cq = [
        10,
        11,
        9,
        13,
        10,
        14,
        11,
        15,
        12,
        16,
        13,
        20,
        30,
        45,
        61,  # <- 이 예에서는 여기가 ddCq 가 4연속 나온 지점이고 이때 cycle 은 15
        100,
        200,
        400,
        480,
        500,
        501,
    ]
    pred = (
        10  # 4번 연속 양의 값을 가지는 ddCq 가 나왔을 때의 최초의 cycle 값 15 - 4 - 1
    )
    for cq_value in obs_cq:
        empty_obj.add_one_cq_value(cq_value)
        empty_obj.pcr_machine_cal_dCq()
        obs = empty_obj.find_baseline_end(consecutive_number=4)
        if not obs:
            continue
        break
    assert pred == obs


def test_find_baseline_end2(empty_obj):
    """negative 의 경우. ddCq 가 4연속 + 인 경우가 없음
    이럴 경우는 max cycle 인 40 에서 4연속의 시작점인 36을 baseline end로 반환"""
    obs_cq = [
        3624,
        3589,
        3600,
        3605,
        3617,
        3632,
        3619,
        3631,
        3644,
        3660,
        3659,
        3669,
        3687,
        3689,
        3694,
        3704,
        3709,
        3707,
        3730,
        3731,
        3740,
        3743,
        3758,
        3763,
        3762,
        3782,
        3783,
        3788,
        3803,
        3799,
        3814,
        3837,
        3845,
        3844,
        3867,
        3890,
        3895,
        3914,
        3920,
        3938,
    ]
    pred = 36
    for cq_value in obs_cq:
        empty_obj.add_one_cq_value(cq_value)
        empty_obj.pcr_machine_cal_dCq()
        obs = empty_obj.find_baseline_end(consecutive_number=4)
        if not obs:
            continue
        break
    assert pred == obs


def test_find_baseline_end3(empty_obj):
    """총 15cycle 을 최대 cycle 수로 세팅 했을 때
    baseline 을 13에서 찾기 때문에 너무 끝부분에서 baseline 을 찾았다고 생각하고 최대값인 15-4인 11을 반환"""
    obs_cq = [
        10,  # 1
        11,  # 2
        9,  # 3
        13,  # 4
        10,  # 5
        14,  # 6
        11,  # 7
        15,  # 8
        12,  # 9 <- 여기가 baseline end
        16,  # 10
        13,  # 11
        20,  # 12
        30,  # 13
        45,  # 14
        61,  # 15 <- 이 예에서는 여기가 ddCq 가 4연속 나온 지점이고 이때 cycle 은 15. 그리고 여기가 마지막 cycle 임
    ]
    pred = 11  # 15 cycle 까지만 사용자가 running 하였다고 하면 4연속이 마지막 cycle 에서 나온 것이어서 15-4
    empty_obj.max_cycle = 15
    for cq_value in obs_cq:
        empty_obj.add_one_cq_value(cq_value)
        empty_obj.pcr_machine_cal_dCq()
        obs = empty_obj.find_baseline_end(consecutive_number=4)
        if not obs:
            continue
        break
    assert pred == obs


def test_find_baseline_end4(empty_obj):
    """처음 값 3개는 무시하고 있는데 baseline_end 를 그 부분에서 찾았을 때에는 무시한다."""
    obs_cq = [
        11075,
        10742,
        11108,  # <- 원래는 여기가 baseline_end 그러나 첫 3개는 무시하므로 건너뜀
        11077,
        10878,
        10832,
        10793,
        10783,
        10788,
        10551,
        10564,
        10532,
        10569,
        10477,
        10371,
        10450,
        10433,
        10267,
        10491,
        10369,
        10327,
        10280,
        10336,
        10369,
        10370,
        10454,
        10531,
        10780,
        11133,
        11716,
        12443,
        13248,
        14058,
        14796,
        15573,
        16349,
        17055,
        17760,
        18467,
        19080,
        19722,
        20316,
        20897,
        21436,
        21925,
    ]
    pred = 25
    empty_obj.max_cycle = len(obs_cq)
    for cq_value in obs_cq:
        empty_obj.add_one_cq_value(cq_value)
        empty_obj.pcr_machine_cal_dCq()
        obs = empty_obj.find_baseline_end(consecutive_number=4, ignore_init_cycles=3)
        if not obs:
            continue
        break
    assert pred == obs


def test_baseline_end_dynamic_thr(empty_obj):
    obs_cq = [
        8510,
        8557,
        8355,
        8309,
        8269,
        8229,
        8180,
        8165,
        8152,
        8141,
        8133,
        8128,
        8117,
        8118,
        8105,
        8104,
        8110,
        8109,
        8112,
        8110,
        8133,
        8166,
        8262,
        8404,
        8656,
        9140,
        10073,
        11632,
        13595,
        15244,
        16332,
        17097,
        17681,
        18191,
        18637,
        19015,
        19351,
        19641,
        19904,
        20134,
        20327,
        20517,
        20685,
        20814,
        20936,
    ]
    pred = 19
    for cq_value in obs_cq:
        empty_obj.add_one_cq_value(cq_value)
        empty_obj.pcr_machine_cal_dCq()
        obs = empty_obj.find_baseline_end(ddcq_threshold_scale=0.01)
        if not obs:
            continue
        break
    assert pred == obs


def test_find_baseline_end_6consecutive(empty_obj):
    obs_cq = [
        9655,
        9719,
        9670,
        9708,
        9668,
        9630,
        9665,
        9667,
        9555,
        9590,
        9515,
        9468,
        9519,
        9468,
        9437,
        9444,
        9408,
        9393,
        9370,
        9404,
        9372,
        9347,
        9310,
        9387,  # <- 24 cycle. 여기가 dynamic T보다 큰 지역 시작
        9382,
        9408,
        9523,
        9553,
        9829,
        10079,  # <- 30 cycle 여기까지 총 7연속 dynamic T 보다 큰 지역
        10518,
        10988,
        11504,
        12030,
        12433,
        12994,
        13485,
        13930,
        14387,
        14847,
        15253,
        15701,
        16099,
        16471,
    ]
    pred = 22
    # 일단 다 넣고
    for cq_value in obs_cq[:29]:  # cycle은 1부터 시작. 29가 30cycle
        empty_obj.add_one_cq_value(cq_value)
        empty_obj.pcr_machine_cal_dCq()
    # n(ddCq > T) 가 4보다 큰 연속된 영역을 찾는다. 이 경우는 총 5개 연속 영역
    obs = empty_obj.find_baseline_end(consecutive_number=4, ddcq_threshold_scale=0.01)
    assert pred == obs


def test_find_consecutive_regions(empty_obj):
    """return 0-based index"""
    obs = [
        math.nan,  # index 0
        math.nan,  # index 1
        -37.66666666666606,  # index 2
        23.33333333333212,
        -34.66666666666606,
        3.6666666666660603,  # index 5
        -0.999999999998181,
        13.999999999998181,
        -24.66666666666606,
        0.0,
        -25.66666666666606,  # index 10
        21.66666666666606,
        5.333333333332121,
        8.000000000001819,
        5.333333333332121,
        -14.66666666666606,  # index 15
        5.0,
        5.33333333333394,
        -10.0,
        23.33333333333212,
        -5.66666666666606,  # index 20
        -0.6666666666660603,
        -23.66666666666788,
        36.33333333333394,  # index 23
        6.66666666666606,
        21.00000000000182,  # index 25
        12.666666666664241,
        11.66666666666788,
        83.33333333333394,
        45.0,
        136.33333333333212,  # index 30
        64.66666666666788,
        88.66666666666606,
        29.0,  # index 33
        -22.33333333333394,
    ]
    pred = [(23, 33)]
    assert empty_obj.find_consecutive_regions(obs, consecutive_number=11, T=0) == pred


def test_find_all_negative_dcqs(empty_obj):
    """입력 list 를 slicing 했을 때 모두 < 0 이면 True 발생"""
    obs_dcq = [np.nan, -1, -2, -3, 4, -5, -6, -7, -8]
    pred = True
    obs = empty_obj.find_all_negative_dcqs(obs_dcq, 1, 3)
    assert obs == pred
    pred = False
    obs = empty_obj.find_all_negative_dcqs(obs_dcq, 1, 4)
    assert obs == pred
    pred = False
    obs = empty_obj.find_all_negative_dcqs(obs_dcq, 0, 3)
    assert obs == pred


def test_find_all_positivie_dcqs(empty_obj):
    """입력 list 를 slicing 했을 때 모두 >=0 이면 True 발생"""
    obs_dcq = [np.nan, 1, 2, 3, -4, 5, 6, 7, 8]
    pred = True
    obs = empty_obj.find_all_positive_dcqs(obs_dcq, 1, 3)
    assert obs == pred
    pred = False
    obs = empty_obj.find_all_positive_dcqs(obs_dcq, 1, 4)
    assert obs == pred
    pred = False
    obs = empty_obj.find_all_positive_dcqs(obs_dcq, 0, 3)
    assert obs == pred


def test_find_baseline_end_prev_consecutive(empty_obj):
    """현재 입력 cycle 에서는 ddCq 가 음수여도 이전에 n(ddCq>T)>=4 인 지역이 있었으면 그 지역을 근거로 baseline_end 를 설정하기"""
    obs_cq = [
        9655,
        9719,
        9670,
        9708,
        9668,
        9630,
        9665,
        9667,
        9555,
        9590,
        9515,
        9468,
        9519,
        9468,
        9437,
        9444,
        9408,
        9393,
        9370,
        9404,
        9372,
        9347,
        9310,
        9387,  # <- 24 cycle. 여기가 dynamic T보다 큰 지역 시작. scale=0.1% (0.001)
        9382,
        9408,
        9523,
        9553,
        9829,
        10079,  # <- 30 cycle 여기까지 총 7연속 dynamic T 보다 큰 지역
        10518,
        10988,
        11504,
        12030,
        12433,  # <- 35 cycle 여기에서 ddCq 가 -22 가 됨.
        12994,
        13485,
        13930,
        14387,
        14847,
        15253,
        15701,
        16099,
        16471,
    ]
    pred = 22
    # 일단 다 넣고
    for cq_value in obs_cq[:35]:
        empty_obj.add_one_cq_value(cq_value)
        empty_obj.pcr_machine_cal_dCq()
    obs = empty_obj.find_baseline_end(consecutive_number=4, ddcq_threshold_scale=0.001)
    assert pred == obs


def test_find_baseline_end_prev_consecutive2(empty_obj):
    """현재 입력 cycle 에서는 ddCq 가 음수여도 이전에 n(ddCq>T)>=4 인 지역이 있었으면 그 지역을 근거로 baseline_end 를 설정하기"""
    obs_cq = [
        9655,
        9719,
        9670,
        9708,
        9668,
        9630,
        9665,
        9667,
        9555,
        9590,
        9515,
        9468,
        9519,
        9468,
        9437,
        9444,
        9408,
        9393,
        9370,
        9404,
        9372,
        9347,
        9310,
        9387,  # <- 24 cycle. 여기가 dynamic T보다 큰 지역 시작. scale=0.1% (0.001)
        9382,
        9408,
        9523,
        9553,
        9829,
        10079,  # <- 30 cycle 여기까지 총 7연속 dynamic T 보다 큰 지역
        10518,
        10988,
        11504,
        12030,
        12433,  # <- 35 cycle 여기에서 ddCq 가 -22 가 됨.
        12994,
        13485,
        13930,
        14387,
        14847,
        15253,
        15701,
        16099,
        16471,
    ]
    pred = 22
    # 일단 다 넣고
    for cq_value in obs_cq:
        empty_obj.add_one_cq_value(cq_value)
        empty_obj.pcr_machine_cal_dCq()
    obs = empty_obj.find_baseline_end(
        consecutive_number=4, ddcq_threshold_scale=0.001, max_cycle=30
    )
    assert pred == obs


def test_find_baseline_exception(empty_obj):
    """issue 252: n(ddCt > T) >= 4이어도 dCt 가 음수일 경우는 무시"""
    obs_cq = [
        1404,
        1401,
        1360,
        1305,  # <-4cycle. ddCq 가 T보다 커지는 시점
        1264,
        1235,
        1232,  # <- 7cycle. ddCq 4연속 지점
        1228,
        1231,
    ]
    pred = False
    for cq_value in obs_cq:
        empty_obj.add_one_cq_value(cq_value)
        empty_obj.pcr_machine_cal_dCq()
    obs = empty_obj.find_baseline_end(consecutive_number=4, ddcq_threshold_scale=0.001)
    assert pred == obs


def test_find_baseline_end_start0(empty_obj):
    """test_find_baseline_end_prev_consecutive2 와 같은 테스트이나 start=0으로 했을 때
    그래도 똑같이 22가 나와야 한다. 첫번째 cycle부터 n번째가 baseline end라는 뜻이기 때문"""
    obs_cq = [
        9655,
        9719,
        9670,
        9708,
        9668,
        9630,
        9665,
        9667,
        9555,
        9590,
        9515,
        9468,
        9519,
        9468,
        9437,
        9444,
        9408,
        9393,
        9370,
        9404,
        9372,
        9347,
        9310,
        9387,  # <- 24 cycle. 여기가 dynamic T보다 큰 지역 시작. scale=0.1% (0.001)
        9382,
        9408,
        9523,
        9553,
        9829,
        10079,  # <- 30 cycle 여기까지 총 7연속 dynamic T 보다 큰 지역
        10518,
        10988,
        11504,
        12030,
        12433,  # <- 35 cycle 여기에서 ddCq 가 -22 가 됨.
        12994,
        13485,
        13930,
        14387,
        14847,
        15253,
        15701,
        16099,
        16471,
    ]
    pred = 22
    # 일단 다 넣고
    for cq_value in obs_cq:
        empty_obj.add_one_cq_value(cq_value, start=0)
        empty_obj.pcr_machine_cal_dCq()
    obs = empty_obj.find_baseline_end(
        consecutive_number=4, ddcq_threshold_scale=0.001, max_cycle=30
    )
    assert pred == obs


def test_baseline_lm(empty_obj):
    obs_cq = [
        10,
        11,
        9,
        13,
        10,
        14,
        11,
        15,
        12,  # <- 여기가 baseline end. cycle 은 9. list index 는 8
        16,
        13,
        20,
        30,
        45,
        61,  # <- 이 예에서는 여기가 ddCq 가 4연속 나온 지점이고 이때 cycle 은 15
        100,
        200,
        400,
        480,
        500,
        501,
    ]
    pred = 0.045
    for cq_value in obs_cq:
        ### PCR 장치를 모사하는 부분. 한개 사이클이 지날 때마다 Cq 값을 입력한다.
        ### 그러다가 ddCq 가 발견되면 baseline end 를 반환받고 루프를 끝낸다.
        empty_obj.add_one_cq_value(cq_value)
        empty_obj.pcr_machine_cal_dCq()
        baseline_end = empty_obj.find_baseline_end(consecutive_number=4)
        if not baseline_end:
            continue
        break
    obs = empty_obj.baseline_lm(baseline_end)
    assert abs(pred - obs) < 0.01


def test_slope_to_angle(obj2):
    pred = 45
    assert pred == obj2.slope_to_angle(1, "degree")
    pred = 0.7854
    assert abs(pred - obj2.slope_to_angle(1, "radian")) < 0.1
    pred = 63.43
    assert abs(pred - obj2.slope_to_angle(2, "degree")) < 0.1


def test_linear_regression(empty_obj):
    obs_x = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    obs_y = [1, 2, 3, 4, 5, 6, 7, 8, 9]  # 완전 일직선
    pred = 0.05
    obs_slope, obs_intercept, obs_pval = empty_obj.lm(obs_x, obs_y)
    assert obs_slope == 1.0
    assert obs_pval < pred
    obs_x = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    obs_y = [0.9, 2.1, 3.5, 4.1, 4.9, 6.1, 7.2, 7.8, 9]  # 거의 일직선
    pred = 0.05
    obs_slope, obs_intercept, obs_pval = empty_obj.lm(obs_x, obs_y)
    assert obs_pval < pred
    obs_x = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    obs_y = [1, 10, 2, -19, 5, 100, 0, 80, 9]  # 중구난방
    pred = 0.05
    obs_slope, obs_intercept, obs_pval = empty_obj.lm(obs_x, obs_y)
    assert obs_pval > pred


@pytest.mark.skip(reason="아직 fine tuning 알고리즘이 확정되지 않음")
def test_fining_baseline_end(empty_obj):
    """일차적으로 찾은 baseline end 를 가지고 좀 더 정확하게 찾아본다. 현재 제안된 알고리즘
    1. linear regression 을 한다.
    2. baseline +-1 에서 linear regression 을 한다
    3. 가장 작은 p-value 를 가지는 baseline 을 최종 결정한다.
    아직 잘 작동하는 것인지 확인 되지 않았음. 일단 skip 할 것
    """
    obs_cq = [
        10,
        11,
        9,
        13,
        10,
        14,
        11,
        15,
        12,
        16,
        13,  # <- 여기가 실질적인 baseline end. 11 cycle
        20,
        30,  # <- 여기가 1차적으로 찾은 baseline end. 13 cycle
        45,
        61,  # <- 이 예에서는 여기가 ddCq 가 4연속 나온 지점이고 이때 cycle 은 15
        100,
        200,
        400,
        480,
        500,
        501,
    ]
    pred = 11
    for cq_value in obs_cq:
        empty_obj.add_one_cq_value(cq_value)
        empty_obj.pcr_machine_cal_dCq()
        crude_baseline_end = empty_obj.find_baseline_end(consecutive_number=4)
        if not crude_baseline_end:
            continue
        break
    fine_baseline_end = empty_obj.fining_baseline_end(crude_baseline_end)
    assert pred == fine_baseline_end


def test_get_baseline_end_cq(empty_obj):
    obs = {
        "cycle": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "Cq": [5, 4, 3, 2, 1, 2, 3, 4, 5, 6],
    }
    obs_cycle = 5
    pred = 1
    assert obs["Cq"][4] == empty_obj.get_baseline_end_cq(obs_cycle, input_cq=obs)
    obs_cycle = 1
    pred = 5
    assert obs["Cq"][0] == empty_obj.get_baseline_end_cq(obs_cycle, input_cq=obs)


def test_rotate_cq(empty_obj):
    """cycle 5일때를 baseline 으로 잡았을 때 좌표 회전
    -45도 회전을 한다. 그러면 baseline end 이후의 점들은 y축에 평행하게 계산이 되는데
    그 값들을 그대로 이용하면서 cycle은 기존 값을 사용하여 표시한다.
    즉 기울기가 -1 인 그래프를 -45도 회전이동하면 baseline end 이후 값들은 baseline 위로 그대로 올라가게 된다.
    +     +                    +                                 +
     +   +    = 회전이동 =>     +  = cycle 을 그대로 이용 =>      +
      + +                      +                               +
       +                    ++++                           ++++
    1234567                 1234                           1234567
    Args:
        empty_obj (_type_): _description_
    """
    obs = {
        "cycle": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "Cq": [5, 4, 3, 2, 1, 2, 3, 4, 5, 6],
        "Cq_smoothing": [500, 400, 300, 200, 100, 200, 300, 400, 500, 600],
    }
    obs_baseline_end = 5
    pred = {
        "cycle": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "Cq": [
            1.00000,
            1.00000,
            1.00000,
            1.0,
            1.0,
            2.4142,
            3.8284,
            5.2426,
            6.6569,
            8.0711,
        ],
        "Cq_smoothing": [
            1.00000 * 500,
            1.00000 * 500,
            1.00000 * 500,
            1.0 * 500,
            1.0 * 500,
            2.4142 * 500,
            3.8284 * 500,
            5.2426 * 500,
            6.6569 * 500,
            8.0711 * 500,
        ],
    }
    assert (
        empty_obj.rotate_cqs(obs_baseline_end, input_cq=obs, ignore_init_cycles=0)
        == pred
    )


def test_rotate_cq_start0(empty_obj):
    """입력 cycle 이 0에서부터 시작할 때 rotate 된 것도 cycle 이 0에서 시작해야 함
    baseline end 는 1-based 로 들어가야 함. 위와 동일한 경우이기 때문에 결과도 동일해야 함"""
    obs = {
        "cycle": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        "Cq": [5, 4, 3, 2, 1, 2, 3, 4, 5, 6],
        "Cq_smoothing": [500, 400, 300, 200, 100, 200, 300, 400, 500, 600],
    }
    obs_baseline_end = 5  # 5번째 cycle 이 baseline end
    pred = {
        "cycle": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        "Cq": [
            1.00000,
            1.00000,
            1.00000,
            1.0,
            1.0,
            2.4142,
            3.8284,
            5.2426,
            6.6569,
            8.0711,
        ],
        "Cq_smoothing": [
            1.00000 * 500,
            1.00000 * 500,
            1.00000 * 500,
            1.0 * 500,
            1.0 * 500,
            2.4142 * 500,
            3.8284 * 500,
            5.2426 * 500,
            6.6569 * 500,
            8.0711 * 500,
        ],
    }
    obs_res = empty_obj.rotate_cqs(obs_baseline_end, input_cq=obs, ignore_init_cycles=0)
    assert obs_res == pred


def test_rotate_cqs_error(empty_obj):
    """baseline_end 가 3 (ignore_init_cycle)보다 1 큰 4일 경우에는 linear regression 에 딱 1개의 vector 만 들어가기 때문에 에러 발생
    그래서 임시적으로 ignore_init_cycle 을 1 줄여서 통과 시킴"""
    obs = {
        "cycle": [
            1,
            2,
            3,
            4,  # <- 여기가 baseline_end 임
            5,
            6,
            7,
            8,
            9,
            10,
        ],
        "Cq": [
            20486 / 20486,
            20701 / 20486,
            20491.66667 / 20486,
            20271.66667 / 20486,
            19833 / 20486,
            19583.66667 / 20486,
            19335.33333 / 20486,
            19119.33333 / 20486,
            18953.33333 / 20486,
            18825.33333 / 20486,
        ],
        "Cq_smoothing": [
            20486,
            20701,
            20491.66667,
            20271.66667,
            19833,
            19583.66667,
            19335.33333,
            19119.33333,
            18953.33333,
            18825.33333,  # <- 여기가 4번째 ddCq > 0 인 지점
        ],
    }
    obs_baseline_end = 4
    pred = {
        "cycle": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "Cq": [
            0.9678,
            0.989,
            0.9895,
            0.9895,
            0.9789,
            0.9774,
            0.976,
            0.9762,
            0.9789,
            0.9834,
        ],
        "Cq_smoothing": [
            19826.3508,
            20260.654,
            20270.897,
            20270.897,
            20053.7454,
            20023.0164,
            19994.336,
            19998.4332,
            20053.7454,
            20145.9324,
        ],
    }
    assert (
        empty_obj.rotate_cqs(obs_baseline_end, input_cq=obs, ignore_init_cycles=3)
        == pred
    )
