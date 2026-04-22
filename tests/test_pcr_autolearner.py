import pytest
from pytest import approx
from src.pcr_autolearner import AutoLearner


def compare_dicts_with_approx(dict1, dict2):
    if dict1.keys() != dict2.keys():
        return False
    for key in dict1:
        if isinstance(dict1[key], tuple):
            if len(dict1[key]) != len(dict2[key]):
                return False
            for i in range(len(dict1[key])):
                if dict1[key][i] != approx(dict2[key][i]):
                    return False
        else:
            if dict1[key] != approx(dict2[key]):
                return False
    return True


@pytest.fixture(scope="function")
def obj():
    res = AutoLearner()
    res.read_dir("tests/pcr_autolearner_inputs/")
    # res.set_target("07_RnaseP2")")
    yield res
    del res


def test_auto_learner_initialization():
    learner = AutoLearner()
    assert learner is not None


def test_file_version_check():
    """파일 버전 체크를 한다. 버전은 현재 2가지를 지원한다.
    version1: ----, ==== 이런 것이 없고 RFU 값만 있음
    version2: ----, ==== 이런 것이 있고 RFU, SUBSTRACTION, NORMALIZATION 값이 있음
    """
    res = AutoLearner()
    assert (
        res.check_file_version(
            "tests/pcr_autolearner_inputs/pcr_autolearner_input_20250525_170000_1.csv"
        )
        == "version1"
    )
    assert (
        res.check_file_version(
            "tests/pcr_autolearner_inputs/250528_pcr_autolearner_input_version2.csv"
        )
        == "version2"
    )
    assert (
        res.check_file_version(
            "tests/pcr_autolearner_inputs/DX10024F2001_20260220_101320_MRMP_edmund_version3.csv"
        )
        == "version3"
    )


def test_read_input():
    """한개의 csv 파일을 읽는다. input cycle 은 무조건 45 cycle 로 맞춘다"""
    res = AutoLearner()
    res.read_input(
        "tests/pcr_autolearner_inputs/pcr_autolearner_input_20250525_170000_1.csv"
    )
    rec = res.cq_container[0]
    assert rec.cycles == list(range(45))
    assert rec.sample_names == {
        1: "FAM_Ch#1",
        2: "FAM_Ch#2",
        3: "FAM_Ch#3",
        4: "FAM_Ch#4",
        5: "Cal560_Ch#1",
        6: "Cal560_Ch#2",
        7: "Cal560_Ch#3",
        8: "Cal560_Ch#4",
        9: "Cal610_Ch#1",
        10: "Cal610_Ch#2",
        11: "Cal610_Ch#3",
        12: "Cal610_Ch#4",
        13: "Quasar670_Ch#1",
        14: "Quasar670_Ch#2",
        15: "Quasar670_Ch#3",
        16: "Quasar670_Ch#4",
    }
    assert rec.meta_info[7] == {
        "Target Name": "07_RnaseP2",
        "Th_nor": "0.5",
        "Th_sub": "1000",
        "Cutoff": "40",
        "Ct Value": "0.0",
        "Pos/Neg": "",
    }
    assert rec.Cq_values[1] == [
        3011.0,
        3061.0,
        3011.0,
        3066.0,
        3026.0,
        3025.0,
        3041.0,
        3022.0,
        3030.0,
        3031.0,
        3082.0,
        3082.0,
        3036.0,
        3083.0,
        3089.0,
        3038.0,
        3085.0,
        3094.0,
        3037.0,
        3090.0,
        3039.0,
        3045.0,
        3098.0,
        3099.0,
        3101.0,
        3065.0,
        3081.0,
        3166.0,
        3241.0,
        3371.0,
        3652.0,
        4169.0,
        4909.0,
        5776.0,
        6530.0,
        7159.0,
        7510.0,
        7775.0,
        8002.0,
        8091.0,  # <- 원래는 여기가 40 cycle 로 마지막인데 45 cycle 로 늘리면서 중복해서 5번 반복
        8091.0,
        8091.0,
        8091.0,
        8091.0,
        8091.0,
    ]


def test_read_input_version2():
    """한개의 csv 파일을 읽는다. input cycle 은 무조건 45 cycle 로 맞춘다"""
    res = AutoLearner()
    res.read_input(
        "tests/pcr_autolearner_inputs/250528_pcr_autolearner_input_version2.csv"
    )
    rec = res.cq_container[0]
    assert (
        rec.fn
        == "tests/pcr_autolearner_inputs/250528_pcr_autolearner_input_version2.csv"
    )
    assert rec.cycles == list(range(45))
    assert rec.sample_names == {
        1: "FAM_Ch#1",
        2: "FAM_Ch#2",
        3: "FAM_Ch#3",
        4: "FAM_Ch#4",
        5: "Cal560_Ch#1",
        6: "Cal560_Ch#2",
        7: "Cal560_Ch#3",
        8: "Cal560_Ch#4",
        9: "Cal610_Ch#1",
        10: "Cal610_Ch#2",
        11: "Cal610_Ch#3",
        12: "Cal610_Ch#4",
        13: "Quasar670_Ch#1",
        14: "Quasar670_Ch#2",
        15: "Quasar670_Ch#3",
        16: "Quasar670_Ch#4",
    }
    assert rec.meta_info[7] == {
        "Target Name": "07_NDM",
        "Th_nor": "0.1",
        "Th_sub": "1000",
        "Cutoff": "40",
        "Ct Value": "0.0",
        "Pos/Neg": "",
    }
    assert rec.Cq_values[1] == [
        1169,
        1165,
        1166,
        1165,
        1165,
        1163,
        1165,
        1166,
        1166,
        1167,
        1167,
        1167,
        1165,
        1167,
        1167,
        1167,
        1165,
        1166,
        1166,
        1168,
        1167,
        1168,
        1168,
        1169,
        1168,
        1169,
        1168,
        1168,
        1168,
        1169,
        1167,
        1169,
        1168,
        1169,
        1168,
        1170,
        1170,
        1169,
        1170,
        1168,
        1168,
        1168,
        1168,
        1168,
        1168,
    ]


def test_read_input_version3():
    """한개의 csv 파일을 읽는다. input cycle 은 무조건 45 cycle 로 맞춘다"""
    res = AutoLearner()
    res.read_input(
        "tests/pcr_autolearner_inputs/DX10024F2001_20260220_101320_MRMP_edmund_version3.csv"
    )
    rec = res.cq_container[0]
    assert (
        rec.fn
        == "tests/pcr_autolearner_inputs/DX10024F2001_20260220_101320_MRMP_edmund_version3.csv"
    )
    assert rec.cycles == list(range(45))
    assert rec.sample_names == {
        1: "FAM_Ch#1",
        2: "FAM_Ch#2",
        3: "FAM_Ch#3",
        4: "FAM_Ch#4",
        5: "Cal560_Ch#1",
        6: "Cal560_Ch#2",
        7: "Cal560_Ch#3",
        8: "Cal560_Ch#4",
        9: "Cal610_Ch#1",
        10: "Cal610_Ch#2",
        11: "Cal610_Ch#3",
        12: "Cal610_Ch#4",
        13: "Quasar670_Ch#1",
        14: "Quasar670_Ch#2",
        15: "Quasar670_Ch#3",
        16: "Quasar670_Ch#4",
    }
    assert rec.meta_info[1] == {
        "Target Name": "01_Blue#01",
        "Pos/Neg": "Abnormal",
        "Ct Value": "8.06",
        "Cutoff": "40",
        "Threshold": "0",
    }
    assert rec.meta_info[16] == {
        "Target Name": "16_Red#04",
        "Pos/Neg": "Positive",
        "Ct Value": "15.19",
        "Cutoff": "40",
        "Threshold": "0",
    }
    assert rec.Cq_values[1][:5] == [
        349.0,
        348.0,
        350.0,
        351.0,
        350.0,
    ]


def test_read_timestamp(obj):
    """read_input() 할 때 fn 에서 timestamp 를 읽어온다.
    pcr_autolearner_input_20250525_170000_1.csv 의 경우 20250525_170000 이다."""
    pred = "20250525_170000"
    obs = obj.parse_timestamp(
        "tests/pcr_autolearner_inputs//pcr_autolearner_input_20250525_170000_1.csv"
    )
    assert obs == pred

    pred = "20260220_101320"
    obs = obj.parse_timestamp(
        "tests/pcr_autolearner_inputs/DX10024F2001_20260220_101320_MRMP_edmund_version3.csv"
    )
    assert obs == pred


def test_read_directory():
    """디렉토리 내의 모든 csv 파일들을 읽는다.
    read_dir() 할 때 CqValue class 에 fn 외에 시간을 파싱해서 넣어 준다.
    그 다음 pcr_autolearner::set_target에서 sample_names, Cq_values, meta_info, Cq_init_value 의 key로 iii_CqValue.gen_time 으로 key 를 잡는다.
    """
    res = AutoLearner()
    res.read_dir("tests/pcr_autolearner_inputs/")
    assert len(res.cq_container) == 5
    assert res.cq_container[0].cycles == list(range(45))
    assert (
        res.cq_container[0].fn
        == "tests/pcr_autolearner_inputs//250528_pcr_autolearner_input_simple_version2.csv"
    )
    assert (
        res.cq_container[1].fn
        == "tests/pcr_autolearner_inputs//250528_pcr_autolearner_input_version2.csv"
    )
    assert res.cq_container[1].gen_timestamp == "250528"
    assert (
        res.cq_container[3].fn
        == "tests/pcr_autolearner_inputs//pcr_autolearner_input_20250525_170000_1.csv"
    )
    assert res.cq_container[3].gen_timestamp == "20250525_170000"
    assert (
        res.cq_container[4].fn
        == "tests/pcr_autolearner_inputs//pcr_autolearner_input_20250525_170800_2.csv"
    )
    assert res.cq_container[4].gen_timestamp == "20250525_170800"
    assert res.abnormal_files == [
        "tests/pcr_autolearner_inputs//empty.csv",
        # "tests/pcr_autolearner_inputs//test.ini", # csv가 아닌 파일은 무시함
    ]


def test_get_target_names(obj):
    """디렉토리 내의 모든 csv 파일들을 읽은 뒤 각 파일내에 있는 item 이름들을 목록화해서 반환받는다."""
    pred = sorted(
        [
            "01_Blue#01",
            "01_InfA",
            "02_Blue#02",
            "02_IMP",
            "02_InfB",
            "03_Blue#03",
            "03_KPC",
            "03_RnaseP1",
            "04_Blue#04",
            "04_GES",
            "04_Sars-Cov2_N",
            "05_Green#01",
            "05_RSVA",
            "06_Green#02",
            "06_RSVB",
            "06_VIM",
            "07_Green#03",
            "07_NDM",
            "07_RnaseP2",
            "08_Entero",
            "08_Green#04",
            "08_OXA-48",
            "09_CoV_E229",
            "09_Orange#01",
            "10_BS",
            "10_CoV_NL63",
            "10_Orange#02",
            "11_BS",
            "11_Orange#03",
            "11_RnaseP3",
            "12_Adeno",
            "12_BS",
            "12_Orange#04",
            "13_Boca",
            "13_Red#01",
            "14_PIV_1_2_3_4",
            "14_Red#02",
            "15_Red#03",
            "15_RnaseP4",
            "16_External control",
            "16_Red#04",
        ]
    )
    assert obj.get_target_names() == pred


def test_set_target(obj):
    """입력 디렉토리에 있는 파일들에서 item 에 대한 target name이 같은 것들로만
    meta_info, Cq_values, sample_names 등을 만든다.
    원본 값들을 복사해 오는 형식이다. 원본값 변형 없음
    """
    obj.set_target("07_RnaseP2")
    assert obj.meta_info == {
        1: {
            "Target Name": "07_RnaseP2",
            "Th_nor": "0.5",
            "Th_sub": "1000",
            "Cutoff": "40",
            "Ct Value": "0.0",
            "Pos/Neg": "",
        },
        2: {
            "Target Name": "07_RnaseP2",
            "Th_nor": "0.5",
            "Th_sub": "1000",
            "Cutoff": "40",
            "Ct Value": "0.0",
            "Pos/Neg": "",
        },
    }
    assert obj.sample_names == {1: "Cal560_Ch#3", 2: "Cal560_Ch#3"}
    assert obj.cycles == list(range(45))
    assert obj.gen_timestamps == {1: "20250525_170000", 2: "20250525_170800"}
    assert obj.file_names == {
        1: "tests/pcr_autolearner_inputs//pcr_autolearner_input_20250525_170000_1.csv",
        2: "tests/pcr_autolearner_inputs//pcr_autolearner_input_20250525_170800_2.csv",
    }


def test_set_target_kpc(obj):
    obj.set_target("03_KPC")
    pred_meta = {
        "Target Name": "03_KPC",
        "Th_nor": "0.1",
        "Th_sub": "1500",
        "Cutoff": "40",
        "Ct Value": "32.2",
        "Pos/Neg": "Positive",
    }
    assert obj.meta_info == dict(zip(range(1, 3), [pred_meta] * 2))
    assert obj.sample_names == {
        1: "FAM_Ch#3",
        2: "FAM_Ch#3",
    }
    assert obj.cycles == list(range(45))
    assert obj.gen_timestamps == {
        1: "250528",
        2: "250528",
    }
    assert obj.Cq_values == {
        1: [
            100.0,
            200.0,
            300.0,
            400.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
            500.0,
        ],
        2: [
            6588.0,
            6243.0,
            6144.0,
            6058.0,
            6038.0,
            6003.0,
            5970.0,
            5943.0,
            5951.0,
            5930.0,
            5906.0,
            5908.0,
            5885.0,
            5869.0,
            5852.0,
            5846.0,
            5837.0,
            5814.0,
            5825.0,
            5823.0,
            5795.0,
            5798.0,
            5781.0,
            5779.0,
            5777.0,
            5781.0,
            5779.0,
            5795.0,
            5815.0,
            5887.0,
            5977.0,
            6180.0,
            6442.0,
            6768.0,
            7102.0,
            7426.0,
            7736.0,
            8027.0,
            8290.0,
            8538.0,
            8538.0,
            8538.0,
            8538.0,
            8538.0,
            8538.0,
        ],
    }


def test_get_expected_value(obj):
    obj.set_target("07_RnaseP2")
    assert obj.get_expected_value(1, key="Pos/Neg") == ""
    assert obj.get_expected_value(2, key="Pos/Neg") == ""
    # assert obj.get_expected_value(7, key="Pos/Neg") == ""


def test_get_sample_size(obj):
    obj.set_target("07_RnaseP2")
    pred = 2
    obs = obj.get_sample_size()
    assert obs == pred
    obj.set_target("02_InfB")
    pred = 2
    obs = obj.get_sample_size()
    assert obs == pred


def test_get_pos_num(obj):
    "지금은 Positive 만 +1, Negative 와 빈칸은 negative 임"
    obj.set_target("07_RnaseP2")
    pred = 0
    obs = obj.get_pos_num()
    assert obs == pred
    obj.set_target("02_InfB")
    pred = 2
    obs = obj.get_pos_num()
    assert obs == pred


def test_baseline_endq(obj):
    """Cq 값들을 모두 입력한 다음 baseline_end 가 허용된 cycle 내에서 있으면 1 아니면 0
    InfB 에 대해서는 28 cycle 에서 baseline end 가 있음"""
    obj.set_target("02_InfB")
    pred = (
        1  # positive case sample number 1: scale factor 0.01일때 26cycle 에서 baseline end.
    )
    obs = obj.baseline_endq(1, cycle_range=(3, 38))
    assert obs == pred
    pred = 1
    obs = obj.baseline_endq(2, cycle_range=(3, 38))
    assert obs == pred
    pred = 0  # negative case. baseline end at 28 cycle.
    obs = obj.baseline_endq(2, cycle_range=(3, 27), scaling_factor=0.0)
    assert obs == pred


def test_cal_loss_function(obj):
    """https://github.com/io2oi/ctlabs-primerwrapper/issues/264#issue-2804358209
    calculating loss function
    L = 1/K * Sum_i (y_i - e_i)
    테스트입력에서 Positive 13개, Negative 0개, no label 3개임.
    그런데 Negative 가 실제로는 Positive
    여기서 Loss function이 0이 되도록 cycle range, scaling factor 등을 조정해 나가야 함
    K=n(Positive),
    y_i = 1 when the ith sample is positive and
    y_i = 0 when the ith sample is negative,
    e_i는 i번째 샘플이 positive 일 때 baseline end 가 허용범위 cycle 내에서 있으면 1 아니면 0
    """
    obj.set_target("02_InfB")
    pred = 1 / 2 * ((1 - 1) + (1 - 1))
    obs = obj.cal_loss_function(cycle_range=(3, 35))
    assert obs == pred


def test_iteration_training_scaling_factor(obj):
    """https://github.com/io2oi/ctlabs-primerwrapper/issues/267
    loss function 이 0이 아니면 답이 틀렸다는 이야기.
    - 일 경우: negative 를 positive 로 인식 했다는 뜻. scaling factor 를 늘려본다.
    + 일 경우: Positive 를 negative 로 인식 했다는 뜻. scaling factor 를 줄여본다.
    """
    obj.set_target("02_InfB")
    ## step size 를 작게하고 iteration 을 늘리면 좀 더 세밀하게 scaling_factor 범위를 찾음
    pred = {
        "scaling_factor": (0.002, 0.093),
        "loss_value": (0.0, 0.0),
    }
    obs = obj.iteration_training_scaling_factor(
        cycle_range=(3, 35), epsilon=1e-4, max_iteration=100, step_size=1e-3
    )
    assert compare_dicts_with_approx(obs, pred)
    ## step size 를 크게 하고 iteration 을 줄이면 왔다갔다 하다가 멈춤 범위를 찾음
    pred = {"scaling_factor": (0.02, 0.090), "loss_value": (0.0, 0.0)}
    obs = obj.iteration_training_scaling_factor(
        cycle_range=(3, 35), epsilon=1e-4, max_iteration=10, step_size=1e-2
    )
    assert compare_dicts_with_approx(obs, pred)


def test_find_max_range(obj):
    """[(cycle_range, loss_value), ...]
    가 있을 때 loss_value 의 절대값이 작은 것부터 큰 것까지
    그리고 최대 범위의 cycle_range 를 반환
    """
    obs_input = [((0, 1), 0.01), ((0, 10), 0.01), ((3, 39), 0.01), ((10, 30), -0.02)]
    pred = ((3, 39), 0.01)
    obs = obj.find_max_range(obs_input)
    assert obs == pred
    obs_input = [((0, 40), -0.01), ((5, 20), 0.01), ((3, 18), 0.1)]
    pred = ((0, 40), -0.01)
    obs = obj.find_max_range(obs_input)
    assert obs == pred
    obs_input = [
        ((3, 33), 0.23076923076923078),
        ((3, 34), 0.07692307692307693),
        ((3, 35), 0.0),
        ((3, 36), -0.23076923076923078),
        ((3, 37), -0.23076923076923078),
        ((3, 38), -0.23076923076923078),
        ((4, 33), 0.23076923076923078),
        ((4, 34), 0.07692307692307693),
        ((4, 35), 0.0),
        ((4, 36), -0.23076923076923078),
        ((4, 37), -0.23076923076923078),
        ((4, 38), -0.23076923076923078),
    ]
    pred = ((3, 35), 0.0)
    obs = obj.find_max_range(obs_input)
    assert obs == pred


def test_iteration_training_cycle_range(obj):
    """After getting scaling_factor"""
    obj.set_target("02_InfB")
    pred = {
        "cycle_range": ((3, 43), (3, 43)),
        "loss_value": (0.0, 0.0),
    }
    obs = obj.iteration_training_cycle_range((0.01, 0.017))
    assert obs == pred


def test_find_horizontal_regions(obj):
    pred = [(4, 10)]
    obs_cq_values = [
        0,
        1,
        2,
        3,
        3,
        3,
        3,
        3,
        3,
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
    ]

    cycles = range(1, 21)  # cycle 이 1부터 시작하므로 1~20 까지의 range 를 사용
    obs = obj.find_horizontal_regions(
        obs_cq_values,
        cycles,
    )
    assert obs == pred


def test_find_horizontal_regions_init_noise(obj):
    """처음 3 cycle (0-based) 부분에 noise 가 있는 경우 무시.
    아래의 경우 3~7 cycle이 수평구간이지만 3cycle 이 포함되어 있으므로 noise 로 취급한다."""
    pred = [(4, 7), (11, 17)]
    obs_cq_values = [
        -1,  # 1 cycle
        2,
        5,  # noise 로 취급하므로 제거. 3 cycle
        5,
        5,
        5,
        5,  # 7 cycle
        6,
        7,
        8,
        9,  # 11 cycle
        9,
        9,
        9,
        9,
        9,
        9,  # 17 cycle
        10,
        11,
        12,
        13,
        14,
    ]

    cycles = range(1, 23)  # cycle 이 1부터 시작하므로 1~20 까지의 range 를 사용
    obs = obj.find_horizontal_regions(obs_cq_values, cycles, min_length=4)
    assert obs == pred


def test_find_horizontal_region_realrfu(obj):
    """원본 RFU plot png 파일을 읽어서 horizontal region 을 찾는다.
    첫 region 의 시작과 끝을 이용하면 될 듯.
    issue 312 참조"""
    pred = [(12, 29), (39, 44)]
    # 아래는 Adeno 2 sample. Negative 로 잘 못 labeling 되어있는 것이라고 회의 시간에 이야기 하였음
    obs_cq_values = [
        2236.0,
        2200.0,
        2388.0,
        2551.0,
        2701.0,
        2862.0,
        2991.0,
        3056.0,
        3138.0,
        3196.0,
        3309.0,
        3389.0,
        3476.0,
        3500.0,
        3542.0,
        3552.0,
        3575.0,
        3551.0,
        3563.0,
        3548.0,
        3552.0,
        3568.0,
        3579.0,
        3574.0,
        3598.0,
        3610.0,
        3593.0,
        3628.0,
        3655.0,
        3721.0,
        3813.0,
        4002.0,
        4312.0,
        4633.0,
        4930.0,
        5195.0,
        5421.0,
        5638.0,
        5834.0,
        6015.0,
        6015.0,
        6015.0,
        6015.0,
        6015.0,
        6015.0,
    ]
    cycles = range(45)  # cycle 이 0부터 시작하므로 0~44 까지의 range 를 사용
    obs = obj.find_horizontal_regions(obs_cq_values, cycles)
    assert obs == pred
    cycles = range(1, 46)  # CYCLE 1 부터 시작하므로 1~45 까지의 range 를 사용
    pred = [(13, 30), (40, 45)]
    obs = obj.find_horizontal_regions(obs_cq_values, cycles)
    assert obs == pred


def test_baseline_substraction(obj):
    """선택된 item target name 에 대한 Cq 값들에 대해서
    1. horizontal region 을 찾는다. 못 찾으면 3~8 cycle을 horizontal region 으로 설정
    2. pcr_analyzer.add_one_cq_value()를 이용해서 값을 넣은 뒤
     'Cq_smooting' 값을 이용해서 baseline 의 평균 값을 구하고
    3. rotation 을  horizontal region 의 끝을 cycle 을 중심으로 해준다
    4. 그런다음 self.Cq_values[item]에 substraction 값을 넣는다.

    아래 테스트에서 5 cycle 이 45cycle에 맞추어 더 들어가 있어서 +5 와 [:40]이 붙어있음
    """
    obj.set_target("02_InfB")
    pred = [
        -40.333514285714045,
        -39.09791428571407,
        -36.62671428571366,
        -27.66861428571383,
        -10.061314285713706,
        1.9857857142860666,
        7.2370857142864224,
        2.603585714286055,
        2.912485714286049,
        3.2213857142860434,
        5.383685714286003,
        5.692585714285997,
        5.692585714285997,
        4.45698571428602,
        1.3679857142860783,
        2.294685714286061,
        2.603585714286055,
        4.148085714286026,
        3.5302857142860375,
        1.0590857142860841,
        -2.0299142857138577,
        -1.412114285714324,
        -0.7943142857143357,
        0.13238571428655632,
        -2.95661428571384,
        -4.810014285713805,
        -5.736714285713788,
        -3.5744142857138286,
        -3.2655142857138344,
        -0.7943142857143357,
        8.781585714285939,
        35.65588571428589,
        88.16888571428626,
        182.3833857142863,
        334.0532857142862,
        547.8120857142862,
        803.8901857142864,
        1071.3975857142868,
        1325.6222857142857,
        1557.9150857142868,
    ]
    obj.baseline_substraction()
    assert len(obj.Cq_values[2]) == len(pred) + 5
    assert obj.Cq_values[1][:40] == pred


def test_normalization(obj):
    """baseline substraction 이후에 normalization 을 해준다.
    normalization 은 첫번째 cycle의 값을 이용해서 모든 baseline substraction 값들을 나누어 주면 된다.
    2번째 csv 에 있는 InfB 의 경우 첫번째 cycle 값이 3089.0 이다.
    마지막 5개 cycle 은 45 cycle 이기 때문에 5번 반복된 값이 들어가 있다. (rotation 때문에 값들이 약간 변경됨)
    """
    obj.set_target("02_InfB")
    pred = [
        value / 3089.0
        for value in [
            -40.333514285714045,
            -39.09791428571407,
            -36.62671428571366,
            -27.66861428571383,
            -10.061314285713706,
            1.9857857142860666,
            7.2370857142864224,
            2.603585714286055,
            2.912485714286049,
            3.2213857142860434,
            5.383685714286003,
            5.692585714285997,
            5.692585714285997,
            4.45698571428602,
            1.3679857142860783,
            2.294685714286061,
            2.603585714286055,
            4.148085714286026,
            3.5302857142860375,
            1.0590857142860841,
            -2.0299142857138577,
            -1.412114285714324,
            -0.7943142857143357,
            0.13238571428655632,
            -2.95661428571384,
            -4.810014285713805,
            -5.736714285713788,
            -3.5744142857138286,
            -3.2655142857138344,
            -0.7943142857143357,
            8.781585714285939,
            35.65588571428589,
            88.16888571428626,
            182.3833857142863,
            334.0532857142862,
            547.8120857142862,
            803.8901857142864,
            1071.3975857142868,
            1325.6222857142857,
            1557.9150857142868,
        ]
    ]
    obj.normalization(smoothing_size=3)
    assert obj.Cq_values[2][:40] == pred


def test_normalization_2(obj):
    """baseline substraction 이후에 normalization 을 해준다.
    normalization primitive 는 baseline substraction 할 때에 rotation 을 하지 않고
    단순히 첫 cycle의 값으로 baseline substraction 만 해준다.
    """
    obj.set_target("02_InfB")
    pred = [
        0.0,
        0.0009711880867594691,
        0.002374015323189862,
        0.005935038307974583,
        0.012301715765619942,
        0.016833926837164132,
        0.019100032372936226,
        0.018236754073594524,
        0.018992122585518457,
        0.01974749109744254,
        0.021042408546455164,
        0.021797777058379243,
        0.022445235782885556,
        0.022661055357720946,
        0.02233732599546779,
        0.02320060429480949,
        0.02395597280673357,
        0.025142980468328427,
        0.025574619617999353,
        0.025358800043163963,
        0.025035070680910806,
        0.02579043919283474,
        0.026653717492176588,
        0.027624905578936057,
        0.027193266429265135,
        0.0273011762166829,
        0.027624905578936057,
        0.028919823027948682,
        0.02967519153987262,
        0.03107801877630301,
        0.03485486133592312,
        0.04413510305384703,
        0.06183230819035287,
        0.09291032696665588,
        0.1426567389662242,
        0.21247437142548833,
        0.29599654688680266,
        0.38329556490773714,
        0.46617028164454516,
        0.5420308621992014,
    ]
    obj.normalization(smoothing_size=3, mode="primitive")
    assert obj.Cq_values[2][:40] == pred


def test_sigmoid_fitting(obj):
    """sigmoid fitting을 한다. #320 참조"""
    obs_x = range(1, 41)
    obs_y = [
        0.00000000e00,
        1.06581410e-14,
        4.37871961e-12,
        3.27679217e-10,
        9.31322575e-09,
        1.43489068e-07,
        1.44884058e-06,
        1.07374067e-05,
        6.28325923e-05,
        3.05166468e-04,
        1.27463248e-03,
        4.69964015e-03,
        1.55963323e-02,
        4.72512867e-02,
        1.31872339e-01,
        3.39885078e-01,
        8.03364879e-01,
        1.70737745e00,
        3.16609039e00,
        5.00000000e00,
        6.75211651e00,
        8.06847196e00,
        8.90555625e00,
        9.39050485e00,
        9.66011492e00,
        9.80837733e00,
        9.89030092e00,
        9.93613003e00,
        9.96217021e00,
        9.97721545e00,
        9.98605500e00,
        9.99133390e00,
        9.99453609e00,
        9.99650768e00,
        9.99773884e00,
        9.99851800e00,
        9.99901739e00,
        9.99934133e00,
        9.99955387e00,
        9.99969483e00,
    ]
    pred = {"y0": 0, "a": 10, "x0": 20, "b": -15}
    obs = obj.sigmoid_fitting(obs_x, obs_y, function="logistic")
    for key in pred:
        assert obs[key] == approx(pred[key], rel=1e-5, abs=1e-5), f"Key {key} failed"


def test_read_input_version1_simple():
    res = AutoLearner()
    res.read_input(
        "tests/pcr_autolearner_inputs/250528_pcr_autolearner_input_simple_version2.csv"
    )
    rec = res.cq_container[0]
    assert rec.cycles == list(range(45))
    assert rec.sample_names == {
        1: "FAM_Ch#1",
        2: "FAM_Ch#2",
        3: "FAM_Ch#3",
        4: "FAM_Ch#4",
    }
    assert rec.meta_info[3] == {
        "Target Name": "03_KPC",
        "Th_nor": "0.1",
        "Th_sub": "1500",
        "Cutoff": "40",
        "Ct Value": "32.2",
        "Pos/Neg": "Positive",
    }
    assert rec.Cq_values[3] == [100, 200, 300, 400, 500] + [500] * 40


def test_get_key_from_label_additional(obj):
    # empty or None -> None
    assert obj.get_key_from_label("") is None
    assert obj.get_key_from_label(None) is None

    # leading spaces then number
    assert obj.get_key_from_label("   42_some") == 42

    # underscore with zero-padded number -> int conversion drops leading zeros
    assert obj.get_key_from_label("label_007") == 7

    # digits embedded: re.search finds the first numeric group
    assert obj.get_key_from_label("abc123def456") == 123

    # plus sign before number: match fails but search finds the digit
    assert obj.get_key_from_label("+3") == 3

    # no digits at all -> None
    assert obj.get_key_from_label("no digits here") is None

    assert obj.get_key_from_label("3_20240529_120000") == 3

    assert obj.get_key_from_label("3 (tr_P:pr_N)") == 3

    assert obj.get_key_from_label("3") == 3
