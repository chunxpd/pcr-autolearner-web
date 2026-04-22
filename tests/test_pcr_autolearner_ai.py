import pytest
from src.pcr_autolearner_ai import PCRPredicter, load_data, load_data_target


@pytest.fixture(scope="module")
def obj():
    """PCRPredicter 객체를 생성하고 초기화합니다."""
    # PCRPredicter 객체를 생성합니다.
    obj = PCRPredicter()
    yield obj
    del obj


def test_load_data(obj):
    """데이터를 로드하는지 테스트합니다."""
    from src.pcr_autolearner import AutoLearner

    runner = AutoLearner()
    obs_inputdata, obs_labels, obs_targets, obs_source_files = load_data(
        runner, "tests/pcr_autolearner_inputs/", normalization_mode="primitive"
    )
    assert len(obs_inputdata) == 65
    assert len(obs_labels) == 65
    assert len(obs_targets) == 41


def test_load_data_target(obj):
    """데이터를 로드하는지 테스트합니다."""
    from src.pcr_autolearner import AutoLearner

    runner = AutoLearner()
    runner.read_dir("tests/pcr_autolearner_inputs/")
    obs_inputdata, obs_labels, obs_target, obs_keys = load_data_target(
        runner, "02_InfB", normalization_mode="primitive"
    )
    pred = [
        [
            0.0,
            0.00097119,
            0.00237402,
            0.00593504,
            0.01230172,
            0.01683393,
            0.01910003,
            0.01823675,
            0.01899212,
            0.01974749,
            0.02104241,
            0.02179778,
            0.02244524,
            0.02266106,
            0.02233733,
            0.0232006,
            0.02395597,
            0.02514298,
            0.02557462,
            0.0253588,
            0.02503507,
            0.02579044,
            0.02665372,
            0.02762491,
            0.02719327,
            0.02730118,
            0.02762491,
            0.02891982,
            0.02967519,
            0.03107802,
            0.03485486,
            0.0441351,
            0.06183231,
            0.09291033,
            0.14265674,
            0.21247437,
            0.29599655,
            0.38329556,
            0.46617028,
            0.54203086,
            0.59015863,
            0.61314341,
            0.61314341,
            0.61314341,
            0.61314341,
        ],
        [
            0.0,
            0.00097119,
            0.00237402,
            0.00593504,
            0.01230172,
            0.01683393,
            0.01910003,
            0.01823675,
            0.01899212,
            0.01974749,
            0.02104241,
            0.02179778,
            0.02244524,
            0.02266106,
            0.02233733,
            0.0232006,
            0.02395597,
            0.02514298,
            0.02557462,
            0.0253588,
            0.02503507,
            0.02579044,
            0.02665372,
            0.02762491,
            0.02719327,
            0.02730118,
            0.02762491,
            0.02891982,
            0.02967519,
            0.03107802,
            0.03485486,
            0.0441351,
            0.06183231,
            0.09291033,
            0.14265674,
            0.21247437,
            0.29599655,
            0.38329556,
            0.46617028,
            0.54203086,
            0.59015863,
            0.61314341,
            0.61314341,
            0.61314341,
            0.61314341,
        ],
    ]
    assert len(obs_inputdata) == 2
    assert len(obs_inputdata[0]) == 45
    assert obs_labels.tolist() == ["Positive", "Positive"]
    assert obs_target == "02_InfB"
    assert obs_keys == [1, 2]  # 순서대로 1,2 이렇게 증가함


def test_ddcq(obj):
    import numpy as np

    obs = np.array([[0, 1, 2, 3, 4, 5], [0, 10, 20, 30, 40, 50]])
    pred = [
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
    ]
    ddcq = np.diff(np.diff(obs, axis=1), axis=1)
    assert ddcq.tolist() == pred


def test_load_data_ddcq(obj):
    """데이터를 로드하는지 테스트합니다."""
    from src.pcr_autolearner import AutoLearner

    runner = AutoLearner()
    obs_inputdata, obs_labels, obs_targets, obs_source_files = load_data(
        runner,
        "tests/pcr_autolearner_inputs/",
        normalization_mode="primitive",
        ddcq=True,
    )
    assert len(obs_inputdata) == 65
    assert len(obs_inputdata[0]) == 43  # ddcq가 되면 45-2 가 되어서 43이 됨
    assert len(obs_labels) == 65
    assert len(obs_targets) == 41


def test_load_data_target_ddcq(obj):
    """데이터를 로드하는지 테스트합니다."""
    from src.pcr_autolearner import AutoLearner

    runner = AutoLearner()
    runner.read_dir("tests/pcr_autolearner_inputs/")
    obs_inputdata, obs_labels, obs_target, obs_keys = load_data_target(
        runner, "02_InfB", normalization_mode="primitive", ddcq=True
    )
    assert len(obs_inputdata) == 2
    assert len(obs_inputdata[0]) == 43  # ddcq가 되면 45-2 가 되어서 43이 됨
    assert obs_labels.tolist() == ["Positive", "Positive"]
    assert obs_target == "02_InfB"
    assert obs_keys == [1, 2]  # 순서대로 1,2 이렇게 증가함
