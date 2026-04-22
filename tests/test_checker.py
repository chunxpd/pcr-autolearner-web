import pytest
from src.check_status import ChkStatus


@pytest.fixture(scope="module")
def obj():
    mock_config = {
        "TEMPDIR": "tests",
        "jobid": "tests_jobiid_00000",
    }
    res = ChkStatus(mock_config, "aaaaaaaaa")
    yield res
    del res


def test_check_num_fasta(obj):
    pred = 50
    obs = obj.check_num_fasta()
    assert obs == pred
