import datetime
import pytest
from src.common import Common


@pytest.fixture
def mock_datetime(monkeypatch):
    specific_date = datetime.datetime(2023, 1, 1)

    class MockDatetime:
        @classmethod
        def now(cls):
            return specific_date

    monkeypatch.setattr(datetime, "datetime", MockDatetime)


def test_suffix_replace():
    obs_fn, obs_suffix, obs_tar_suff = ("foo.fa.fa", ".fa", ".log")
    obs = Common.suffix_replace(obs_fn, obs_suffix, obs_tar_suff)
    pred = "foo.fa.log"
    assert pred == obs
    obs_fn, obs_suffix, obs_tar_suff = ("foo.fa.fasta", ".fasta", ".log")
    obs = Common.suffix_replace(obs_fn, obs_suffix, obs_tar_suff)
    pred = "foo.fa.log"
    assert pred == obs
    obs_fn, obs_suffix, obs_tar_suff = ("foo.txt", ".fa", ".log")
    obs = Common.suffix_replace(obs_fn, obs_suffix, obs_tar_suff)
    pred = ""
    assert pred == obs


def test_get_date_string(mock_datetime):
    pred = "20230101"
    obs = Common.get_date_string()
    assert pred == obs


def test_get_rmgap_fa():
    obs = "test.fa"
    pred = "test_rmgap.fa"
    assert pred == Common.get_rmgap_fa_fn(obs)
    obs = "testfa.fa"
    pred = "testfa_rmgap.fa"
    assert pred == Common.get_rmgap_fa_fn(obs)


def test_change_fa_suffix_with():
    obs = "test.fa"
    pred = "test.helo"
    assert pred == Common.change_fa_suffix_with(obs, ".helo")
    obs = "test.fasta"
    pred = "test.helo"
    assert pred == Common.change_fa_suffix_with(obs, ".helo")
    obs = "test.ff"
    with pytest.raises(
        ValueError,
        match=f"input fasta file MUST end with .fa or .fasta. But your file is {obs}",
    ):
        Common.change_fa_suffix_with(obs, ".helo")
