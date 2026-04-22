import pytest
from src.download_2genebank import DwGeneBank2

# from Bio.Seq import UndefinedSequenceError


@pytest.fixture(scope="module")
def obj():
    config = {
        "email": "io2oi@yahoo.co.kr",
        "jobid": "tests_jobiid_00000",
        "TEMPDIR": "tests",
        "GB_STORAGE": "tests/GBfiles",
        "FA_STORAGE": "tests/FASTA",
    }
    res = DwGeneBank2(config, "aaaaaaaaa")
    yield res
    del res


def test_gen_hsp_rec(obj):
    # pred = [
    #     {"Hsp_hit-from": 4875425, "Hsp_hit-strand": 1, "Hsp_hit-to": 4875445, "gbid": "CP101406.1"},
    #     {"Hsp_hit-from": 4877573, "Hsp_hit-strand": 1, "Hsp_hit-to": 4877593, "gbid": "CP101406.1"},
    #     {"Hsp_hit-from": 3834276, "Hsp_hit-strand": 1, "Hsp_hit-to": 3834296, "gbid": "CP102346.1"},
    #     {"Hsp_hit-from": 5040025, "Hsp_hit-strand": 1, "Hsp_hit-to": 5040045, "gbid": "CP102346.1"},
    #     {"Hsp_hit-from": 3834276, "Hsp_hit-strand": 1, "Hsp_hit-to": 3834296, "gbid": "CP102344.1"},
    #     {"Hsp_hit-from": 5040025, "Hsp_hit-strand": 1, "Hsp_hit-to": 5040045, "gbid": "CP102344.1"},
    # ]
    # obs = obj.gen_hsp_rec()
    # assert pred == obs[:6]
    with pytest.raises(ValueError):
        obj.gen_hsp_rec()




def test_merge_two_hsps(obj):
    """서로 겹치는 Hsp 들은 합친다. strand 의 방향이 다르면 에러 발생한다.
    CP101406.1 첫번째
        HSP1: 100 -------- 500
    CP101406.1 두번째
        HSP1:        400 --------------  1000
    결과물
        HSP1: 100 ---------------------  1000
    """
    obs1 = {"Hsp_hit-from": 100, "Hsp_hit-strand": 1, "Hsp_hit-to": 500, "gbid": "CP101406.1"}
    obs2 = {"Hsp_hit-from": 400, "Hsp_hit-strand": 1, "Hsp_hit-to": 1000, "gbid": "CP101406.1"}
    pred = {"Hsp_hit-from": 100, "Hsp_hit-strand": 1, "Hsp_hit-to": 1000, "gbid": "CP101406.1"}
    assert pred == obj.merge_two_hsps(obs1, obs2)
    obs1 = {"Hsp_hit-from": 100, "Hsp_hit-strand": 1, "Hsp_hit-to": 500, "gbid": "CP101406.1"}
    obs2 = {"Hsp_hit-from": 600, "Hsp_hit-strand": 1, "Hsp_hit-to": 1000, "gbid": "CP101406.1"}
    pred = False
    assert pred == obj.merge_two_hsps(obs1, obs2)
    obs1 = {"Hsp_hit-from": 1001, "Hsp_hit-strand": 1, "Hsp_hit-to": 1500, "gbid": "CP101406.1"}
    obs2 = {"Hsp_hit-from": 600, "Hsp_hit-strand": 1, "Hsp_hit-to": 1000, "gbid": "CP101406.1"}
    pred = False
    assert pred == obj.merge_two_hsps(obs1, obs2)


def test_combine_hsp_recs(obj):
    """서로 겹치는 Hsp 들은 합친다. strand 의 방향이 다르면 에러 발생한다.
    CP101406.1 첫번째
        HSP1: 100 -------- 500
        HSP2:                                        1100 --------------- 1500
    CP101406.1 두번째
        HSP1:        400 --------------  1000
        HSP2:                                                                      2100 ---------- 2600
    결과물
        HSP1: 100 ---------------------  1000
        HSP2:                                        1100 --------------- 1500
        HSP3:                                                                      2100 ---------- 2600
    """
    obs1 = [
        {"Hsp_hit-from": 100, "Hsp_hit-strand": 1, "Hsp_hit-to": 500, "gbid": "CP101406.1"},
        {"Hsp_hit-from": 1100, "Hsp_hit-strand": 1, "Hsp_hit-to": 1500, "gbid": "CP101406.1"},
        {"Hsp_hit-from": 500, "Hsp_hit-strand": 1, "Hsp_hit-to": 1000, "gbid": "CP102346.1"},
    ]
    obs2 = [
        {"Hsp_hit-from": 400, "Hsp_hit-strand": 1, "Hsp_hit-to": 1000, "gbid": "CP101406.1"},
        {"Hsp_hit-from": 2100, "Hsp_hit-strand": 1, "Hsp_hit-to": 2600, "gbid": "CP101406.1"},
        {"Hsp_hit-from": 10000, "Hsp_hit-strand": 1, "Hsp_hit-to": 11000, "gbid": "CP102346.1"},
    ]
    pred = [
        {"Hsp_hit-from": 100, "Hsp_hit-strand": 1, "Hsp_hit-to": 1000, "gbid": "CP101406.1"},
        {"Hsp_hit-from": 2100, "Hsp_hit-strand": 1, "Hsp_hit-to": 2600, "gbid": "CP101406.1"},
        {"Hsp_hit-from": 10000, "Hsp_hit-strand": 1, "Hsp_hit-to": 11000, "gbid": "CP102346.1"},
        {"Hsp_hit-from": 1100, "Hsp_hit-strand": 1, "Hsp_hit-to": 1500, "gbid": "CP101406.1"},
        {"Hsp_hit-from": 500, "Hsp_hit-strand": 1, "Hsp_hit-to": 1000, "gbid": "CP102346.1"},
    ]
    assert pred == obj.combine_hsp_recs(obs1, obs2)
