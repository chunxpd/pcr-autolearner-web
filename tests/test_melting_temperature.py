# pylint: disable=missing-module-docstring
# pylint: disable=missing-function-docstring
# pylint: disable=redefined-outer-name
import pytest
from src.cal_melting_temperature import MeltingTemp


@pytest.fixture(scope="module")
def obj():
    res = MeltingTemp(isdeoxy=True)
    yield res
    del res


def test_cal_gc(obj):
    obs_seq = "AAAGGGGCCCGGAAGGGGACGGTTG"
    obj.set_sequence(obs_seq)
    pred = 100 * 17 / 25
    assert pred, pred == obj.cal_gc(obs_seq)


def test_set_sequence(obj):
    obs_seq = "AAAGGGGCCCGGAAGGGGACGGTTG"
    obj.set_sequence(obs_seq)
    assert "min" == obj.choice
    assert 6 == obj.a_count
    assert 2 == obj.t_count
    assert 13 == obj.g_count
    assert 4 == obj.c_count
    assert 0 == obj.s_count
    assert 17 == obj.e_gc_min
    assert len(obs_seq) - 8 == obj.e_gc_max  # 17
    assert 100 * 17 / 25 == obj.gc_val_min


def test_adj_tm(obj):
    """refer to http://biotools.nubic.northwestern.edu/OligoCalc.html
    I've check the concordance of the results of the website and the equation below
    Tm= 100.5 + (41 * (yG+zC)/(wA+xT+yG+zC)) - (820/(wA+xT+yG+zC)) + 16.6*log10([Na+])
    nA = 6
    nT = 2
    nG = 13
    nC = 4
    len = 25
    Na+ = 0.050 # 50mM -> 0.050 microM
    Tm = 100.5 + (41 * 17 / 25) - (820/25) + 16.6*log10(0.050)
    salt concentration: Na+ 50 mM
    """
    pred = 73.98
    obs_seq = "AAAGGGGCCCGGAAGGGGACGGTTG"
    obj.set_sequence(obs_seq)
    assert pred == obj.cal_adj_melting_temp(salt_concentration=50)


def test_count_neighbors(obj):
    obs = "ATATGGG"
    assert 2 == obj.count_neighbors("GG", obs)
    assert 2 == obj.count_neighbors("AT", obs)
    assert 1 == obj.count_neighbors("TA", obs)


def test_rlnk(obj):
    """primer_concentration: nM"""
    obs_seq = "AAAGGGGCCCGGAAGGGGACGGTTG"
    obj.set_sequence(obs_seq)
    pred = 33.404
    assert pred == obj.cal_RlnK(primer_concentration=50)
    pred = 34.419
    assert pred == obj.cal_RlnK(primer_concentration=30)


def test_deltag(obj):
    obs_seq = "AAAGGGGCCCGGAAGGGGACGGTTG"
    obj.set_sequence(obs_seq)
    pred = 39.4
    assert pred == obj.DeltaG()


def test_deltah(obj):
    obs_seq = "AAAGGGGCCCGGAAGGGGACGGTTG"
    obj.set_sequence(obs_seq)
    pred = 232.9
    assert pred == obj.DeltaH()


def test_deltas(obj):
    obs_seq = "AAAGGGGCCCGGAAGGGGACGGTTG"
    obj.set_sequence(obs_seq)
    pred = 607.8
    assert pred == obj.DeltaS()


def test_nn_tm(obj):
    """refer to http://biotools.nubic.northwestern.edu/OligoCalc.html
    I've check the concordance of the results of the website and the equation below
    """
    pred = 63.42
    obs_seq = "AAAGGGGCCCGGAAGGGGACGGTTG"
    obj.set_sequence(obs_seq)
    assert pred == obj.NearestNeighborTM(primer_concentration=50, salt_concentration=50)
