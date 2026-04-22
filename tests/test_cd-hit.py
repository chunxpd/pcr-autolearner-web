import pytest
from src.cd_hit import CDHIT


@pytest.fixture(scope="module")
def obj():
    res = CDHIT()
    yield res
    del res


def test_run_cdhit(obj):
    pred = [
        "tests/cd_hit_example.fa",
        "tests/cd_hit_example_100.fa.clstr",
        "tests/cd_hit_example_100.fa",
        "tests/cd_hit_example_lb.fa.clstr",
        "tests/cd_hit_example_lb.fa",
    ]
    obs = obj.run_cdhit("tests/cd_hit_example.fa")
    assert sorted(pred) == sorted(obs)


def test_sub_parse_lb_clstr(obj):
    pred = 3
    obs = obj.sub_parse_lb_clstr("tests/cd_hit_example_lb.fa.clstr")
    assert pred == len(obs)


def test_parse_lb_clstr(obj):
    pred = ["seq0:This", "seq30:This"]
    obs = obj.parse_lb_clstr("tests/cd_hit_example_lb.fa.clstr")
    assert pred == obs


def test_read_fa_file(obj):
    pred = [
        "seq99:This",
        "seq98:This",
        "seq70:This",
        "seq50:This",
        "seq30:This",
        "seq0:This",
        "seq100-1:This",
    ]
    obs = obj.read_fa_file("tests/cd_hit_example_100.fa")
    assert pred == obs


def test_exclude_fa_file(obj):
    pred = ["seq99:This", "seq98:This", "seq70:This", "seq50:This", "seq100-1:This"]
    obs_recs = [
        "seq99:This",
        "seq98:This",
        "seq70:This",
        "seq50:This",
        "seq30:This",
        "seq0:This",
        "seq100-1:This",
    ]
    obs_exclist = ["seq0", "seq30"]
    obs = obj.exclude_fa_file(obs_recs, obs_exclist)
    assert pred == obs


def test_write_fa_file(obj):
    pred = 5
    obs = obj.write_fa_file(
        "tests/cd_hit_example100_lb_exclude.fa",
        in_fa="tests/cd_hit_example_100.fa",
        exc_list=["seq0", "seq30"],
    )
    assert pred == obs


def test_fa_replace(obj):
    pred = "hello_lb.fa"
    obs_fa = "hello.fa"
    assert pred == obj._fa_replace(obs_fa, "_lb.fa")
    pred = "hello.fa_lb.fa"
    obs_fa = "hello.fa.fa"
    assert pred == obj._fa_replace(obs_fa, "_lb.fa")
