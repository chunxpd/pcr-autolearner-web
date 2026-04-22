import pytest
import os


from src.hsp_pair_score import Score


@pytest.fixture(scope="module")
def obj():
    res = Score()
    yield res
    del res


def test_preprocessing(obj):
    obs = [
        {
            "E": 0.000183945,
            "length": 6993871,
            "title": (
                "gi|1858710525|gb|CP054795.1| Mycolicibacterium smegmatis strain "
                "FDAARGOS_679 chromosome"
            ),
            "query_len": 21,
            "aln_len": 21,
            "identity": 21,
            "aln_pct": 100.0,
            "id_pct": 100.0,
        },
        {
            "E": 0.000183945,
            "length": 6993871,
            "title": (
                "gi|1858710525|gb|CP054795.1| Mycolicibacterium smegmatis strain "
                "FDAARGOS_679 chromosome"
            ),
            "query_len": 21,
            "aln_len": 21,
            "identity": 21,
            "aln_pct": 100.0,
            "id_pct": 100.0,
        },
    ]
    pred = {
        (
            "gi|1858710525|gb|CP054795.1| Mycolicibacterium smegmatis strain "
            "FDAARGOS_679 chromosome"
        ): {
            "E": 0.000183945,
            "length": 6993871,
            "title": (
                "gi|1858710525|gb|CP054795.1| Mycolicibacterium smegmatis strain "
                "FDAARGOS_679 chromosome"
            ),
            "query_len": 21,
            "aln_len": 21,
            "identity": 21,
            "aln_pct": 100.0,
            "id_pct": 100.0,
        }
    }
    assert pred == obj.preprocessing(obs)


def test_preprocessing2(obj):
    """인위적으로 E 값들을 바꾸어보았음. 제일 처음 것만 나올 것"""
    obs = [
        {
            "E": 0.001,
            "length": 6993871,
            "title": (
                "gi|1858710525|gb|CP054795.1| Mycolicibacterium smegmatis strain "
                "FDAARGOS_679 chromosome"
            ),
            "query_len": 21,
            "aln_len": 21,
            "identity": 21,
            "aln_pct": 100.0,
            "id_pct": 100.0,
        },
        {
            "E": 0.002,
            "length": 6993871,
            "title": (
                "gi|1858710525|gb|CP054795.1| Mycolicibacterium smegmatis strain "
                "FDAARGOS_679 chromosome"
            ),
            "query_len": 21,
            "aln_len": 21,
            "identity": 21,
            "aln_pct": 100.0,
            "id_pct": 100.0,
        },
    ]
    pred = {
        (
            "gi|1858710525|gb|CP054795.1| Mycolicibacterium smegmatis strain "
            "FDAARGOS_679 chromosome"
        ): {
            "E": 0.001,
            "length": 6993871,
            "title": (
                "gi|1858710525|gb|CP054795.1| Mycolicibacterium smegmatis strain "
                "FDAARGOS_679 chromosome"
            ),
            "query_len": 21,
            "aln_len": 21,
            "identity": 21,
            "aln_pct": 100.0,
            "id_pct": 100.0,
        }
    }
    assert pred == obj.preprocessing(obs)


def test_preprocessing3(obj):
    """인위적 예: 두번째 것이 반환되어 나올 것"""
    obs = [
        {
            "E": 0.003,
            "length": 6993871,
            "title": (
                "gi|1858710525|gb|CP054795.1| Mycolicibacterium smegmatis strain "
                "FDAARGOS_679 chromosome"
            ),
            "query_len": 21,
            "aln_len": 20,
            "identity": 20,
            "aln_pct": 20 / 21,
            "id_pct": 20 / 21,
        },
        {
            "E": 0.002,
            "length": 6993871,
            "title": (
                "gi|1858710525|gb|CP054795.1| Mycolicibacterium smegmatis strain "
                "FDAARGOS_679 chromosome"
            ),
            "query_len": 21,
            "aln_len": 21,
            "identity": 21,
            "aln_pct": 100.0,
            "id_pct": 100.0,
        },
    ]
    pred = {
        (
            "gi|1858710525|gb|CP054795.1| Mycolicibacterium smegmatis strain "
            "FDAARGOS_679 chromosome"
        ): {
            "E": 0.002,
            "length": 6993871,
            "title": (
                "gi|1858710525|gb|CP054795.1| Mycolicibacterium smegmatis strain "
                "FDAARGOS_679 chromosome"
            ),
            "query_len": 21,
            "aln_len": 21,
            "identity": 21,
            "aln_pct": 100.0,
            "id_pct": 100.0,
        }
    }
    assert pred == obj.preprocessing(obs)


def test_combine_hsp_sets(obj):
    obs1 = {
        "id1": {"aln_pct": 100.0, "id_pct": 100.0, "E": 0.001},
        "id2": {"aln_pct": 100.0, "id_pct": 100.0, "E": 0.001},
    }
    obs2 = {
        "id1": {"aln_pct": 100.0, "id_pct": 100.0, "E": 0.001},
        "id3": {"aln_pct": 100.0, "id_pct": 100.0, "E": 0.001},
    }
    obs3 = {
        "id2": {"aln_pct": 100.0, "id_pct": 100.0, "E": 0.001},
        "id3": {"aln_pct": 100.0, "id_pct": 100.0, "E": 0.001},
    }
    pred = {
        "id1": [
            {"aln_pct": 100.0, "id_pct": 100.0, "E": 0.001},
            {"aln_pct": 100.0, "id_pct": 100.0, "E": 0.001},
            {},
        ],
        "id2": [
            {"aln_pct": 100.0, "id_pct": 100.0, "E": 0.001},
            {},
            {"aln_pct": 100.0, "id_pct": 100.0, "E": 0.001},
        ],
        "id3": [
            {},
            {"aln_pct": 100.0, "id_pct": 100.0, "E": 0.001},
            {"aln_pct": 100.0, "id_pct": 100.0, "E": 0.001},
        ],
    }
    assert pred == obj.combine_hsp_sets(obs1, obs2, obs3)
    assert pred == obj.combine_hsp_sets(*[obs1, obs2, obs3])


def test_mk_csv_file(obj):
    """gb id 에 , 가 있는 경우가 있다. csv 의 deliminator 를 tab 으로 바꾼다.
    예) gi|1320742873|gb|CP018897.1| Bordetella parapertussis isolate H889 chromosome, complete genome
    """
    obs = {
        "id1": [
            {"aln_pct": 100.0, "id_pct": 100.0, "E": 0.001},
            {"aln_pct": 99.0, "id_pct": 100.0, "E": 0.001},
            {},
        ],
        "id2": [
            {"aln_pct": 100.0, "id_pct": 90.0, "E": 0.001},
            {},
            {"aln_pct": 100.0, "id_pct": 100.0, "E": 0.001},
        ],
        "id3": [
            {},
            {"aln_pct": 100.0, "id_pct": 100.0, "E": 0.001},
            {"aln_pct": 100.0, "id_pct": 80.0, "E": 0.001},
        ],
    }
    pred = (
        "id\taln_pct1\tid_pct1\tE1\taln_pct2\tid_pct2\tE2\taln_pct3\tid_pct3\tE3",
        "id1\t100.0\t100.0\t0.001\t99.0\t100.0\t0.001\t\t\t",
        "id2\t100.0\t90.0\t0.001\t\t\t\t100.0\t100.0\t0.001",
        "id3\t\t\t\t100.0\t100.0\t0.001\t100.0\t80.0\t0.001",
    )
    assert "\n".join(pred) == obj.mk_csv_file(obs, delim="\t")
    pred = (
        "id,aln_pct1,id_pct1,E1,aln_pct2,id_pct2,E2,aln_pct3,id_pct3,E3",
        "id1,100.0,100.0,0.001,99.0,100.0,0.001,,,",
        "id2,100.0,90.0,0.001,,,,100.0,100.0,0.001",
        "id3,,,,100.0,100.0,0.001,100.0,80.0,0.001",
    )
    assert "\n".join(pred) == obj.mk_csv_file(obs, delim=",")
