import pytest
from src.oligoanalyzer import OligoAnalyzer


@pytest.fixture(scope="module")
def obj():
    res = OligoAnalyzer()
    res.set_workbook("tests/OligoInput/testinput_all.xlsx")
    yield res
    del res


def test_get_ws_primer_range(obj):
    pred = ("A3:F23", "A26:B32")
    obj.set_sheet("candida")
    obs = obj._read_sheet_get_primer_range()
    assert pred == obs
    pred = ("A3:F5", "A8:B8")
    obj.set_sheet("pan candida")
    obs = obj._read_sheet_get_primer_range()
    assert pred == obs


def test_read_input(obj):
    pred = {
        "pan candida": {
            "forward": {"seq": "CGCGCTACACTGACGGA", "length": 17, "GC%": 64.7, "Tm": 57.7},
            "probe": {"seq": "TGCGCTTACTAGGAATTCCTCGTTGAAG", "length": 28, "GC%": 46.4, "Tm": 60.2},
            "reverse": {"seq": "CAGGGACGTAATCAACGCA", "length": 19, "GC%": 52.6, "Tm": 55.1},
            "template": "CGCGCTACACTGACGGAGCCAGCGAGTATTAACCTTGGCCGAGAGGTCTGGGAAATCTTGTGAAACTCCGTCGTGCTGGGGATAGAGCATTGTAATTATTGCTCTTCAACGAGGAATTCCTAGTAAGCGCAAGTCATCAGCTTGCGTTGATTACGTCCCTG",
        }
    }
    obj.set_workbook("tests/OligoInput/testinput_all.xlsx")
    assert pred == obj.read_sheet("pan candida")


def test_read_input2(obj):
    pred = {
        "C. auris": {
            "forward": {"GC%": 23.1, "Tm": 49.8, "length": 26, "seq": "ATCATTATTGATATTTTGCATACACA"},
            "probe": {
                "GC%": 31.0,
                "Tm": 55.1,
                "length": 29,
                "seq": "GATTTGGATTTTAAAACTAACCCAACGTT",
            },
            "reverse": {"GC%": 23.1, "Tm": 48.9, "length": 26, "seq": "GTTTTATGTTTTTGTTTTAGTTGAAC"},
            "template": "ATCATTATTGATATTTTGCATACACACTGATTTGGATTTTAAAACTAACCCAACGTTAAGTTCAACTAAAACAAAAACATAAAAC",
        },
        "C. glabrata": {
            "forward": {"GC%": 60.0, "Tm": 57.9, "length": 20, "seq": "TCTGGGGAGGACCAGTGTAG"},
            "probe": {"GC%": 52.0, "Tm": 61.4, "length": 25, "seq": "GAATGCTATTTCTCCTGCCTGCGCT"},
            "reverse": {"GC%": 55.0, "Tm": 56.9, "length": 20, "seq": "ACACTCCCAGGTCTTTGTCG"},
            "template": "TCTGGGGAGGACCAGTGTAGACACTCAGGAGGTTCCTAAAATATTTTCTCTGCTGTGAATGCTATTTCTCCTGCCTGCGCTTAAGTGCGCGGTTGGTGGGTGTTCTGCAGTGGGGGGAGGGAGCCGACAAAGACCTGGGAGTGT",
        },
        "C. guilliemondii": {
            "forward": {"GC%": 58.8, "Tm": 56.2, "length": 17, "seq": "CTTAACTGCGCGGCGAA"},
            "probe": {"GC%": 39.3, "Tm": 56.3, "length": 28, "seq": "CAGTGTCTTTTTGATACAGAACTCTTGC"},
            "reverse": {"GC%": 57.9, "Tm": 54.6, "length": 19, "seq": "GCCCAACCTATCTCTAGGC"},
            "template": "CTTAACTGCGCGGCGAAAAACCTTACACACAGTGTCTTTTTGATACAGAACTCTTGCTTTGGTTTGGCCTAGAGATAGGTTGGGC",
        },
        "C. krusei": {
            "forward": {"GC%": 34.8, "Tm": 50.2, "length": 23, "seq": "GGAAGGATCATTACTGTGATTTA"},
            "probe": {"GC%": 57.1, "Tm": 61.0, "length": 21, "seq": "TACACTGCGTGAGCGGAACGA"},
            "reverse": {"GC%": 32.0, "Tm": 52.7, "length": 25, "seq": "TGCTATATTCCACATTTTAGGTGTT"},
            "template": "GGAAGGATCATTACTGTGATTTAGTACTACACTGCGTGAGCGGAACGAAAACAAAAACACCTAAAATGTGGAATATAGCA",
        },
        "c. albicans": {
            "forward": {"GC%": 39.1, "Tm": 54.9, "length": 23, "seq": "TGCACCACATGTGTTTTTCTTTG"},
            "probe": {"GC%": 82.4, "Tm": 65.8, "length": 17, "seq": "CGGTGGGCCCAGCCTGC"},
            "reverse": {"GC%": 23.1, "Tm": 48.8, "length": 26, "seq": "TGATAAAAAATTGGTTGTAAGTTTAG"},
            "template": "TGCACCACATGTGTTTTTCTTTGAAACAAACTTGCTTTGGCGGTGGGCCCAGCCTGCCGCCAGAGGTCTAAACTTACAACCAATTTTTTATCA",
        },
        "c. parapsilosis": {
            "forward": {"GC%": 33.3, "Tm": 51.4, "length": 24, "seq": "GTGCTTAACTGCATTTTTTCTTAC"},
            "probe": {"GC%": 56.5, "Tm": 60.5, "length": 23, "seq": "TAGGCCTTCTATATGGGGCCTGC"},
            "reverse": {
                "GC%": 22.2,
                "Tm": 50.7,
                "length": 27,
                "seq": "TGACATTAAATAAAATTTGGTTGAGTT",
            },
            "template": "GTGCTTAACTGCATTTTTTCTTACACATGTGTTTTTCTTTTTTTGAAAACTTTGCTTTGGTAGGCCTTCTATATGGGGCCTGCCAGAGATTAAACTCAACCAAATTTTATTTAATGTCA",
        },
        "c. tropicalis": {
            "forward": {"GC%": 34.8, "Tm": 52.8, "length": 23, "seq": "TGCACCACATGTGTTTTTTATTG"},
            "probe": {"GC%": 60.0, "Tm": 63.7, "length": 25, "seq": "CTTTGGTGGCGGGAGCAATCCTACC"},
            "reverse": {"GC%": 25.0, "Tm": 48.4, "length": 24, "seq": "TTTGACTGTAAATAAAAAGTTTGG"},
            "template": "TGCACCACATGTGTTTTTTATTGAACAAATTTCTTTGGTGGCGGGAGCAATCCTACCGCCAGAGGTTATAACTAAACCAAACTTTTTATTTACAGTCAAA",
        },
    }
    obj.set_workbook("tests/OligoInput/testinput_all.xlsx")
    assert pred == obj.read_sheet("candida")


def test_read_input_to_list(obj):
    pred = [
        ["pan candida", "Forward", "CGCGCTACACTGACGGA", 17, 64.7, 57.7],
        ["pan candida", "probe", "TGCGCTTACTAGGAATTCCTCGTTGAAG", 28, 46.4, 60.2],
        ["pan candida", "reverse", "CAGGGACGTAATCAACGCA", 19, 52.6, 55.1],
        [
            "pan candida",
            "template",
            "CGCGCTACACTGACGGAGCCAGCGAGTATTAACCTTGGCCGAGAGGTCTGGGAAATCTTGTGAAACTCCGTCGTGCTGGGGATAGAGCATTGTAATTATTGCTCTTCAACGAGGAATTCCTAGTAAGCGCAAGTCATCAGCTTGCGTTGATTACGTCCCTG",
        ],
    ]
    obs = {
        "pan candida": {
            "Forward": {"seq": "CGCGCTACACTGACGGA", "length": 17, "GC%": 64.7, "Tm": 57.7},
            "probe": {"seq": "TGCGCTTACTAGGAATTCCTCGTTGAAG", "length": 28, "GC%": 46.4, "Tm": 60.2},
            "reverse": {"seq": "CAGGGACGTAATCAACGCA", "length": 19, "GC%": 52.6, "Tm": 55.1},
            "template": "CGCGCTACACTGACGGAGCCAGCGAGTATTAACCTTGGCCGAGAGGTCTGGGAAATCTTGTGAAACTCCGTCGTGCTGGGGATAGAGCATTGTAATTATTGCTCTTCAACGAGGAATTCCTAGTAAGCGCAAGTCATCAGCTTGCGTTGATTACGTCCCTG",
        }
    }
    assert pred == obj.read_sheet_to_list(obs)


def test_step1_tm_filter_1(obj):
    """1단계:
    입력된 primer/probe를 기반으로 melting temperature 필터.
    * primer/probe 의 길이를 늘렸다 줄였다, 위치를 이동하면서 Tm 을 계산하고 사전 정의된 범위안의 primer/probe 를 남긴다.
        fw/rv primer 의 경우: min(25bp, len(primer)+5bp) 까지만 늘린다.
        primer 당 입력 길이 + 5bp 까지 5’과 3’으로 늘리고 최대 25bp까지 늘리면서
        단 Nbp 를 늘렸을 경우는 5' 과 3' 방향의 합이 Nbp 이다. 제각각 Nbp 가 아님. 대략 아래와 같은 코드이다. 일종의 convolution 임 (아래 pseudo code 참조)
    * probe의 경우는 최대 길이 30bp까지
    * 해당 부위의 amplicon 에 degenerate base 가 있으면 가능한 모든 base 를 반영하여 만든다
    """
    obs = {
        "forward": {
            "GC%": 34.8,
            "Tm": 52.8,
            "length": 23,
            "seq": "TGCACCACATGTGTTTTTTATTG",
        },  # 23bp
        "probe": {
            "GC%": 60.0,
            "Tm": 63.7,
            "length": 25,
            "seq": "CTTTGGTGGCGGGAGCAATCCTACC",
        },  # 25bp
        "reverse": {
            "GC%": 25.0,
            "Tm": 48.4,
            "length": 24,
            "seq": "TTTGACTGTAAATAAAAAGTTTGG",
        },  # 21bp
        "template": "TGCACCACATGTGTTTTTTATTGAACAAATTTCTTTGGTGGCGGGAGCAATCCTACCGCCAGAGGTTATAACTAAACCAAACTTTTTATTTACAGTCAAA",
    }
    assert 25 == obj.step1_1_get_max_length(obs, "forward")
    assert 25 == obj.step1_1_get_max_length(obs, "reverse")
    assert 30 == obj.step1_1_get_max_length(obs, "probe")
    obs = {
        "forward": {
            "GC%": 34.8,
            "Tm": 52.8,
            "length": 23,
            "seq": "TGCACCACATGTGT",  # 14bp
        },  # 23bp
        "probe": {
            "GC%": 60.0,
            "Tm": 63.7,
            "length": 25,
            "seq": "CTTTGGTGGCGGGAGCAATCCTACC",
        },  # 25bp
        "reverse": {
            "GC%": 25.0,
            "Tm": 48.4,
            "length": 24,
            "seq": "TTTGACTGTAAATAAAAAGTTTGG",
        },  # 21bp
        "template": "TGCACCACATGTGTTTTTTATTGAACAAATTTCTTTGGTGGCGGGAGCAATCCTACCGCCAGAGGTTATAACTAAACCAAACTTTTTATTTACAGTCAAA",
    }
    assert 19 == obj.step1_1_get_max_length(obs, "forward")


def test_step1_tm_filter_2_0(obj):
    """primer 가 template 에서 어디에서 matching 되는지 찾는다.
    gap 을 허용하지 않는다. degenerate base 도 반영한다.
    """
    pred = 0
    obs_primer = "TGCACCACATGTGTTTTTTATTG"
    obs_template = "TGCACCACATGTGTTTTTTATTGAACAAATTTCTTTGGTGGCGGGAGCAATCCTACCGCCAGAGGTTATAACTAAACCAAACTTTTTATTTACAGTCAAA"
    assert pred == obj.step1_2_find_align_index(obs_primer, obs_template)
    obs_primer = "TGCWCCACATGTGTTTTTTATTG"  # 4bp 위치에 W 를 도입
    obs_template = "TGCACCACATGTGTTTTTTATTGAACAAATTTCTTTGGTGGCGGGAGCAATCCTACCGCCAGAGGTTATAACTAAACCAAACTTTTTATTTACAGTCAAA"
    assert pred == obj.step1_2_find_align_index(obs_primer, obs_template)
    obs_primer = "GCWCCACATGTGTTTTTTATTG"  # 위 서열에서 5' 첫번째 bp 제거.
    obs_template = "TGCACCACATGTGTTTTTTATTGAACAAATTTCTTTGGTGGCGGGAGCAATCCTACCGCCAGAGGTTATAACTAAACCAAACTTTTTATTTACAGTCAAA"
    assert 1 == obj.step1_2_find_align_index(obs_primer, obs_template)
    pred = 0
    obs_primer = "TGCACCACATGTGTTTTTTATTGAA"
    obs_template = "TGCACCACATGTGTTTTTTATTGAACAAATTTCTTTGGTGGCGGGAGCAATCCTACCGCCAGAGGTTATAACTAAACCAAACTTTTTATTTACAGTCAAA"
    assert pred == obj.step1_2_find_align_index(obs_primer, obs_template)


def test_step1_tm_filter_2(obj):
    """step1.2 : template sequence 에서 max len 만큼 하나씩 5' 또는 3' 방향으로 서열을 늘려간다.
    reverse primer 의 경우에는 template 에서 가져오기 때문에 입력 서열만 reverse_complement 하고 template 는 그대로 가지고 온다.
    """
    obs = {
        "forward": {
            "GC%": 34.8,
            "Tm": 52.8,
            "length": 23,
            "seq": "TGCACCACATGTGTTTTTTATTG",
        },  # 23bp
        "probe": {
            "GC%": 60.0,
            "Tm": 63.7,
            "length": 25,
            "seq": "CTTTGGTGGCGGGAGCAATCCTACC",
        },  # 25bp
        "reverse": {
            "GC%": 25.0,
            "Tm": 48.4,
            "length": 24,
            "seq": "TTTGACTGTAAATAAAAAGTTTGG",
        },  # 24bp, reverse complement
        "template": "TGCACCACATGTGTTTTTTATTGAACAAATTTCTTTGGTGGCGGGAGCAATCCTACCGCCAGAGGTTATAACTAAACCAAACTTTTTATTTACAGTCAAA",
    }
    pred = [
        "TGCACCACATGTGTTTTTTATTG",
        "TGCACCACATGTGTTTTTTATTGA",
        "TGCACCACATGTGTTTTTTATTGAA",
    ]
    assert set(pred) == set(obj.step1_2_extend_seq(obs, "forward"))
    pred = [
        "CCAAACTTTTTATTTACAGTCAAA",
        "ACCAAACTTTTTATTTACAGTCAAA",
    ]
    assert set(pred) == set(obj.step1_2_extend_seq(obs, "reverse"))
    pred = [
        "CTTTGGTGGCGGGAGCAATCCTACC",
        "CTTTGGTGGCGGGAGCAATCCTACCG",
        "CTTTGGTGGCGGGAGCAATCCTACCGC",
        "CTTTGGTGGCGGGAGCAATCCTACCGCC",
        "CTTTGGTGGCGGGAGCAATCCTACCGCCA",
        "CTTTGGTGGCGGGAGCAATCCTACCGCCAG",
        "TCTTTGGTGGCGGGAGCAATCCTACC",
        "TCTTTGGTGGCGGGAGCAATCCTACCG",
        "TCTTTGGTGGCGGGAGCAATCCTACCGC",
        "TCTTTGGTGGCGGGAGCAATCCTACCGCC",
        "TCTTTGGTGGCGGGAGCAATCCTACCGCCA",
        "TTCTTTGGTGGCGGGAGCAATCCTACC",
        "TTCTTTGGTGGCGGGAGCAATCCTACCG",
        "TTCTTTGGTGGCGGGAGCAATCCTACCGC",
        "TTCTTTGGTGGCGGGAGCAATCCTACCGCC",
        "TTTCTTTGGTGGCGGGAGCAATCCTACC",
        "TTTCTTTGGTGGCGGGAGCAATCCTACCG",
        "TTTCTTTGGTGGCGGGAGCAATCCTACCGC",
        "ATTTCTTTGGTGGCGGGAGCAATCCTACC",
        "ATTTCTTTGGTGGCGGGAGCAATCCTACCG",
        "AATTTCTTTGGTGGCGGGAGCAATCCTACC",
    ]
    assert set(pred) == set(obj.step1_2_extend_seq(obs, "probe"))


def test_step1_tm_filter_3(obj):
    """생성된 서열들에 대해서 melting temperature 계산하기
    온도 계산은 salt (Na+) adjusted  방법과 nearest neighbor 방법 두 가지가 있다.
    http://biotools.nubic.northwestern.edu/OligoCalc.html
    """
    pred = {
        "TGCACCACATGTGT": 40.83,
        "TGCACCACATGTGTT": 43.37,
        "TGCACCACATGTGTTT": 45.59,
        "TGCACCACATGTGTTTT": 47.55,
        "TGCACCACATGTGTTTTT": 49.29,
        "TGCACCACATGTGTTTTTT": 50.85,
    }
    obs = {
        "forward": {
            "GC%": 34.8,
            "Tm": 52.8,
            "length": 23,
            "seq": "TGCACCACATGTGT",  # 14bp
        },  # 23bp
        "probe": {
            "GC%": 60.0,
            "Tm": 63.7,
            "length": 25,
            "seq": "CTTTGGTGGCGGGAGCAATCCTACC",
        },  # 25bp
        "reverse": {
            "GC%": 25.0,
            "Tm": 48.4,
            "length": 24,
            "seq": "TTTGACTGTAAATAAAAAGTTTGG",
        },  # 21bp
        "template": "TGCACCACATGTGTTTTTTATTGAACAAATTTCTTTGGTGGCGGGAGCAATCCTACCGCCAGAGGTTATAACTAAACCAAACTTTTTATTTACAGTCAAA",
    }
    assert pred == obj.step1_3_cal_tm(obs, "forward", tm_method="adj")
    pred = {
        "TGCACCACATGTGT": 38.63,
        "TGCACCACATGTGTT": 40.62,
        "TGCACCACATGTGTTT": 42.38,
        "TGCACCACATGTGTTTT": 43.95,
        "TGCACCACATGTGTTTTT": 45.35,
        "TGCACCACATGTGTTTTTT": 46.62,
    }
    assert pred == obj.step1_3_cal_tm(obs, "forward", tm_method="nn")


def test_step1_tm_filter_4(obj):
    """filtering the generated sequeces by a given Tm range"""
    pred = {"TGCACCACATGTGTTT": 45.59}
    obs = obj.step1_4_filtering_tm(
        {
            "TGCACCACATGTGT": 40.83,
            "TGCACCACATGTGTT": 43.37,
            "TGCACCACATGTGTTT": 45.59,
            "TGCACCACATGTGTTTT": 47.55,
            "TGCACCACATGTGTTTTT": 49.29,
            "TGCACCACATGTGTTTTTT": 50.85,
        },
        [45, 46],
    )
    assert pred == obs


def test_step1_tm_filter_5(obj):
    """최종 결과를 반환한다.
    reverse 의 경우에는 계산 과정 중에서 5'->3' 상보염기로 변환하였으므로 다시 원래대로 3'->5' 으로 변환한다.
    """
    pred = {
        "forward": {
            "TGCACCACATGTGTTTTTTATTG": 57.51,
            "TGCACCACATGTGTTTTTTATTGA": 58.4,
            "TGCACCACATGTGTTTTTTATTGAA": 59.22,
        },
        "probe": {
            "AATTTCTTTGGTGGCGGGAGCAATCCTACC": 72.07,
            "ATTTCTTTGGTGGCGGGAGCAATCCTACC": 71.83,
            "ATTTCTTTGGTGGCGGGAGCAATCCTACCG": 73.43,
            "CTTTGGTGGCGGGAGCAATCCTACC": 70.7,
            "CTTTGGTGGCGGGAGCAATCCTACCG": 72.59,
            "CTTTGGTGGCGGGAGCAATCCTACCGC": 74.34,
            "CTTTGGTGGCGGGAGCAATCCTACCGCC": 75.97,
            "CTTTGGTGGCGGGAGCAATCCTACCGCCA": 76.07,
            "CTTTGGTGGCGGGAGCAATCCTACCGCCAG": 77.53,
            "TCTTTGGTGGCGGGAGCAATCCTACC": 71.02,
            "TCTTTGGTGGCGGGAGCAATCCTACCG": 72.83,
            "TCTTTGGTGGCGGGAGCAATCCTACCGC": 74.51,
            "TCTTTGGTGGCGGGAGCAATCCTACCGCC": 76.07,
            "TCTTTGGTGGCGGGAGCAATCCTACCGCCA": 76.17,
            "TTCTTTGGTGGCGGGAGCAATCCTACC": 71.31,
            "TTCTTTGGTGGCGGGAGCAATCCTACCG": 73.04,
            "TTCTTTGGTGGCGGGAGCAATCCTACCGC": 74.66,
            "TTCTTTGGTGGCGGGAGCAATCCTACCGCC": 76.17,
            "TTTCTTTGGTGGCGGGAGCAATCCTACC": 71.58,
            "TTTCTTTGGTGGCGGGAGCAATCCTACCG": 73.24,
            "TTTCTTTGGTGGCGGGAGCAATCCTACCGC": 74.8,
        },
        "reverse": {"TTTGACTGTAAATAAAAAGTTTGG": 54.98, "TTTGACTGTAAATAAAAAGTTTGGT": 55.94},
    }
    obs = {
        "forward": {
            "GC%": 34.8,
            "Tm": 52.8,
            "length": 23,
            "seq": "TGCACCACATGTGTTTTTTATTG",
        },  # 23bp
        "probe": {
            "GC%": 60.0,
            "Tm": 63.7,
            "length": 25,
            "seq": "CTTTGGTGGCGGGAGCAATCCTACC",
        },  # 25bp
        "reverse": {
            "GC%": 25.0,
            "Tm": 48.4,
            "length": 24,
            "seq": "TTTGACTGTAAATAAAAAGTTTGG",
        },  # 21bp
        "template": "TGCACCACATGTGTTTTTTATTGAACAAATTTCTTTGGTGGCGGGAGCAATCCTACCGCCAGAGGTTATAACTAAACCAAACTTTTTATTTACAGTCAAA",
    }
    assert pred == obj.step1_tm_filter_main(obs, [0, 100, 60, 80], tm_method="adj")

    pred = {
        "forward": {
            "TGCACCACATGTGTTTTTTATTG": 51.07,
            "TGCACCACATGTGTTTTTTATTGA": 52.33,
            "TGCACCACATGTGTTTTTTATTGAA": 53.05,
        },
        "probe": {
            "AATTTCTTTGGTGGCGGGAGCAATCCTACC": 63.18,
            "ATTTCTTTGGTGGCGGGAGCAATCCTACC": 62.94,
            "ATTTCTTTGGTGGCGGGAGCAATCCTACCG": 64.97,
            "CTTTGGTGGCGGGAGCAATCCTACC": 61.4,
            "CTTTGGTGGCGGGAGCAATCCTACCG": 63.77,
            "CTTTGGTGGCGGGAGCAATCCTACCGC": 65.38,
            "CTTTGGTGGCGGGAGCAATCCTACCGCC": 66.38,
            "CTTTGGTGGCGGGAGCAATCCTACCGCCA": 67.26,
            "CTTTGGTGGCGGGAGCAATCCTACCGCCAG": 68.2,
            "TCTTTGGTGGCGGGAGCAATCCTACC": 62.11,
            "TCTTTGGTGGCGGGAGCAATCCTACCG": 64.36,
            "TCTTTGGTGGCGGGAGCAATCCTACCGC": 65.88,
            "TCTTTGGTGGCGGGAGCAATCCTACCGCC": 66.83,
            "TCTTTGGTGGCGGGAGCAATCCTACCGCCA": 67.68,
            "TTCTTTGGTGGCGGGAGCAATCCTACC": 62.41,
            "TTCTTTGGTGGCGGGAGCAATCCTACCG": 64.57,
            "TTCTTTGGTGGCGGGAGCAATCCTACCGC": 66.04,
            "TTCTTTGGTGGCGGGAGCAATCCTACCGCC": 66.95,
            "TTTCTTTGGTGGCGGGAGCAATCCTACC": 62.68,
            "TTTCTTTGGTGGCGGGAGCAATCCTACCG": 64.77,
            "TTTCTTTGGTGGCGGGAGCAATCCTACCGC": 66.18,
        },
        "reverse": {"TTTGACTGTAAATAAAAAGTTTGG": 48.1, "TTTGACTGTAAATAAAAAGTTTGGT": 49.3},
    }
    assert pred == obj.step1_tm_filter_main(obs, [0, 100, 60, 80], tm_method="nn")


def test_step1_tm_filter_5_different_conc(obj):
    """최종 결과를 반환한다.
    reverse 의 경우에는 계산 과정 중에서 5'->3' 상보염기로 변환하였으므로 다시 원래대로 3'->5' 으로 변환한다.
    concentration 변화를 줄 수 있다.
    """
    pred = {
        "forward": {
            "TGCACCACATGTGTTTTTTATTG": 49.46,
            "TGCACCACATGTGTTTTTTATTGA": 50.72,
            "TGCACCACATGTGTTTTTTATTGAA": 51.44,
        },
        "probe": {
            "AATTTCTTTGGTGGCGGGAGCAATCCTACC": 61.57,
            "ATTTCTTTGGTGGCGGGAGCAATCCTACC": 61.33,
            "ATTTCTTTGGTGGCGGGAGCAATCCTACCG": 63.36,
            "CTTTGGTGGCGGGAGCAATCCTACC": 59.8,
            "CTTTGGTGGCGGGAGCAATCCTACCG": 62.17,
            "CTTTGGTGGCGGGAGCAATCCTACCGC": 63.77,
            "CTTTGGTGGCGGGAGCAATCCTACCGCC": 64.77,
            "CTTTGGTGGCGGGAGCAATCCTACCGCCA": 65.66,
            "CTTTGGTGGCGGGAGCAATCCTACCGCCAG": 66.59,
            "TCTTTGGTGGCGGGAGCAATCCTACC": 60.5,
            "TCTTTGGTGGCGGGAGCAATCCTACCG": 62.75,
            "TCTTTGGTGGCGGGAGCAATCCTACCGC": 64.28,
            "TCTTTGGTGGCGGGAGCAATCCTACCGCC": 65.22,
            "TCTTTGGTGGCGGGAGCAATCCTACCGCCA": 66.07,
            "TTCTTTGGTGGCGGGAGCAATCCTACC": 60.8,
            "TTCTTTGGTGGCGGGAGCAATCCTACCG": 62.96,
            "TTCTTTGGTGGCGGGAGCAATCCTACCGC": 64.43,
            "TTCTTTGGTGGCGGGAGCAATCCTACCGCC": 65.34,
            "TTTCTTTGGTGGCGGGAGCAATCCTACC": 61.08,
            "TTTCTTTGGTGGCGGGAGCAATCCTACCG": 63.16,
            "TTTCTTTGGTGGCGGGAGCAATCCTACCGC": 64.58,
        },
        "reverse": {"TTTGACTGTAAATAAAAAGTTTGG": 46.49, "TTTGACTGTAAATAAAAAGTTTGGT": 47.69},
    }
    obs = {
        "forward": {
            "GC%": 34.8,
            "Tm": 52.8,
            "length": 23,
            "seq": "TGCACCACATGTGTTTTTTATTG",
        },  # 23bp
        "probe": {
            "GC%": 60.0,
            "Tm": 63.7,
            "length": 25,
            "seq": "CTTTGGTGGCGGGAGCAATCCTACC",
        },  # 25bp
        "reverse": {
            "GC%": 25.0,
            "Tm": 48.4,
            "length": 24,
            "seq": "TTTGACTGTAAATAAAAAGTTTGG",
        },  # 21bp
        "template": "TGCACCACATGTGTTTTTTATTGAACAAATTTCTTTGGTGGCGGGAGCAATCCTACCGCCAGAGGTTATAACTAAACCAAACTTTTTATTTACAGTCAAA",
    }
    assert pred == obj.step1_tm_filter_main(
        obs, [0, 100, 0, 100], tm_method="nn", salt_concentration=40
    )
    pred = {
        "forward": {
            "TGCACCACATGTGTTTTTTATTG": 50.78,
            "TGCACCACATGTGTTTTTTATTGA": 52.04,
            "TGCACCACATGTGTTTTTTATTGAA": 52.77,
        },
        "probe": {
            "AATTTCTTTGGTGGCGGGAGCAATCCTACC": 62.96,
            "ATTTCTTTGGTGGCGGGAGCAATCCTACC": 62.7,
            "ATTTCTTTGGTGGCGGGAGCAATCCTACCG": 64.74,
            "CTTTGGTGGCGGGAGCAATCCTACC": 61.14,
            "CTTTGGTGGCGGGAGCAATCCTACCG": 63.52,
            "CTTTGGTGGCGGGAGCAATCCTACCGC": 65.13,
            "CTTTGGTGGCGGGAGCAATCCTACCGCC": 66.14,
            "CTTTGGTGGCGGGAGCAATCCTACCGCCA": 67.04,
            "CTTTGGTGGCGGGAGCAATCCTACCGCCAG": 67.97,
            "TCTTTGGTGGCGGGAGCAATCCTACC": 61.85,
            "TCTTTGGTGGCGGGAGCAATCCTACCG": 64.12,
            "TCTTTGGTGGCGGGAGCAATCCTACCGC": 65.65,
            "TCTTTGGTGGCGGGAGCAATCCTACCGCC": 66.6,
            "TCTTTGGTGGCGGGAGCAATCCTACCGCCA": 67.46,
            "TTCTTTGGTGGCGGGAGCAATCCTACC": 62.16,
            "TTCTTTGGTGGCGGGAGCAATCCTACCG": 64.33,
            "TTCTTTGGTGGCGGGAGCAATCCTACCGC": 65.81,
            "TTCTTTGGTGGCGGGAGCAATCCTACCGCC": 66.73,
            "TTTCTTTGGTGGCGGGAGCAATCCTACC": 62.45,
            "TTTCTTTGGTGGCGGGAGCAATCCTACCG": 64.54,
            "TTTCTTTGGTGGCGGGAGCAATCCTACCGC": 65.96,
        },
        "reverse": {"TTTGACTGTAAATAAAAAGTTTGG": 47.81, "TTTGACTGTAAATAAAAAGTTTGGT": 49.02},
    }
    assert pred == obj.step1_tm_filter_main(
        obs, [0, 100, 0, 100], tm_method="nn", primer_concentration=40
    )


def test_step2_hairpin_filter_1_run_mfold(obj):
    """step1 의 결과물에서부터 시작."""
    obs = {
        "forward": {
            "TGCACCACATGTGTTTTTTATTG": 49.46,
            "TGCACCACATGTGTTTTTTATTGA": 50.72,
            "TGCACCACATGTGTTTTTTATTGAA": 51.44,
        },
        "probe": {
            "AATTTCTTTGGTGGCGGGAGCAATCCTACC": 61.57,
            "ATTTCTTTGGTGGCGGGAGCAATCCTACC": 61.33,
            "ATTTCTTTGGTGGCGGGAGCAATCCTACCG": 63.36,
            "CTTTGGTGGCGGGAGCAATCCTACC": 59.8,
            "CTTTGGTGGCGGGAGCAATCCTACCG": 62.17,
            "CTTTGGTGGCGGGAGCAATCCTACCGC": 63.77,
            "CTTTGGTGGCGGGAGCAATCCTACCGCC": 64.77,
            "CTTTGGTGGCGGGAGCAATCCTACCGCCA": 65.66,
            "CTTTGGTGGCGGGAGCAATCCTACCGCCAG": 66.59,
            "TCTTTGGTGGCGGGAGCAATCCTACC": 60.5,
            "TCTTTGGTGGCGGGAGCAATCCTACCG": 62.75,
            "TCTTTGGTGGCGGGAGCAATCCTACCGC": 64.28,
            "TCTTTGGTGGCGGGAGCAATCCTACCGCC": 65.22,
            "TCTTTGGTGGCGGGAGCAATCCTACCGCCA": 66.07,
            "TTCTTTGGTGGCGGGAGCAATCCTACC": 60.8,
            "TTCTTTGGTGGCGGGAGCAATCCTACCG": 62.96,
            "TTCTTTGGTGGCGGGAGCAATCCTACCGC": 64.43,
            "TTCTTTGGTGGCGGGAGCAATCCTACCGCC": 65.34,
            "TTTCTTTGGTGGCGGGAGCAATCCTACC": 61.08,
            "TTTCTTTGGTGGCGGGAGCAATCCTACCG": 63.16,
            "TTTCTTTGGTGGCGGGAGCAATCCTACCGC": 64.58,
        },
        "reverse": {"TTTGACTGTAAATAAAAAGTTTGG": 46.49, "TTTGACTGTAAATAAAAAGTTTGGT": 47.69},
    }
    pred = {
        "forward": {
            "TGCACCACATGTGTTTTTTATTG": {
                "Tm": 49.46,
                "mfold": {
                    "1": {
                        "Tm(Celsius)": 85.0,
                        "dG": -2.92,
                        "dH": -21.8,
                        "dS": -60.87,
                        "fold": (
                            "\n"
                            " \n"
                            "T--------|   CA\n"
                            "         GCAC  \\\n"
                            "         TGTG  C\n"
                            "GTTATTTTT^   TA\n"
                            "  20        10\n"
                            " \n"
                            " \n"
                            " "
                        ),
                    },
                    "Foldingat(Celsius)": 37.0,
                    "Mg++": 0.0,
                    "Na+": 50.0,
                },
            },
            "TGCACCACATGTGTTTTTTATTGA": {"Tm": 50.72, "mfold": {}},
            "TGCACCACATGTGTTTTTTATTGAA": {"Tm": 51.44, "mfold": {}},
        }
    }
    assert (
        pred["forward"]["TGCACCACATGTGTTTTTTATTG"]
        == obj.step2_run_mfold(obs)["forward"]["TGCACCACATGTGTTTTTTATTG"]
    )


def test_step2_hairpin_filter_2_mfold_parser(obj):
    obs = (
        "<html><head><title>Folding "
        "energy details</title><META "
        'HTTP-EQUIV="Content-Type" '
        'CONTENT="text/html"; '
        'charset=UTF-8"></head>\n'
        "Linear DNA folding at 37 "
        "&#176;C. [Na<SUP>+</SUP>] = "
        "50, [Mg<SUP>++</SUP>] = "
        "0.0.\n"
        "<html><head><title>Folding "
        "energy details</title><META "
        'HTTP-EQUIV="Content-Type" '
        'CONTENT="text/html"; '
        'charset=UTF-8"></head>\n'
        '<body bgcolor="#dfc890" '
        'text="#000000" '
        'link="#8e2323" '
        'vlink="#2f4f4f" '
        'alink="#00ffff"><b>\n'
        '<a name="STRUCTURE_1"> '
        "Structure 1</a>\n"
        "Folding bases      1 to     "
        "25 of seq\n"
        " <br>\n"
        " dG =     -6.43 &nbsp; dH "
        "=    -40.20 &nbsp; dS =   "
        "-108.88 &nbsp; "
        "T<sub>m</sub> =   96.1 ℃\n"
        " \n"
        " \n"
        " <pre>\n"
        "        10\n"
        "CTTT|  GGC    G\n"
        "    GGT   GGGA C\n"
        "    CCA   TCCT A\n"
        "----^  ---    A\n"
        "           20\n"
        " \n"
        " \n"
        " </pre>\n"
        "<hr size=3 width=80% "
        "noshade><p><!--endof   1 "
        "-->\n"
        '<a name="STRUCTURE_2"> '
        "Structure 2</a>\n"
        "Folding bases      1 to     "
        "25 of seq\n"
        " <br>\n"
        " dG =     -5.73 &nbsp; dH "
        "=    -39.70 &nbsp; dS =   "
        "-109.53 &nbsp; "
        "T<sub>m</sub> =   89.3 ℃\n"
        " \n"
        " \n"
        " <pre>\n"
        "        10\n"
        "CTTT   |GCG   G\n"
        "    GGTG   GGA C\n"
        "    CCAT   CCT A\n"
        "----   ^---   A\n"
        "           20\n"
        " \n"
        " \n"
        " </pre>\n"
        "<hr size=3 width=80% "
        "noshade><p><!--endof   2 "
        "-->\n"
        "<address><code>mfold</code> "
        "version 3.6<br>M. Zuker, "
        "Rensselaer Polytechnic "
        "Institute</address></body></html>\n"
    )
    pred = {
        "Foldingat(Celsius)": 37.0,
        "Na+": 50.0,
        "Mg++": 0.0,
        "1": {
            "dG": -6.43,
            "dH": -40.20,
            "dS": -108.88,
            "Tm(Celsius)": 96.1,
            "fold": "\n        10\nCTTT|  GGC    G\n    GGT   GGGA C\n    CCA   TCCT A\n----^  ---    A\n           20\n \n \n ",
        },
        "2": {
            "dG": -5.73,
            "dH": -39.70,
            "dS": -109.53,
            "Tm(Celsius)": 89.3,
            "fold": "\n        10\nCTTT   |GCG   G\n    GGTG   GGA C\n    CCAT   CCT A\n----   ^---   A\n           20\n \n \n ",
        },
    }
    assert pred == obj.step2_mfold_parser(obs)


def test_step2_hairpin_filter_2_1_html_parser1(obj):
    pred = {"Foldingat(Celsius)": 37, "Na+": 50, "Mg++": 0.0}
    obs = "Linear DNA folding at 37 &#176;C. [Na<SUP>+</SUP>] = 50, [Mg<SUP>++</SUP>] = 0.0."
    assert pred == obj.step2_html_parser1(obs)


def test_step2_hairpin_filter_2_2_html_parser2(obj):
    pred = {"dG": -6.43, "dH": -40.20, "dS": -108.88, "Tm(Celsius)": 96.1}
    obs = " dG =     -6.43 &nbsp; dH =    -40.20 &nbsp; dS =   -108.88 &nbsp; T<sub>m</sub> =   96.1 ℃"
    assert pred == obj.step2_html_parser2(obs)


def test_step2_homodimer_get_reverseseq(obj):
    pred = "ATGCTTT"
    obs = "TTTCGTA"
    assert pred == obj.step2_get_reverseseq(obs)


def test_step2_homodimer_get_selfmaxmatching(obj):
    """
       GTCAGCTAGTAGCTAGCGCA
       |  ||||| |||||  |
    ACGCGATCGATGATCGACTG
    """
    pred = ["AGCTA", "TAGCT"]
    obs = "GTCAGCTAGTAGCTAGCGCA"
    assert pred == obj.step2_get_selfmaxmatching(obs)


def test_step2_homodimer_cal_tm(obj):
    """oligo calc 사이트에서는 7bp 이상일 때만 NN 방법으로 계산이 가능하다. deltaG, deltaS, deltaH 도 마찬가지
    일단은 넘어간다.
    어쨌든 현재 구현된 것 중에서 가장 Tm 높은 값과 deltaG 가장 낮은 값을 출력한다"""
    obs = ["AGCTA", "TAGCT"]
    pred = (14.0, -1.2)
    assert pred == obj.step2_homodimer_cal_tm(obs)


def test_step2_hoomodimer_main(obj):
    """input 은 step2_run_mfold 의 결과물"""
    obs = {
        "forward": {
            "TGCACCACATGTGTTTTTTATTG": {
                "Tm": 49.46,
                "mfold": {
                    "1": {
                        "Tm(Celsius)": 85.0,
                        "dG": -2.92,
                        "dH": -21.8,
                        "dS": -60.87,
                        "fold": (
                            "\n"
                            " \n"
                            "T--------|   CA\n"
                            "         GCAC  \\\n"
                            "         TGTG  C\n"
                            "GTTATTTTT^   TA\n"
                            "  20        10\n"
                            " \n"
                            " \n"
                            " "
                        ),
                    },
                    "Foldingat(Celsius)": 37.0,
                    "Mg++": 0.0,
                    "Na+": 50.0,
                },
            },
            "TGCACCACATGTGTTTTTTATTGA": {"Tm": 50.72, "mfold": {}},
            "TGCACCACATGTGTTTTTTATTGAA": {"Tm": 51.44, "mfold": {}},
        }
    }
    pred = {
        "forward": {
            "TGCACCACATGTGTTTTTTATTG": {
                "Tm": 49.46,
                "homodimer": {"Tm": 24.0, "deltaG": -5.7},  # obs 결과에 이 부분이 더해 짐
                "mfold": {
                    "1": {
                        "Tm(Celsius)": 85.0,
                        "dG": -2.92,
                        "dH": -21.8,
                        "dS": -60.87,
                        "fold": (
                            "\n"
                            " \n"
                            "T--------|   CA\n"
                            "         GCAC  \\\n"
                            "         TGTG  C\n"
                            "GTTATTTTT^   TA\n"
                            "  20        10\n"
                            " \n"
                            " \n"
                            " "
                        ),
                    },
                    "Foldingat(Celsius)": 37.0,
                    "Mg++": 0.0,
                    "Na+": 50.0,
                },
            },
            "TGCACCACATGTGTTTTTTATTGA": {
                "Tm": 50.72,
                "mfold": {},
                "homodimer": {"Tm": 24.0, "deltaG": -5.7},
            },
            "TGCACCACATGTGTTTTTTATTGAA": {
                "Tm": 51.44,
                "mfold": {},
                "homodimer": {"Tm": 24.0, "deltaG": -5.7},
            },
        }
    }
    assert pred == obj.step2_homodimer_main(obs)


def test_step2_selection_done(obj):
    """step2_homodimer_main 의 결과물에서 남기고 싶은 것만 남기는 부분

    Args:
        obj (_type_): _description_
    """
    obs = {
        "forward": {
            "TGCACCACATGTGTTTTTTATTG": {
                "Tm": 49.46,
                "homodimer": {"Tm": 24.0, "deltaG": -5.7},  # obs 결과에 이 부분이 더해 짐
                "mfold": {
                    "1": {
                        "Tm(Celsius)": 85.0,
                        "dG": -2.92,
                        "dH": -21.8,
                        "dS": -60.87,
                        "fold": (
                            "\n"
                            " \n"
                            "T--------|   CA\n"
                            "         GCAC  \\\n"
                            "         TGTG  C\n"
                            "GTTATTTTT^   TA\n"
                            "  20        10\n"
                            " \n"
                            " \n"
                            " "
                        ),
                    },
                    "Foldingat(Celsius)": 37.0,
                    "Mg++": 0.0,
                    "Na+": 50.0,
                },
            },
            "TGCACCACATGTGTTTTTTATTGA": {
                "Tm": 50.72,
                "mfold": {},
                "homodimer": {"Tm": 24.0, "deltaG": -5.7},
            },
            "TGCACCACATGTGTTTTTTATTGAA": {
                "Tm": 51.44,
                "mfold": {},
                "homodimer": {"Tm": 24.0, "deltaG": -5.7},
            },
        }
    }
    pred = {
        "forward": {
            "TGCACCACATGTGTTTTTTATTG": {
                "Tm": 49.46,
                "homodimer": {"Tm": 24.0, "deltaG": -5.7},  # obs 결과에 이 부분이 더해 짐
                "mfold": {
                    "1": {
                        "Tm(Celsius)": 85.0,
                        "dG": -2.92,
                        "dH": -21.8,
                        "dS": -60.87,
                        "fold": (
                            "\n"
                            " \n"
                            "T--------|   CA\n"
                            "         GCAC  \\\n"
                            "         TGTG  C\n"
                            "GTTATTTTT^   TA\n"
                            "  20        10\n"
                            " \n"
                            " \n"
                            " "
                        ),
                    },
                    "Foldingat(Celsius)": 37.0,
                    "Mg++": 0.0,
                    "Na+": 50.0,
                },
            }
        }
    }
    obs_selection_dict = {"forward": ["TGCACCACATGTGTTTTTTATTG"]}
    assert "TGCACCACATGTGTTTTTTATTG" in pred["forward"]
    assert "TGCACCACATGTGTTTTTTATTGA" not in pred["forward"]
    assert "TGCACCACATGTGTTTTTTATTGAA" not in pred["forward"]
    assert pred == obj.step2_selection_done(obs, obs_selection_dict)


def test_step3_heterodimer(obj):
    """2단계의 output 에서부터 시작한다.
    2단계 output: step2_homodimer_main()
    {
        "forward": {
            seq1: {},
            seq2: {}, ...
        },
        "reverse": {},
        "probe": {}
    }
    각 seq 는 다음의 데이터를 가진다.
    seq: {
        "Tm": float, #step1 결과
        "homodimer": {"Tm": float, "deltaG": float}, #step2 homodimer 연산 결과
        "mfold": { #step2 hairpin 계산 결과
            "1": {
                "Tm(Celsius): float,
                "dG": float,
                "dH": float,
                "dS": float,
                "fold": str
            }, ...,
            "Foldingat(Celsius)": float,
            "Mg++": float,
            "Na+": flaot
        },
        }
    }
    step3 에서는 일단 입력이 들어온 seq 들에 대해서 선택한 서열만 이용하고
    forward, reverse, probe 의 모든 서열을 대상으로 하여서 (homodimer 를 제외) heteroprimer 를 계산한다.
    계산 방법은 homodimer 와 동일 하다. 다만 output 의 구조체가 input 을 유지할 수는 없다. 따라서 output 은 별도로 받는 것으로 한다.
    """
    obs = {
        "forward": {
            "TGCACCACATGTGTTTTTTATTG": {
                "Tm": 50.0,
                "homodimer": {"Tm": 18.0, "deltaG": -10.0},
                "mfold": {"1": {}, "Foldingat(Celsius)": 36.0, "Mg++": 0.0, "Na+": 50},
            },
        },
        "probe": {
            "AATTTCTTTGGTGGCGGGAGCAATCCTACC": {
                "Tm": 61.57,
                "homodimer": {},
                "mfold": {},
            },
        },
        "reverse": {
            "TTTGACTGTAAATAAAAAGTTTGG": {"Tm": 46.49, "homodimer": {}, "mfold": {}},
        },
    }
    pred = {
        ("forward", "TGCACCACATGTGTTTTTTATTG", "probe", "AATTTCTTTGGTGGCGGGAGCAATCCTACC"): {
            "Tm": 16.0,
            "deltaG": -1.8,
        },
        ("forward", "TGCACCACATGTGTTTTTTATTG", "reverse", "TTTGACTGTAAATAAAAAGTTTGG"): {
            "Tm": 8.0,
            "deltaG": 1.4,
        },
        ("probe", "AATTTCTTTGGTGGCGGGAGCAATCCTACC", "reverse", "TTTGACTGTAAATAAAAAGTTTGG"): {
            "Tm": 10.0,
            "deltaG": 1.1,
        },
    }
    assert pred == obj.step3_heterodimer_main(obs)


def test_step3_heterodimer_all_strains_1(obj):
    """이 함수는 입력을 한 개의 strain 에 대한 것을 하는 것이 아니라 모든 strain 에 대해서 계산하는 것이다.
    다음의 요구사항이 추가된다.

    * primer-primer binding 에서 end biding 갯수 3' end 의 binding max 값, total binding, max binding 갯수
    * end binding 에서 GC 갯수
    * probe-probe binding 은 신경 안 써도 됨

    {
        "strain1": {
            "forward": {
                seq1: {},
                seq2: {}, ...
            },
            "reverse": {},
            "probe": {}
        },
        "strain2": {...},
        "strain3": {...},...
    }
    """
    obs = {
        "strain1": {
            "forward": {
                "TGCACCACATGTGTTTTTTATTG": {
                    "Tm": 50.0,
                    "homodimer": {"Tm": 18.0, "deltaG": -10.0},
                    "mfold": {"1": {}, "Foldingat(Celsius)": 36.0, "Mg++": 0.0, "Na+": 50},
                },
            },
            "probe": {
                "AATTTCTTTGGTGGCGGGAGCAATCCTACC": {
                    "Tm": 61.57,
                    "homodimer": {},
                    "mfold": {},
                },
            },
            "reverse": {
                "TTTGACTGTAAATAAAAAGTTTGG": {"Tm": 46.49, "homodimer": {}, "mfold": {}},
            },
        },
        "strain2": {
            "forward": {
                "GGCACCACATGTGTTTTTTATTG": {
                    "Tm": 50.0,
                    "homodimer": {"Tm": 18.0, "deltaG": -10.0},
                    "mfold": {"1": {}, "Foldingat(Celsius)": 36.0, "Mg++": 0.0, "Na+": 50},
                },
            },
            "probe": {
                "GATTTCTTTGGTGGCGGGAGCAATCCTACC": {
                    "Tm": 61.57,
                    "homodimer": {},
                    "mfold": {},
                },
            },
            "reverse": {
                "GTTGACTGTAAATAAAAAGTTTGG": {"Tm": 46.49, "homodimer": {}, "mfold": {}},
            },
        },
    }
    pred = {
        "forward": {
            "TGCACCACATGTGTTTTTTATTG": {
                "strain": "strain1",
                "Tm": 50.0,
                "homodimer": {"Tm": 18.0, "deltaG": -10.0},
                "mfold": {"1": {}, "Foldingat(Celsius)": 36.0, "Mg++": 0.0, "Na+": 50},
            },
            "GGCACCACATGTGTTTTTTATTG": {
                "strain": "strain2",
                "Tm": 50.0,
                "homodimer": {"Tm": 18.0, "deltaG": -10.0},
                "mfold": {"1": {}, "Foldingat(Celsius)": 36.0, "Mg++": 0.0, "Na+": 50},
            },
        },
        "probe": {
            "AATTTCTTTGGTGGCGGGAGCAATCCTACC": {
                "strain": "strain1",
                "Tm": 61.57,
                "homodimer": {},
                "mfold": {},
            },
            "GATTTCTTTGGTGGCGGGAGCAATCCTACC": {
                "strain": "strain2",
                "Tm": 61.57,
                "homodimer": {},
                "mfold": {},
            },
        },
        "reverse": {
            "TTTGACTGTAAATAAAAAGTTTGG": {
                "strain": "strain1",
                "Tm": 46.49,
                "homodimer": {},
                "mfold": {},
            },
            "GTTGACTGTAAATAAAAAGTTTGG": {
                "strain": "strain2",
                "Tm": 46.49,
                "homodimer": {},
                "mfold": {},
            },
        },
    }
    assert pred == obj.step3_combine_all_strains(obs)


def test_step3_heterodimer_all_strains_main(obj):
    obs = {
        "forward": {
            "TGCACCACATGTGTTTTTTATTG": {
                "strain": "strain1",
                "Tm": 50.0,
                "homodimer": {"Tm": 18.0, "deltaG": -10.0},
                "mfold": {"1": {}, "Foldingat(Celsius)": 36.0, "Mg++": 0.0, "Na+": 50},
            },
            "GGCACCACATGTGTTTTTTATTG": {
                "strain": "strain2",
                "Tm": 50.0,
                "homodimer": {"Tm": 18.0, "deltaG": -10.0},
                "mfold": {"1": {}, "Foldingat(Celsius)": 36.0, "Mg++": 0.0, "Na+": 50},
            },
        },
        "probe": {
            "AATTTCTTTGGTGGCGGGAGCAATCCTACC": {
                "strain": "strain1",
                "Tm": 61.57,
                "homodimer": {},
                "mfold": {},
            },
            "GATTTCTTTGGTGGCGGGAGCAATCCTACC": {
                "strain": "strain2",
                "Tm": 61.57,
                "homodimer": {},
                "mfold": {},
            },
        },
        "reverse": {
            "TTTGACTGTAAATAAAAAGTTTGG": {
                "strain": "strain1",
                "Tm": 46.49,
                "homodimer": {},
                "mfold": {},
            },
            "GTTGACTGTAAATAAAAAGTTTGG": {
                "strain": "strain2",
                "Tm": 46.49,
                "homodimer": {},
                "mfold": {},
            },
        },
    }
    pred = {
        (
            "strain1:forward",
            "TGCACCACATGTGTTTTTTATTG",
            "strain2:forward",
            "GGCACCACATGTGTTTTTTATTG",
        ): {
            "Tm": 12.0,
            "deltaG": 0.1,
            "tt_binding": 7,
            "pairwise": "        TGCACCACATGTGTTTTTTATTG\n        |.|||.|.|..|\nGGCACCACATGTGTTTTTTATTG\n",
        },
        (
            "strain1:forward",
            "TGCACCACATGTGTTTTTTATTG",
            "strain1:probe",
            "AATTTCTTTGGTGGCGGGAGCAATCCTACC",
        ): {
            "Tm": 16.0,
            "deltaG": -1.8,
            "tt_binding": 8,
            "pairwise": "        TGCACCACATGTGTTTTTTATTG\n          ||||.|....|||\nAATTTCTTTGGTGGCGGGAGCAATCCTACC\n",
        },
        (
            "strain1:forward",
            "TGCACCACATGTGTTTTTTATTG",
            "strain2:probe",
            "GATTTCTTTGGTGGCGGGAGCAATCCTACC",
        ): {
            "Tm": 16.0,
            "deltaG": -1.8,
            "tt_binding": 8,
            "pairwise": "        TGCACCACATGTGTTTTTTATTG\n          ||||.|....|||\nGATTTCTTTGGTGGCGGGAGCAATCCTACC\n",
        },
        (
            "strain1:forward",
            "TGCACCACATGTGTTTTTTATTG",
            "strain1:reverse",
            "TTTGACTGTAAATAAAAAGTTTGG",
        ): {
            "Tm": 8.0,
            "deltaG": 1.4,
            "tt_binding": 11,
            "pairwise": "TGCACCACATGTGTTTTTTATTG\n      ||||.|.|||||.|\nTTTGACTGTAAATAAAAAGTTTGG\n",
        },
        (
            "strain1:forward",
            "TGCACCACATGTGTTTTTTATTG",
            "strain2:reverse",
            "GTTGACTGTAAATAAAAAGTTTGG",
        ): {
            "Tm": 8.0,
            "deltaG": 1.4,
            "tt_binding": 11,
            "pairwise": "TGCACCACATGTGTTTTTTATTG\n      ||||.|.|||||.|\nGTTGACTGTAAATAAAAAGTTTGG\n",
        },
        (
            "strain2:forward",
            "GGCACCACATGTGTTTTTTATTG",
            "strain1:probe",
            "AATTTCTTTGGTGGCGGGAGCAATCCTACC",
        ): {
            "Tm": 16.0,
            "deltaG": -1.8,
            "tt_binding": 8,
            "pairwise": "        GGCACCACATGTGTTTTTTATTG\n          ||||.|....|||\nAATTTCTTTGGTGGCGGGAGCAATCCTACC\n",
        },
        (
            "strain2:forward",
            "GGCACCACATGTGTTTTTTATTG",
            "strain2:probe",
            "GATTTCTTTGGTGGCGGGAGCAATCCTACC",
        ): {
            "Tm": 16.0,
            "deltaG": -1.8,
            "tt_binding": 8,
            "pairwise": "        GGCACCACATGTGTTTTTTATTG\n          ||||.|....|||\nGATTTCTTTGGTGGCGGGAGCAATCCTACC\n",
        },
        (
            "strain2:forward",
            "GGCACCACATGTGTTTTTTATTG",
            "strain1:reverse",
            "TTTGACTGTAAATAAAAAGTTTGG",
        ): {
            "Tm": 8.0,
            "deltaG": 1.4,
            "tt_binding": 11,
            "pairwise": "GGCACCACATGTGTTTTTTATTG\n      ||||.|.|||||.|\nTTTGACTGTAAATAAAAAGTTTGG\n",
        },
        (
            "strain2:forward",
            "GGCACCACATGTGTTTTTTATTG",
            "strain2:reverse",
            "GTTGACTGTAAATAAAAAGTTTGG",
        ): {
            "Tm": 8.0,
            "deltaG": 1.4,
            "tt_binding": 11,
            "pairwise": "GGCACCACATGTGTTTTTTATTG\n      ||||.|.|||||.|\nGTTGACTGTAAATAAAAAGTTTGG\n",
        },
        (
            "strain1:probe",
            "AATTTCTTTGGTGGCGGGAGCAATCCTACC",
            "strain1:reverse",
            "TTTGACTGTAAATAAAAAGTTTGG",
        ): {
            "Tm": 10.0,
            "deltaG": 1.1,
            "tt_binding": 7,
            "pairwise": "       AATTTCTTTGGTGGCGGGAGCAATCCTACC\n        ||||.|||\nTTTGACTGTAAATAAAAAGTTTGG\n",
        },
        (
            "strain1:probe",
            "AATTTCTTTGGTGGCGGGAGCAATCCTACC",
            "strain2:reverse",
            "GTTGACTGTAAATAAAAAGTTTGG",
        ): {
            "Tm": 10.0,
            "deltaG": 1.1,
            "tt_binding": 7,
            "pairwise": "       AATTTCTTTGGTGGCGGGAGCAATCCTACC\n        ||||.|||\nGTTGACTGTAAATAAAAAGTTTGG\n",
        },
        (
            "strain2:probe",
            "GATTTCTTTGGTGGCGGGAGCAATCCTACC",
            "strain1:reverse",
            "TTTGACTGTAAATAAAAAGTTTGG",
        ): {
            "Tm": 10.0,
            "deltaG": 1.1,
            "tt_binding": 7,
            "pairwise": "       GATTTCTTTGGTGGCGGGAGCAATCCTACC\n        ||||.|||\nTTTGACTGTAAATAAAAAGTTTGG\n",
        },
        (
            "strain2:probe",
            "GATTTCTTTGGTGGCGGGAGCAATCCTACC",
            "strain2:reverse",
            "GTTGACTGTAAATAAAAAGTTTGG",
        ): {
            "Tm": 10.0,
            "deltaG": 1.1,
            "tt_binding": 7,
            "pairwise": "       GATTTCTTTGGTGGCGGGAGCAATCCTACC\n        ||||.|||\nGTTGACTGTAAATAAAAAGTTTGG\n",
        },
        (
            "strain1:reverse",
            "TTTGACTGTAAATAAAAAGTTTGG",
            "strain2:reverse",
            "GTTGACTGTAAATAAAAAGTTTGG",
        ): {
            "Tm": 12.0,
            "deltaG": 0.3,
            "tt_binding": 7,
            "pairwise": "  TTTGACTGTAAATAAAAAGTTTGG\n    ||||..|.|.|\nGTTGACTGTAAATAAAAAGTTTGG\n",
        },
    }
    assert pred == obj.step3_heterodimer_all_strains_main(obs)


def test_step3_heterodimer_all_strains_main_probe_pair(obj):
    obs = {
        "probe": {
            "AATTTCTTTGGTGGCGGGAGCAATCCTACC": {
                "strain": "strain1",
                "Tm": 61.57,
                "homodimer": {},
                "mfold": {},
            },
            "GATTTCTTTGGTGGCGGGAGCAATCCTACC": {
                "strain": "strain2",
                "Tm": 61.57,
                "homodimer": {},
                "mfold": {},
            },
        },
    }
    pred = {}
    assert pred == obj.step3_heterodimer_all_strains_main(obs)


def test_step3_get_crossmaxmatching(obj):
    """
      TTTGACTGTAAATAAAAAGTTTGG
        ||||  | | |    |
    GTTGACTGTAAATAAAAAGTTTGG
    결국 gap 없는 smith-watermann 에서 local alignment 찾는 문제임.
    """
    obs1 = "TTTGACTGTAAATAAAAAGTTTGG"
    obs2 = "GTTGACTGTAAATAAAAAGTTTGG"
    pred = ["TGAC"]
    assert pred == obj.step3_get_crossmaxmatching(obs1, obs2)


def test_step3_smith_waterman(obj):
    """
      TTTGACTGTAAATAAAAAGTTTGG
        ||||  | | |    |
    GTTGACTGTAAATAAAAAGTTTGG
    결국 gap 없는 smith-watermann 에서 local alignment 찾는 문제임.
    다만 match 가 동일한 nt 에 대해서 주는 것이 아니라 complement base 에 주어야 하는 것이 문제.
    """
    obs1 = "TTTGACTGTAAATAAAAAGTTTGG"
    obs2 = "GTTGACTGTAAATAAAAAGTTTGG"
    pred = "  TTTGACTGTAAATAAAAAGTTTGG\n    ||||..|.|.|\nGTTGACTGTAAATAAAAAGTTTGG\n"
    assert pred == obj.step3_smith_waterman(obs1, obs2)
