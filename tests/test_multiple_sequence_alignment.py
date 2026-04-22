#!/usr/bin/env python

from collections import OrderedDict
from io import StringIO
import pytest
from src.multiple_sequence_alignment import MSA
from src.config import Config

MockFasta = """>MT108470.1:1-521:Aspergillus_fumigatus_strain_CNRMA17.179_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds HSP:3-21
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC
>MT108469.1:1-521:Aspergillus_fumigatus_strain_CNRMA16.56_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds HSP:3-21
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC
>MT108468.1:1-521:Aspergillus_fumigatus_strain_CNRMA15.166_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds HSP:3-21
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC
>MT108467.1:1-521:Aspergillus_fumigatus_strain_CNRMA17.580_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds HSP:3-21
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC"""

MockALN = """CLUSTAL O(1.2.4) multiple sequence alignment


gb|CP025153.1|/1100;>gb|CP025170.1|/1100;>gb|CP017626.1|/1100;>gb|AF247190.1|/1626;>ref|XM_712140.2|/15      TGAAAAAGCA
gb|CP025161.1|/1100;>gb|CP025178.1|/1100;>gb|AF247189.1|/1626;>ref|XM_712142.1|/13                           TGAAAAAGCA
emb|FM992691.1|/123;>ref|XM_002420030.1|                                                                     ----------
gb|AF184907.1|/1233                                                                                          ----------
gb|M90812.1|YSAPHR1                                                                                          ----------
gb|CP032014.1|/1100                                                                                          TGAAAAAGCA
gb|CP115801.1|/1882                                                                                          ----------

"""

MockALN2 = """CLUSTAL O(1.2.4) multiple sequence alignment


gb|CP025153.1|/1100;>gb|CP025170.1|/1100;>gb|CP017626.1|/1100;>gb|AF247190.1|/1626;>ref|XM_712140.2|/15      TGAAAAAGNA
gb|CP025161.1|/1100;>gb|CP025178.1|/1100;>gb|AF247189.1|/1626;>ref|XM_712142.1|/13                           TGAAAAAGCA
emb|FM992691.1|/123;>ref|XM_002420030.1|                                                                     ----------
gb|AF184907.1|/1233                                                                                          ----------
gb|M90812.1|YSAPHR1                                                                                          ----------
gb|CP032014.1|/1100                                                                                          XGAAAAAGCA
gb|CP115801.1|/1882                                                                                          ----------

"""
MockALN3 = """CLUSTAL format alignment by MAFFT FFT-NS-2 (v7.490)

FJ969509;MC:578                                                                                      ataca
GQ117112;MC:147                                                                                      ataca
FJ998208;MC:113                                                                                      ataca
GQ132144;MC:92                                                                                       ataca

"""


@pytest.fixture(scope="module")
def obj():
    config = Config()
    res = MSA(config)
    yield res
    del res


def test_rm_sequences(obj):
    """특정 id 들을 제거한다"""
    obs = [
        "MT108470.1:1-521:Aspergillus_fumigatus_strain_CNRMA17.179_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds/1-100"
    ]
    pred = """>MT108469.1:1-521:Aspergillus_fumigatus_strain_CNRMA16.56_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds HSP:3-21
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC
>MT108468.1:1-521:Aspergillus_fumigatus_strain_CNRMA15.166_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds HSP:3-21
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC
>MT108467.1:1-521:Aspergillus_fumigatus_strain_CNRMA17.580_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds HSP:3-21
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC
"""
    assert pred == obj.rm_sequences(StringIO(MockFasta), obs)
    obs = [
        "MT108470.1:1-521:Aspergillus_fumigatus_strain_CNRMA17.179_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds/1-100",
        "MT108469.1:1-521:Aspergillus_fumigatus_strain_CNRMA16.56_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds/111-123",
    ]
    pred = """>MT108468.1:1-521:Aspergillus_fumigatus_strain_CNRMA15.166_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds HSP:3-21
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC
>MT108467.1:1-521:Aspergillus_fumigatus_strain_CNRMA17.580_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds HSP:3-21
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC
"""
    assert pred == obj.rm_sequences(StringIO(MockFasta), obs)


def test_rc_sequences(obj):
    """특정 id 들을 reverse complement 한다"""
    obs = [
        "MT108470.1:1-521:Aspergillus_fumigatus_strain_CNRMA17.179_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds/1-100"
    ]
    pred = """>MT108470.1:1-521:Aspergillus_fumigatus_strain_CNRMA17.179_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds HSP:3-21
GGCTCAGCACACATCCGGACCGCGTGATTCTAGACAACTCTGAAGTGGTGCTGCGATTAT
>MT108469.1:1-521:Aspergillus_fumigatus_strain_CNRMA16.56_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds HSP:3-21
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC
>MT108468.1:1-521:Aspergillus_fumigatus_strain_CNRMA15.166_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds HSP:3-21
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC
>MT108467.1:1-521:Aspergillus_fumigatus_strain_CNRMA17.580_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds HSP:3-21
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC
"""
    assert pred == obj.rc_sequences(StringIO(MockFasta), obs)
    obs = [
        "MT108470.1:1-521:Aspergillus_fumigatus_strain_CNRMA17.179_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds/1-100",
        "MT108467.1:1-521:Aspergillus_fumigatus_strain_CNRMA17.580_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds/111-130",
    ]
    pred = """>MT108470.1:1-521:Aspergillus_fumigatus_strain_CNRMA17.179_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds HSP:3-21
GGCTCAGCACACATCCGGACCGCGTGATTCTAGACAACTCTGAAGTGGTGCTGCGATTAT
>MT108469.1:1-521:Aspergillus_fumigatus_strain_CNRMA16.56_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds HSP:3-21
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC
>MT108468.1:1-521:Aspergillus_fumigatus_strain_CNRMA15.166_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds HSP:3-21
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC
>MT108467.1:1-521:Aspergillus_fumigatus_strain_CNRMA17.580_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds HSP:3-21
GGCTCAGCACACATCCGGACCGCGTGATTCTAGACAACTCTGAAGTGGTGCTGCGATTAT
"""
    assert pred == obj.rc_sequences(StringIO(MockFasta), obs)


# def test_atgc_counter(obj):
# """clustalo 의 output 을 읽고 각 position 별 ATGC 및 total count 를 count 한다
# TGAAAAAGCA
# TGAAAAAGCA
# ----------
# ----------
# ----------
# TGAAAAAGCA
# ----------
# MSA 서열상에 N 이나 X 와 같은 ATGC 가 없다면 모두 count 하지 않고 TT 에만 표기한다.
# """
# pred = {
#     "A": [0, 0, 3, 3, 3, 3, 3, 0, 0, 3],
#     "G": [0, 3, 0, 0, 0, 0, 0, 3, 0, 0],
#     "T": [3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#     "C": [0, 0, 0, 0, 0, 0, 0, 0, 3, 0],
#     "TT": [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
# }
# obs = obj.aln_count_nucleotides(StringIO(MockALN))
# assert pred == obs
# ### sequence 상에 N 이 있을 때는 TT 에는 넣고 AGTC 에는 넣지 않는다.
# pred = {
#     "A": [0, 0, 3, 3, 3, 3, 3, 0, 0, 3],
#     "G": [0, 3, 0, 0, 0, 0, 0, 3, 0, 0],
#     "T": [2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#     "C": [0, 0, 0, 0, 0, 0, 0, 0, 2, 0],
#     "TT": [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
# }
# obs = obj.aln_count_nucleotides(StringIO(MockALN2))
# assert pred == obs


def test_allnt_counter(obj):
    """aln_count_nucleotides 에 대한 test.
    1. degeneracy 를 반영한다.
    2. 출력 순서를 https://en.wikiversity.org/wiki/Gene_transcriptions/Degenerate_nucleotides 에 있는 Symbol 순서 (두번째 것)에 따라 출력하는지 본다.
    Symbol[15]	Description	Bases represented
    A	adenosine
    C	cytidine
    G	guanosine
    T	thymidine
    U	uridine
    W	weak
    S	strong
    M	[Amine] amino
    K	[Ketone] keto
    R	[Purine] purine
    Y	[Pyrimidine] pyrimidine
    B	not A (B comes after A)
    D	not C (D comes after C)
    H	not G (H comes after G)
    V	not T (V comes after T and U)
    N or X	any base (not a gap)
    """
    pred = {
        "A": [0, 0, 3, 3, 3, 3, 3, 0, 0, 3],
        "C": [0, 0, 0, 0, 0, 0, 0, 0, 2, 0],
        "G": [0, 3, 0, 0, 0, 0, 0, 3, 0, 0],
        "T": [2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "U": [0] * 10,
        "W": [0] * 10,
        "S": [0] * 10,
        "M": [0] * 10,
        "K": [0] * 10,
        "R": [0] * 10,
        "Y": [0] * 10,
        "B": [0] * 10,
        "D": [0] * 10,
        "H": [0] * 10,
        "V": [0] * 10,
        "TT": [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    }
    obs = obj.aln_count_nucleotides(StringIO(MockALN2))
    assert OrderedDict(pred) == OrderedDict(obs)


def test_allcnt_counter_mafft(obj):
    """MAFFT test"""
    pred = {
        "A": [4, 0, 4, 0, 4],
        "C": [0, 0, 0, 4, 0],
        "G": [0, 0, 0, 0, 0],
        "T": [0, 4, 0, 0, 0],
        "U": [0, 0, 0, 0, 0],
        "W": [0, 0, 0, 0, 0],
        "S": [0, 0, 0, 0, 0],
        "M": [0, 0, 0, 0, 0],
        "K": [0, 0, 0, 0, 0],
        "R": [0, 0, 0, 0, 0],
        "Y": [0, 0, 0, 0, 0],
        "B": [0, 0, 0, 0, 0],
        "D": [0, 0, 0, 0, 0],
        "H": [0, 0, 0, 0, 0],
        "V": [0, 0, 0, 0, 0],
        "TT": [4, 4, 4, 4, 4],
    }
    obs = obj.aln_count_nucleotides(StringIO(MockALN3))
    assert OrderedDict(pred) == OrderedDict(obs)


def test_mk_consensus_seq(obj):
    """consensus sequence generator
    when the ratios nt1 and nt2 are same, regular expression [nt1nt2] is depicted
    ex: A: 50%, C: 50% -> [AC] where the order is alphabetical ascending.
    """
    pred = "T,G,A,A,A,A,A,G,C,A"
    obs = {
        "A": [0, 0, 3, 3, 3, 3, 3, 0, 0, 3],
        "C": [0, 0, 0, 0, 0, 0, 0, 0, 2, 0],
        "G": [0, 3, 0, 0, 0, 0, 0, 3, 0, 0],
        "T": [2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "U": [0] * 10,
        "W": [0] * 10,
        "S": [0] * 10,
        "M": [0] * 10,
        "K": [0] * 10,
        "R": [0] * 10,
        "Y": [0] * 10,
        "B": [0] * 10,
        "D": [0] * 10,
        "H": [0] * 10,
        "V": [0] * 10,
        "TT": [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    }
    assert pred == obj.mk_consensus_seq(obs)
    pred = "[AT],[CG]"
    obs = {
        "A": [2, 0],
        "C": [0, 2],
        "G": [0, 2],
        "T": [2, 0],
        "U": [0] * 2,
        "W": [0] * 2,
        "S": [0] * 2,
        "M": [0] * 2,
        "K": [0] * 2,
        "R": [0] * 2,
        "Y": [0] * 2,
        "B": [0] * 2,
        "D": [0] * 2,
        "H": [0] * 2,
        "V": [0] * 2,
        "TT": [2, 2],
    }
    assert pred == obj.mk_consensus_seq(obs)


def test_print_out_aln_count(obj):
    pred = (
        "nt,1,2,3,4,5,6,7,8,9,10\n"
        "A,0,0,3,3,3,3,3,0,0,3\n"
        "C,0,0,0,0,0,0,0,0,3,0\n"
        "G,0,3,0,0,0,0,0,3,0,0\n"
        "T,3,0,0,0,0,0,0,0,0,0\n"
        "U,0,0,0,0,0,0,0,0,0,0\n"
        "W,0,0,0,0,0,0,0,0,0,0\n"
        "S,0,0,0,0,0,0,0,0,0,0\n"
        "M,0,0,0,0,0,0,0,0,0,0\n"
        "K,0,0,0,0,0,0,0,0,0,0\n"
        "R,0,0,0,0,0,0,0,0,0,0\n"
        "Y,0,0,0,0,0,0,0,0,0,0\n"
        "B,0,0,0,0,0,0,0,0,0,0\n"
        "D,0,0,0,0,0,0,0,0,0,0\n"
        "H,0,0,0,0,0,0,0,0,0,0\n"
        "V,0,0,0,0,0,0,0,0,0,0\n"
        "TT,3,3,3,3,3,3,3,3,3,3\n"
        "con,T,G,A,A,A,A,A,G,C,A"
    )
    obs = {
        "A": [0, 0, 3, 3, 3, 3, 3, 0, 0, 3],
        "C": [0, 0, 0, 0, 0, 0, 0, 0, 3, 0],
        "G": [0, 3, 0, 0, 0, 0, 0, 3, 0, 0],
        "T": [3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "U": [0] * 10,
        "W": [0] * 10,
        "S": [0] * 10,
        "M": [0] * 10,
        "K": [0] * 10,
        "R": [0] * 10,
        "Y": [0] * 10,
        "B": [0] * 10,
        "D": [0] * 10,
        "H": [0] * 10,
        "V": [0] * 10,
        "TT": [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    }
    assert pred == obj.print_out_aln_count(obs)
