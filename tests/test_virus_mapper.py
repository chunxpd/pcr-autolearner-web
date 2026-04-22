import pysam
import pytest
from src.virus_mapper import Mapper
from io import StringIO
from Bio import SeqIO


@pytest.fixture(scope="module")
def mock_ref_seq():
    sequence = """>sequence1
ATG---CATGC---ATGC
>sequence2
ATG---CGCA---TGCGC
>sequence3
ATG---CGC--TTT--GC
"""
    fasta_records = list(SeqIO.parse(StringIO(sequence), "fasta"))
    yield fasta_records
    del fasta_records


@pytest.fixture(scope="module")
def mock_input_seq():
    sequence = """>in_seq1
ACGGT[AG]CCT
>in_seq2
TTAC[AT]TTAG[GC]AA
"""
    fasta_records = list(SeqIO.parse(StringIO(sequence), "fasta"))
    yield fasta_records
    del fasta_records


@pytest.fixture(scope="module")
def mock_sam():
    sam = "tests/test_virus_query.sam"
    res = pysam.AlignmentFile(sam, "r")
    yield res
    del res


@pytest.fixture(scope="module")
def obj():
    res = Mapper()
    yield res
    del res


@pytest.fixture(scope="module")
def gap_dic(obj, mock_ref_seq):
    res = obj.mk_gap_dic(mock_ref_seq)
    yield res
    del res


def test_mk_gap_dic(obj, mock_ref_seq):
    """sequence1
    ori. pos : 1 3   4   8   9  12
    sequence1: ATG---CATGC---ATGC
    gap. ps  : 1 3   7       15
    >sequence3
    ori. pos : 1 3   4 6  7 9  10
    sequence3: ATG---CGC--TTT--GC
    gap. ps  :
    """
    pred = {
        "sequence1": [(3, 6), (8, 11)],
        "sequence2": [(3, 6), (7, 10)],
        "sequence3": [(3, 6), (6, 8), (9, 11)],
    }
    assert pred == obj.mk_gap_dic(mock_ref_seq)


def test_adj_map_pos(obj, gap_dic):
    """SAM/BAM file: the 4th column: POS 1-based leftmost mapping POSition from SAM specification
    ori. pos : 1 3   4   8   9  12
    sequence1: ATG---CATGC---ATGC
    gap. ps  : 1 3   7       15
    query1   : ATG---C
    => 1
    sequence1: ATG---CATGC---ATGC
    query1   :       CATGC
    => adjust from 4 to 7
    sequence1: ATG---CATGC---ATGC
    query1   :               ATGC
    => adjust from 9 to 15
    """
    obs = 1  # ex: the short sequence is ATGC
    pred = 1
    assert pred == obj.adj_map_pos(obs, "sequence1", gap_dic)
    obs = 4  # ex: the query seq is CATGC
    pred = 7
    assert pred == obj.adj_map_pos(obs, "sequence1", gap_dic)
    obs = 9
    pred = 15
    assert pred == obj.adj_map_pos(obs, "sequence1", gap_dic)


def test_rm_gaps(obj, mock_ref_seq):
    pred = """>sequence1
ATGCATGCATGC
>sequence2
ATGCGCATGCGC
>sequence3
ATGCGCTTTGC"""
    assert pred == seqrec_to_str(obj.rm_gaps(mock_ref_seq))


def test_gen_seqs_from_degenerancy(obj, mock_input_seq):
    obs = "ACGGT[AG]CCT"
    pred = ["ACGGTACCT", "ACGGTGCCT"]
    assert pred == obj.gen_seqs_from_degerancy(obs)
    obs = "TTAC[AT]TTAG[GC]AA"
    pred = sorted(
        [
            "TTAC A TTAG G AA".replace(" ", ""),
            "TTAC A TTAG C AA".replace(" ", ""),
            "TTAC T TTAG G AA".replace(" ", ""),
            "TTAC T TTAG C AA".replace(" ", ""),
        ]
    )
    assert pred == obj.gen_seqs_from_degerancy(obs)
    obs = "AAAA[GC]TTT[ACG]TTT"
    pred = [
        "AAAA C TTT A TTT".replace(" ", ""),
        "AAAA C TTT C TTT".replace(" ", ""),
        "AAAA C TTT G TTT".replace(" ", ""),
        "AAAA G TTT A TTT".replace(" ", ""),
        "AAAA G TTT C TTT".replace(" ", ""),
        "AAAA G TTT G TTT".replace(" ", ""),
    ]
    assert pred == obj.gen_seqs_from_degerancy(obs)


def test_gen_fasta_from_degenerancy(obj, mock_input_seq):
    """ """
    pred_seq = (
        ">in_seq1_1\n"
        "ACGGTACCT\n"
        ">in_seq1_2\n"
        "ACGGTGCCT\n"
        ">in_seq2_1\n"
        "TTACATTAGCAA\n"
        ">in_seq2_2\n"
        "TTACATTAGGAA\n"
        ">in_seq2_3\n"
        "TTACTTTAGCAA\n"
        ">in_seq2_4\n"
        "TTACTTTAGGAA"
    )
    assert pred_seq == seqrec_to_str(obj.gen_fasta_from_degenerancy(mock_input_seq))


def seqrec_to_str(fa_rec):
    return "\n".join(["\n".join([f">{v.id}", f"{v.seq}"]) for v in fa_rec])


def test_parse_sam_file(obj, mock_sam):
    """SAM format: 1-based mapping position. Only mapped queries are returned"""
    pred = {
        "NC_002205.1": {"query1_1": 1, "query1_2": 1},
        "NC_002207.1": {"query3_1": 87, "query3_2": 87},
        "NC_002208.1": {"query4_1": 1, "query4_2": 1},
        "NC_002209.1": {"query5_1": 227},
    }
    assert pred == obj.parse_sam_file(mock_sam)


def test_transform_data_to_matrix(obj):
    obs = {
        "NC_002205.1": {"query1_1": 1, "query1_2": 1},
        "NC_002207.1": {"query3_1": 87, "query3_2": 87},
        "NC_002208.1": {"query4_1": 1, "query4_2": 1},
        "NC_002209.1": {"query5_1": 227},
    }
    pred = (
        [
            "query1_1",
            "query1_2",
            "query3_1",
            "query3_2",
            "query4_1",
            "query4_2",
            "query5_1",
        ],
        ["NC_002205.1", "NC_002207.1", "NC_002208.1", "NC_002209.1"],
        {
            "NC_002205.1": {
                "query1_1": 1,
                "query1_2": 1,
                "query3_1": None,
                "query3_2": None,
                "query4_1": None,
                "query4_2": None,
                "query5_1": None,
            },
            "NC_002207.1": {
                "query1_1": None,
                "query1_2": None,
                "query3_1": 87,
                "query3_2": 87,
                "query4_1": None,
                "query4_2": None,
                "query5_1": None,
            },
            "NC_002208.1": {
                "query1_1": None,
                "query1_2": None,
                "query3_1": None,
                "query3_2": None,
                "query4_1": 1,
                "query4_2": 1,
                "query5_1": None,
            },
            "NC_002209.1": {
                "query1_1": None,
                "query1_2": None,
                "query3_1": None,
                "query3_2": None,
                "query4_1": None,
                "query4_2": None,
                "query5_1": 227,
            },
        },
    )
    assert pred == obj.transform_data_to_matrix(obs)
