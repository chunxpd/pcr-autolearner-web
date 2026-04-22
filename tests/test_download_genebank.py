import pytest
from src.download_genebank import DwGeneBank
from Bio.Seq import UndefinedSequenceError


class MockHSP:
    """
    |  Members:
    |      - score           BLAST score of hit.  (float)
    |      - bits            Number of bits for that score.  (float)
    |      - expect          Expect value.  (float)
    |      - num_alignments  Number of alignments for same subject.  (int)
    |      - identities      Number of identities (int) if using the XML parser.
    |        Tuple of number of identities/total aligned (int, int)
    |        if using the (obsolete) plain text parser.
    |      - positives       Number of positives (int) if using the XML parser.
    |        Tuple of number of positives/total aligned (int, int)
    |        if using the (obsolete) plain text parser.
    |      - gaps            Number of gaps (int) if using the XML parser.
    |        Tuple of number of gaps/total aligned (int, int) if
    |        using the (obsolete) plain text parser.
    |      - align_length    Length of the alignment. (int)
    |      - strand          Tuple of (query, target) strand.
    |      - frame           Tuple of 1 or 2 frame shifts, depending on the flavor.
    |
    |      - query           The query sequence.
    |      - query_start     The start residue for the query sequence.  (1-based)
    |      - query_end       The end residue for the query sequence.  (1-based)
    |      - match           The match sequence.
    |      - sbjct           The sbjct sequence.
    |      - sbjct_start     The start residue for the sbjct sequence.  (1-based)
    |      - sbjct_end       The end residue for the sbjct sequence.  (1-based)

    """

    def __init__(self):
        """
        <Hsp_num>5</Hsp_num>
        <Hsp_bit-score>26.2635</Hsp_bit-score>
        <Hsp_score>13</Hsp_score>
        <Hsp_evalue>10.9313</Hsp_evalue>
        <Hsp_query-from>3</Hsp_query-from>
        <Hsp_query-to>15</Hsp_query-to>
        <Hsp_hit-from>954515</Hsp_hit-from>
        <Hsp_hit-to>954503</Hsp_hit-to>
        <Hsp_query-frame>1</Hsp_query-frame>
        <Hsp_hit-frame>-1</Hsp_hit-frame>
        <Hsp_identity>13</Hsp_identity>
        <Hsp_positive>13</Hsp_positive>
        <Hsp_gaps>0</Hsp_gaps>
        <Hsp_align-len>13</Hsp_align-len>
        <Hsp_qseq>TGCTCGCAACCAC</Hsp_qseq>
        <Hsp_hseq>TGCTCGCAACCAC</Hsp_hseq>
        <Hsp_midline>|||||||||||||</Hsp_midline>
        """
        self.score = 13
        self.bits = 26.2635
        self.expec = 10.9313
        self.num_alignments = 13
        self.identities = 13
        self.positives = 13
        self.gaps = 0
        self.align_length = 13
        self.strand = ("Plus", "Plus")
        self.frame = (1, -1)
        self.query = "TGCTCGCAACCAC"
        self.query_start = 15
        self.query_end = 15
        self.match = "|||||||||||||"
        self.sbjct = "TGCTCGCAACCAC"
        self.sbjct_start = 954515
        self.sbjct_end = 954503


@pytest.fixture(scope="module")
def obj():
    config = {
        "email": "io2oi@yahoo.co.kr",
        "jobid": "tests_jobiid_00000",
        "TEMPDIR": "tests",
        "GB_STORAGE": "tests/GBfiles",
        "FA_STORAGE": "tests/FASTA",
    }
    res = DwGeneBank(config, "aaaaaaaaa")
    yield res
    del res


def test_setting_check(obj):
    pred = "tests/tests_jobiid_00000-aaaaaaaaa-genebank.fasta"
    obs = obj.outfn
    assert pred == obs


@pytest.mark.skip("한 번 성공하면 그 뒤로는 gb 파일이 존재함. test 성공 못함")
def test_download_all_gbfiles(obj):
    ## 새로운 gbid 가 하나 이상 있을 때
    pred = "tests/tests_jobiid_00000/tests_jobiid_00000_tests_jobiid_00000_aaaaaaaaa_all.gb"
    obs_hsps = [{"gbid": "AP022596.1"}, {"gbid": "CP023583.1"}, {"gbid": "BX571966.1"}]
    obs = obj.download_all_gbfiles(obs_hsps)
    assert pred == obs
    ## 새로운 gbid 가 하나도 없을 때
    pred = ""
    obs_hsps = [{"gbid": "AP022596.1"}, {"gbid": "CP023583.1"}]
    obs = obj.download_all_gbfiles(obs_hsps)
    assert pred == obs


def test_gbid_list_iterator(obj):
    pred = ["AP022596.1,CP023583.1", "BX571966.1"]
    obs_list = ["AP022596.1", "CP023583.1", "BX571966.1"]
    obs = obj.gbid_list_iterator(obs_list, 2)
    assert pred == obs
    pred = ["AP022596.1", "CP023583.1", "BX571966.1"]
    obs_list = ["AP022596.1", "CP023583.1", "BX571966.1"]
    obs = obj.gbid_list_iterator(obs_list, 3)
    assert pred == obs
    pred = ["AP022596.1,CP023583.1,BX571966.1"]
    obs_list = ["AP022596.1", "CP023583.1", "BX571966.1"]
    obs = obj.gbid_list_iterator(obs_list, 4)
    assert pred == obs


def test_gbid_list_iterator2(obj):
    pred = [
        [{"gbid": "AP022596.1"}],
        [{"gbid": "CP023583.1"}],
        [{"gbid": "BX571966.1"}],
    ]
    obs_hsps = [{"gbid": "AP022596.1"}, {"gbid": "CP023583.1"}, {"gbid": "BX571966.1"}]
    obs = obj.gbid_list_iterator(obs_hsps, 3, joinq=False)
    assert pred == obs
    pred = [
        [{"gbid": "AP022596.1"}, {"gbid": "CP023583.1"}],
        [{"gbid": "BX571966.1"}],
    ]
    obs_hsps = [{"gbid": "AP022596.1"}, {"gbid": "CP023583.1"}, {"gbid": "BX571966.1"}]
    obs = obj.gbid_list_iterator(obs_hsps, 2, joinq=False)
    assert pred == obs


@pytest.mark.skip("눈으로 확인 하였음. test 하기 불가능")
def test_split_each_gb(obj):
    pred = ""
    obs_gb_all = "tests/tests_jobiid_00000/tests_jobiid_00000_tests_jobiid_00000_all.gb"
    obs = obj.split_each_gb(obs_gb_all)
    assert pred == obs


def test_get_entrez_sequence_info(obj):
    """gene feature 에서 join 이 들어있거나 max_updown_len 보다 gene 이 크거나 할 경우는 포함하지 않음.
    이유는 아직 정확하게 parsing 할 수 있지가 않기 때문.
    return type:
    {
        "seq": SeqIO.Seqfeature.seq,
        "organism": <str>,
        "description": <str>,
        "features": [{
            "start": SeqIO.Seqfeature.location.start, # 0-based
            "end": SeqIO.Seqfeature.location.end,     # 1-based
            "strand": SeqIO.Seqfeature.location.strand, # -1: reverse complement
        },...]
    }
    """
    obs = obj.get_entrez_sequence_info("CP054795.1")
    pred = 6993871
    assert pred == len(obs["seq"])
    assert obs["seq"][:10] == "GGCGAACGCA"
    pred = 6673  # gene 의 갯수
    assert pred == len(obs["features"])
    assert obs["organism"] == "Mycolicibacterium smegmatis"
    assert obs["description"] == "Mycolicibacterium smegmatis strain FDAARGOS_679 chromosome"


def test_get_hsp_hit_genomic_position(obj):
    """
    Hsp_hit-frame 이 -1 이면 Hsp_hit-from 과 Hsp_hit-to 를 서로 바꾸어 5'->3'으로 반환받는다.
    <Hsp>
      <Hsp_num>5</Hsp_num>
      <Hsp_bit-score>26.2635</Hsp_bit-score>
      <Hsp_score>13</Hsp_score>
      <Hsp_evalue>10.9313</Hsp_evalue>
      <Hsp_query-from>3</Hsp_query-from>
      <Hsp_query-to>15</Hsp_query-to>
      <Hsp_hit-from>954515</Hsp_hit-from>
      <Hsp_hit-to>954503</Hsp_hit-to>
      <Hsp_query-frame>1</Hsp_query-frame>
      <Hsp_hit-frame>-1</Hsp_hit-frame>
      <Hsp_identity>13</Hsp_identity>
      <Hsp_positive>13</Hsp_positive>
      <Hsp_gaps>0</Hsp_gaps>
      <Hsp_align-len>13</Hsp_align-len>
      <Hsp_qseq>TGCTCGCAACCAC</Hsp_qseq>
      <Hsp_hseq>TGCTCGCAACCAC</Hsp_hseq>
      <Hsp_midline>|||||||||||||</Hsp_midline>
    </Hsp>
    """
    mock_hsp = MockHSP()
    mock_hsp.strand = ("Plus", "Minus")
    mock_hsp.sbjct_start = 954515
    mock_hsp.sbjct_end = 954503
    pred = {
        "Hsp_hit-from": 954503,
        "Hsp_hit-to": 954515,
        "Hsp_hit-strand": -1,
    }
    obs = obj.get_hsp_hit_genomic_position(mock_hsp)
    assert pred == obs


def test_hsp_select_criteria(obj):
    """hsp 를 가지고 올 기준 수립
    2022년 10월 12일 영상회의 결과: https://github.com/io2oi/ctlabs-primerwrapper/issues/16 참조
    query_length 에 대해서 alignment 가 80% 이상 되고 identity는 90% 이상일 때 가지고 옴
    query_length 를 21bp 라고 할 때
    80% 는 21*0.8 = 16.8 즉 17bp 이상이어야 하고
    identity 는 17bp * 0.9 = 15.3 즉 16bp 이상이어야 한다."""
    obs_query_length = 21
    mock_hsp = MockHSP()
    mock_hsp.align_length = 17
    mock_hsp.identities = 16
    obs = obj.hsp_select_criteria(obs_query_length, mock_hsp)
    assert True == obs
    mock_hsp.align_length = 16
    mock_hsp.identities = 16
    obs = obj.hsp_select_criteria(obs_query_length, mock_hsp)
    assert False == obs
    mock_hsp.align_length = 17
    mock_hsp.identities = 15
    obs = obj.hsp_select_criteria(obs_query_length, mock_hsp)
    assert False == obs
    mock_hsp.align_length = 17
    mock_hsp.identities = 15
    obs = obj.hsp_select_criteria(obs_query_length, mock_hsp, al_thr_coeff=0.6, id_thr_coeff=0.6)
    assert True == obs


def test_find_overlab_entrez_genes(obj):
    pred = (
        {"end": 2714465, "start": 2711298, "strand": -1},
        {"end": 2716351, "start": 2714824, "strand": -1},
    )
    obs_hsp = {"gbid": "CP054795.1", "Hsp_hit-from": 2714562, "Hsp_hit-to": 2714582}
    obs = obj.find_overlab_entrez_genes(obs_hsp)
    assert pred == obs


def test_adj_overlab_entrez_genes(obj):
    # 변화 없을 때. 즉 길이가 < 5KB
    obs_reg = (
        {"end": 2714465, "start": 2711298, "strand": -1},
        {"end": 2716351, "start": 2714824, "strand": -1},
    )
    obs_gbdic = {"Hsp_hit-strand": -1}
    pred = (
        {"end": 2714465, "start": 2711298, "strand": -1},
        {"end": 2716351, "start": 2714824, "strand": -1},
    )
    obs = obj._adj_overlab_entrez_genes(2714562, 2714582, 10000000, obs_reg, obs_gbdic)
    assert pred == obs
    # 변화 있을 때. 즉 길이가 > 5KB
    obs_reg = (
        {"end": 2714465, "start": 0, "strand": -1},
        {"end": 10000000, "start": 2714824, "strand": -1},
    )
    pred = (
        {"end": 2714562, "start": 2714562 - 5000, "strand": -1},
        {"end": 2714582 + 5000, "start": 2714582, "strand": -1},
    )
    obs = obj._adj_overlab_entrez_genes(2714562, 2714582, 10000000, obs_reg, obs_gbdic)
    assert pred == obs


def test_find_overlab_entrez_genes_errorcase1(obj):
    """gb file 에 gene feature 가 없는 경우
    HSP 의 Hsp_hit-from - 50KB ~ Hsp_hit-to + 50KB 를 반환한다.
    """
    pred = (
        {"end": 3834276, "start": 3834276 - 5000, "strand": 1},
        {"end": 3834296 + 5000, "start": 3834296, "strand": 1},
    )
    obs_hsp = {
        "gbid": "CP102342.1",
        "Hsp_hit-from": 3834276,
        "Hsp_hit-to": 3834296,
        "Hsp_hit-strand": 1,
    }
    obs = obj.find_overlab_entrez_genes(obs_hsp)
    assert pred == obs


def test_find_overlab_entrez_genes_errorcase2(obj):
    """gb file 에 gene feature 가 1개만 있는 경우.
    위와 동일하게 처리하자.
    """
    pred = (
        {"start": 5616 - 5000, "end": 5616, "strand": 1},
        {"start": 5627, "end": 5627 + 5000, "strand": 1},
    )
    obs_hsp = {
        "gbid": "AJ238027.1",
        "Hsp_hit-from": 5616,
        "Hsp_hit-to": 5627,
        "Hsp_hit-strand": 1,
    }
    obs = obj.find_overlab_entrez_genes(obs_hsp)
    assert pred == obs


def test_find_overlab_entrez_gnenes_errorcase3(obj):
    """issue #61 case
    HSP 의 영역이 gb file 에 정의된 gene 영역의 제일 마지막보다 더 뒤쪽일 때
    바로 앞의 gene 영역과 맨 마지막까지의 서열을 가지고 온다.
    """
    pred = (
        {"end": 7466, "start": 6074, "strand": 1},
        {"end": 7642, "start": 7579, "strand": -1},
    )
    obs_hsp = {
        "gbid": "AF377339.1",
        "Hsp_hit-from": 7563,
        "Hsp_hit-to": 7579,
        "Hsp_hit-strand": -1,
    }
    obs = obj.find_overlab_entrez_genes(obs_hsp)
    assert pred == obs


def test_find_overlab_entrez_gnenes_errorcase4(obj):
    """issue #61 case
    HSP 의 영역이 gb file 에 정의된 gene 영역의 제일 앞보다 더 앞일 때
    0 부터 첫번째 gene 영역을 가지고 온다.
    """
    pred = (
        {"end": 100, "start": 0, "strand": -1},
        {"end": 1708, "start": 445, "strand": -1},
    )
    obs_hsp = {
        "gbid": "AF377339.1",
        "Hsp_hit-from": 100,
        "Hsp_hit-to": 200,
        "Hsp_hit-strand": -1,
    }
    obs = obj.find_overlab_entrez_genes(obs_hsp)
    assert pred == obs


def test_find_overlab_entrez_genes_realcase1(obj):
    """??? 왜 길이가 거의 10Kbp 인가?
    >HG917972.2:1156347-1166341 HSP:Mycobacterium marinum E11 main chromosome genome:1158680-1158696
    통과조건: obj.max_updown_len *2 + query primer seq length 보다 pred 의 길이가 더 짧아야 한다.
    """
    pred = (
        {"end": 1161456, "start": 1157925, "strand": 1},
        {"end": 1161456, "start": 1157925, "strand": 1},
    )
    obs_hsp = {
        "gbid": "HG917972.2",
        "Hsp_hit-from": 1158680,
        "Hsp_hit-to": 1158696,
        "Hsp_hit-strand": 1,
    }
    obs = obj.find_overlab_entrez_genes(obs_hsp)
    assert pred == obs
    qlen = 1158696 - 1158680
    assert abs(obs[1]["end"] - obs[0]["start"]) < qlen + obj.max_updown_len * 2


def test_find_overlab_entrez_genes_realcase2(obj):
    """issue102: AF035716.1
        Entrez genes : s1--------e1  s2-------e2  s3-----------e3
        Query HSP    :                               !--------$
    위와 같이 query HSP 가 gene들 중 제일 마지막에 있을 경우 downstream 을 찾으면
        s - $ 의 부호는
        Entrez genes : (-)--------e1 (-)-------e2 (-)-----------e3
        Query HSP    :                               !--------$
    가 된다.
    따라서 이 경우는 전체 gene 들을 모두 돌면서 그 안에 있는 것을 확인하는 부분이 있어 줘야 한다.
    gbdic:
    {
        'gbid': 'AF035716.1',
        'seq': Seq('CCTTATACAGCTCAATTAGTGATTAGTGAATAGATAAAAATATGTTAAATATAT...AGT'),
        'features': [
            {'start': ExactPosition(0), 'end': ExactPosition(33), 'strand': 1},
            {'start': ExactPosition(154), 'end': ExactPosition(652), 'strand': 1},
            {'start': ExactPosition(2485), 'end': ExactPosition(6076), 'strand': 1}
            ],
        'organism': 'Clostridioides difficile',
        'description': 'Clostridium difficile toxin B (toxB) gene, partial cds; putative regulatory protein (txe2) gene, complete cds; and toxin A (toxA) gene, partial rearranged sequence'
    }
    """
    pred = ({"end": 6076, "start": 2485, "strand": -1}, {"end": 6076, "start": 2485, "strand": -1})
    obs_hsp = {"gbid": "AF035716.1", "Hsp_hit-from": 3253, "Hsp_hit-to": 3282, "Hsp_hit-strand": -1}
    obs = obj.find_overlab_entrez_genes(obs_hsp)
    assert pred == obs
    pred = ({"end": 33, "start": 0, "strand": -1}, {"end": 33, "start": 0, "strand": -1})
    obs_hsp = {"gbid": "AF035716.1", "Hsp_hit-from": 1, "Hsp_hit-to": 20, "Hsp_hit-strand": -1}
    obs = obj.find_overlab_entrez_genes(obs_hsp)
    assert pred == obs


def test_find_overlab_entrez_genes_realcase3(obj):
    """길이가 너무 긴 결과는 self.max_updown_len * 2 + querylen 보다 짧아야 한다.
    아래의 GB 파일에서 아래의 hsp 영역에 맞는 결과는
    (
        {"start": 13304669, "end": 13305264, "strand": -1},
        {"start": 13848781, "end": 13849738, "strand": -1},
    ) 으로 545069bp 이다. 하지만 그러면 self.max_updown_len 5000 bp 보다 너무 길게 되므로
    Hsp_hit-from - 5000 ~ Hsp_hit-to + 5000 bp 까지만 가지고 온다.
    """
    pred = (
        {
            "start": 13481458 - obj.max_updown_len,
            "end": 13481458,
            "strand": 1,
        },
        {"start": 13481478, "end": 13481478 + obj.max_updown_len, "strand": 1},
    )
    obs_hsp = {
        "gbid": "LR778263.1",
        "Hsp_hit-from": 13481458,
        "Hsp_hit-to": 13481478,
        "Hsp_hit-strand": 1,
    }
    obs = obj.find_overlab_entrez_genes(obs_hsp)
    assert pred == obs


def test_get_ntsequence(obj):
    ### + strand
    pred = (
        "tcgaggaccgcgtcgtcatggtcgcgttgaccctggccgcaggcgacaccgaactgccgagaaactcgtcg"
        "acgagatcatcgacggccggttccagccggccacccccacgttcctcaactcgggcaagaagcagcgcggc"
        "gagccggtgtcctgcttcctgttgcgcatcgaggacaa"
    ).upper()
    obs_overlab = (
        {"end": 6993060, "start": 6993000, "strand": 1},
        {"end": 6993180, "start": 6993120, "strand": 1},
    )
    obs = obj.get_ntseqeunce("CP054795.1", obs_overlab)
    assert len(obs) == obs_overlab[1]["end"] - obs_overlab[0]["start"]
    assert len(pred) == obs_overlab[1]["end"] - obs_overlab[0]["start"]
    assert pred == str(obs)
    ### - strand
    pred = (
        "TTGTCCTCGATGCGCAACAGGAAGCAGGACACCGGCTCGCCGCGCTGCTTCTTGCCCGAGTTGAGGAACGT"
        "GGGGGTGGCCGGCTGGAACCGGCCGTCGATGATCTCGTCGACGAGTTTCTCGGCAGTTCGGTGTCGCCTGC"
        "GGCCAGGGTCAACGCGACCATGACGACGCGGTCCTCGA"
    )
    obs_overlab = (
        {"end": 6993060, "start": 6993000, "strand": -1},
        {"end": 6993180, "start": 6993120, "strand": -1},
    )
    obs = obj.get_ntseqeunce("CP054795.1", obs_overlab)
    assert pred == str(obs)


def test_get_ntsequence2(obj):
    """full sequence 의 길이에 대한 threshold 보다 작을 경우는 전체 서열을 가지고 온다"""
    pred = 6993871  # CP054795.1 의 서열의 전체 길이
    obs_overlab = (
        {"end": 6993060, "start": 6993000, "strand": 1},
        {"end": 6993180, "start": 6993120, "strand": 1},
    )
    obs = obj.get_ntseqeunce("CP054795.1", obs_overlab, lenthr=10**100)
    assert pred == len(obs)


def test_get_ntsequence_errorcase1(obj):
    """gb file 에 gene feature 가 없고 sequence 도 없는 경우.
    Entrez 에서 fasta 파일을 config['FA_STORAGE'] 의 경로에 저장하고 실행한다.
    """
    pred = "TGATTTACGATCTGGAACAAGGCCAAGCTCATTGGAACAA"
    obs_overlab = (
        {"end": 10, "start": 0, "strand": 1},
        {"end": 40, "start": 30, "strand": 1},
    )
    obs = obj.get_ntseqeunce("LR594554.2", obs_overlab)
    assert pred == obs
    # with pytest.raises(UndefinedSequenceError):
    #     obs = obj.get_ntseqeunce("LR594554.2", obs_overlab)


def test_get_ntsequence_realcase1(obj):
    """아래의 경우 서열 길이가 9868bp 로 오고 있음.
    lengthr 보다 gbid 의 서열이 짧으면 overlab 정보 상관 없이 전체서열을 가지고 온다.
    lengthr 보다 overlab 이 크면 3220-1206 인 2014bp 을 가지고 온다.
    >JF271804.1:1206-3220 HSP:Mycobacterium sp. strain MOTT-01 rpoBC operon sequence:2206-2220
    """
    pred = 2014
    obs_overlab = (
        {"start": 1206, "end": 2205, "strand": 1},
        {"start": 2220, "end": 3220, "starnd": 1},
    )
    obs = obj.get_ntseqeunce("JF271804.1", obs_overlab, lenthr=500)
    assert pred == len(obs)
    pred = 9868
    obs_overlab = (
        {"start": 1206, "end": 2205, "strand": 1},
        {"start": 2220, "end": 3220, "starnd": 1},
    )
    obs = obj.get_ntseqeunce("JF271804.1", obs_overlab, lenthr=10000)
    assert pred == len(obs)


# def test_get_ntsequence_realcase2(obj):
#     """issue102 case: gene 사이의 길이가 매우 긴 경우.
#     이런 경우에는 애초에 overlab 구간에서 matching region 앞 뒤로 up/down stream 의 길이를 조정해서 왔어야 한다."""
#     pred = 13463557 - 13380019
#     obs_overlab = (
#         {"start": 13380019, "end": 13402141, "strand": 1},
#         {"start": 13430850, "end": 13463557, "strand": 1},
#     )
#     obs = obj.get_ntseqeunce("LR778263.1", obs_overlab, lenthr=10000)
#     assert pred == len(obs)


def test_convert_seq_to_fasta(obj):
    """clustalo 에서 종명이 잘리는 문제때문에 id 부분에 gbid:seqrange:종명 으로 넣고
    description 으로 HSP:Hsp_from:Hsp_to 정보를 넣는다.
    """
    pred = (
        ">CP054795.1:complement(6993001-6993180):Mycolicibacterium_smegmatis_strain_FDAARGOS_679_chromosome HSP:2714563-2714582\n"
        "TTGTCCTCGATGCGCAACAGGAAGCAGGACACCGGCTCGCCGCGCTGCTTCTTGCCCGAG\n"
        "TTGAGGAACGTGGGGGTGGCCGGCTGGAACCGGCCGTCGATGATCTCGTCGACGAGTTTC\n"
        "TCGGCAGTTCGGTGTCGCCTGCGGCCAGGGTCAACGCGACCATGACGACGCGGTCCTCGA\n"
    )
    obs_seq = (
        "TTGTCCTCGATGCGCAACAGGAAGCAGGACACCGGCTCGCCGCGCTGCTTCTTGCCCGAGTTGAGGAACGT"
        "GGGGGTGGCCGGCTGGAACCGGCCGTCGATGATCTCGTCGACGAGTTTCTCGGCAGTTCGGTGTCGCCTGC"
        "GGCCAGGGTCAACGCGACCATGACGACGCGGTCCTCGA"
    )
    obs_overlab = (
        {"end": 6993060, "start": 6993000, "strand": -1},
        {"end": 6993180, "start": 6993120, "strand": -1},
    )
    obs_meta = {
        "name": "Mycolicibacterium smegmatis strain FDAARGOS_679 chromosome",
        "hsp": {"gbid": "CP054795.1", "Hsp_hit-from": 2714562, "Hsp_hit-to": 2714582},
    }
    obs = obj.convert_seq_to_fasta(obs_meta, obs_seq, obs_overlab)
    assert pred == obs


def test_get_entrez_fasta(obj):
    pred = "AGTTGACTCTACCCTAGTTGACTCTACCCTAGTTGACTCT"
    obs = obj.get_entrez_fasta(
        "LR812045.1",
        (
            {"start": 0, "end": 10, "strand": 1},
            {"start": 30, "end": 40, "strand": 1},
        ),
    )
    assert pred == obs


def test_get_entrez_fasta2(obj):
    # reverse seq of ATAACCAATCAGATTTTCACTTTGAGGATGGAAGGCTGCGAGACTATAAAAATCATTTTCCAGTAGATTT
    pred = "AAATCTACTGGAAAATGATTTTTATAGTCTCGCAGCCTTCCATCCTCAAAGTGAAAATCTGATTGGTTAT"
    obs = obj.get_entrez_fasta(
        "OX104105.1",
        (
            {"start": 0, "end": 30, "strand": -1},
            {"start": 40, "end": 70, "strand": -1},
        ),
    )
    assert pred == obs
