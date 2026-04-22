from io import StringIO
import pytest
from src.merge_fasta import MERGE_FASTA

MockFasta = """>MT108470.1:1-521:Aspergillus_fumigatus_strain_CNRMA17.179_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds HSP:3-21
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC
>MT108469.1:1-521:Aspergillus_fumigatus_strain_CNRMA16.56_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds HSP:3-21
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC
>MT108468.1:1-521:Aspergillus_fumigatus_strain_CNRMA15.166_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds HSP:3-21
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC
>MT108467.1:1-521:Aspergillus_fumigatus_strain_CNRMA17.580_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds HSP:3-21
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC"""

MockFasta2 = """>seq1
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC
>seq2
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC
>seq3
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC
>seq4
ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCC"""


@pytest.fixture(scope="module")
def obj():
    res = MERGE_FASTA()
    yield res
    del res


def test_rmgap_infasta(obj):
    """입력된 fasta 파일에서 gap 을 제거한다."""
    obs = "tests/merge_fasta_example.fa"
    pred = "tests/merge_fasta_example_rmgap.fa"
    assert pred == obj.rmgap_infasta(obs)


def test_run_cdhit100(obj):
    pred = [
        "tests/merge_fasta_example_rmgap.fa",
        "tests/merge_fasta_example_rmgap_100.fa.clstr",
    ]
    obs = obj.run_cdhit100("tests/merge_fasta_example_rmgap.fa")
    assert pred == obs


def test_sub_parse_100_clstr(obj):
    pred = [
        [
            ">Cluster 0",
            " >gb|CP025153.1|/1100",
            " >gb|CP025170.1|/1100",
            " >gb|CP017626.1|/1100",
            " >gb|AF247190.1|/1626",
            " >ref|XM_712140.2|/15",
        ],
        [
            ">Cluster 1",
            " >gb|CP025161.1|/1100",
            " >gb|CP025178.1|/1100",
            " >gb|AF247189.1|/1626",
            " >ref|XM_712142.1|/13",
        ],
        [">Cluster 2", " >emb|FM992691.1|/123", " >ref|XM_002420030.1|"],
        [">Cluster 3", " >gb|AF184907.1|/1233"],
        [">Cluster 4", " >gb|M90812.1|YSAPHR1"],
        [">Cluster 5", " >gb|CP032014.1|/1100"],
        [">Cluster 6", " >gb|CP115801.1|/1882"],
    ]
    obs = obj.sub_parse_100_clstr("tests/merge_fasta_example_rmgap_100.fa.clstr")
    assert pred == obs


def test_get_fa_recs(obj):
    """partial match 식으로 가지고 오게 하니까 중복된 것들을 너무 많이 가지고 온다.
    그래서 exact match 식으로 바꾸었음

    Args:
        obj (_type_): _description_
    """
    obs = [
        ">Cluster 0",
        " >gb|CP025153.1|/11001",
        " >gb|CP025170.1|/11001",
    ]  # , ' >gb|CP017626.1|/1100', ' >gb|AF247190.1|/1626', ' >ref|XM_712140.2|/15']
    pred = [
        (
            "TGAAAAAGCATAAGAAGGAAAAGAAAGACAAAAAGGAGAAGAAGGAAAAGAAGAAAAAGAAGGAAAAGAAGG"
            "ACAAGAAGGACAAGAAGGACAAGAAAAACAAGAAAGACAAGAAAGACAAGAAGGACAAGAAGGACAAGAAGG"
            "ACAAAATACGAACTGGTTCAGATGAAACACTAGTATCAAAAGAATCTAGTGCTACTCCGCCCCCAATTGCCA"
            "CTCGACTTTCAGCTAGATCAAAGTGGATTAAGCAAAAAAGGGCATCAGTGATGGATTCAAAAGCATTAAACG"
            "AAATTTTCATGATCACCAATTGAGAAACTTAGGGATGGTATAAATAAACTAGGACAATTTGTAATGCATGTA"
            "CTACAATAGATAGAGAATAGATATGTGTACACAAGTTGATTTCAAATTGGAACTTTTAAATAATGAACTCTT"
            "ATTAAACTATATATATTCGTCCTATACAACAGAACCCTTCATATCGGTTAGTCTTCATTGATTTAAAAAACA"
            "ACGGACATACCACCAACAAATGCAGTAACAATAGTAATAATTGAAACCAATTTGACCATAGACATTTGTTGA"
            "GTTGCTTTAACTCCAGAGCTTGAGCTGGACCCAGAGCTGGTGCTGCTGCTTGATGATCCAGAAGCAGATGCA"
            "GAGGAAGATCCAGATTTGGAGCTTCCACCAGAGGAAGATGCCTTACCAGTAGCTACTCCACTAACAGCTTTG"
            "CAGCTGCCACTAGCACTAGCCTTGCTATTGATTGAAGCACTGCCTTTGAAATCACAAGCACTGGATTTCTTG"
            "TCTTGGTCAAGGTAATACTGGTTCAACACATATGACAAACGATCTTTATCAGAACAGAAGGAAGCAACACCA"
            "TATTCACCTTTGCTACCATCGGCAGAAATAGCCGAGCAATCGATATAACCACAAACTTGACCAAACAAGTCA"
            "CCATAGTCTTCTTTATCAACATCGTCAGCAACAACACATTCTAATGACTTGGAAATACAATCACA"
        ),
        (
            "TGAAAAAGCATAAGAAGGAAAAGAAAGACAAAAAGGAGAAGAAGGAAAAGAAGAAAAAGAAGGAAAAGAAGG"
            "ACAAGAAGGACAAGAAGGACAAGAAAAACAAGAAAGACAAGAAAGACAAGAAGGACAAGAAGGACAAGAAGG"
            "ACAAAATACGAACTGGTTCAGATGAAACACTAGTATCAAAAGAATCTAGTGCTACTCCGCCCCCAATTGCCA"
            "CTCGACTTTCAGCTAGATCAAAGTGGATTAAGCAAAAAAGGGCATCAGTGATGGATTCAAAAGCATTAAACG"
            "AAATTTTCATGATCACCAATTGAGAAACTTAGGGATGGTATAAATAAACTAGGACAATTTGTAATGCATGTA"
            "CTACAATAGATAGAGAATAGATATGTGTACACAAGTTGATTTCAAATTGGAACTTTTAAATAATGAACTCTT"
            "ATTAAACTATATATATTCGTCCTATACAACAGAACCCTTCATATCGGTTAGTCTTCATTGATTTAAAAAACA"
            "ACGGACATACCACCAACAAATGCAGTAACAATAGTAATAATTGAAACCAATTTGACCATAGACATTTGTTGA"
            "GTTGCTTTAACTCCAGAGCTTGAGCTGGACCCAGAGCTGGTGCTGCTGCTTGATGATCCAGAAGCAGATGCA"
            "GAGGAAGATCCAGATTTGGAGCTTCCACCAGAGGAAGATGCCTTACCAGTAGCTACTCCACTAACAGCTTTG"
            "CAGCTGCCACTAGCACTAGCCTTGCTATTGATTGAAGCACTGCCTTTGAAATCACAAGCACTGGATTTCTTG"
            "TCTTGGTCAAGGTAATACTGGTTCAACACATATGACAAACGATCTTTATCAGAACAGAAGGAAGCAACACCA"
            "TATTCACCTTTGCTACCATCGGCAGAAATAGCCGAGCAATCGATATAACCACAAACTTGACCAAACAAGTCA"
            "CCATAGTCTTCTTTATCAACATCGTCAGCAACAACACATTCTAATGACTTGGAAATACAATCACA"
        ),
    ]
    assert pred == obj.get_fa_recs("tests/merge_fasta_example_rmgap.fa", obs)


def test_get_faid_recs(obj):
    """입력된 cd-hit 의 cluster 파일에 해당하는 fasta 파일의 id 들을 다 가지고 온다.
    출력 형식은 csv 파일 형식으로
    <cluster 번호>,<id들의 갯수>,<id1>,<id2>,...
    """
    obs = [
        ">Cluster 0",
        " >gb|CP025153.1|/1100",
        " >gb|CP025170.1|/1100",
    ]  # , ' >gb|CP017626.1|/1100', ' >gb|AF247190.1|/1626', ' >ref|XM_712140.2|/15']
    pred = ",".join(["Cluster 0", "2", "gb|CP025153.1|/11001", "gb|CP025170.1|/11001"])
    assert pred == obj.get_faid_recs("tests/merge_fasta_example_rmgap.fa", obs)
    obs = [">Cluster 0", " >MT108470.1:1-521", " >MT108469.1:1-521"]
    pred = ",".join(
        [
            "Cluster 0",
            "2",
            "MT108470.1:1-521:Aspergillus_fumigatus_strain_CNRMA17.179_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds",
            "MT108469.1:1-521:Aspergillus_fumigatus_strain_CNRMA16.56_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds",
        ]
    )
    assert pred == obj.get_faid_recs(StringIO(MockFasta), obs)


def test_recover_faid_recs(obj):
    """unique id dict 에 따라서 id 들을 바꾼다
    fasta 파일은 mk_faheader_unique 에 의해서 이미 seqNNN 으로 바뀐 상태일 때 작동한다.
    """
    obs = [">Cluster 0", " >seq1", " >seq2"]
    obs_uniq_dict = {
        "seq1": "MT108470.1:1-521:Aspergillus_fumigatus_strain_CNRMA17.179_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds",
        "seq2": "MT108469.1:1-521:Aspergillus_fumigatus_strain_CNRMA16.56_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds",
    }
    pred = ",".join(
        [
            "Cluster 0",
            "2",
            "MT108470.1:1-521:Aspergillus_fumigatus_strain_CNRMA17.179_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds",
            "MT108469.1:1-521:Aspergillus_fumigatus_strain_CNRMA16.56_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds",
        ]
    )
    assert pred == obj.recover_faid_recs(StringIO(MockFasta2), obs, obs_uniq_dict)


def test_save_faid_recs(obj):
    """get_faid_recs 의 결과물을 저장한다."""
    obs_in = [
        [
            ">Cluster 0",
            " >gb|CP025153.1|/1100",
            " >gb|CP025170.1|/1100",
            " >gb|CP017626.1|/1100",
            " >gb|AF247190.1|/1626",
            " >ref|XM_712140.2|/15",
        ],
        [
            ">Cluster 1",
            " >gb|CP025161.1|/1100",
            " >gb|CP025178.1|/1100",
            " >gb|AF247189.1|/1626",
            " >ref|XM_712142.1|/13",
        ],
        [">Cluster 2", " >emb|FM992691.1|/123", " >ref|XM_002420030.1|"],
        [">Cluster 3", " >gb|AF184907.1|/1233"],
        [">Cluster 4", " >gb|M90812.1|YSAPHR1"],
        [">Cluster 5", " >gb|CP032014.1|/1100"],
        [">Cluster 6", " >gb|CP115801.1|/1882"],
    ]
    pred = "tests/merge_fasta_example_rmgap_id.csv"
    obs = obj.save_faid_recs("tests/merge_fasta_example_rmgap.fa", obs_in)
    assert pred == obs


def test_local_alignment(obj):
    """두 개의 100% 동일한 서열들을 하나로 합친다.
    ATGCGGG
     TGCGGGT
    ---->
    ATGCGGGT
    의 첫번째 단계: local pairwise alignment
    """
    obs = ["ATGCGGG", "TGCGGGT"]
    pred = ["ATGCGGG", ".TGCGGGT"]
    assert pred == obj.local_alignment(obs)


def test_get_consensus_seq(obj):
    """두 개의 100% 동일한 서열들을 하나로 합친다.
    ATGCGGG
     TGCGGGT
    ---->
    ATGCGGGT
    의 두번째 단계: consensus sequence
    """
    obs = ["ATGCGGG", ".TGCGGGT"]
    pred = "ATGCGGGT"
    assert pred == obj.get_consensus_seq(obs)


def test_get_consensus_seq_main(obj):
    """세 개 이상의 동일한 서열들을 하나로 합하고 consensus seq 를 구한다.
    TATGCGGG
    .ATGCGGG
    ..TGCGGGT
    .....GGGTAAA
    ---->
    TATGCGGGTAAA
    """
    obs = ["TATGCGGG", "ATGCGGG", "TGCGGGT", "GGGTAAA"]
    pred = "TATGCGGGTAAA"
    assert pred == obj.get_consensus_seq_main(obs)


def test_mk_faheader_unique(obj, mocker):
    # mocker.patch("uuid.uuid4", return_value="mocked_uuid")
    obs = "tests/merge_fasta_example.fa"
    pred = (
        obs,
        {
            "seq1": "gb|CP025161.1|/11001",
            "seq10": "gb|M90812.1|YSAPHR1A/1625",
            "seq11": "ref|XM_712140.2|/1507",
            "seq12": "ref|XM_712142.1|/1311",
            "seq13": "gb|AF184907.1|/1233",
            "seq14": "emb|FM992691.1|/1233",
            "seq15": "ref|XM_002420030.1|/1146",
            "seq2": "gb|CP025178.1|/11001",
            "seq3": "gb|CP025153.1|/11001",
            "seq4": "gb|CP025170.1|/11001",
            "seq5": "gb|CP017626.1|/11001",
            "seq6": "gb|CP032014.1|/11001",
            "seq7": "gb|CP115801.1|/1882",
            "seq8": "gb|AF247189.1|/1626",
            "seq9": "gb|AF247190.1|/1626",
        },
    )
    assert pred == obj.mK_faheader_unique(obs)


def test_recover_clusterid_to_orifaid(obj):
    obs = [[">Cluster 0", " >seq1", " >seq2"], [">Cluster 1", " >seq3", " >seq4"]]
    obs_uniq_dict = {
        "seq1": "MT108470.1:1-521:Aspergillus_fumigatus_strain_CNRMA17.179_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds",
        "seq2": "MT108469.1:1-521:Aspergillus_fumigatus_strain_CNRMA16.56_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds",
        "seq3": "ABC12345.1:1-100:hello",
        "seq4": "ABC6789.1:101-200:world",
    }
    pred = [
        [
            ">Cluster 0",
            " >MT108470.1:1-521:Aspergillus_fumigatus_strain_CNRMA17.179_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds",
            " >MT108469.1:1-521:Aspergillus_fumigatus_strain_CNRMA16.56_14-alpha_sterol_demethylase_(cyp51A)_gene,_complete_cds",
        ],
        [">Cluster 1", " >ABC12345.1:1-100:hello", " >ABC6789.1:101-200:world"],
    ]
    assert pred == obj.recover_clusterid_to_orifaid(obs, obs_uniq_dict)


def test_recover_faheader(obj):
    obs = (
        "tests/merge_fasta_example.fa",
        {
            "seq1": "gb|CP025161.1|/11001",
            "seq10": "gb|M90812.1|YSAPHR1A/1625",
            "seq11": "ref|XM_712140.2|/1507",
            "seq12": "ref|XM_712142.1|/1311",
            "seq13": "gb|AF184907.1|/1233",
            "seq14": "emb|FM992691.1|/1233",
            "seq15": "ref|XM_002420030.1|/1146",
            "seq2": "gb|CP025178.1|/11001",
            "seq3": "gb|CP025153.1|/11001",
            "seq4": "gb|CP025170.1|/11001",
            "seq5": "gb|CP017626.1|/11001",
            "seq6": "gb|CP032014.1|/11001",
            "seq7": "gb|CP115801.1|/1882",
            "seq8": "gb|AF247189.1|/1626",
            "seq9": "gb|AF247190.1|/1626",
        },
    )
    pred = "tests/merge_fasta_example.fa"
    assert pred == obj.recover_faheader(obs)


def test_set_merge_count(obj):
    """각 faid 뒤에 클러스터내의 merge 숫자들을 표시한다."""
    obs = [
        [
            ">Cluster 0",
            " >gb|CP025153.1|/1100",
            " >gb|CP025170.1|/1100",
            " >gb|CP017626.1|/1100",
            " >gb|AF247190.1|/1626",
            " >ref|XM_712140.2|/15",
        ],
        [
            ">Cluster 1",
            " >gb|CP025161.1|/1100",
            " >gb|CP025178.1|/1100",
            " >gb|AF247189.1|/1626",
            " >ref|XM_712142.1|/13",
        ],
        [">Cluster 2", " >emb|FM992691.1|/123", " >ref|XM_002420030.1|"],
        [">Cluster 3", " >gb|AF184907.1|/1233"],
        [">Cluster 4", " >gb|M90812.1|YSAPHR1"],
        [">Cluster 5", " >gb|CP032014.1|/1100"],
        [">Cluster 6", " >gb|CP115801.1|/1882"],
    ]
    pred = {
        "emb|FM992691.1|/123": "emb|FM992691.1|/123;MC:2",
        "gb|AF184907.1|/1233": "gb|AF184907.1|/1233;MC:1",
        "gb|AF247189.1|/1626": "gb|AF247189.1|/1626;MC:4",
        "gb|AF247190.1|/1626": "gb|AF247190.1|/1626;MC:5",
        "gb|CP017626.1|/1100": "gb|CP017626.1|/1100;MC:5",
        "gb|CP025153.1|/1100": "gb|CP025153.1|/1100;MC:5",
        "gb|CP025161.1|/1100": "gb|CP025161.1|/1100;MC:4",
        "gb|CP025170.1|/1100": "gb|CP025170.1|/1100;MC:5",
        "gb|CP025178.1|/1100": "gb|CP025178.1|/1100;MC:4",
        "gb|CP032014.1|/1100": "gb|CP032014.1|/1100;MC:1",
        "gb|CP115801.1|/1882": "gb|CP115801.1|/1882;MC:1",
        "gb|M90812.1|YSAPHR1": "gb|M90812.1|YSAPHR1;MC:1",
        "ref|XM_002420030.1|": "ref|XM_002420030.1|;MC:2",
        "ref|XM_712140.2|/15": "ref|XM_712140.2|/15;MC:5",
        "ref|XM_712142.1|/13": "ref|XM_712142.1|/13;MC:4",
    }
    assert pred == obj.set_merge_count(obs)


def test_set_merge_count2(obj):
    obs = [
        [
            ">Cluster 0",
            " >gb|CP025153.1|/1100",
            " >gb|CP025170.1|/1100;MC:2",
            " >gb|CP017626.1|/1100",
            " >gb|AF247190.1|/1626",
            " >ref|XM_712140.2|/15",
        ],
        [
            ">Cluster 1",
            " >gb|CP025161.1|/1100",
            " >gb|CP025178.1|/1100",
            " >gb|AF247189.1|/1626;MC:4",
            " >ref|XM_712142.1|/13",
        ],
    ]
    pred = {
        "gb|AF247189.1|/1626;MC:4": "gb|AF247189.1|/1626;MC:7",
        "gb|AF247190.1|/1626": "gb|AF247190.1|/1626;MC:6",
        "gb|CP017626.1|/1100": "gb|CP017626.1|/1100;MC:6",
        "gb|CP025153.1|/1100": "gb|CP025153.1|/1100;MC:6",
        "gb|CP025161.1|/1100": "gb|CP025161.1|/1100;MC:7",
        "gb|CP025170.1|/1100;MC:2": "gb|CP025170.1|/1100;MC:6",
        "gb|CP025178.1|/1100": "gb|CP025178.1|/1100;MC:7",
        "ref|XM_712140.2|/15": "ref|XM_712140.2|/15;MC:6",
        "ref|XM_712142.1|/13": "ref|XM_712142.1|/13;MC:7",
    }
    assert pred == obj.set_merge_count(obs)


def test_parse_mc_info(obj):
    obs = " >ref|XM_712142.1|/13;MC:7"
    pred = 7
    assert pred == obj.parse_mc_info(obs)
    obs = ">refaasfasdfasfasfa"
    pred = 0
    assert pred == obj.parse_mc_info(obs)


def test_change_mc_info(obj):
    obs = " >ref|XM_712142.1|/13;MC:7"
    obs_cnt = 9
    pred = " >ref|XM_712142.1|/13;MC:9"
    assert pred == obj.change_mc_info(obs, obs_cnt)
    obs = ">refaasfasdfasfasfa"
    obs_cnt = 1
    pred = ">refaasfasdfasfasfa;MC:1"
    assert pred == obj.change_mc_info(obs, obs_cnt)
    obs = " >gb|CP025153.1|/1100"
    obs_cnt = 6
    pred = " >gb|CP025153.1|/1100;MC:6"
    assert pred == obj.change_mc_info(obs, obs_cnt)


def test_clusterone_merge_count(obj):
    """한 개 클러스터 안에서 ;MC:nnn 이 있을 때 ;MC:nnn 이 없는 faid 들의 갯수를 더한다.
    MC:nnn 이 없을 때는 그냥 count 한다."""
    obs = [
        ">Cluster 0",
        " >gb|CP025153.1|/1100",
        " >gb|CP025170.1|/1100;MC:2",
        " >gb|CP017626.1|/1100",
        " >gb|AF247190.1|/1626",
        " >ref|XM_712140.2|/15",
    ]
    pred = [
        ">Cluster 0",
        " >gb|CP025153.1|/1100;MC:6",
        " >gb|CP025170.1|/1100;MC:6",
        " >gb|CP017626.1|/1100;MC:6",
        " >gb|AF247190.1|/1626;MC:6",
        " >ref|XM_712140.2|/15;MC:6",
    ]
    assert pred == obj.clusterone_merge_count(obs)


def test_clusterone_merge_count0(obj):
    obs = [
        ">Cluster 0",
        " >gb|CP025153.1|/1100",
        " >gb|CP025170.1|/1100",
        " >gb|CP017626.1|/1100",
        " >gb|AF247190.1|/1626",
        " >ref|XM_712140.2|/15",
    ]
    pred = [
        ">Cluster 0",
        " >gb|CP025153.1|/1100;MC:5",
        " >gb|CP025170.1|/1100;MC:5",
        " >gb|CP017626.1|/1100;MC:5",
        " >gb|AF247190.1|/1626;MC:5",
        " >ref|XM_712140.2|/15;MC:5",
    ]
    assert pred == obj.clusterone_merge_count(obs)
