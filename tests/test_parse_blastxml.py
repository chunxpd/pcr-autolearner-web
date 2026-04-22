import pytest
import os
from src.parse_blastxml import Parser


@pytest.fixture(scope="module")
def obj():
    cur_dir = os.path.dirname(os.path.abspath(__file__))
    test_xmlfn = f"{cur_dir}/blast_include_example.xml"
    res = Parser(test_xmlfn)
    yield res
    del res


@pytest.fixture(scope="module")
def obj_err():
    cur_dir = os.path.dirname(os.path.abspath(__file__))
    test_xmlfn = f"{cur_dir}/blast_parseerror_example.xml"
    res = Parser(test_xmlfn)
    yield res
    del res


def test_get_hsps(obj):
    pred = [
        {
            "E": 0.000183945,
            "length": 6993871,
            "match": "|||||||||||||||||||||",
            "query": "GATGCTCGCAACCACTATCCA",
            "sbjct": "GATGCTCGCAACCACTATCCA",
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
            "match": "|||||||||||||||||||||",
            "query": "GATGCTCGCAACCACTATCCA",
            "sbjct": "GATGCTCGCAACCACTATCCA",
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
    obs = obj.get_hsps()
    assert obs[:2] == pred


def test_gen_text(obj):
    pred = """<table class='table table-fixed table-striped'>
<thead><tr>
<th>title</th><th>query_len</th><th>aln_pct</th><th>id_pct</th><th>E</th><th>match</th>
</tr>
</thead>
<tbody>
<tr>
<td>gi|1858710525|gb|CP054795.1| Mycolicibacterium smegmatis strain FDAARGOS_679 chromosome</td>
<td>21</td>
<td>100.0</td>
<td>100.0</td>
<td>0.000183945</td>
<td>GATGCTCGCAACCACTATCCA<br>|||||||||||||||||||||<br>GATGCTCGCAACCACTATCCA</td>
</tr>
</tbody>
</table>
"""
    obs_rec = [
        {
            "E": 0.000183945,
            "length": 6993871,
            "match": "|||||||||||||||||||||",
            "query": "GATGCTCGCAACCACTATCCA",
            "sbjct": "GATGCTCGCAACCACTATCCA",
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
    obs = obj.gen_text(obs_rec)
    assert obs == pred


def test_parseerror(obj_err):
    pred = [
        {
            "E": 0.731715,
            "aln_len": 22,
            "aln_pct": 100.0,
            "id_pct": 100.0,
            "identity": 22,
            "length": 10972864,
            "match": "||||||||||||||||||||||",
            "query": "GAAGAACAGGTTCTCCAGCAGG",
            "query_len": 22,
            "sbjct": "GAAGAACAGGTTCTCCAGCAGG",
            "title": "gi|1860121958|gb|CP054925.1| Amycolatopsis sp. Hca4 chromosome, "
            "complete genome",
        },
    ]
    obs = obj_err.get_hsps()
    assert obs[:1] == pred
