from xml.etree.ElementPath import prepare_descendant
import pytest
from src.call_blastn import CallBlast
import os.path


pytestmark = pytest.mark.skip(reason="Development suspended")


@pytest.fixture(scope="module")
def obj():
    home_path = os.path.expanduser("~")
    config = {
        "jobid": "tests_jobiid_00000",
        "TEMPDIR": "tests",
        "HITLISTSIZE": 100,
        "BLASTN": f"{home_path}/.local/bin/blastn",
        "NT": f"{home_path}/nt/nt",
        "SC_TAXSEARCH": f"{home_path}/.local/bin/get_species_taxids.sh",
    }
    res = CallBlast(config)
    yield res
    del res


def test_reverse_complement(obj):
    """strand 가 -1 이면 reverse complement 하여 반환.
    strand 가 +1 이면 아무것도 하지 않고 그대로 반환"""
    pred = "TTAAACCCCGGGG"
    obs_in_sequence = "CCCCGGGGTTTAA"
    obs = obj.reverse_complement(obs_in_sequence, strand=-1)
    assert obs == pred
    pred = "CCCCGGGGTTTAA"
    obs = obj.reverse_complement(obs_in_sequence, strand=1)
    assert obs == pred


def test_reverse_complement_multiple_nt(obj):
    pred = "GYTTBCGGTAGATGTCSAGCA"
    obs = obj.reverse_complement("TGCTSGACATCTACCGVAARC", strand=-1)
    assert obs == pred


def test_parse_taxid(obj):
    """all[filter] NOT txid1772[ORGN] NOT txid1773[ORGN]"""
    obs = "Mycobacterium (taxid: 1763), Canis lupus (taxid: 9612),"
    pred = "all[filter] NOT txid1763[ORGN] NOT txid9612[ORGN]"
    assert pred == obj.parse_taxid(obs)


def test_parse_taxid_comma(obj):
    obs = "Mycobacterium (taxid: 1763), Canis lupus (taxid: 9612),"
    pred = "1763,9612"
    assert pred == obj.parse_taxid(obs, sep=",")


def test_query_blast_standalone(obj):
    pred = "tests/tests_jobiid_00000/tests_jobiid_00000-CAGGGCATCGCGCACGCAT.xml"
    obs = obj.query_blast_standalone("CAGGGCATCGCGCACGCAT", strand=1, megablastq=False)
    assert pred == obs


def test_query_blast_exclude_taxa_standalone(obj):
    pred = "tests/tests_jobiid_00000/tests_jobiid_00000-CAGGGCATCGCGCACGCAT_excludetaxa.xml"
    obs = obj.query_blast_exclude_taxa_standalone(
        "CAGGGCATCGCGCACGCAT",
        "Mycobacterium (taxid: 1763), Canis lupus (taxid: 9612),",
        strand=1,
        megablastq=False,
    )
    assert pred == obs


def test_get_species_taxids(obj):
    pred = "tests/taxid.list"
    obs = obj.get_species_taxids(
        "Mycobacterium (taxid: 1763), Canis lupus (taxid: 9612),", "tests/taxid.list"
    )
    assert pred == obs
