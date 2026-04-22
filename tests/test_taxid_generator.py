import pytest
from src.taxid_generator import TxGenerator


@pytest.fixture(scope="module")
def obj():
    config = {
        "TAXORIPATH": "data/ori/fullnamelineage.dmp",
        "TAXJSPATH": "data/fullnamelineage.js",
        "TAXDBPATH": "data/fullnamelineage.db",
        "AUTOCOMPLETESIZE": 5,
    }
    res = TxGenerator(config)
    yield res
    del res


def test_get_taxinfo(obj):
    pred = "root (taxid: 1)"
    obs = obj.get_taxinfo()
    assert pred == obs[0]
    pred = "cellular organisms (taxid: 131567)"
    assert pred == obs[1]


def test_get_taxinfo_query(obj):
    pred = ["Mycobacterium phage Keshu (taxid: 1567471)"]
    obs = obj.get_taxinfo_query("mycobacterium phage Keshu")
    assert pred == obs
    ### 결과 값은 딱 5개만 가지고 온다.
    pred = [
        "Mycobacterium phage Pinkman (taxid: 1897426)",
        "Mycobacterium phage PinheadLarry (taxid: 2283299)",
        "Mycobacterium phage Pinkcreek (taxid: 2891199)",
        "Mycobacterium phage PinkPlastic (taxid: 2502451)",
        "Mycobacterium phage Pinto (taxid: 1498204)",
    ]
    ### 대소문자 상관 없이 같은 결과를 가지고 온다
    for taxname in [
        "mycobacterium phage Pin",
        "mycobacterium phage pin",
        "Mycobacterium phage Pin",
    ]:
        obs = obj.get_taxinfo_query(taxname)
        assert pred == obs
