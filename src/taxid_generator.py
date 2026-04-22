import os.path
import json
import sqlite3
from src.config import Config


class TxGenerator:
    def __init__(self, config: dict):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        self.orifn = f"{dir_path}/../static/{config['TAXORIPATH']}"
        self.outfn = f"{dir_path}/../static/{config['TAXJSPATH']}"
        self.dbfn = f"{dir_path}/../static/{config['TAXDBPATH']}"
        self.con = sqlite3.connect(self.dbfn, check_same_thread=False)
        self.cursor = self.con.cursor()
        self.fetch_size = config["AUTOCOMPLETESIZE"]

    def gen_js(self):
        res = []
        with open(self.orifn, "r") as handle:
            for lll in handle:
                recs = lll.split("|")
                taxid = recs[0].replace("\t", "")
                name = recs[1].replace("\t", "")
                res.append(f"{name} (taxid: {taxid})")
        with open(self.outfn, "w") as handle:
            json.dump(res, handle, indent=5)

    def gen_sqlite3(self):
        exe_str = "DROP TABLE IF EXISTS lineage "
        self.cursor.execute(exe_str)
        exe_str = "CREATE TABLE lineage(Name text, UNIQUE(Name))"
        self.cursor.execute(exe_str)
        with open(self.orifn, "r") as handle:
            for lll in handle:
                recs = lll.split("|")
                taxid = recs[0].replace("\t", "")
                name = recs[1].replace("\t", "").replace("'", "''")
                exec_str = f"INSERT INTO lineage VALUES('{name} (taxid: {taxid})')"
                self.cursor.execute(exec_str)
        self.con.commit()
        self.con.close()

    def get_taxinfo(self):
        with open(self.outfn) as handle:
            res = json.load(handle)
        return res

    def get_taxinfo_query(self, taxname: str) -> list:
        """get_taxinfo_query(self, taxname: str)-> list
        sqlite3 에 taxname 을 던져서 유사한 결과를 가지고 온다.
        대소문자를 구별하지 않는다.

        Args:
            taxname (str): taxname 부분 문자열

        Returns:
            list: sqlite3.fetchall()
        """
        exe_str = rf"SELECT * FROM lineage WHERE Name LIKE '%{taxname}%'"
        res = self.cursor.execute(exe_str)
        res = res.fetchmany(self.fetch_size)
        res = [v[0] for v in res]
        return res


if __name__ == "__main__":
    config = Config()
    txgen = TxGenerator(config.config)
    txgen.gen_js()
    txgen.gen_sqlite3()
