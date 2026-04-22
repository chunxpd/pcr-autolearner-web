from src.parse_blastxml import Parser


class Score:
    def __init__(self):
        pass

    def preprocessing(self, hsps: list) -> dict:
        """입력된 HSP 입력 값에서 같은 title 을 가지고 있는 것들의 경우 가장
        낮은 E-value 를 보이는 것 하나만 남기고 지우기.
        같은 E-value 이면 제일 처음 것만 남기기

        Args:
            hsps (list): parse_blastxml.get_hsps() 의 반환 값
                Ex::
                [
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
                    },...
                ]

        Returns:
            dict: preprocessing 결과. title 이 key, hsp record 가 value 인 dictionary
                Ex::
                {
                    (
                        "gi|1858710525|gb|CP054795.1| Mycolicibacterium smegmatis strain "
                        "FDAARGOS_679 chromosome"
                    ): {
                        "E": 0.000183945,
                        "query_len": 21,
                        "aln_len": 21,
                        "identity": 21,
                        "aln_pct": 100.0,
                        "id_pct": 100.0,
                    },
                    ...
                }
        """
        res = {}
        for rec in hsps:
            if rec["title"] in res:
                res[rec["title"]].append(rec)
            else:
                res[rec["title"]] = [rec]
        for kkk, vvv in res.items():
            res[kkk] = sorted(vvv, key=lambda x: x["E"])[0]
        return res

    def combine_hsp_sets(self, *hspses) -> dict:
        """불특정 갯수의 preprocess 된 hsps dict 를 입력 받아서
        gb id 아래에 필요한 정보를 반환하는 기능

        Aarguments:
            *hspses (dict): 1개 이상의 preprocess 된 hsp dict 들
        Returns:
            dict: _description_
        """
        res = dict()
        keys = []
        for hsps_dic in hspses:
            keys.extend(hsps_dic.keys())
        ukeys = list(set(keys))
        for kkk in ukeys:
            res[kkk] = [] if kkk not in res else res[kkk]
            for hsps_dic in hspses:
                if kkk in hsps_dic:
                    res[kkk].append(hsps_dic[kkk])
                else:
                    res[kkk].append({})
        return res

    def mk_csv_file(self, combine_set: dict, delim="\t") -> str:
        """self.combine_hsp_sets 의 결과를 csv 파일로 만든다.
        헤더를 붙인다.

        Args:
            combine_set (dict): 각 GB id 별 입력 blastn xml 파일의 HSP 들의 결과값들

        Returns:
            str: cvs 문자열
        """
        res = list()
        out_keys = ["aln_pct", "id_pct", "E"]  # for print out Hsp attributes
        for gbid, recs in combine_set.items():
            tmpcont = [gbid]
            for rec in recs:
                if rec == {}:
                    tmpcont.extend([""] * len(out_keys))
                    continue
                for kkk in out_keys:
                    tmpcont.append(f"{rec[kkk]}")
            res.append(delim.join(tmpcont))
        res.insert(0, self._mk_header(len(recs), delim))
        return "\n".join(res)

    def _mk_header(self, file_num: int, delim: str) -> str:
        header = ["id"]
        header.extend(
            [
                delim.join([f"aln_pct{i}", f"id_pct{i}", f"E{i}"])
                for i in range(1, file_num + 1)
            ]
        )
        header = delim.join(header)
        return header


if __name__ == "__main__":
    import sys

    runner = Score()
    hspses = []
    for fn in sys.argv[1:]:
        parser = Parser(fn)
        prc_hsps = runner.preprocessing(parser.get_hsps())
        hspses.append(prc_hsps)
    cmb_hspes = runner.combine_hsp_sets(*hspses)
    outstr = runner.mk_csv_file(cmb_hspes)
    sys.stdout.write(outstr)
