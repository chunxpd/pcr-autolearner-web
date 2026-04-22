import Bio.Blast.NCBIXML as bxml


class Parser:
    def __init__(self, xmlfn):
        self.xmlfn = xmlfn
        self.blast_records = self.get_blastrecords()

    def get_blastrecords(self):
        """get_blastrecords()
        XML 파일에 xml tag 가 아닌 다른 문자열이 들어있으면 에러 발생.
        가끔 CREATE_VIEW 가 들어있을 때가 있음.
        이러한 에러를 발생시키는 문자열을 제거하기 위해 임시로 `_correction`이라는 파일을 만들고
        다시 xml 파일을 출력함
        """
        tmp_o_fn = f"{self.xmlfn}_correction"
        with open(self.xmlfn, "r") as ihandle, open(tmp_o_fn, "w") as ohandle:
            for line in ihandle:
                if line.strip() in ["CREATE_VIEW"]:
                    continue
                ohandle.write(line)
        return bxml.parse(open(tmp_o_fn, "r"))

    def get_hsps(self):
        res = []
        for blast_record in self.blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    hspcont = {
                        "title": alignment.title,
                        "length": alignment.length,
                        "E": hsp.expect,
                        "query": hsp.query,
                        "match": hsp.match,
                        "sbjct": hsp.sbjct,
                        "query_len": blast_record.query_length,
                        "aln_len": hsp.align_length,
                        "identity": hsp.identities,
                        "aln_pct": round(
                            100 * hsp.align_length / blast_record.query_length, 2
                        ),
                        "id_pct": round(
                            100 * hsp.identities / blast_record.query_length, 2
                        ),
                    }
                    res.append(hspcont)
        return res

    def gen_text(self, hsps: list) -> str:
        res = "<table class='table table-fixed table-striped'>\n"
        res += "<thead><tr>\n"
        for kkk in ["title", "query_len", "aln_pct", "id_pct", "E", "match"]:
            res += f"<th>{kkk}</th>"
        # res += "".join(f"<th>{v}</th>" for v in hsps[0].keys())
        res += "\n</tr>\n</thead>\n"
        res += "<tbody>\n"
        for hsp in hsps:
            res += "<tr>\n"
            for ele in ["title", "query_len", "aln_pct", "id_pct", "E"]:
                res += f"<td>{hsp[ele]}</td>\n"
            res += f"<td>{hsp['query']}<br>{hsp['match']}<br>{hsp['sbjct']}</td>\n"
            res += "</tr>\n"
        res += "</tbody>\n</table>\n"
        return res


if __name__ == "__main__":
    pass
