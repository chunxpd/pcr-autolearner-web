from src.download_genebank import DwGeneBank


class StrandError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class DwGeneBank2(DwGeneBank):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def merge_two_hsps(self, hsp1: dict, hsp2: dict) -> dict:
        """서로 겹치는 Hsp 들은 합친다. strand 의 방향이 다르면 에러 발생한다.
        첫 번째
            HSP1: 100 -------- 500
        두 번째
            HSP1:        400 --------------  1000
        결과물
            HSP1: 100 ---------------------  1000

        Args:
            hsp1 (dict): gen_hsp_rec() 의 반환
            hsp2 (dict): gen_hsp_rec() 의 반환

        Returns:
            dict: 예:: {"Hsp_hit-from": 100, "Hsp_hit-strand": 1, "Hsp_hit-to": 500, "gbid": "CP101406.1"}
        """
        res = False
        start1 = hsp1["Hsp_hit-from"]
        start2 = hsp2["Hsp_hit-from"]
        end1 = hsp1["Hsp_hit-to"]
        end2 = hsp2["Hsp_hit-to"]
        index1 = end1 < start2
        index2 = start1 > end2
        if hsp1["gbid"] != hsp2["gbid"]:  # gbid 는 동일해야 함
            return res
        if index1 | index2:  # overlab 이 없을 때
            return res
        if hsp1["Hsp_hit-strand"] * hsp2["Hsp_hit-strand"] < 0:  # strand 의 방향은 일치해야 함
            raise StrandError("Strand orientations are not same")
        res = {
            "Hsp_hit-from": min(start1, start2),
            "Hsp_hit-strand": hsp1["Hsp_hit-strand"],
            "Hsp_hit-to": max(end1, end2),
            "gbid": hsp1["gbid"],
        }
        return res

    def combine_hsp_recs(self, hsp1: list, hsp2: list) -> list:
        """서로 겹치는 Hsp 들은 합치고 그렇지 않은 것들은 보존하고 리턴하는 함수

        Args:
            hsp1 (list): _description_
            hsp2 (list): _description_

        Returns:
            list: _description_
        """
        res = []
        already_used = []
        for one_hsp in hsp1:
            for two_hsp in hsp2:
                tmp_res = self.merge_two_hsps(one_hsp, two_hsp)
                if not tmp_res:
                    if not one_hsp in res and not one_hsp in already_used:
                        res.append(one_hsp)
                    if not two_hsp in res and not two_hsp in already_used:
                        res.append(two_hsp)
                    continue
                res.append(tmp_res)
                already_used.append(one_hsp)
                already_used.append(two_hsp)
        return res
