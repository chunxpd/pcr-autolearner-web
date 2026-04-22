import glob
import logging
import os
import subprocess, shlex

from Bio import SeqIO

from src.common import Common


class CDHIT:
    """
    진행 순서는 다음과 같다.
        1. cd-hit -c 1.0 을 해서 100% 를 만든다
        2. 100% 를 대상으로 -c 0.4 를 한다.
        3. 2의 clstr 파일에서 가장 큰 클러스터를 제외한 나머지 클러스터의 record 들을 1의 결과에서 제거한다.
        예)::
        cd-hit -i foo.fa -o foo100.fa -c 1.0 -n 5
        cd-hit -i foo100.fa -o foo40.fa -c 0.4 -n 2 -sc 1 -sf 1
    """

    def __init__(self):
        pass

    def _fa_replace(self, in_fa: str, tar_suffix: str) -> str:
        """_summary_

        Args:
            in_fa (str): fasta 파일이름.
            tar_suffix (str): 교체하려는 확장자

        Returns:
            str: 결과물
        """
        if in_fa.rfind(".fasta") >= 0:
            res = Common.suffix_replace(in_fa, ".fasta", tar_suffix)
        elif in_fa.rfind(".fa") >= 0:
            res = Common.suffix_replace(in_fa, ".fa", tar_suffix)
        return res

    def get_cd_hit100_fn(self, in_fa: str, type="cluster"):
        if type == "cluster":
            return self._fa_replace(in_fa, "_100.fa.clstr")
        else:
            return self._fa_replace(in_fa, "_100.fa")

    def get_cd_hitlb_fn(self, in_fa: str, type="cluster"):
        if type == "cluster":
            return self._fa_replace(in_fa, "_lb.fa.clstr")
        else:
            return self._fa_replace(in_fa, "_lb.fa")

    def get_cd_hit_filtered_fn(self, in_fa: str) -> str:
        """cd-hit 100 ~ cd-hit lower bound 까지를 포함한다는 이름"""
        return self._fa_replace(in_fa, "_100_lb.fa")

    def run_cdhit(self, in_fa: str, lower_bound: float = 0.4) -> bool:
        """run_cdhit(self, in_fa: str, lower_boudn: float = 0.4) -> bool:

        word size 는 threshold 에 따라 다르다. 표준설명서에는 다음과 같다.
        -n 5 for thresholds 0.7 ~ 1.0
        -n 4 for thresholds 0.6 ~ 0.7
        -n 3 for thresholds 0.5 ~ 0.6

        lower_bound 는 0.4 보다 작으면 시스템 폴트가 난다. 0.4 가 최소인 듯 함
        Args:
            in_fa (str): _description_
            lower_boudn (float, optional): _description_. Defaults to 0.4.
        """
        in_fa_100 = self.get_cd_hit100_fn(in_fa, type="fasta")
        in_fa_lb = self.get_cd_hitlb_fn(in_fa, type="fasta")

        lower_bound = max(lower_bound, 0.4)
        word_size = self.get_word_size(lower_bound)

        exec_str100 = f"cd-hit -i {in_fa} -o {in_fa_100} -c 1.0 -n 5 -sc 1 -sf 1 -M 0"
        exec_str_lb = f"cd-hit -i {in_fa_100} -o {in_fa_lb} -c {lower_bound} -n {word_size} -sc 1 -sf 1 -M 0"
        exec_str100_res = subprocess.run(
            shlex.split(exec_str100),
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
        )
        exec_str_lb_res = subprocess.run(
            shlex.split(exec_str_lb),
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
        )
        # exec_str100_res = os.system(exec_str100)
        # exec_str_lb_res = os.system(exec_str_lb)
        if exec_str100_res.returncode > 0 or exec_str_lb_res.returncode > 0:
            raise Exception("cd-hit running error")
        res = [
            in_fa,
            in_fa_100,
            self.get_cd_hit100_fn(in_fa, type="cluster"),
            in_fa_lb,
            self.get_cd_hitlb_fn(in_fa, type="cluster"),
        ]
        return res

    def get_word_size(self, threshold: float) -> int:
        word_size = 5
        if threshold >= 0.7:
            word_size = 5
        elif threshold >= 0.6:
            word_size = 4
        elif threshold >= 0.5:
            word_size = 3
        elif threshold >= 0.4:
            word_size = 2
        return word_size

    def get_cd_hit_res_lists(self, in_fa: str) -> list:
        basename = in_fa.replace(".fa", "")
        res = glob.glob(f"{basename}*.fa*")
        return res

    def parse_lb_clstr(self, lb_clstr: str) -> list:
        """parse_lb_clstr(self, lb_clstr: str) -> list
        lower bound cluster file 을 파싱해서 fasta 파일의 헤드들을 가지고 온다.
        가장 큰 규모의 클러스터는 빼고 나머지에서만 가지고 온다.

        Args:
            lb_clstr (str): cd-hit 결과물인 cluster 파일

        Returns:
            list: fa 의 레코드의 헤드들
        """
        res = []
        clstr_recs = self.sub_parse_lb_clstr(lb_clstr)
        for iii in range(1, len(clstr_recs)):
            cluster = clstr_recs[iii]
            if cluster[0] == ">Cluster 0":
                continue
            for rec in cluster:
                if rec[0] == ">":
                    continue
                rec_head = rec.split(",")[1].split("...")[0]
                res.append(rec_head.replace(" ", "")[1:])
        return res

    def sub_parse_lb_clstr(self, lb_clstr: str) -> list:
        res = []
        with open(lb_clstr, "r") as handle:
            for lll in handle:
                if lll[0] == ">":
                    try:
                        res.append(tmpcont)
                    except:
                        pass
                    tmpcont = [lll]
                    continue
                tmpcont.append(lll)
            res.append(tmpcont)
        return res

    def exclude_fa_file(self, total_heads: list, exc_heads: list) -> list:
        res = []
        for rrr in total_heads:
            if self.partial_matchq(exc_heads, rrr):  # rrr in exc_heads:
                continue
            res.append(rrr)
        return res

    def partial_matchq(self, ori_list: list, tar_str: str):
        """tar_str 이 ori_list 의 element 들과 하나씩 비교하여 일부이면 True 반환"""
        for o_str in ori_list:
            if tar_str in o_str or o_str in tar_str:
                return True
        return False

    def exact_matchq(self, ori_list: list, tar_str: str):
        """tar_str 이 ori_list 의 element 들과 하나씩 비교하여 정확히 일치하는 것이 있으면 True 반환"""
        for o_str in ori_list:
            if tar_str == o_str:
                return True
        return False

    def read_fa_file(self, in_fa: str):
        res = []
        for record in SeqIO.parse(in_fa, "fasta"):
            res.append(record.id)
        return res

    def write_fa_file(self, out_fa: str, in_fa: str = "", exc_list: list = []) -> int:
        all_heads = self.read_fa_file(in_fa)
        filtered_heads = self.exclude_fa_file(all_heads, exc_list)
        with open(out_fa, "w") as out_handle, open(in_fa, "r") as inhandle:
            cont = []
            for record in SeqIO.parse(inhandle, "fasta"):
                if record.id not in filtered_heads:
                    continue
                cont.append(record)
            res = SeqIO.write(cont, out_handle, "fasta")
        return res

    def concat_files(self, config, jobid, file_list) -> str:
        """입력된 fasta file 들을 모두 concat 하여 jobid.fa 로 만들고
        아웃풋 파일 이름을 반환한다.

        Args:
            config (_type_): _description_
            jobid (_type_): _description_
            file_list (_type_): _description_

        Returns:
            str: _description_
        """
        tempdir = Common.get_upload_path(config, jobid)
        fff = [f"{Common.get_upload_file_path(config,jobid,v)}" for v in file_list]
        ofa = Common.get_upload_file_path(config, jobid, jobid + ".fa")
        exec_str = f"cat {' '.join(fff)} > {ofa}"
        os.system(exec_str)
        return ofa

    def main(self, in_fa: str, lower_bound: float):
        """
        생성되는 파일
        in_fasta 에서 부터 run_cdhit 으로부터 아래 파일들이 자동 생성
        {in_fasta%.fa}100.fa
        {in_fasta%.fa}100.fa.clstr
        {in_fasta%.fa}lb.fa
        {in_fasta%.fa}lb.fa.clstr
        최종 파일은
        {in_fasta%.fa}100_lb.fa # cd-hit 100 ~ lower bound 까지만 남겼다는 뜻
        """
        res = {"res": None, "count": None}
        _ = self.run_cdhit(in_fa, lower_bound=lower_bound)
        exc_list = self.parse_lb_clstr(self.get_cd_hitlb_fn(in_fa, type="cluster"))
        final_file = self.get_cd_hit_filtered_fn(in_fa)
        cdhit100fn = self.get_cd_hit100_fn(in_fa, type="fasta")
        count = self.write_fa_file(final_file, in_fa=cdhit100fn, exc_list=exc_list)
        res["res"] = final_file
        res["count"] = count
        return res


if __name__ == "__main__":
    import sys

    logging.basicConfig(
        level=logging.INFO,
        datefmt="%m/%d/%Y %I:%M:%S %p",
        format="%(asctime)s:%(levelname)s:%(message)s",
    )

    runner = CDHIT()
    runner.main(sys.argv[1], 0.4)
