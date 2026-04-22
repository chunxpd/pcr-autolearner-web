import os
import re
import uuid
from src.cd_hit import CDHIT
from src.common import Common
import subprocess, shlex
from Bio.Seq import Seq
from Bio import Align
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import argparse


class MERGE_FASTA(CDHIT):
    def __init__(self):
        self.GAP = "."
        self.pattern = re.compile(r";MC:(\d+)")  # merge count pattern

    def rmgap_infasta(self, fa_name: str) -> str:
        """입력된 fasta 파일에서 gap 을 제거한다.
        반환은 <fata_prefix>_rmgap.fa 의 파일명이다.

        Args:
            fa_name (str): 입력 fasta file name

        Returns:
            str: 출력 fasta file name
        """
        res = ""
        outfn = Common.get_rmgap_fa_fn(fa_name)
        output_cont = []
        for record in SeqIO.parse(fa_name, "fasta"):
            record.seq = record.seq.replace("-", "")
            output_cont.append(record)
        SeqIO.write(output_cont, outfn, "fasta")
        return outfn

    def run_cdhit100(self, in_fa: str) -> list:
        """cd-hit 을 돌린다. threshold 는 100

        Args:
            in_fa (str): _description_

        Returns:
            bool: _description_
        """
        in_fa_100 = self.get_cd_hit100_fn(in_fa, type="fasta")
        word_size = self.get_word_size(1.0)
        exec_str100 = f"cd-hit -i {in_fa} -o {in_fa_100} -c 1.0 -n 5 -sc 1 -sf 1 -M 0"
        exec_str100_res = subprocess.run(
            shlex.split(exec_str100),
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
        )
        if exec_str100_res.returncode > 0:
            raise Exception("cd-hit running error")
        res = [in_fa, self.get_cd_hit100_fn(in_fa, type="cluster")]
        return res

    def sub_parse_100_clstr(self, clstr: str) -> list:
        """100 clstr 파일 결과를 읽어온다.

        Args:
            clstr (str): _description_

        Returns:
            list: _description_
        """
        res = []
        with open(clstr, "r") as handle:
            for lll in handle:
                if lll[0] == ">":
                    try:
                        res.append(tmpcont)
                    except:
                        pass
                    tmpcont = [lll.strip()]
                    continue
                tmpcont.append(lll.split(",")[1].split("...")[0])
            res.append(tmpcont)
        return res

    def parse_100_clstr(self, clstr: str) -> list:
        """_summary_

        Args:
            clstr (str): _description_

        Returns:
            list: _description_
        """
        res = []
        clstr_recs = self.sub_parse_100_clstr(clstr)
        for cluster in clstr_recs:
            for rec in cluster:
                if rec[0] == ">":
                    try:
                        res.append(tmpcont)
                    except:
                        pass
                    tmpcont = [rec]
                    continue
                rec_head = rec.split(",")[1].split("...")[0]
                res.append(rec_head.replace(" ", "")[1:])
        return res

    def get_fa_recs(self, fa_name: str, id_list: list) -> list:
        """fasta 파일을 열어서 id 에 있는 것들의 서열만 가지고 온다.

        Args:
            fa_name (str): fasta file name
            id_list (list): cd-hit cluster file's id list

        Returns:
            list: _description_
        """

        res = []
        id_list = [v.strip().replace(">", "") for v in id_list]
        for record in SeqIO.parse(fa_name, "fasta"):
            # if self.partial_matchq(id_list, record.id.replace(">", "")):
            if self.exact_matchq(id_list, record.id.replace(">", "")):
                res.append(str(record.seq))
        return res

    def get_faid_recs(self, fa_name: str, id_list: list) -> str:
        """입력된 cd-hit 의 cluster 파일에 해당하는 fasta 파일의 id 들을 다 가지고 온다.
        출력 형식은 csv 파일 형식으로
        <cluster 번호>,<id들의 갯수>,<id1>,<id2>,...

        Args:
            fa_name (str): fasta file name
            id_list (list): cd-hit cluster file's id list

        Returns:
            str: _description_
        """
        container = []
        id_list = [v.strip().replace(">", "") for v in id_list]
        for record in SeqIO.parse(fa_name, "fasta"):
            if self.partial_matchq(id_list, record.id.replace(">", "")):
                container.append(str(record.id))
        n_recs = len(container)
        res = [id_list[0], str(n_recs)] + container
        return ",".join(res)

    def recover_faid_recs(
        self, fa_name: str, id_list: list, unique_id_dict: dict
    ) -> str:
        """입력된 cd-hit 의 cluster 파일에 해당하는 fasta 파일의 id 들을 다 가지고 온다.
        출력 형식은 csv 파일 형식으로
        <cluster 번호>,<id들의 갯수>,<id1>,<id2>,...

        Args:
            fa_name (str): fasta file name
            id_list (list): cd-hit cluster file's id list

        Returns:
            str: _description_
        """
        container = []
        id_list = [v.strip().replace(">", "") for v in id_list]
        for record in SeqIO.parse(fa_name, "fasta"):
            # if self.partial_matchq(id_list, record.id.replace(">", "")):
            if self.exact_matchq(id_list, record.id.replace(">", "")):
                container.append(str(unique_id_dict[record.id]))
        n_recs = len(container)
        res = [id_list[0], str(n_recs)] + container
        return ",".join(res)

    def recover_clusterid_to_orifaid(self, clusters: list, uniq_id_dict) -> list:
        """cd-hit100전에 uniq id 로 변환한 상태에서의 cluster 결과를 원래의 faid 로 변환하여 반환한다.
        recover_faid_recs() 와 근본적으로는 같은 작업이다.

        Args:
            clusterid_list (list): 예) [[">Cluster 0", " >seq1", " >seq2"], [">Cluster 1",...],...]
            uniq_id_dict (_type_): 예) {
                                        "seq1": "MT108470.1:1-521:Aspergillus_fum...",
                                        "seq2": "MT108469.1:1-521:Aspergillus_fum...",
                                    }

        Returns:
            list: 예) [">Cluster 0", " >MT108470.1:...", " >MT108469.1:..."]
        """
        res = []
        for clusterid_list in clusters:
            id_list = [v.strip().replace(">", "") for v in clusterid_list[1:]]
            res.append([clusterid_list[0]] + [f" >{uniq_id_dict[v]}" for v in id_list])
        return res

    def save_faid_recs(
        self, in_fasta: str, clusterid_list: list, uniq_id_dict: dict = None
    ) -> str:
        """get_faid_recs 의 결과물을 저장한다.

        Args:
            in_fasta (str): _description_
            clusterid_list (list): _description_

        Returns:
            str: _description_
        """
        if in_fasta.endswith(".fa"):
            outfn = in_fasta.replace(".fa", "_id.csv")
        elif in_fasta.endswith(".fasta"):
            outfn = in_fasta.replace(".fasta", "_id.csv")
        with open(outfn, "w") as ohandle:
            for id_list in clusterid_list:
                if uniq_id_dict is None:
                    ohandle.write(self.get_faid_recs(in_fasta, id_list) + "\n")
                else:
                    ohandle.write(
                        self.recover_faid_recs(in_fasta, id_list, uniq_id_dict) + "\n"
                    )
        return outfn

    def local_alignment(self, in_seq_list: list) -> str:
        """두 개의 100% 동일한 서열들을 하나로 합친다.
        를 위한 1단계: local pairwise alignment

        Args:
            seq_list (list): _description_

        Returns:
            str: _description_
        """
        res = []
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -1000
        alignments = aligner.align(in_seq_list[0], in_seq_list[1])
        aln = alignments[0]
        t_segments, q_segments = aln.aligned
        if len(t_segments) == 0:
            return [in_seq_list[0], in_seq_list[1]]

        t_start, t_end = t_segments[0][0], t_segments[-1][1]
        q_start, q_end = q_segments[0][0], q_segments[-1][1]

        prefix1 = in_seq_list[0][:t_start]
        prefix2 = in_seq_list[1][:q_start]
        max_prefix_len = max(len(prefix1), len(prefix2))

        mid1 = aln[0].replace("-", self.GAP)
        mid2 = aln[1].replace("-", self.GAP)

        seq1 = (
            (self.GAP * (max_prefix_len - len(prefix1)))
            + prefix1
            + mid1
            + in_seq_list[0][t_end:]
        )
        seq2 = (
            (self.GAP * (max_prefix_len - len(prefix2)))
            + prefix2
            + mid2
            + in_seq_list[1][q_end:]
        )
        return [seq1, seq2]

    def get_consensus_seq(self, in_aln_list: list) -> str:
        """
        참조: https://biopython.org/docs/1.75/api/Bio.Align.AlignInfo.html

        Args:
            in_aln_list (list): _description_

        Returns:
            str: _description_
        """

        seq_list = []
        max_seq_len = max([len(seq) for seq in in_aln_list])  # 가장 긴 서열 길이
        in_aln_list = [
            seq + self.GAP * (max_seq_len - len(seq)) for seq in in_aln_list
        ]  # 3'에 빈칸 패딩하고 길이 동일하게 맞춤
        for seq in in_aln_list:
            seq_list.append(SeqRecord(Seq(seq)))
        alignment = MultipleSeqAlignment(seq_list)
        consensus = SummaryInfo(alignment)
        res = str(consensus.dumb_consensus(threshold=-1))
        return res

    def get_consensus_seq_main(self, in_seq_list: list) -> str:
        """세 개 이상의 서열을 입력 받아서 consensus seq 를 만든다.
        원래는 MSA 를 하고 바로 consensus seq 를 만들어야 하는데 여기서는 pairwise alinment 를 반복해서 한다.

        Args:
            in_seq_list (list): _description_

        Returns:
            str: _description_
        """
        seq0 = in_seq_list[0]
        for iii in range(1, len(in_seq_list)):
            tmp_in = seq0, in_seq_list[iii]
            tmp_aln = self.local_alignment(tmp_in)
            seq0 = self.get_consensus_seq(tmp_aln)
        return seq0

    def mK_faheader_unique(self, fa_name: str) -> tuple:
        """fasta file 의 header 가 11글자를 넘어가면 cd-hit 에서는 구별하지를 못한다.
        그리고 마지막 id 저장할 때 id 에 , 가 있으면 csv 에서 구별자이기 때문에 _ 로 바꾸어준다.
        이에 따라서
        1. 입력 fasta 파일의 id 들을 seqNNN 으로 바꾸고
        2. id 에 , 가 있으면 _ 로 바꾸어서 반환한다.

        Args:
            fa_name (str): _description_

        Returns:
            tuple: _description_
        """
        res = {}
        # outfn = fa_name + "_" + str(uuid.uuid4())
        output_cont = []
        s_num = 1
        for record in SeqIO.parse(fa_name, "fasta"):
            newid = f"seq{s_num}"
            res[newid] = record.id.replace(",", "_")
            record.id = newid
            output_cont.append(record)
            s_num += 1
        SeqIO.write(output_cont, fa_name, "fasta")
        return (fa_name, res)

    def recover_faheader(self, header_tup: tuple):
        """_summary_

        Args:
            header_tup (tuple): _description_
        """
        output_cont = []
        for record in SeqIO.parse(header_tup[0], "fasta"):
            record.id = header_tup[1][record.id]
            record.description = ""
            output_cont.append(record)
        SeqIO.write(output_cont, header_tup[0], "fasta")
        return header_tup[0]

    def set_merge_count(self, clusters: list) -> dict:
        """cluster 파일에서 각 클러스터별 갯수를 아이디에 표시
        ;MC:nnn 이렇게 표시한다.
        만약 한 클러스터 안에 MC:nnn 이 있을 경우에는 거기에 더한다.
        그리고 원래의 id 에 대응 시키고 dict 를 반환한다.

        Args:
            cluster (list): _description_

        Returns:
            dict: _description_
        """
        res = {}
        for cluster in clusters:
            tmp_clst = self.clusterone_merge_count(cluster)
            for kkk, vvv in zip(cluster[1:], tmp_clst[1:]):
                res[kkk.strip().replace(">", "")] = vvv.strip().replace(">", "")
        return res

    def parse_mc_info(self, faid: str):
        """

        Args:
            faid (str): _description_
        """
        # pattern = re.compile(r";MC:(\d+)")
        pattern = self.pattern
        matches = pattern.findall(faid)
        if len(matches) == 0:
            return 0
        else:
            return int(matches[0])

    def change_mc_info(self, faid, newcnt: int) -> str:
        """_summary_

        Args:
            faid (_type_): _description_
            newcnt (int): _description_

        Returns:
            str: _description_
        """
        # pattern = re.compile(r";MC:(\d+)")
        pattern = self.pattern
        replacement = f";MC:{newcnt}"
        res = re.sub(pattern, replacement, faid)
        if res == faid:
            res = f"{faid};MC:{newcnt}"
        return res

    def clusterone_merge_count(self, one_cluster: list):
        """한 개 클러스터 안에서 ;MC:nnn 이 있을 때 ;MC:nnn 이 없는 faid 들의 갯수를 더한다

        Args:
            one_cluster (list): ex:: ['>Cluster 0', ' >faid', ' >faid2',...]
        """
        cont = []
        zero_cnt = 0
        for faid in one_cluster[1:]:
            if self.parse_mc_info(faid) == 0:
                zero_cnt += 1
            else:
                cont.append(self.parse_mc_info(faid))
        tar_mc = max(cont) + zero_cnt if len(cont) > 0 else zero_cnt
        res = [one_cluster[0]]
        for faid in one_cluster[1:]:
            res.append(self.change_mc_info(faid, tar_mc))
        return res

    def main(self):
        import sys

        parser = argparse.ArgumentParser(
            description="From fasta, 100% identical sequences are merged to a consensus sequence"
        )
        parser.add_argument("filename", help="Input FASTA file")
        parser.add_argument(
            "-s",
            "--sort",
            dest="sort_way",
            choices=["size", "ascending", "descending"],
            help="sort the output fasta ordered by cluster size, ascending, or descending ids",
            default="false",
        )
        parser.add_argument(
            "--debugq",
            default=False,
            type=bool,
            help="if True, all temperary files are not removed.",
        )
        args = parser.parse_args(sys.argv[1:])

        fa_name = args.filename
        ## 입력 fasta 서열에서 - 제거
        fa_rmgap_name = self.rmgap_infasta(fa_name)
        out_fasta_name = Common.change_fa_suffix_with(fa_rmgap_name, "_merge.fasta")
        ## cd-hit 을 실행하여 cluster 를 생성하기 전에 blast 에서 나온 서열 id 들이 너무 길기 때문에 seqNNN형태로 변환
        fa_rmgap_name, uniq_id_dict = self.mK_faheader_unique(fa_rmgap_name)
        ## 동일한 서열들을 cluster 생성
        self.run_cdhit100(fa_rmgap_name)
        ## 동일한 서열 cluster 내의 id 들을 추출
        clusters = self.sub_parse_100_clstr(
            self.get_cd_hit100_fn(fa_rmgap_name, type="cluster")
        )
        ## 같은 cluster 내의 서열들을 하나로 합치고 출력
        out_fa_cont = []
        for cls in clusters:
            tmp_fa_seqs = self.get_fa_recs(fa_rmgap_name, cls)
            tmp_consensus = self.get_consensus_seq_main(tmp_fa_seqs)
            tmp_seq_record = SeqRecord(Seq(tmp_consensus))
            tmp_seq_record.id = cls[1].replace(" >", "")
            tmp_seq_record.description = ""
            out_fa_cont.append(tmp_seq_record)
            # tmp_header = cls[1]
            # print(f">{tmp_header}\n{tmp_consensus}")
        if args.sort_way == "ascending":
            out_fa_cont = sorted(out_fa_cont, key=lambda x: x.id, reverse=False)
        elif args.sort_way == "descending":
            out_fa_cont = sorted(out_fa_cont, key=lambda x: x.id, reverse=True)
        SeqIO.write(out_fa_cont, out_fasta_name, "fasta")
        ## 같은 cluster 내의 아이디 들을 파일에 저장
        fa_id_fn = self.save_faid_recs(
            fa_rmgap_name, clusters, uniq_id_dict=uniq_id_dict
        )
        ## 각 클러스터에서 MC:nnn 을 설정하고 원래 {faid: faid;MC:nnn} 으로 반환
        faid_clusters = self.recover_clusterid_to_orifaid(clusters, uniq_id_dict)
        uniq_id_dict2 = self.set_merge_count(faid_clusters)
        ## 다시 id 들을 원복
        self.recover_faheader((out_fasta_name, uniq_id_dict))  # , uniq_id_dict2))
        ## ID 들에 faid;MC:nnn 으로 변경
        self.recover_faheader((out_fasta_name, uniq_id_dict2))
        print(f"{out_fasta_name} and {fa_id_fn} are created")
        ## debug 모드일 때는 중간 산물들을 남겨둔다
        if not args.debugq:
            os.remove(fa_rmgap_name)
            os.remove(self.get_cd_hit100_fn(fa_rmgap_name))
            os.remove(self.get_cd_hit100_fn(fa_rmgap_name, type="fasta"))


if __name__ == "__main__":
    runner = MERGE_FASTA()
    runner.main()
