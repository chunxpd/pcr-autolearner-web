import logging
import os
import shlex
import subprocess

from Bio import SeqIO, AlignIO
from io import StringIO
from src.common import Common


class MSA:
    """
    1. 주어진 jobid 에 해당하는 fasta 파일의 경로를 받는다.
    2. clustalo 를 실행하고
    3. high homology region (HHR)을 찾고
    4. aln 파일을 내려받을 수 있게 하고
    5. 어떤 부위가 HHR 인지를 표 형식으로 화면에 출력한다.
    """

    def __init__(self, config):
        """일단은 jobid 가 정해지면 fasta 파일의 경로와 이름이 정해진다.
        현재는 하드코딩으로 진행하고 있다. 나중에 모든 클래스에서 해당 성격의 파일에 대한
        경로와 이름을 반환하는 것을 만들어야 한다.

        Args:
            config (_type_): _description_
        """
        self.config = config
        self.n_cores = os.cpu_count()
        self.nt_list = [
            "A",
            "C",
            "G",
            "T",
            "U",
            "W",
            "S",
            "M",
            "K",
            "R",
            "Y",
            "B",
            "D",
            "H",
            "V",
        ]

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

    def run_msa(self, in_fa: str, method="clustalo") -> str:
        """
        clustalo 는 다음과 같이 실행한다.
        clustalo -i <in fasta> -o <out.aln> --outfmt=clu -l <log>
                --force -v --auto

        Args:
            method (str, optional): _description_. Defaults to "clustalo". another is mafft

        Returns:
            out_fn (str) : clustalo 결과인 clustalw 형식의 aln 파일 이름 리턴
        """
        # in_fa = self.config["genebankfn"]
        in_fa = in_fa
        out_fn = self._fa_replace(in_fa, f"-{method}.aln")
        log_fn = self._fa_replace(in_fa, f"-{method}.log")
        logging.info(f"MSA method: {method} out: {out_fn}, log: {log_fn}")
        if method == "clustalo":
            # if os.path.isfile(out_fn) and os.path.getsize(out_fn) > 0:
            #     return out_fn
            execstr = (
                f"clustalo -i {in_fa} -o {out_fn} --outfmt=clu "
                f"-l {log_fn} --force -v --auto --threads={self.n_cores}"
            )
            logging.info(f"..... in {__name__}:: msa cmd : {execstr}")
            # _ = os.popen(execstr)  # 이 명령어는 command 를 bg 로 돌리고 바로 반환된다.
            execrt = subprocess.run(shlex.split(execstr))
        elif method == "mafft":
            self.sub_run_msa_mafft(in_fa, out_fn, log_fn)
        return out_fn

    def sub_run_msa_mafft(
        self, in_fa: str, out_fn: str, log_fn: str, aln_option: str = "auto"
    ):
        """refer to
        https://github.com/io2oi/ctlabs-primerwrapper/issues/110
        and https://github.com/io2oi/ctlabs-primerwrapper/issues/137

        Args:
            in_fa (str): _description_
            aln_option (str): one of fasta_optiona, accurate_option, auto_option (default)
        """
        input_file = in_fa
        output_file = out_fn
        log_file = log_fn
        fast_option = "--retree 2 --maxiterate 1000"  # FFT-NS-i
        accurate_option = (
            "--ep 0 --genafpair --maxiterate 1000 --maxiterate 1000"  # E-INS-i
        )
        auto_option = "--auto"
        if aln_option == "auto":
            msa_option = auto_option
        elif aln_option == "fast":
            msa_option = fast_option
        elif aln_option == "accurate":
            msa_option == accurate_option
        else:
            raise (f"There is no option like {aln_option} for MAFFT")
        command_string = f"mafft {msa_option} --clustalout --thread {self.n_cores} --namelength 100 {input_file}"
        logging.info(f"..... in {__name__} :: msa cmd : {command_string}")
        command_args = shlex.split(command_string)
        try:
            with open(output_file, "w") as output_handle, open(
                log_file, "w"
            ) as log_handle:
                subprocess.run(
                    command_args, stdout=output_handle, stderr=log_handle, check=True
                )

        except subprocess.CalledProcessError as e:
            logging.error("Error:", e.stderr.decode("utf-8"))

    def find_hrr(self, aln_fn: str) -> list:
        """msa 결과물에서 * 를 보고 인덱스를 리턴한다.

        Args:
            aln_fn (str): self.run_msa 의 결과물인 aln 파일의 경로

        Returns:
            list: bed 와 유사하게 (0-based, 1-based) 의 index pair 의 list
                [(10, 20), (30, 50), ...]
        """
        ## 현재는 구현하지 않고 Mock output 만 제공
        res = [(10, 20), (30, 50)]
        return res

    def hrr_list_to_html(self, hrr_regs: list) -> str:
        """self.run_hrr 의 결과를 html table 로 변환한다"""
        # res = "<table class='tablp'>"
        # res += "<tr>"
        res = "<tr>"
        for kkk in ["start", "end"]:
            res += f"<th>{kkk}</th>"
        res += "</tr>"
        for hrr in hrr_regs:
            res += "<tr>"
            res += f"<td>{hrr[0]}</td>"
            res += f"<td>{hrr[1]}</td>"
            res += "</tr>"
        # res += "</table>"
        return res

    def _existq_(self, iid: str, container):
        res = False
        for each_id in container:
            if iid in each_id or each_id in iid:
                res = True
                break
        return res

    def rm_sequences(self, in_fasta: StringIO, rm_ids: list) -> str:
        """in_fasta 파일을 열고 rm_ids 에 있는 것들은 제거 한다.
        in_fasta 파일 내의 id 와 rm_ids 의 아이디 들은 어느 한 쪽으로는 완전 일치해야 한다.

        Args:
            in_fasta (str): fasta file path
            rm_ids (list): id list

        Returns:
            str: fasta 형식의 string 또는 list
        """
        res = []
        for rec in SeqIO.parse(in_fasta, "fasta"):
            identifier = rec.id
            desc = rec.description
            seq = rec.seq
            if self._existq_(identifier, rm_ids):
                continue
                # for each_id in rm_ids:
                #     if identifier in each_id or each_id in identifier:
                #         continue
            res.append(SeqIO.SeqRecord(seq=seq, id=identifier, description=desc))
        return "".join([r.format("fasta") for r in res])

    def rc_sequences(self, in_fasta: StringIO, rc_ids: list) -> str:
        """in_fasta 파일을 열고 rc_ids 에 있는 것들은 reverse complement한다.
        in_fasta 파일 내의 id 와 rc_ids 의 아이디 들은 어느 한 쪽으로는 완전 일치해야 한다.

        Args:
            in_fasta (str): fasta file path
            rc_ids (list): id list

        Returns:
            str: fasta 형식의 string 또는 list
        """
        res = []
        for rec in SeqIO.parse(in_fasta, "fasta"):
            identifier = rec.id
            desc = rec.description
            seq = rec.seq
            if self._existq_(identifier, rc_ids):
                seq = seq.reverse_complement()
            # for each_id in rc_ids:
            #     if identifier in each_id or each_id in identifier:
            #         seq = seq.reverse_complement()
            res.append(SeqIO.SeqRecord(seq=seq, id=identifier, description=desc))
        return "".join([r.format("fasta") for r in res])

    def aln_count_nucleotides(self, in_aln: StringIO) -> dict:
        """clustalo 의 output 을 읽고 각 position 별 ATGC 및 total count 를 count 한다.
        서열에 N 이 있을 경우 count 하지 않고 Total count 에는 더한다.
        dict 의 출력 순서는
        출력 순서를 https://en.wikiversity.org/wiki/Gene_transcriptions/Degenerate_nucleotides 에 있는 Symbol 순서 (두번째 것)에 따라 출력하는지 본다.
        Symbol[15]	Description	Bases represented
        A	adenosine
        C	cytidine
        G	guanosine
        T	thymidine
        U	uridine
        W	weak
        S	strong
        M	[Amine] amino
        K	[Ketone] keto
        R	[Purine] purine
        Y	[Pyrimidine] pyrimidine
        B	not A (B comes after A)
        D	not C (D comes after C)
        H	not G (H comes after G)
        V	not T (V comes after T and U)
        N or X	any base (not a gap)

        Args:
            in_aln (StringIO): clustalo output 파일

        Returns:
            dict: 예:: {"A": [], "C": [], "G": [], "T": [], ,"U": [], "W": [], "S": [], "M":[],
                        "K":[], "R":[], "Y": [], "B":[], "D": [], "H": [], "V":[], TT": []}
        """
        alignment = AlignIO.read(in_aln, "clustal")
        total_len = alignment.get_alignment_length()
        nt_list = self.nt_list
        # [
        #     "A",
        #     "C",
        #     "G",
        #     "T",
        #     "U",
        #     "W",
        #     "S",
        #     "M",
        #     "K",
        #     "R",
        #     "Y",
        #     "B",
        #     "D",
        #     "H",
        #     "V",
        # ]
        res = dict()
        for nt in nt_list:
            res[nt] = [0] * total_len
        res["TT"] = [0] * total_len
        for rec in alignment:
            for iii in range(total_len):
                nt = rec.seq[iii].upper()
                if nt == "-":
                    continue
                # if nt in ["A", "T", "G", "C"]:
                if nt in nt_list:
                    res[nt][iii] += 1
                res["TT"][iii] += 1
        return res

    def mk_consensus_seq(self, cnt_dict: dict) -> str:
        """from the output aln_count_nucleotides(), to make the consensus sequence.
        If all ratio are same of nucleotides, regular expression should be used such as [AC]

        Args:
            cnt_dict (dict): _description_

        Returns:
            str: _description_
        """
        res = []
        aln_length = len(cnt_dict["A"])
        for iii in range(aln_length):
            tmpcont = []
            for nt in self.nt_list:
                tmpcont.append(cnt_dict[nt][iii])
            maxcnt = max(tmpcont)
            max_indices = [
                index for index, item in enumerate(tmpcont) if item == maxcnt
            ]
            if len(max_indices) == 1:
                tmpres = self.nt_list[max_indices[0]]
            else:
                tmpres = ""
                for jjj in max_indices:
                    tmpres += self.nt_list[jjj]
                tmpres = "[" + tmpres + "]"
            res.append(tmpres)
        return ",".join(res)

    def print_out_aln_count(self, cnt_dict: dict, file_out: str = None) -> str:
        """_summary_

        Args:
            cnt_dict (dict): aln_count_nucleotides() 의 반환값
            file_outq (str): output file name. .aln changed to -ntcnt.csv, None => just return the string

        Returns:
            str: string output
        """
        ttlen = len(cnt_dict["A"])
        res = []
        header = ["nt"] + [str(v) for v in range(1, ttlen + 1)]
        header = ",".join(header)
        res.append(header)
        for base, cnt_list in cnt_dict.items():
            tmpout = [base] + [str(v) for v in cnt_list]
            res.append(",".join(tmpout))
        consensus_seq = self.mk_consensus_seq(cnt_dict)
        res.append("con," + consensus_seq)
        if file_out is not None:
            cnt_fn = Common.suffix_replace(file_out, ".aln", "-ntcnt.csv")
            with open(cnt_fn, "w") as handle:
                handle.write("\n".join(res))
            return cnt_fn
        return "\n".join(res)


def fa_replace(in_fa: str, tar_suffix: str) -> str:
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


def gen_id_list(in_file: str) -> list:
    """file 을 읽고 \n 을 제거하고 리스트로 반환한다.
    idlist 는 fasta id 에서 > 를 제외한 문자열이다."""
    res = []
    with open(in_file, "r") as handle:
        for line in handle:
            res.append(line.strip())
    return res


def wr_fasta_str(in_fa_str: str, config: dict, ori_fn: str, tar_suffix: str) -> str:
    """fasta 형식의 문자열을 uuid 로 저장하고 파일이름을 반환한다.

    Args:
        in_fa_str (str): _description_

    Returns:
        str: _description_
    """
    # adj_fn = uuid.uuid4()
    # final_in_fasta = f"{config['TEMPDIR']}/{str(adj_fn)}.fasta"
    final_in_fasta = fa_replace(ori_fn, tar_suffix)
    with open(final_in_fasta, "w") as handle:
        handle.write(adj_fasta)
    return final_in_fasta


if __name__ == "__main__":
    import sys
    from src.config import Config
    import argparse
    import uuid

    logging.basicConfig(
        level=logging.INFO,
        datefmt="%m/%d/%Y %I:%M:%S %p",
        format="%(asctime)s:%(levelname)s:%(message)s",
    )

    parser = argparse.ArgumentParser(
        prog="multiple_sequence_alignment.py", description="multiple sequence alignment"
    )
    parser.add_argument("fasta", type=str, help="input fasta file")
    parser.add_argument(
        "--rmlist", type=str, help="A file containg the fasta id lists of removal"
    )
    parser.add_argument(
        "--rclist",
        type=str,
        help="A file containg the fasta id lists of reverse complement",
    )
    parser.add_argument(
        "--method",
        choices=["clustalo", "mafft"],
        default="clustalo",
        help="MSA method. Select one",
    )

    args = parser.parse_args(sys.argv[1:])
    config = Config()
    msa = MSA(config)
    final_in_fasta = args.fasta

    if args.rmlist:
        rm_list = gen_id_list(args.rmlist)
        adj_fasta = msa.rm_sequences(args.fasta, rm_list)
        final_in_fasta = wr_fasta_str(
            adj_fasta, config.config, final_in_fasta, "-rmids.fasta"
        )
    if args.rclist:
        rc_list = gen_id_list(args.rclist)
        adj_fasta = msa.rc_sequences(args.fasta, rc_list)
        final_in_fasta = wr_fasta_str(
            adj_fasta, config.config, final_in_fasta, "-rcids.fasta"
        )
    aln_fn = msa.run_msa(final_in_fasta, method=args.method)
    cnt_dict = msa.aln_count_nucleotides(aln_fn)
    cnt_fn = msa.print_out_aln_count(cnt_dict, file_out=aln_fn)
    logging.info(f"MSA result file: {aln_fn}")
    logging.info(f"MSA nt count file: {cnt_fn}")
