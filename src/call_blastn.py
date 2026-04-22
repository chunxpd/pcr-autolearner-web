import logging
import os
import os.path
import shlex
import subprocess
from time import sleep
import multiprocessing

import Bio.Blast.NCBIWWW as blast
import Bio.Blast.NCBIXML as bxml
from Bio.Seq import Seq

from src.common import Common


class CallBlast:
    semaphore = 0

    def __init__(self, config):
        self.uuid = config["jobid"]
        self.tempdir = config["TEMPDIR"]
        self.hitlist_size = int(config["HITLISTSIZE"])
        self.blastn = config["BLASTN"]
        self.nt = config["NT"]
        self.sh_get_species_taxids = config["SC_TAXSEARCH"]
        self.testq = config["TESTQ"]
        self.cpu_num = self.get_cpunum()

    def get_cpunum(self):
        return multiprocessing.cpu_count()

    def query_blast(
        self, in_sequence: str, strand: int = 1, megablastq: bool = False
    ) -> str:
        """Blast web 에 in_sequence 를 던진다.
            query 결과는 tempdir/<uuid>-<in_sequence>.xml 로 저장한다.
        Args:
            in_sequence (str): a primer sequence
            strand (int) : (optional) +1 : forward primer, -1: reverse primer
            megablastq (bool): True/False

        Returns:
            str: 저장된 파일 이름
        """
        outfn = Common.get_ncbiall_xml_fn(
            {"TEMPDIR": self.tempdir, "jobid": self.uuid}, in_sequence
        )
        if self.testq and os.path.isfile(outfn):
            return outfn
        # outfn = f"{self.tempdir}/{self.uuid}-{in_sequence}.xml"
        # if os.path.isfile(outfn):
        #     return outfn

        while True:
            if CallBlast.semaphore > 0:
                print(f"Wait: semaphore ({CallBlast.semaphore}), blastn:{in_sequence}")
                sleep(10)
                continue
            else:
                CallBlast.semaphore += 1
                print(f"Start...blastn:{in_sequence}")
                try:
                    res = blast.qblast(
                        "blastn",
                        "nt",
                        self.reverse_complement(in_sequence, strand),
                        short_query=True,
                        megablast=megablastq,
                        hitlist_size=self.hitlist_size,
                    )
                except Exception:
                    return "qblast Error"
                with open(outfn, "w") as handle:
                    handle.write(res.read())
                res.close()
                CallBlast.semaphore -= 1
                break
        return outfn

    def query_blast_exclude_taxa(
        self, in_sequence: str, exc_taxa: str, strand: int = 1, megablastq: bool = False
    ) -> str:
        """query_blast_exclude_taxa(
            self, in_sequence: str, exc_taxa: str, strand: int = 1, megablastq: bool = False
        ) -> str:
            Blast web 에 in_sequence 를 던진다.
                query 결과는 tempdir/<uuid>-<in_sequence>.xml 로 저장한다.
                exc_taxa 는 ; 구분자로 되어있는 txid 이다.

            Args:
                in_sequence (str): a primer sequence
                exc_taxa (str): txid;txid;...
                strand (int): (optional) +1: forward strand, -1: reverse strand
                megablastq (bool): True/False

            Returns:
                str: 저장 파일 이름
        """
        outfn = Common.get_ncbiall_xml_exc_fn(
            {"TEMPDIR": self.tempdir, "jobid": self.uuid}, in_sequence
        )
        if self.testq and os.path.isfile(outfn):
            return outfn
        # outfn = f"{self.tempdir}/{self.uuid}-{in_sequence}_excludetaxa.xml"
        # if os.path.isfile(outfn):
        #     return outfn

        while True:
            if CallBlast.semaphore > 0:
                print(
                    f"Wait exclude: semaphore ({CallBlast.semaphore}), "
                    + f"blastn:{in_sequence}, exclude: {self.parse_taxid(exc_taxa)}"
                )
                sleep(10)
                continue
            else:
                print(f"Start...exclude:{in_sequence}")
                CallBlast.semaphore += 1
                exclude_str = self.parse_taxid(exc_taxa)
                try:
                    res = blast.qblast(
                        "blastn",
                        "nt",
                        self.reverse_complement(in_sequence, strand),
                        short_query=True,
                        entrez_query=exclude_str,
                        megablast=megablastq,
                        hitlist_size=self.hitlist_size,
                    )
                except Exception:
                    return "qblast Error"

                with open(outfn, "w") as handle:
                    handle.write(res.read())
                res.close()
                CallBlast.semaphore -= 1
                break
        return outfn

    def reverse_complement(self, in_sequence: str, strand: int) -> str:
        """reverse_complement(self, in_sequence: str) -> str
        DNA sequence in_sequence 를 reverse complement 하여 반환한다.

        Args:
            in_sequence (str): DNA 서열
            strand (int): +1 이면 forward 이므로 입력 값 그대로 반환. -1 이면 r.c.

        Returns:
            str: reverse complement DNA 서열
        """
        res = ""
        if strand == 1:
            res = in_sequence
        else:
            ori_seq = Seq(in_sequence)
            tar_seq = ori_seq.reverse_complement()
            res = str(tar_seq)
        return res

    def parse_taxid(self, exc_str: str, sep: str = "entrez") -> str:
        """parse_taxid(self, exc_str: str) -> str:

        Args:
            exc_str (str): excluding string
                EX:: Mycobacterium (taxid: 1763), Canis lupus (taxid: 9612),
            sep (str): return format
                entrez: see return example
                ,: comma seperate
        Returns:
            str: entrez query string
                ex1:: all[filter] NOT txid1763[ORGN] NOT txid9612[ORGN]
                ex2:: 1763,9612
        """
        res = []
        parse_res = exc_str.split(",")
        for rec in parse_res:
            if "taxid:" not in rec:
                continue
            tmplist = rec.split("taxid:")
            res.append(tmplist[1].replace(" ", "").replace(")", ""))
        if sep == "entrez":
            res = ["NOT txid" + v + "[ORGN]" for v in res]
            res = " ".join(res)
            res = "all[filter] " + res
        elif sep == ",":
            res = ",".join(res)
        return res

    def query_blast_standalone(
        self, in_sequence: str, strand: int = 1, megablastq: bool = False
    ) -> str:
        """Blast standalone 에 in_sequence 를 던진다.
            query 결과는 tempdir/<uuid>-<in_sequence>.xml 로 저장한다.
        Args:
            in_sequence (str): a primer sequence
            strand (int) : (optional) +1 : forward primer, -1: reverse primer
            megablastq (bool): True/False

        Returns:
            str: 저장된 파일 이름.
                만약 config['TESTQ']가 True 이고 기존에 돌려 놓은 결과가 있으면 바로 리턴
        """
        os.makedirs(f"{self.tempdir}/{self.uuid}", exist_ok=True)
        outfn = Common.get_ncbiall_xml_fn(
            {"TEMPDIR": self.tempdir, "jobid": self.uuid}, in_sequence
        )
        if self.testq and os.path.isfile(outfn):
            return outfn
        in_fa = Common.get_ncbiall_fa_fn(
            {"TEMPDIR": self.tempdir, "jobid": self.uuid}, in_sequence
        )
        task = "megablast" if megablastq else "blastn-short"
        with open(in_fa, "w") as handle:
            final_seq = self.reverse_complement(in_sequence, strand)
            handle.write(f">primer\n{final_seq}")

        while True:
            semaphore_thr = int(self.cpu_num / 2)
            if CallBlast.semaphore > semaphore_thr:
                logging.info(
                    f"Wait: semaphore ({CallBlast.semaphore}), blastn:{in_sequence}"
                )
                sleep(60)
                continue
            else:
                CallBlast.semaphore += 1
                logging.info(f"Start...blastn:{in_sequence}")
                execstr = (
                    f"{self.blastn} -query {in_fa} -db {self.nt} "
                    f"-outfmt 5 -out {outfn} -num_threads {self.cpu_num} -task {task} "
                    f"-max_target_seqs {self.hitlist_size} "
                    f"-evalue 1000 -word_size 7 -gapopen 5 -gapextend 2 -reward 1 -penalty -3"
                )
                logging.info(f"blastn: {execstr}")
                res = subprocess.run(shlex.split(execstr))
                try:
                    res.check_returncode()
                except subprocess.CalledProcessError:
                    return "blastn Error"
                CallBlast.semaphore -= 1
                break
        return outfn

    def query_blast_exclude_taxa_standalone(
        self, in_sequence: str, exc_taxa: str, strand: int = 1, megablastq: bool = False
    ) -> str:
        """query_blast_exclude_taxa(
            self, in_sequence: str, exc_taxa: str, strand: int = 1, megablastq: bool = False
        ) -> str:
            Blast standalone 에 in_sequence 를 던진다.
                query 결과는 tempdir/<uuid>-<in_sequence>.xml 로 저장한다.
                exc_taxa 는 ; 구분자로 되어있는 txid 이다.

            Args:
                in_sequence (str): a primer sequence
                exc_taxa (str): txid;txid;...
                strand (int): (optional) +1: forward strand, -1: reverse strand
                megablastq (bool): True/False

            Returns:
                str: 저장 파일 이름
                만약 config['TESTQ']가 True 이고 기존에 돌려 놓은 결과가 있으면 바로 리턴
        """
        os.makedirs(f"{self.tempdir}/{self.uuid}", exist_ok=True)
        outfn = Common.get_ncbiall_xml_exc_fn(
            {"TEMPDIR": self.tempdir, "jobid": self.uuid}, in_sequence
        )
        if self.testq and os.path.isfile(outfn):
            return outfn

        in_fa = Common.get_ncbiall_fa_exc_fn(
            {"TEMPDIR": self.tempdir, "jobid": self.uuid}, in_sequence
        )
        taxidlist = Common.get_ncbiall_taxidlist_exc_fn(
            {"TEMPDIR": self.tempdir, "jobid": self.uuid}, in_sequence
        )
        task = "megablast" if megablastq else "blastn-short"
        with open(in_fa, "w") as handle:
            final_seq = self.reverse_complement(in_sequence, strand)
            handle.write(f">primer\n{final_seq}")

        while True:
            semaphore_thr = int(self.cpu_num / 2)
            if CallBlast.semaphore > semaphore_thr:
                logging.info(
                    f"Wait exclude: semaphore ({CallBlast.semaphore}), "
                    + f"blastn:{in_sequence}, exclude: {self.parse_taxid(exc_taxa, sep=',')}"
                )
                sleep(60)
                continue
            else:
                logging.info(f"Start...exclude:{in_sequence}")
                CallBlast.semaphore += 1
                exclude_list = self.get_species_taxids(exc_taxa, taxidlist)
                execstr = (
                    f"{self.blastn} -query {in_fa} -db {self.nt} "
                    f"-outfmt 5 -out {outfn} -num_threads {self.cpu_num} -task {task} "
                    f"-negative_taxidlist {exclude_list} "
                    f"-max_target_seqs {self.hitlist_size} "
                    f"-evalue 1000 -word_size 7 -gapopen 5 -gapextend 2 -reward 1 -penalty -3"
                )
                logging.info(f"blastn: {execstr}")
                res = subprocess.run(shlex.split(execstr))
                try:
                    res.check_returncode()
                except subprocess.CalledProcessError:
                    return f"blastn exclude taxs Error: {in_sequence}"
                CallBlast.semaphore -= 1
                break
        return outfn

    def get_species_taxids(self, exc_taxa: str, taxidfn: str):
        """
        Mycobacterium (taxid: 1763), Canis lupus (taxid: 9612),
        와 같은 문자열을 받아들이고나서 입력된 taxanomy id 의 하위 taxa 들을 모두 가지고 온다.

        Args:
            exc_taxa (str): _description_
        """
        res = []
        if os.path.exists(taxidfn):
            os.remove(taxidfn)
        logging.info(f"START:: make exc taxid lists of {exc_taxa} : {taxidfn}")
        taxids = self.parse_taxid(exc_taxa, sep=",")
        for iid in taxids.split(","):
            if iid == "":
                continue
            execstr = f"{self.sh_get_species_taxids} -t {iid}"
            tmpout = subprocess.run(shlex.split(execstr), capture_output=True)
            res.extend(str(tmpout.stdout, "utf-8").split("\n"))
        res = list(set(res))
        with open(taxidfn, "w") as handle:
            handle.write("\n".join(res))
        logging.info(f"END:: make exc taxid lists of {exc_taxa} : {taxidfn}")
        return taxidfn


if __name__ == "__main__":
    pass
