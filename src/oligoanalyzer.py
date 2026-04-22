import math
import os
import re
import subprocess
from openpyxl import load_workbook
from Bio.Seq import Seq
from Bio import Align
from src.common import Common
from src.cal_melting_temperature import MeltingTemp
from src.config import Config
import src.smith_waterman as smwa
import uuid
from bs4 import BeautifulSoup


class OligoAnalyzer:
    def __init__(self):
        self.workbook = None
        self.melting_temp_obj = MeltingTemp(isdeoxy=True)
        self.salt_concentration = 50
        self.primer_conc = 50
        self.config = Config().config

    def set_workbook(self, wbname: str):
        """_summary_

        Args:
            wbname (str): _description_
        """
        res = load_workbook(filename=wbname)
        self.workbook = res
        return res

    def set_sheet(self, wsname: str):
        """set an active sheet of self.workbook

        Args:
            wsname (str): _description_
        """
        self.worksheet = self.workbook[wsname]

    def set_concentration(self, salt_conc: float = 50, primer_conc: float = 50) -> bool:
        """setter of self.salt_concentration, primer_concentration
        concentration unit: refer to http://biotools.nubic.northwestern.edu/OligoCalc.html
            salt (Na+) : mM
            primer : nM

        Args:
            salt_conc (float): default 50mM
            primer_conc (float): default 50nM

        Returns:
            bool: _description_
        """
        self.salt_concentration = salt_conc
        self.primer_conc = primer_conc
        return True

    def read_sheet(self, sheetname: str) -> dict:
        """input excel 파일에서 종별 fw, probe, rv 서열 및 파라미터를 읽고
        또 template sequence 를 읽어 온다.
        이때 빈칸이나 문자열 끝의 unicode 는 삭제한다.

        Args:
            sheetname (str): 엑셀 파일내의 sheet 이름.

        Return:
            dict: example::
            {"strain name":{
                "Froward": {"seq": <str>, "length": <int>, "GC%": <float>, "Tm": <float>},
                "probe": {상동}
                "reverse": {상동}
                "template": <str>
                },
                {...},...
            }
        """
        res = {}
        if self.workbook is None:
            raise Exception("no workbook is set")
        self.set_sheet(sheetname)
        primer_range, template_range = self._read_sheet_get_primer_range()
        for iii, row in enumerate(self.worksheet[primer_range]):
            if iii % 3 == 0:
                strain_name = row[0].value
                res[strain_name] = {}
            res[strain_name][row[1].value.lower()] = {
                "seq": row[2].value.replace(" ", "").replace("\xa0", ""),
                "length": int(row[3].value),
                r"GC%": float(row[4].value),
                "Tm": float(row[5].value),
            }
        for row in self.worksheet[template_range]:
            strain_name = row[0].value
            res[strain_name]["template"] = row[1].value.replace(" ", "").replace("\xa0", "")
        return res

    def read_sheet_to_list(self, sheet_dict: dict) -> list:
        """read_sheet() 의 ouput 을 list 로 바꾼다

        Args:
            sheet_dict (dict): _description_

        Returns:
            list: _description_
        """
        res = []
        for species, seqdict in sheet_dict.items():
            for seqtype, seqrec in seqdict.items():
                if seqtype != "template":
                    tmplist = [
                        species,
                        seqtype,
                        seqrec["seq"],
                        seqrec["length"],
                        seqrec["GC%"],
                        seqrec["Tm"],
                    ]
                else:
                    tmplist = [species, seqtype, seqrec]
                res.append(tmplist)
        return res

    def _read_sheet_get_primer_range(self) -> tuple:
        """엑셀파일은 1행부터 시작이다. 0행이 아님.

        Args:
            worksheet (_type_): _description_

        Return:
            tuple : primer range, template range
        """
        primer_start_rn = 0
        primer_end_rn = 0
        template_start_rn = 0
        for iii, row in enumerate(self.worksheet.iter_rows()):
            if row[0].value == "종":
                primer_start_rn = iii + 1
            if row[0].value == "target":
                primer_end_rn = iii - 1
                template_start_rn = iii + 1
        strain_num = int((primer_end_rn - primer_start_rn + 1) / 3)
        return (
            f"A{primer_start_rn+1}:F{primer_end_rn}",
            f"A{template_start_rn+1}:B{template_start_rn+strain_num}",
        )

    def step1_1_get_max_length(self, seq_rec: dict, seq_key: str) -> int:
        """입력된 seq_rec 에서 서열의 길이와 seq_key 종류를 바탕으로 탐구할 최대 길이를 설정한다.
        primer 의 경우는 최대 25bp 까지 입력된 서열길이 + 5bp 까지
        probe 의 경우는 최대 30bp 까지 입력된 서열길이 + 5bp 까지

        Args:
            seq_rec (dict): self.read_sheet 의 결과물 중 한 개의 리코드
            seq_key (str): forward, probe, reverse 중 하나

        Returns:
            int: max length
        """
        res = None
        tar_seq = seq_rec[seq_key]["seq"]
        if seq_key in ["forward", "reverse"]:
            res = min(25, len(tar_seq) + 5)
        else:
            res = min(30, len(tar_seq) + 5)
        return res

    def step1_2_extend_seq(self, seq_rec: dict, seq_key: str) -> list:
        """step1.2 : template sequence 에서 max len 만큼 하나씩 5', 3' 방향으로 서열을 늘려간다.
        ----
        pseudo code: 25bp 를 최대로 할 경우

        for lll in range(len(primer_seq)+1, min(25, len(primer_seq)+5)):
            for five_prime_len in range(lll-len(primer_seq)):
                three_prime_len = lll-five_prime_len

        Args:
            seq_rec (dict): _description_
            seq_key (str): _description_

        Returns:
            list: _description_
        """
        res = []
        tar_seq = seq_rec[seq_key]["seq"]
        if seq_key == "reverse":
            tar_seq = str(Seq(tar_seq).reverse_complement())
        tar_seq_len = len(tar_seq)
        template_seq = seq_rec["template"]
        max_seq_len = self.step1_1_get_max_length(seq_rec, seq_key)
        align_pos = self.step1_2_find_align_index(tar_seq, template_seq)
        for cur_len in range(max_seq_len, tar_seq_len - 1, -1):
            aug_len = cur_len - tar_seq_len
            for five_prime_pos in range(max(0, align_pos - aug_len), align_pos + 1):
                three_prime_pos = five_prime_pos + cur_len
                temperary_seq = (
                    template_seq[five_prime_pos:align_pos]
                    + tar_seq
                    + template_seq[align_pos + tar_seq_len : three_prime_pos]
                )
                res.append(temperary_seq)
        return res

    def step1_2_find_align_index(self, tar_seq: str, template_seq: str) -> int:
        """primer 가 template 에서 어디에서 matching 되는지 찾는다.
        gap 을 허용하지 않는다. degenerate base 도 반영한다.
        Args:
            tar_seq (str): _description_
            template_seq (str): _description_

        Returns:
            int: _description_
        """
        res = None
        res = template_seq.find(tar_seq)
        max_degenerate_base = 1
        if res < 0:  # not find case
            tar_seq_len = len(tar_seq)
            template_seq_len = len(template_seq)
            for iii in range(template_seq_len - tar_seq_len):
                tmp_seq = template_seq[iii : iii + len(tar_seq)]
                same_cnt = tar_seq_len
                for a, b in zip(tar_seq, tmp_seq):
                    if a != b:
                        same_cnt -= 1
                if (tar_seq_len - same_cnt) <= max_degenerate_base:
                    res = iii
                    break
        return res

    def step1_3_cal_tm(self, seq_rec: dict, seq_key: str, tm_method: str = "nn") -> dict:
        """self.melting_temp_obj 에서는 두개의 온도 계산이 가능
        if args.type == "adj": salt (NA+) adjusted
            tmval = runner.cal_adj_melting_temp(salt_concentration=args.na)
        elif args.type == "nn": nearest neighbor
            tmval = runner.NearestNeighborTM(
                primer_concentration=args.primer_conc, salt_concentration=args.na
            )

        Args:
            seq_rec (dict): _description_
            seq_key (str): _description_
            tm_method (str): adj or nn. default adj

        Returns:
            dict: _description_
        """
        res = {}
        generated_seqs = self.step1_2_extend_seq(seq_rec, seq_key)
        for seq in generated_seqs:
            self.melting_temp_obj.set_sequence(seq)
            if tm_method == "adj":
                res[seq] = self.melting_temp_obj.cal_adj_melting_temp(
                    salt_concentration=self.salt_concentration
                )
            elif tm_method == "nn":
                res[seq] = self.melting_temp_obj.NearestNeighborTM(
                    primer_concentration=self.primer_conc,
                    salt_concentration=self.salt_concentration,
                )
        return res

    def step1_4_filtering_tm(self, seq_tm: dict, tm_range: list) -> dict:
        """주어진 tm_range 에 있는 서열들을 반환한다.

        Args:
            seq_tm (dict): _description_
            tm_range (list): _description_

        Returns:
            dict: _description_
        """
        res = {}
        for iseq, itm in seq_tm.items():
            if not (itm > tm_range[0] and itm < tm_range[1]):
                continue
            res[iseq] = itm
        return res

    def step1_tm_filter_main(
        self,
        seq_rec: dict,
        filter_tm_range: list,
        tm_method: str = "nn",
        salt_concentration: float = 50,
        primer_concentration: float = 50,
    ) -> dict:
        """_summary_

        Args:
            seq_rec (dict): _description_
            filter_tm_range (list): [primer_from, primer_to, probe_from, probe_to]

        Returns:
            dict: {"forward|probe|reverse": {"seq1": Tm1, ...},...}. reverse 의 경우 원래 input 인 3'->5' 서열로 반환
                ex::
                {
                    "forward": {
                        "TGCACCACATGTGTTTTTTATTG": 57.51, ...
                    },
                    "probe": {
                        "AATTTCTTTGGTGGCGGGAGCAATCCTACC": 72.07, ...
                    },
                    "reverse": {
                        "TTTGACTGTAAATAAAAAGTTTGG": 54.98, ...
                    }
                }
        """
        res = {}
        self.set_concentration(salt_conc=salt_concentration, primer_conc=primer_concentration)
        for seq_key in ["forward", "reverse", "probe"]:
            generated_seq_tms = self.step1_3_cal_tm(seq_rec, seq_key, tm_method=tm_method)
            if seq_key == "probe":
                filterd_seq_tms = self.step1_4_filtering_tm(generated_seq_tms, filter_tm_range[2:])
            else:
                filterd_seq_tms = self.step1_4_filtering_tm(generated_seq_tms, filter_tm_range[:2])
            if seq_key == "reverse":
                tmpdic = {}
                for iseq, itm in filterd_seq_tms.items():
                    tmpseq = Seq(iseq)
                    tmpdic[str(tmpseq.reverse_complement())] = itm
                filterd_seq_tms = tmpdic
            res[seq_key] = filterd_seq_tms

        return res

    def step2_run_mfold(self, step1dic: dict, qt_label=None):
        """local 에 설치되어있는 mfold 를 실행한다.
        Args:
            step1dic (dict): self.ste1_tm_filter_main 의 결과물
            qt_label () : oligoanalyzer_gui 에서 호출할 때 사용하는 변수
        """
        res = step1dic
        jobid = uuid.uuid1()
        cur_dir = os.getcwd()
        working_dir = f"{self.config['TEMPDIR']}/{jobid}"
        os.makedirs(working_dir)
        os.chdir(working_dir)
        if qt_label is not None:
            _tmp_cnt = 1
        for seqtype, seqdic in step1dic.items():
            for eachseq, melttemp in seqdic.items():
                mfold_out = self._run_mfold_local(seqtype, eachseq, jobid)
                res[seqtype][eachseq] = {
                    "Tm": melttemp,
                    "mfold": self.step2_mfold_parser(mfold_out),
                }
                if qt_label is not None:
                    qt_label[1].setText(f"status: {qt_label[0]}'s {_tmp_cnt} sequence(s) done")
                    qt_label[2].processEvents()
                    _tmp_cnt += 1
        os.chdir(cur_dir)
        return res

    def _run_mfold_local(self, seqtype: str, sequence: str, jobid: str) -> dict:
        """서열을 받은 다음 tempdir/uuid/ 를 생성하고 거기 아래에서 mfold 를 실행 한 다음 결과를 리턴한다.
        mfold 명령어는 다음과 같다.
        참고: issue 93
        mfold 사용시에 세팅해야 할 옵션
            NA=RNA (default) or DNA
            LC=sequence type (default = linear) or circular
            SEQ=file_name 여기서 확장자 아래는 제외한 filename 을 이용해서 subdirectory 를 만든다고 적혀있음
            RUN_TYPE=text or html
        ```
        mfold SEQ=mfoldtest.fasta NA=DNA RUN_TYPE=html
        ```
        Args:
            seqtype (str): forward, reverse, probe
            sequence (str): primer, probe 의 서열. 보통 50bp 미만이다.
            jobid (str): uuid

        Returns:
            dict: _description_
        """
        # fasta_full_path = Common.get_mfold_file_path(self.config, jobid, sequence)
        os.makedirs(f"{seqtype}-{sequence}")
        curdir = os.getcwd()
        os.chdir(f"{seqtype}-{sequence}")
        fastain = f"{sequence}.fasta"
        with open(fastain, "w") as handle:
            handle.write(f">seq\n{sequence}")
        ret = subprocess.run(["mfold", f"SEQ={fastain}", "NA=DNA", "RUN_TYPE=html", "NA_CONC=50"])
        try:
            with open(fastain.replace(".fasta", ".out.html"), "r") as handle:
                res = handle.read()
        except FileNotFoundError as e:
            res = {"error": e}
        os.chdir(curdir)
        return res

    def step2_mfold_parser(self, htmlcont: str):
        res = {
            "Foldingat(Celsius)": None,
            "Na+": None,
            "Mg++": None,
        }
        if not isinstance(htmlcont, str):
            return res
        try:
            soup = BeautifulSoup(htmlcont, "html.parser")
            structures = soup.find_all("pre")
            for iii, structure in enumerate(structures):
                res[str(iii + 1)] = {"fold": structure.text}
            iii = 1
            for line in htmlcont.split("\n"):
                if "Linear DNA" in line:
                    res["Foldingat(Celsius)"] = line
                    res.update(self.step2_html_parser1(line))
                if line[:3] == " dG":
                    # res[str(iii)].update({"dG": line})
                    res[str(iii)].update(self.step2_html_parser2(line))
                    iii += 1
        except TypeError:
            pass
        return res

    def step2_html_parser1(self, htmlcont) -> dict:
        """input: 'Linear DNA folding at 37 &#176;C. [Na<SUP>+</SUP>] = 50, [Mg<SUP>++</SUP>] = 0.0.'

        Args:
            htmlcont (_type_): _description_
        """
        res = {"Foldingat(Celsius)": None, "Na+": None, "Mg++": None}
        temp, cont = htmlcont.split(";")
        tmp1 = temp.replace("Linear DNA folding at ", "")
        tmp1 = tmp1.replace("&#176", "")
        res["Foldingat(Celsius)"] = float(tmp1)
        nacont, mgcont = cont.split(",")
        res["Na+"] = float(nacont.split("=")[1])
        res["Mg++"] = float(mgcont[:-1].split("=")[1])  # 0.0. 이기 때문에 마지막 . 을 제거해야 함.
        return res

    def step2_html_parser2(self, htmlcont) -> dict:
        """input: ' dG =     -6.43 &nbsp; dH =    -40.20 &nbsp; dS =   -108.88 &nbsp; T<sub>m</sub> =   96.1 ℃'

        Args:
            htmlcont (_type_): _description_

        Returns:
            dict: _description_
        """
        # Define regular expressions for each variable
        dG_pattern = r"dG\s*=\s*([-\d\.]+)"
        dH_pattern = r"dH\s*=\s*([-\d\.]+)"
        dS_pattern = r"dS\s*=\s*([-\d\.]+)"
        Tm_pattern = r"T<sub>m<\/sub>\s*=\s*([-\d\.]+)"

        # Use regular expressions to extract values
        dG_match = re.search(dG_pattern, htmlcont)
        dH_match = re.search(dH_pattern, htmlcont)
        dS_match = re.search(dS_pattern, htmlcont)
        Tm_match = re.search(Tm_pattern, htmlcont)

        # Create dictionary with extracted values
        res = {
            "dG": float(dG_match.group(1)),
            "dH": float(dH_match.group(1)),
            "dS": float(dS_match.group(1)),
            "Tm(Celsius)": float(Tm_match.group(1)),
        }
        return res

    def step2_get_reverseseq(self, inseq: str) -> str:
        """_summary_

        Args:
            inseq (str): _description_

        Returns:
            str: _description_
        """
        res = inseq[::-1]
        return res

    def _complement_matchq(self, abase: str, bbase: str) -> bool:
        """check whether a and b bases are complement bases
        A-T, G-C are complement
        Args:
            abase (str): A base
            bbase (str): B base

        Returns:
            bool: True or False. complement is True
        """
        comp_list = [["A", "T"], ["C", "G"]]
        in_bp = abase + bbase
        if sorted(in_bp) in comp_list:
            return True
        else:
            return False

    def step2_get_selfmaxmatching(self, inseq: str) -> str:
        """_summary_

        Args:
            inseq (str): _description_

        Returns:
            str: _description_
        """
        res = []
        cont = []
        cntcont = []
        rvseq = self.step2_get_reverseseq(inseq)
        ttlen = len(inseq)
        for iii in range(ttlen):
            mtcnt = 0  # match count
            mtstr = ""  # match string
            for a, b in zip(inseq[:-iii], rvseq[iii:]):
                if self._complement_matchq(a, b):
                    mtcnt += 1
                    mtstr += a
                    continue
                cont.append(mtstr)
                cntcont.append(mtcnt)
                mtcnt = 0
                mtstr = ""
        maxmtlen = max(cntcont)
        for sss in cont:
            if len(sss) == maxmtlen:
                res.append(sss)
        return res

    def step2_homodimer_cal_tm(
        self,
        seqlist: list,
    ) -> list:
        """우선 nearest neighbor 방법으로 계산 한 다음 (원래는 7bp 이상이어야 계산 가능)
        영하 온도가 계산되어서 나오면 adjusted melting temperature 계산

        Args:
            seqlist (list): [matching_seq1, ...]

        Returns:
            list: 가장 높은 Tm 값을 반환 : (Tm, deltaG) 형태
        """
        res = []
        for matseq in seqlist:
            self.melting_temp_obj.set_sequence(matseq)
            tmp_tm = self.melting_temp_obj.NearestNeighborTM(
                primer_concentration=self.primer_conc,
                salt_concentration=self.salt_concentration,
                ignore_seq_len=True,
            )
            delta_g = self.melting_temp_obj.DeltaG(ignore_seq_len=True)
            if tmp_tm < 0:
                tmp_tm = self.melting_temp_obj.cal_adj_melting_temp(
                    salt_concentration=self.salt_concentration
                )
            res.append((tmp_tm, -delta_g))
        sorted(res)
        return res[-1]

    def step2_homodimer_main(self, hairpindic: dict) -> dict:
        """_summary_

        Args:
            hairpindic (dict): _description_

        Returns:
            dict: _description_
        """
        for seq_type in hairpindic:
            for in_seq in hairpindic[seq_type]:
                mat_list = self.step2_get_selfmaxmatching(in_seq)
                maxtm, mindeltg = self.step2_homodimer_cal_tm(mat_list)
                hairpindic[seq_type][in_seq]["homodimer"] = {"Tm": maxtm, "deltaG": mindeltg}
        return hairpindic

    def step2_selection_done(self, step2_final: dict, selection_dict: dict) -> dict:
        """step2_homodimer_main의 결과물에서 selection_dict 에 있는 것들만 남긴다.

        Args:
            step2_final (dict): Ex::
                {
                    "forward": {
                        <seq> : {}, ...
                    },
                    "reverse": {
                        <seq> : {}, ...
                    },
                    "probe": {
                        <seq> : {}, ...
                    }
                }
            selection_dict (dict): Ex::
                {
                    "forward": [<seq>, ...],
                    "reverse": [],
                    "probe": []
                }

        Returns:
            dict: step2_final (dict) 와 동일한 형태. 다만 selection_dict 에 있는 것만 남긴다
        """
        res = {}
        for seq_type in selection_dict:
            if seq_type not in res:
                res[seq_type] = {}
            for seq in selection_dict[seq_type]:
                res[seq_type][seq] = step2_final[seq_type][seq]
        return res

    def step3_get_crossmaxmatching(self, inseq1: str, inseq2: str) -> str:
        """_summary_

        Args:
            inseq (str): _description_

        Returns:
            str: _description_
        """
        res = []
        cont = []
        cntcont = []
        match_pos = []
        ttlen = min(len(inseq1), len(inseq2))
        for iii in range(ttlen):
            mtcnt = 0  # match count
            mtstr = ""  # match string
            for a, b in zip(inseq1[:-iii], inseq2[iii:]):
                if self._complement_matchq(a, b):
                    mtcnt += 1
                    mtstr += a
                    continue
                cont.append(mtstr)
                cntcont.append(mtcnt)
                match_pos.append(iii)
                mtcnt = 0
                mtstr = ""
        maxmtlen = max(cntcont)
        for sss, pos in zip(cont, match_pos):
            if len(sss) == maxmtlen:
                res.append(sss)
        return res

    def step3_heterodimer_main(self, step2dic: dict) -> dict:
        """step2 의 결과물을 가지고 모든 서열 pair 간의 Tm 과 deltaG 를 계산한다.

        Args:
            step2dic (dict): _description_

        Returns:
            dict: _description_
        """
        res = {}
        seqtype_seqlist = []
        for seqtype in step2dic:
            for seq in step2dic[seqtype]:
                seqtype_seqlist.append((seqtype, seq))
        for iii in range(len(seqtype_seqlist) - 1):
            for jjj in range(iii + 1, len(seqtype_seqlist)):
                maxmatch = self.step3_get_crossmaxmatching(
                    seqtype_seqlist[iii][1], seqtype_seqlist[jjj][1]
                )
                maxtm, mindeltg = self.step2_homodimer_cal_tm(maxmatch)
                res[seqtype_seqlist[iii] + seqtype_seqlist[jjj]] = {"Tm": maxtm, "deltaG": mindeltg}
        return res

    def step3_combine_all_strains(self, step2dic_all: dict) -> dict:
        """_summary_

        Args:
            step2dic_all (dict): _description_

        Returns:
            dict: _description_
        """
        res = {"forward": {}, "reverse": {}, "probe": {}}
        for strain_name, seq_dict in step2dic_all.items():
            for p_type, sub_dict in seq_dict.items():
                for sequence, seq_dict in sub_dict.items():
                    seq_dict["strain"] = strain_name
                res[p_type].update(sub_dict)
        return res

    def step3_heterodimer_all_strains_main(self, combin_dict: dict) -> dict:
        """step3_combine_all_strains의 결과물부터 시작
        probe-probe pair 는 무시

        Args:
            combin_dict (dict): _description_

        Returns:
            dict: _description_
        """
        res = {}
        DELIM = ":"
        seqtype_seqlist = []
        for seqtype in combin_dict:
            for seq in combin_dict[seqtype]:
                strain_name = combin_dict[seqtype][seq]["strain"]
                seqtype_seqlist.append((f"{strain_name}{DELIM}{seqtype}", seq))
        for iii in range(len(seqtype_seqlist) - 1):
            for jjj in range(iii + 1, len(seqtype_seqlist)):
                if ":probe" in seqtype_seqlist[iii][0] and ":probe" in seqtype_seqlist[jjj][0]:
                    continue
                maxmatch = self.step3_get_crossmaxmatching(
                    seqtype_seqlist[iii][1], seqtype_seqlist[jjj][1]
                )
                maxtm, mindeltg = self.step2_homodimer_cal_tm(maxmatch)
                pairwise_alignment = self.step3_smith_waterman(
                    seqtype_seqlist[iii][1], seqtype_seqlist[jjj][1]
                )
                res[seqtype_seqlist[iii] + seqtype_seqlist[jjj]] = {
                    "Tm": maxtm,
                    "deltaG": mindeltg,
                    "tt_binding": pairwise_alignment.count("|"),
                    "pairwise": pairwise_alignment,
                }
        print(res)
        return res

    def step3_smith_waterman(self, seq1: str, seq2: str):
        """_summary_

        Args:
            seq1 (str): _description_
            seq2 (str): _description_
        """
        # start, end = smwa.smith_waterman(seq1, seq2)
        # return seq1[start:end]
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.match_score = 3
        aligner.mismatch_score = -1
        aligner.open_gap_score = -10000
        # aligner.extend_gap_score = -10000
        comp_seq2 = str(Seq(seq2).complement())
        alignments = aligner.align(seq1, comp_seq2)
        aln = sorted(alignments)[0]

        t_segments, q_segments = aln.aligned
        if len(t_segments) == 0:
            return f"{seq1}\n\n{seq2}\n"

        t_start = t_segments[0][0]
        q_start = q_segments[0][0]

        mid1 = str(aln[0])
        mid2 = str(aln[1])

        markers = ""
        for a, b in zip(mid1, mid2):
            if a == b and a != "-":
                markers += "|"
            else:
                markers += "."

        S1 = max(0, q_start - t_start)
        S2 = max(0, t_start - q_start)
        SM = max(t_start, q_start)

        res = " " * S1 + seq1 + "\n"
        res += " " * SM + markers + "\n"
        res += " " * S2 + seq2 + "\n"
        return res


if __name__ == "__main__":
    pass
