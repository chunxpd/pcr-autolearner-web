from io import StringIO
import logging
import math
import os
import os.path
import yaml

import Bio.Blast.NCBIXML as bxml
import numpy as np
from Bio import Entrez, SeqIO
from Bio.Seq import Seq, UndefinedSequenceError
from Bio.SeqFeature import CompoundLocation
from tqdm import tqdm

from src.common import Common


class DwGeneBank:
    def __init__(self, config, primer: str, max_updown_len: int = 5000):
        self.email = config["email"]
        self.uuid = config["jobid"]
        self.primer = primer
        self.ncbiall = Common.get_ncbiall_xml_fn(config, self.primer)
        self.outfn = Common.get_genbank_fn(config, primer)
        self.config = config
        self.entrez_gb_container = {
            "gbid": None,
            "organism": None,
            "seq": None,
            "features": None,
        }
        self.progress_log = dict()
        self.max_updown_len = max_updown_len

    def wr_progress_log(self, data: dict):
        """
        download() 의 진행상태를 저장한다.
        progress_log.yaml:
        {
            primer1: {
                "N_GBID": <int>,  # download 진행할 gbid 갯수
                "CUR_STEP": <str> # 현재 상태
                "EST_TIME(s)": <int> # 25초*N_GBID
            },
            primer2: {...}, ...
        }
        한 개의 DwGeneBank object 는 한 개의 primer 서열을 갖는다.
        """
        if self.primer in self.progress_log:
            self.progress_log[self.primer].update(data)
        else:
            self.progress_log.update({self.primer: data})
        yaml_out = Common.get_dw_genebank_progress_log_fn(self.config, self.uuid)
        if not os.path.isfile(yaml_out):
            handle = open(yaml_out, "w", encoding="utf-8")
            handle.write("")
            handle.close()
        with open(yaml_out, "r", encoding="utf-8") as handle:
            log_cont = yaml.load(handle, Loader=yaml.FullLoader)
        try:
            log_cont.update(self.progress_log)
        except:
            log_cont = self.progress_log
        with open(yaml_out, "w", encoding="utf-8") as handle:
            yaml.dump(log_cont, handle, default_flow_style=False)

    def download_full_sequence(self):
        """genebank 에 gb 아이디들은 , 분류 리스트가 들어간다
        예 : "CP101406.1,CP097320.1,CP076378.1,CP014860.2"
        Returns:
            _type_: _description_
        """
        Entrez.email = self.email  # Always tell NCBI who you are

        iids = self.gen_idlist(rttype="list")
        if os.path.isfile(self.outfn) and os.path.getsize(self.outfn) > 0:
            return self.outfn
        for iid in iids:
            handle = Entrez.efetch(db="nucleotide", id=iid, rettype="fasta", retmode="text")
            with open(self.outfn, "a") as outfp:
                outfp.write(handle.read())
            handle.close()
        return self.outfn

    def download(
        self, al_thr_coeff: float = 0.8, id_thr_coeff: float = 0.9, final_strand: int = 1
    ) -> str:
        """
        config = {
            "email": "<your email>",
            "jobid": "<uuid>",
            "ncbiallfn": "tests/blast_include_example.xml",
            "ncbiexcfn": "tests/blast_exclude_example.xml",
            "TEMPDIR": "<tempdirpath>",
        }\n
        genebank = DwGeneBank(config)\n
        genebank.download()

        ncbiall 에 있는 HSP 들을 모두 읽어서 해당하는 nucleotide 서열들을
        Entrez 에서 내려받은 이후 (gb file) gb file 을 파싱하여 HSP 의 up/downstream
        유전자들을 outfn 에 저장하는 메써드.

        최초 BLASTn 을 할 때 input 이 reverse 인 상태로 들어갔다면 여기서 output 서열의 방향을 reverse complement 로 바꾼다

        테스트 결과: 한개 GB file 내려받는데 대략 25초 걸림.

        Arguments:
            al_thr_coeff (float) : for gen_hsp_rec() and hsp_select_criteria()
            id_thr_coeff (float) : for gen_hsp_rec() and hsp_select_criteria()
            final_strand (int) : 최종 서열의 방향. 1이면 그대로 -1 이면 reverse complement

        Returns:
            str: outfn
        """
        ### step 0: 이전에 작업하다 만 fasta 파일이 있으면 지운다.
        if os.path.isfile(self.outfn) and os.path.getsize(self.outfn) > 0:
            os.remove(self.outfn)
        self.wr_progress_log({"CUR_STEP": "step0"})
        ### step 1: blastn 결과에서 hsp 들에 대한 정보를 가지고 온다.
        hsps = self.gen_hsp_rec(al_thr_coeff=al_thr_coeff, id_thr_coeff=id_thr_coeff)
        n_hsps = len(hsps)
        logging.info(f"# of hsps: {n_hsps}")
        self.wr_progress_log({"CUR_STEP": "step1: parse blastn"})
        hsps_iterator = self.gbid_list_iterator(hsps, 50, joinq=False)
        for _hsps in hsps_iterator:
            ### step 2: gbid 에 따른 entrez GB file 들을 미리 내려 받는다.
            #### step 2.1: Entrez 요청을 한 번으로 끝내기 위해서 한개 파일에 모든 gbid 들을 요청
            try:
                gb_all_fn = self.download_all_gbfiles(_hsps)
                if not gb_all_fn:
                    logging.info("No new gb ids")
                    continue
                    # raise Exception(
                    #     f"ERROR:: in {__name__}: downloading GB files error"
                    # )
            except Exception as err:
                logging.error(err)
                continue
            #### step 2.2: 모든 gbid 들이 한개 파일에 저장되어있으므로 모두 나누어 저장
            self.wr_progress_log({"CUR_STEP": "step2.2: split the downloaded GB file"})
            gb_splitq = self.split_each_gb(gb_all_fn)
            if not gb_splitq:
                raise Exception(f"ERROR:: in {self.download.__name__}: splitting GB all file error")
        ### step 3: HSP 의 up/downstrem 의 유전자 서열들을 저장한다
        iii = 1
        for hsp in tqdm(hsps):
            gbid = hsp["gbid"]
            overlabinfo = self.find_overlab_entrez_genes(hsp)
            gb_container = self.get_entrez_sequence_info(gbid)
            fa_meta = {"name": gb_container["description"], "hsp": hsp}
            try:
                final_ntseq = self.get_ntseqeunce(gbid, overlabinfo, self.max_updown_len)
                if final_strand == -1:
                    _ori_seq = Seq(final_ntseq)
                    _tar_seq = _ori_seq.reverse_complement()
                    final_ntseq = str(_tar_seq)
                fa_str = self.convert_seq_to_fasta(
                    fa_meta,
                    final_ntseq,
                    # self.get_ntseqeunce(gbid, overlabinfo, self.max_updown_len),
                    overlabinfo,
                )
            except TypeError as e:
                print(f"{e}: {gbid}, {gb_container}, {hsp}")
                raise e
            with open(self.outfn, "a") as outfp:
                outfp.write(fa_str)
            self.wr_progress_log(
                {"CUR_STEP": f"step3: generate fasta of up/downstream of HSPs ({iii}/{n_hsps}"}
            )
            iii += 1
        return self.outfn

    def download_all_gbfiles(self, hsps) -> str:
        """Entrez 로 부터 GB 파일들을 내려받는다.
        config.yaml 의 GB_STORAGE 아래에 내려받으면 됨

        Entrez 가 하루 전체 쿼리 갯수를 1000개로 제한하기 때문에
        현재 저장하고 있는 GB 파일이 없는 것들을 한번의 쿼리로 모두 내려받는다.
        이미 다 로컬에 저장되어 있다면 '' 를 반환한다.

        Args:
            hsps (_type_): _description_

        Return:
            gb_all_fn (str) : hsps 에서 한개 이상의 신규 다운로드가 있을 때. 없으면 ''
        """
        tar_idlist = []
        for hsp in hsps:
            gbid = hsp["gbid"]
            gb_fn = Common.get_gb_fn(self.config, gbid)
            if os.path.isfile(gb_fn) and os.path.getsize(gb_fn) > 0:
                continue
            tar_idlist.append(gbid)
        tar_idlist = list(set(tar_idlist))
        logging.info(f"in {self.download_all_gbfiles.__name__} GB id list: {tar_idlist}")
        self.wr_progress_log(
            {
                "CUR_STEP": "step2.1: start to download GB records",
                "N_GBID": len(tar_idlist),
                "EST_TIME(s)": 25 * len(tar_idlist),
            }
        )
        if len(tar_idlist) == 0:
            return ""
        gb_ids = ",".join(tar_idlist)
        logging.info(f"START:: GB file ({len(tar_idlist)}) downloading")
        gb_all_fn = (
            f"{self.config['TEMPDIR']}/{self.uuid}/" f"{self.uuid}_{self.uuid}_{self.primer}_all.gb"
        )
        with open(gb_all_fn, "w") as handle:
            Entrez.email = self.email
            enz_handle = Entrez.efetch(db="nucleotide", id=gb_ids, rettype="gb", retmode="text")
            handle.write(enz_handle.read())
            enz_handle.close()
        logging.info(f"END:: GB file ({len(tar_idlist)}) downloading")
        return gb_all_fn

    def gbid_list_iterator(self, tar_idlist: list, niter: int, joinq=True) -> list:
        """gbid_list_iterator(self, tar_idlist: list) -> list
        tar_idlist 를 niter 개로 쪼개고 ','.join() 해서 반환한다.

        Args:
            tar_idlist(list): 예)  ["AP022596.1", "CP023583.1", "BX571966.1"]
            niter(int) : 반환하는 iteration의 갯수

        Returns:
            list: 예) niter 가 2 였다면 반환은 2번 iteration 을 할 수 있어야 한다.
        """
        res = []
        n_ele = round(len(tar_idlist) / niter)
        if niter > len(tar_idlist):
            if joinq:
                res.append(",".join(tar_idlist))
            else:
                res.append(tar_idlist)
            return res
        prev_iii = 0
        for iii in range(0, len(tar_idlist), n_ele):
            if iii == 0:
                continue
            if joinq:
                res.append(",".join(tar_idlist[prev_iii:iii]))
            else:
                res.append(tar_idlist[prev_iii:iii])

            prev_iii = iii
        if joinq:
            res.append(",".join(tar_idlist[prev_iii : len(tar_idlist)]))
        else:
            res.append(tar_idlist[prev_iii : len(tar_idlist)])
        return res

    def split_each_gb(self, gb_all_fn: str) -> bool:
        """
        download_all_gbfiles() 의 결과물을 분리하여 개별 gb 파일로 저장한다.

        Args:
            gb_all_fn (str): 애

        Returns:
            bool: _description_
        """
        with open(gb_all_fn, "r") as handle:
            os.makedirs(os.path.dirname(gb_all_fn), exist_ok=True)
            tmpcont = []
            for rec in handle:
                tmpcont.append(rec)
                if rec[:7] == "VERSION":
                    wfn = Common.get_gb_fn(self.config, rec.split()[1])
                if rec[:2] == "//":
                    os.makedirs(os.path.dirname(wfn), exist_ok=True)
                    with open(wfn, "w") as ohandle:
                        ohandle.write("".join(tmpcont))
                    tmpcont = []
        return True

    def gen_idlist(self, rttype="text"):
        """alignment.title 의 형식
        gi|2289144095|gb|CP092364.2| Mycobacterium goodii strain ATCC 700504 chromosome, complete genome
        """
        res = []
        with open(f"{self.ncbiall}") as handle:
            blast_records = bxml.parse(handle)
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    res.append(alignment.title.split("|")[3])
        res = list(set(res))
        res.sort()
        if rttype == "text":
            return ",".join(res)
        else:
            return res

    def gen_hsp_rec(self, nrec=None, al_thr_coeff: float = 0.8, id_thr_coeff: float = 0.9) -> list:
        """blastn 결과인 xml 파일에서 각 HSP 별로 gbid 와 hit 의 index 를 가지고 온다.
        alignment.title 에 있는 형식은 다음과 같으며 GI id 는 항상 있지만 GB id 는 있다없다가 한다.
        gi|85720936|gb|CP000252.1|
        gi|295095013|emb|FP929040.1|
        gi|156095607|ref|XM_001613789.1|
        하지만 상관없이 잘 가지고 온다.

        Args:
            nrec (int, optional): Defaults to "all".
                xml 파일에서 최초 nrec 개를 가지고 온다.

        Return:
            list : [
                {
                "gbid": <genebank id>,
                "Hsp_hit-from": <hit genomic position from>, # 1-based
                "Hsp_hit_to": <hit genomic position to>,      # 1-based
                "Hsp_hit-strand": <+1 or -1>
                }, ...
            ]

        참고: HPS 의 내부 변수
            |  Members:
            |      - score           BLAST score of hit.  (float)
            |      - bits            Number of bits for that score.  (float)
            |      - expect          Expect value.  (float)
            |      - num_alignments  Number of alignments for same subject.  (int)
            |      - identities      Number of identities (int) if using the XML parser.
            |        Tuple of number of identities/total aligned (int, int)
            |        if using the (obsolete) plain text parser.
            |      - positives       Number of positives (int) if using the XML parser.
            |        Tuple of number of positives/total aligned (int, int)
            |        if using the (obsolete) plain text parser.
            |      - gaps            Number of gaps (int) if using the XML parser.
            |        Tuple of number of gaps/total aligned (int, int) if
            |        using the (obsolete) plain text parser.
            |      - align_length    Length of the alignment. (int)
            |      - strand          Tuple of (query, target) strand.
            |      - frame           Tuple of 1 or 2 frame shifts, depending on the flavor.
            |
            |      - query           The query sequence.
            |      - query_start     The start residue for the query sequence.  (1-based)
            |      - query_end       The end residue for the query sequence.  (1-based)
            |      - match           The match sequence.
            |      - sbjct           The sbjct sequence.
            |      - sbjct_start     The start residue for the sbjct sequence.  (1-based)
            |      - sbjct_end       The end residue for the sbjct sequence.  (1-based)
        """
        res = []
        logging.info(f"{self.gen_hsp_rec.__name__} input ncbi file: {self.ncbiall}")
        self.wr_progress_log({"IN_BLASTN": self.ncbiall})
        log_tt_hitnum = 0
        log_tt_hspnum = 0
        log_filtered_hspnum = 0
        with open(f"{self.ncbiall}") as handle:
            blast_records = bxml.parse(handle)
            for blast_record in blast_records:
                query_length = blast_record.query_length
                for alignment in blast_record.alignments:
                    log_tt_hitnum += 1
                    for hsp in alignment.hsps:
                        log_tt_hspnum += 1
                        if not self.hsp_select_criteria(
                            query_length,
                            hsp,
                            al_thr_coeff=al_thr_coeff,
                            id_thr_coeff=id_thr_coeff,
                        ):
                            log_filtered_hspnum += 1
                            continue
                        tmpdic = {}
                        tmpdic["gbid"] = alignment.title.split("|")[3]
                        tmpdic.update(self.get_hsp_hit_genomic_position(hsp))
                        res.append(tmpdic)
        self.wr_progress_log(
            {
                "HSP": {
                    "total_hit": log_tt_hitnum,
                    "total_hsp": log_tt_hspnum,
                    "qc_pass_hsp": log_tt_hspnum - log_filtered_hspnum,
                }
            }
        )
        logging.info(
            (
                f"total {log_tt_hitnum} Hits ({log_tt_hspnum} Hsps) -> "
                f"filtered {log_filtered_hspnum} Hsps, "
                f"remains {log_tt_hspnum-log_filtered_hspnum} Hsps"
            )
        )
        if nrec is not None:
            return res[:nrec]
        else:
            return res

    def get_hsp_hit_genomic_position(self, hsp) -> dict:
        """get_hsp_hit_genomic_position(self, hsp) -> tuple:
        Hsp_hit-frame 이 -1 이면 Hsp_hit-from 과 Hsp_hit-to 를 서로 바꾸어 5'->3'으로 반환받는다.

        Args:
            hsp (Bio.Blast.Record.HSP): HSP record

        Returns:
            tuple: {
                "Hsp_hit-from": int,
                "Hsp_hit-to": int,
                "Hsp_hit-strand": int
                }
        """
        res = {"Hsp_hit-from": None, "Hsp_hit-to": None, "Hsp_hit-strand": None}
        if hsp.strand[1] == "Plus":
            res["Hsp_hit-from"] = hsp.sbjct_start
            res["Hsp_hit-to"] = hsp.sbjct_end
            res["Hsp_hit-strand"] = 1
        else:
            res["Hsp_hit-to"] = hsp.sbjct_start
            res["Hsp_hit-from"] = hsp.sbjct_end
            res["Hsp_hit-strand"] = -1
        return res

    def hsp_select_criteria(
        self,
        query_length: int,
        hsp,
        al_thr_coeff: float = 0.8,
        id_thr_coeff: float = 0.9,
    ) -> bool:
        """hsp_select_criteria(self, query_length: int, hsp: Bio.Blast.Record.HSP) -> bool:
        hsp 를 가지고 올 기준 수립
        2022년 10월 12일 영상회의 결과: https://github.com/io2oi/ctlabs-primerwrapper/issues/16 참조
        query_length 에 대해서 alignment 가 80% 이상 되고 identity는 90% 이상일 때 가지고 옴
        query_length 를 21bp 라고 할 때
        80% 는 21*0.8 = 16.8 즉 17bp 이상이어야 하고
        identity 는 17bp * 0.9 = 15.3 즉 16bp 이상이어야 한다.
        이 비율은 사용자가 조절할 수 있다.

        Args:
            query_length (_type_): blast query 의 길이
            hsp (Bio.Blast.Record.HSP): HSP record
            al_thr_coeff (float): default 0.8. alignment threshold 80%
            id_thr_coeff (float): default 0.9. identity threshold 90%

        Returns:
            bool: True:: 기준 통과. Pass, False:: 기준 통과 못함. Fail
        """
        res = False
        al_thr = math.ceil(query_length * al_thr_coeff)
        id_thr = math.ceil(al_thr * id_thr_coeff)
        if hsp.align_length >= al_thr and hsp.identities >= id_thr:
            res = True
        return res

    def get_entrez_sequence_info(self, gbid: str) -> dict:
        """
        1. Entrez 에서 query 에 해당하는 레코드를 gb 형식으로 저장한다.
        2. 저장된 gb 파일에서 sequence 와 gene 의 location, strand 정보를 가지고 온다.

        Args:
            gbid (str): genebank id (accession number)
                ex:: CP54795.1

        Returns:
            dict: gb 파일 형식의 데이터 레코드에서 필요한 부분 추출
                {
                    "seq": SeqIO.Seqfeature.seq,
                    "organism": <str>
                    "features": [{
                        "start": SeqIO.Seqfeature.location.start, # 0-based
                        "end": SeqIO.Seqfeature.location.end,     # 1-based
                        "strand": SeqIO.Seqfeature.location.strand, # -1: reverse complement
                    }, ... ]
                }
        """
        ### 속도 향상을 위해 한번 내려받고 파싱까지 한 것을 직전에 사용하였으면
        ### 그대로 반환한다.
        if self.entrez_gb_container["gbid"] == gbid:
            res = self.entrez_gb_container
            return res
        res = {"gbid": None, "seq": None, "features": None}
        ### Entrez sequence 저장
        config = self.config
        # config["gbid"] = gbid
        gb_out_fn = Common.get_gb_fn(config, gbid)
        if (
            not os.path.isfile(gb_out_fn) or os.path.getsize(gb_out_fn) == 0
        ):  # gb_out_fn 이 존재하면 실행하지 않음
            os.makedirs(os.path.dirname(gb_out_fn), exist_ok=True)
            with open(gb_out_fn, "w") as handle:
                Entrez.email = self.email
                enz_handle = Entrez.efetch(db="nucleotide", id=gbid, rettype="gb", retmode="text")
                handle.write(enz_handle.read())
                enz_handle.close()
        ### 저장된 sequence 를 읽어오기
        ### Features 에서 gene 부분만 저장한다.
        ### Origin 에 있는 서열만 저장한다.
        try:
            record = SeqIO.read(gb_out_fn, "genbank")
        except ValueError:
            print(f"{gb_out_fn}: ValueError")
            raise ValueError
        res["gbid"] = gbid
        res["seq"] = record.seq
        res["features"] = []
        for feature in record.features:
            if feature.type != "gene":
                continue
            if isinstance(feature.location, CompoundLocation):
                continue
            if abs(feature.location.end - feature.location.start) > self.max_updown_len:
                continue
            res["features"].append(
                {
                    "start": feature.location.start,
                    "end": feature.location.end,
                    "strand": feature.location.strand,
                }
            )
        res["organism"] = record.annotations["organism"]
        res["description"] = record.description
        self.entrez_gb_container = res
        return res

    def find_overlab_entrez_genes(self, gbdic: dict) -> tuple:
        """find_overlab_entrez_genes(self, gbdic: dict) -> tuple
        gbdic 의 up/down stream 의 유전자 영역을 찾는다.
        Entrez genes : s1--------e1  s2-------e2  s3-------e3
        Query HSP    :             !--------$
        와 같을 때
        Upstream 찾기:
            e - ! 의 부호를 계산하면
            Entrez genes : s1--------(-)  s2-------(+)  s3-------(+)
            Query HSP    :             !--------$
            가 되므로 이때는 + 가 나온 gene 위치의 바로 앞 gene 을 가지고 오면 된다.
        Downstream 찾기:
            s - $ 의 부호를 계산하면
            Entrez genes : (-)--------e1  (-)-------e2  (+)-------e3
            Query HSP    :             !--------$
            가 되므로 이때는 + 가 나온 gene 위치를 가지고 오면 된다.
        이 연산을 위해서는 numpy.array 와 numpy.where 를 사용하면 빠르게 계산할 수 있다.

        예외의 경우:
            1.  gb file 에 gene feature 가 없는 경우
                    HSP 의 Hsp_hit-from - 5KB ~ Hsp_hit-to + 5KB 를 반환한다.
            2. 위와 같이 gene 영역을 찾았는데 길이가 10KB를 넘으면 10KB를 반환한다.

        Args:
            gbdic (dict): self.gen_hsp_rec() 의 반환값
                example:: {
                    "gbid": "CP054795.1",
                    "Hsp_hit-from": 2714562, # 0-based
                    "Hsp_hit-to": 2714582,    # 1-based
                    "Hsp_hit-strand": 1 or -1
                }

        Returns:
            tuple: upstream 의 위치, downtream 의 위치 짝. strand 는 gbdic 의 Hsp_hit-strand 와 일치시킨다.
            (
                {"end": 2714465, "start": 2711298, "strand": -1},
                {"end": 2716351, "start": 2714824, "strand": -1},
            )
        """
        res = None
        entrez_data = self.get_entrez_sequence_info(gbdic["gbid"])
        query_start = gbdic["Hsp_hit-from"]  # 0-based index
        query_end = gbdic["Hsp_hit-to"]
        ori_start_list = []
        ori_end_list = []

        for feature in entrez_data["features"]:
            ori_start_list.append(feature["start"])
            ori_end_list.append(feature["end"])
        ori_start_array = np.array(ori_start_list)
        ori_end_array = np.array(ori_end_list)

        ### gene feature 가 3개 이상일 때만 아래 작업 진행. 그 이하는 HSP
        ### 기준 +- self.max_updown_len (default: 5000bp) 을 한다.
        if len(entrez_data["features"]) > 2:
            if ori_end_array.max() < query_start:  # query 가 GENE 영역보다 뒤
                upstream_idx = -1
                res = (
                    entrez_data["features"][upstream_idx],
                    {
                        "start": query_end,
                        "end": len(entrez_data["seq"]),
                        "strand": gbdic["Hsp_hit-strand"],
                    },
                )
            elif query_end < ori_start_array.min():  # query 가 gene 영역보다 앞
                downstream_idx = 0
                res = (
                    {"start": 0, "end": query_start, "strand": gbdic["Hsp_hit-strand"]},
                    entrez_data["features"][downstream_idx],
                )
            else:
                ### hsp 가 gene 안에 있을 경우를 먼저 검사
                res = self._find_overlab_entrez_genes_innerset(
                    ori_start_list, ori_end_list, query_start, query_end, gbdic
                )
                if res is None:
                    ### upstream index 는 최소 0이어야 한다.
                    ### downstream index 는 최대 entrez gene 의 레코드 수여야 한다.
                    try:
                        upstream_idx = max(np.where(ori_end_array - query_start > 0)[0][0] - 1, 0)
                        downstream_idx = min(
                            np.where(ori_start_array - query_end > 0)[0][0],
                            len(entrez_data["features"]),
                        )
                        res = (
                            entrez_data["features"][upstream_idx],
                            entrez_data["features"][downstream_idx],
                        )
                    except IndexError:
                        print(f"Not overlap case: {gbdic}")
        else:
            upstream_pos = max(query_start - self.max_updown_len, 0)
            downstream_pos = min(query_end + self.max_updown_len, len(entrez_data["seq"]))
            res = (
                {
                    "start": upstream_pos,
                    "end": query_start,
                    "strand": gbdic["Hsp_hit-strand"],
                },
                {
                    "start": query_end,
                    "end": downstream_pos,
                    "strand": gbdic["Hsp_hit-strand"],
                },
            )
        res = self._adj_overlab_entrez_genes(
            query_start, query_end, len(entrez_data["seq"]), res, gbdic
        )
        return res

    def _find_overlab_entrez_genes_innerset(
        self, ori_start_list, ori_end_list, query_start, query_end, gbdic
    ):
        """issue102: AF035716.1
            Entrez genes : s1--------e1  s2-------e2  s3-----------e3
            Query HSP    :                               !--------$
        위와 같이 query HSP 가 gene들 중 제일 마지막에 있을 경우 downstream 을 찾으면
            s - $ 의 부호는
            Entrez genes : (-)--------e1 (-)-------e2 (-)-----------e3
            Query HSP    :                               !--------$
        가 된다.
        따라서 이 경우는 전체 gene 들을 모두 돌면서 그 안에 있는 것을 확인하는 부분이 있어 줘야 한다.
        """
        res = None
        for tmpi, tmpj in zip(ori_start_list, ori_end_list):
            if not (tmpi <= query_start and query_end <= tmpj):
                continue
            res = (
                {
                    "start": tmpi,
                    "end": tmpj,
                    "strand": gbdic["Hsp_hit-strand"],
                },
                {
                    "start": tmpi,
                    "end": tmpj,
                    "strand": gbdic["Hsp_hit-strand"],
                },
            )
            break
        return res

    def _adj_overlab_entrez_genes(
        self, query_start: int, query_end: int, seq_len: int, reg: tuple, gbdic: dict
    ) -> tuple:
        """find_overlab_entrez_genes 에서 결과물의 길이가 self.max_updown_len * 2 (defualt: 10Kbp)를
        넘어가면 self.max_updown_len * 2+query seq len 으로 조정
        strand 도 gdbid 의 Hsp_hit-strand 와 일치시킨다.

        Args:
            query_start (int): _description_
            query_end (int): _description_
            seq_len (int): len(entrez_data["seq"])
            reg (tuple): _description_


        Returns:
            tuple: _description_
        """
        if reg[1]["end"] - reg[0]["start"] < self.max_updown_len * 2:
            return reg
        upstream_pos = max(query_start - self.max_updown_len, 0)
        downstream_pos = min(query_end + self.max_updown_len, seq_len)
        reg[0]["start"] = upstream_pos
        reg[0]["end"] = query_start
        reg[0]["strand"] = gbdic["Hsp_hit-strand"]
        reg[1]["start"] = query_end
        reg[1]["end"] = downstream_pos
        reg[1]["strand"] = gbdic["Hsp_hit-strand"]
        return reg

    def get_ntseqeunce(self, gbid: str, overlab_idx: tuple, lenthr: int = 10000) -> Seq:
        """genebank 파일의 sequence 에서 서열을 가지고 온다.
        strand 가 모두 - (reverse) 가 아니면 5'->3' 으로 반환하고
        모두 - 이면 reverse complement 로 반환한다.
        whole gene sequence 의 길이가 lengthr(default: 10Kbp) 이하이면 전체 서열모두 반환.
        이상이면 up/downstream 의 영역에 해당하는 것 반환

        Args:
            gbid (str): genebank accession id
            overlab_idx (tuple): self.find_overlab_entrez_genes() 의 반환
                예:: (
                    {"end": 6993060, "start": 6993000, "strand": 1},
                    {"end": 6993180, "start": 6993120, "strand": 1},
                )

        Returns:
            str: SeqIO.Seq object

        Error:
            만약 genbank 파일에 sequence 가 정의되어있지 않다면
            Bio.Seq.UndefinedSequenceError 가 발생한다.
        """
        entrez_data = self.get_entrez_sequence_info(gbid)
        if len(entrez_data["seq"]) < lenthr:
            res = entrez_data["seq"]
        else:
            start = overlab_idx[0]["start"]
            end = overlab_idx[1]["end"]
            res = entrez_data["seq"][start:end]
        if overlab_idx[0]["strand"] == -1 and overlab_idx[1]["strand"] == -1:
            res = res.reverse_complement()
        try:
            str(res)
        except UndefinedSequenceError:
            res = self.get_entrez_fasta(gbid, overlab_idx)
            logging.warning(
                f"No sequence in GB file:: {gbid} {overlab_idx}, so retrieve a fasta from Entrez again"
            )
        return res

    def convert_seq_to_fasta(self, metainfo: dict, seq: str, overlabinfo: tuple) -> SeqIO.SeqRecord:
        """입력된 서열을 fasta 형식의 문자열로 바꾼다.

        Args:
            metainfo (dict): metainfo
                ex:: {
                    "name": "Mycolicibacterium smegmatis strain FDAARGOS_679 chromosome.",
                    "hsp": {"gbid": "CP054795.1", "Hsp_hit-from": 2714562, "Hsp_hit-to": 2714582},
                }
            seq (str): sequence
            overlabinfo (tuple): HSP 의 up/downstream 위치와 방향 tuple

        Returns:
            str: fasta 형식
        """
        start = overlabinfo[0]["start"] + 1  # 1-base 로 변환
        end = overlabinfo[1]["end"]
        if overlabinfo[0]["strand"] == -1 and overlabinfo[1]["strand"] == -1:
            compq = f"complement({start}-{end})"
        else:
            compq = f"{start}-{end}"
        record = SeqIO.SeqRecord(
            seq=Seq(seq),
            id=f"{metainfo['hsp']['gbid']}:{compq}:{metainfo['name'].replace(' ','_')}",
            name=metainfo["hsp"]["gbid"],
            description=f"HSP:{metainfo['hsp']['Hsp_hit-from']+1}-{metainfo['hsp']['Hsp_hit-to']}",
        )
        return record.format("fasta")

    def get_entrez_fasta(self, gbid: str, overlapinfo: tuple) -> str:
        """Entrez 에 fasta 파일을 하나 요청하여 overlapinfo 에 해당하는 영역을 반환

        Args:
            gbid (str): genbank id
                ex:: CP0001.1
            overlapinfo (tuple): self.find_overlab_entrez_genes() 의 결고

        Returns:
            str: nucleotide sequence
        """
        Entrez.email = self.email
        fa_storage_path = self.config["FA_STORAGE"]
        fa_path = f"{fa_storage_path}/{gbid}.fa"
        if not os.path.exists(fa_path):
            ret = Entrez.efetch(db="nucleotide", id=gbid, rettype="fasta", retmode="text")
            ret_fa = ret.read()
            with open(fa_path, "w") as handle:
                handle.write(ret_fa)
            res_obj = SeqIO.read(StringIO(ret_fa), "fasta")
        else:
            with open(fa_path, "r") as handle:
                res_obj = SeqIO.read(handle, format="fasta")
        try:
            res = res_obj.seq[overlapinfo[0]["start"] : overlapinfo[1]["end"]]
            if overlapinfo[0]["strand"] == -1 and overlapinfo[1]["strand"] == -1:
                res = res.reverse_complement()
        except AttributeError as err:
            raise AttributeError(f"{gbid}:{overlapinfo}=>{err}")
        return str(res)


if __name__ == "__main__":
    import sys
    from src.config import Config
    import argparse

    parser = argparse.ArgumentParser(
        prog="download_genebank.py",
        description=(
            "From blastn xml file, a fasta file is generated "
            "after filtering Hsp's alignment and identity ratio "
            "to the query sequence."
        ),
    )
    parser.add_argument("jobid", type=str, help="jobid for the blastn xml file")
    parser.add_argument("primer", type=str, help="primer or probe sequence of query")
    parser.add_argument(
        "--al_ratio",
        dest="al_thr_coeff",
        type=float,
        default=0.8,
        help="Criteria for filtering HSP. alignment ratio >= 0.8(default)",
    )
    parser.add_argument(
        "--id_ratio",
        dest="id_thr_coeff",
        type=float,
        default=0.9,
        help="Criteria for filtering HSP. identity ratio >= 0.9(default)",
    )
    parser.add_argument(
        "--max_updown_len",
        dest="max_updown_len",
        type=int,
        default=5000,
        help="Maxium base pairs of uptream or downstream gene length (default 5KB)",
    )
    parser.add_argument(
        "--final_strand",
        dest="final_strand",
        type=int,
        default=1,
        help="if not checked reverse when run blastn query, set -1. Then get rev.comp sequences",
    )
    if len(sys.argv[1:]) < 2:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(sys.argv[1:])
    logging.basicConfig(level=logging.INFO)
    config = Config()
    config.config["jobid"] = args.jobid  # sys.argv[1]
    primer = args.primer  # sys.argv[2]
    genebank = DwGeneBank(config.config, primer, max_updown_len=args.max_updown_len)
    # all_gb_path = f"tempdir/{jobid}_{jobid}_{primer}_all.gb"
    # genebank.split_each_gb(all_gb_path)
    genebank.download(
        al_thr_coeff=args.al_thr_coeff,
        id_thr_coeff=args.id_thr_coeff,
        final_strand=args.final_strand,
    )
