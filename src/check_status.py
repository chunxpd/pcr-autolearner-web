import os
from src.common import Common
import yaml


class ChkStatus:
    def __init__(self, config, uuid, primer_seq=""):
        # self.tarfn = (
        #     f"{config['TEMPDIR']}/{config['jobid']}-{primer_seq}-genebank.fasta"
        # )
        self.tarfn = Common.get_genbank_fn(config, primer_seq)
        self.progress_log_fn = Common.get_dw_genebank_progress_log_fn(config, uuid)

    def check_num_fasta(self):
        """fasta 파일에서 현재 내려받은 레코드의 수를 확인한다.
        self.tarfn 파일은 DwGeneBank.download() 에서 만드는 파일이다.
        fasta 파일의 각 레코드의 헤더는 > 으로 시작한다.
        따라서 > 의 갯수를 세어서 반환하는 메쏘드이다.
        """
        try:
            res = os.popen(f"grep -c '^>' {self.tarfn}")
            res = int(res.read())
        except ValueError:
            res = 0
        return res

    def check_progress_log(self):
        with open(self.progress_log_fn) as handle:
            log_data = yaml.load(handle, Loader=yaml.FullLoader)
        return log_data


if __name__ == "__main__":
    pass
