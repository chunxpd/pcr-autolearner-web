import datetime


class Common:
    def __init__(self, config):
        self.config = config

    @classmethod
    def get_genbank_fn(self, config, primer):
        res = f"{config['TEMPDIR']}/{config['jobid']}-{primer}-genebank.fasta"
        return res

    @classmethod
    def get_msa_fn(self, config):
        res = f"{config['TEMPDIR']}/{config['jobid']}-genebank.aln"
        return res

    @classmethod
    def get_gb_fn(self, config, gbid):
        res = f"{config['GB_STORAGE']}/{gbid}.gb"
        return res

    @classmethod
    def get_ncbiall_xml_fn(self, config, primer_seq):
        res = (
            f"{config['TEMPDIR']}/{config['jobid']}/{config['jobid']}-{primer_seq}.xml"
        )
        return res

    @classmethod
    def get_ncbiall_fa_fn(self, config, primer_seq):
        res = f"{config['TEMPDIR']}/{config['jobid']}/{config['jobid']}-{primer_seq}.fa"
        return res

    @classmethod
    def get_ncbiall_xml_exc_fn(self, config, primer_seq):
        res = f"{config['TEMPDIR']}/{config['jobid']}/{config['jobid']}-{primer_seq}_excludetaxa.xml"
        return res

    @classmethod
    def get_ncbiall_fa_exc_fn(self, config, primer_seq):
        res = f"{config['TEMPDIR']}/{config['jobid']}/{config['jobid']}-{primer_seq}_excludetaxa.fa"
        return res

    @classmethod
    def get_ncbiall_taxidlist_exc_fn(self, config, primer_seq):
        res = f"{config['TEMPDIR']}/{config['jobid']}/{config['jobid']}-{primer_seq}_excludetaxa.list"
        return res

    @classmethod
    def get_upload_path(self, config, uuid):
        res = f"{config['TEMPDIR']}/{uuid}"
        return res

    @classmethod
    def get_upload_file_path(self, config, uuid, filename):
        res = f"{config['TEMPDIR']}/{uuid}/{filename}"
        return res

    @classmethod
    def get_dw_genebank_progress_log_fn(self, config, uuid):
        res = f"{config['TEMPDIR']}/{uuid}/{uuid}_dw_gbk.log"
        return res

    @classmethod
    def color_wrapping(self, in_text, color):
        col_dic = {
            "grey": "\x1b[38;20m",
            "yellow": "\x1b[33;20m",
            "red": "\x1b[31;20m",
            "bold_red": "\x1b[31;1m",
            "reset": "\x1b[0m",
        }
        res = (
            f"{col_dic[color]}{in_text}{col_dic['reset']}"
            if color in col_dic
            else in_text
        )
        return res

    @classmethod
    def suffix_replace(self, ori_fn: str, ori_suffix: str, tar_suffix: str):
        """ori_fn에서 ori_suffix 를 찾아서 tar_suffix 로 바꾸어줌"""
        res = ""
        ridx = ori_fn.rfind(ori_suffix)
        if ridx >= 0:
            res = ori_fn[:ridx] + tar_suffix
        return res

    @classmethod
    def get_blastn_pair_score_fn(self, config, uuid, fn):
        res = f"{config['TEMPDIR']}/{uuid}/{fn}"
        return res

    @classmethod
    def get_mfold_file_path(self, config, jobid, sequence):
        res = f"{config['TEMPDIR']}/{jobid}/{sequence}.fasta"
        return res

    @classmethod
    def get_date_string(self):
        # Get the current date in the specified format
        current_date = datetime.datetime.now().strftime("%Y%m%d")
        # Construct the JOBID
        res = f"{current_date}"
        return res

    @classmethod
    def get_rmgap_fa_fn(self, in_fasta: str) -> str:
        res = in_fasta.replace(".fa", "_rmgap.fa")
        return res

    @classmethod
    def change_fa_suffix_with(self, in_fa: str, new_suffix: str) -> str:
        if in_fa.endswith(".fa"):
            res = in_fa.replace(".fa", new_suffix)
        elif in_fa.endswith(".fasta"):
            res = in_fa.replace(".fasta", new_suffix)
        else:
            raise ValueError(
                f"input fasta file MUST end with .fa or .fasta. But your file is {in_fa}"
            )
        return res


if __name__ == "__main__":
    pass
