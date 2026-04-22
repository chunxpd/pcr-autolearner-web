from typing import List
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import product
import re

import pysam


class Mapper:
    def __init__(self):
        pass

    def mk_gap_dic(self, fasta_records: list) -> dict:
        """From the input fasta records of which the format is list(SeqIO.parse),
        output is {seqid: [(gap_init, gap_end)], ...}
        The indices are identical to the python index system (0-base, 1-base).
            ori. pos : 1 3   4   8   9  12
            sequence1: ATG---CATGC---ATGC
            gap. ps  : 1 3   7   11  15
        (3~7), (8~12)
        Args:
            fasta_records (list): list(SeqIO.parse(open(fasta)))

        Returns:
            dict: {seqid: [(gap_init, gap_end)], ...}
        """
        gap_dict = {}
        for record in fasta_records:
            seq_id = record.id
            sequence = str(record.seq)
            gaps = []
            start = None
            for i, base in enumerate(sequence):
                if base == "-":
                    if start is None:
                        start = i
                else:
                    if start is not None:
                        gaps.append((start, i))
                        start = None
            # If the sequence ends with a gap
            if start is not None:
                gaps.append((start, len(sequence)))
            gap_dict[seq_id] = gaps
            for iii in range(0, len(gaps) - 1):
                iii_len = gaps[iii][1] - gaps[iii][0]
                for jjj in range(iii + 1, len(gaps)):
                    gaps[jjj] = (gaps[jjj][0] - iii_len, gaps[jjj][1] - iii_len)
        return gap_dict

    def adj_map_pos(self, sam_pos: int, fa_id: str, gap_dic: dict) -> int:
        """_summary_

        Args:
            sam_pos (in): POS of SAM/BAM file (1-based)
            gap_dic (dict): the return of self.mk_gap_dic()
            fa_id (str): fasta record id which is the key of gap_dic

        Returns:
            int: _description_
        """
        gap_list = gap_dic[fa_id]
        sam_pos_0based = sam_pos - 1
        gap_len_list = [t[1] - t[0] for t in gap_list if t[0] <= sam_pos_0based]
        return sum(gap_len_list) + sam_pos

    def gen_seqs_from_degerancy(self, input_reg_expression):
        """_summary_

        Args:
            fasta_records (_type_): _description_
        """
        # Find all variable parts in the input_reg_expression
        pattern = re.compile(r"\[(.*?)\]")
        matches = pattern.findall(input_reg_expression)

        # Replace each variable part with its options
        fixed_parts = pattern.split(input_reg_expression)
        fixed_parts = fixed_parts[::2]

        # Generate all combinations for variable parts
        sequences = []
        for combination in product(*[list(options) for options in matches]):
            seq = "".join(
                fixed_parts[i] + combination[i] for i in range(len(combination))
            )
            seq += fixed_parts[-1]
            sequences.append(seq)

        # Remove duplicates and sort
        sequences = sorted(set(sequences))

        return sequences

    def gen_fasta_from_degenerancy(
        self, fasta_recs: List[SeqRecord]
    ) -> List[SeqRecord]:
        """From a fasta records whose sequences have degenerate nucleoties as a format of regular expression,
        generate nucleotie sequences

        Args:
            fasta_recs (List[SeqRecord]): fasta sequence records with degenerate nucleotides

        Returns:
            List[SeqRecord]: a list of normal nucleotide sequences
        """
        res = []
        for record in fasta_recs:
            seq_id = record.id
            degenerate_sequence = str(record.seq)
            sequences = self.gen_seqs_from_degerancy(degenerate_sequence)
            for iii, seq in enumerate(sequences, start=1):
                res.append(SeqRecord(Seq(seq), id=f"{seq_id}_{iii}", description=""))
        return res

    def rm_gaps(self, fasta_recs: List[SeqRecord]) -> List[SeqRecord]:
        """from input fasta sequences, remove gaps and store it"""
        res = []
        for recored in fasta_recs:
            seq_id = recored.id
            ori_seq = str(recored.seq)
            rmgap_seq = ori_seq.replace("-", "")
            res.append(
                SeqRecord(Seq(rmgap_seq), id=seq_id, description=recored.description)
            )
        return res

    def parse_sam_file(self, samfile: pysam.AlignmentFile) -> dict:
        """_summary_

        Args:
            samfile (pysam.AlignmentFile): _description_

        Returns:
            dict: _description_
        """
        mapping_dict = {}

        for read in samfile:
            if not read.is_unmapped:
                ref_id = samfile.get_reference_name(read.reference_id)
                query_id = read.query_name
                mapped_pos = read.reference_start + 1  # SAM format is 1-based

                if ref_id not in mapping_dict:
                    mapping_dict[ref_id] = {}
                mapping_dict[ref_id][query_id] = mapped_pos

        return mapping_dict

    def transform_data_to_matrix(self, data: dict):
        # Extract unique reference IDs and query IDs
        ref_ids = sorted(set(data.keys()))  # Reference IDs
        query_ids = sorted(
            set(query_id for ref_map in data.values() for query_id in ref_map)
        )

        # Initialize matrix with None (or a different placeholder for missing values)
        matrix = {
            ref_id: {query_id: None for query_id in query_ids} for ref_id in ref_ids
        }

        # Fill the matrix with mapped positions
        for ref_id, query_map in data.items():
            for query_id, mapped_pos in query_map.items():
                if query_id in matrix[ref_id]:
                    matrix[ref_id][query_id] = mapped_pos

        return query_ids, ref_ids, matrix


if __name__ == "__main__":
    runner = Mapper()
