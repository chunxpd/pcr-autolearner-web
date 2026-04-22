import itertools


class NTSeqGenerator:
    def __init__(self):
        self.tt_gen_seq = 0
        self.nucleotides = ["A", "T", "G", "C"]

    def gen_ntsequences(self, length):
        res = []
        self.tt_gen_seq = 0
        nucleotides = self.nucleotides
        for sequence in itertools.product(nucleotides, repeat=length):
            res.append("".join(sequence))
            self.tt_gen_seq += 1
        return res

    def gen_ntsequences_generator(self, length, gcrange=False):
        def is_within_gcrange(combination):
            gc_content = (
                100
                * (combination.count("G") + combination.count("C"))
                / len(combination)
            )
            return gcrange[0] <= gc_content <= gcrange[1] if gcrange else True

        return (
            "".join(combination)
            for combination in itertools.product(self.nucleotides, repeat=length)
            if is_within_gcrange(combination)
        )

    def get_gen_seq_num(self):
        return self.tt_gen_seq

    def print_out_fasta(self, cont: list):
        """입력된 서열 list 를 fasta 형식으로 저장한다.

        Args:
            cont (list): _description_
        """
        num = 0
        for short_seq in cont:
            print(f">seq{num}\n{short_seq}")
            num += 1

    def print_out_fq(self, mer_length: int, cont):
        """

        Args:
            mer_length (int): 13mer 이면 13
            cont (_type_): cont 는 generator 이거나 list
        """
        num = 0
        for short_seq in cont:
            print(f"@seq{mer_length}mer-{num}\n{short_seq}\n+\n{'#'*mer_length}")
            num += 1

    def gen_ntsequences_from_seed(self, seed: tuple) -> list:
        """입력된 하나의 fasta 서열에서 "A,T,G,C"를 붙인 결과를 출력한다.
        서열 아이디에는 _숫자 를 붙인다.

        Args:
            seed (tuple): _description_

        Returns:
            list: _description_
        """
        res = []
        seq_id, short_seq = seed
        nucleotides = self.nucleotides
        num = 0
        for nt in nucleotides:
            res.append((f"{seq_id}_{num}", f"{short_seq}{nt}"))
            num += 1
        return res


if __name__ == "__main__":
    import sys, argparse

    parser = argparse.ArgumentParser(
        description="generate random sequences with the fixed length."
    )

    parser.add_argument(
        "-gccont",
        nargs=2,
        type=int,
        metavar=("min_value", "max_value"),
        default=False,
        help="Specify the GC content range. 0~100",
    )

    parser.add_argument("mer_length", type=int, help="mer length")

    args = parser.parse_args()

    # Accessing the values
    if args.gccont is not False:
        gccont_min, gccont_max = args.gccont
    mer_length = args.mer_length

    runner = NTSeqGenerator()
    # mer_length = int(sys.argv[1])
    if args.gccont is not False:
        sequences = runner.gen_ntsequences_generator(
            mer_length, gcrange=[gccont_min, gccont_max]
        )
    else:
        sequences = runner.gen_ntsequences_generator(mer_length)
    runner.print_out_fq(mer_length, sequences)
