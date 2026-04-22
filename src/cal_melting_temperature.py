"""This is coded following the "Adjusted Melting Temperature (Tm) Calculation
depicted on http://biotools.nubic.northwestern.edu/OligoCalc.html#helpthermo
"""

import math


class MeltingTemp:
    """
    The basic calculation is
    for < 14bp sequences
    Tm= (wA+xT)*2 + (yG+zC)*4 - 16.6*log10(0.050) + 16.6*log10([Na+])
    where w,x,y,z are the number of the bases A,T,G,C in the sequence, respectively.
    for longer sequences
    Tm= 100.5 + (41 * (yG+zC)/(wA+xT+yG+zC)) - (820/(wA+xT+yG+zC)) + 16.6*log10([Na+])
    """

    def __init__(self, isdeoxy=True):
        """If DNA, isdeoxy is True else RNA False"""
        self.a_count = None
        self.c_count = None
        self.g_count = None
        self.t_count = None
        self.u_count = None
        self.i_count = None
        self.m_count = None
        self.r_count = None
        self.w_count = None
        self.s_count = None
        self.y_count = None
        self.k_count = None
        self.v_count = None
        self.h_count = None
        self.d_count = None
        self.b_count = None
        self.n_count = None
        self.choice = None
        self.e_gc_min = None
        self.e_gc_max = None
        self.gc_val_min = None
        self.gc_val_max = None
        self.sequence = None
        self.isDeoxy = isdeoxy

    def set_sequence(self, sequence):
        """set sequence

        Args:
            sequence (str): a given primer sequence
        """
        self.a_count = self.count_char("A", sequence)
        self.c_count = self.count_char("C", sequence)
        self.g_count = self.count_char("G", sequence)
        self.t_count = self.count_char("T", sequence)
        self.u_count = self.count_char("U", sequence)
        self.i_count = self.count_char("I", sequence)
        self.m_count = self.count_char("M", sequence)
        self.r_count = self.count_char("R", sequence)
        self.w_count = self.count_char("W", sequence)
        self.s_count = self.count_char("S", sequence)
        self.y_count = self.count_char("Y", sequence)
        self.k_count = self.count_char("K", sequence)
        self.v_count = self.count_char("V", sequence)
        self.h_count = self.count_char("H", sequence)
        self.d_count = self.count_char("D", sequence)
        self.b_count = self.count_char("B", sequence)
        self.n_count = self.count_char("N", sequence)
        self.choice = self.rt_choice(sequence)
        self.e_gc_min = self.g_count + self.c_count + self.s_count
        self.e_gc_max = len(sequence) - self.a_count - self.t_count - self.w_count - self.u_count
        self.gc_val_min, self.gc_val_max = self.cal_gc(sequence)
        self.sequence = sequence

        # count Nearest Neighbors
        self.aa_count = (
            self.count_neighbors("AA", self.sequence)
            + self.count_neighbors("TT", self.sequence)
            + self.count_neighbors("UU", self.sequence)
        )
        self.at_count = self.count_neighbors("AT", self.sequence) + self.count_neighbors(
            "AU", self.sequence
        )
        self.ta_count = self.count_neighbors("TA", self.sequence) + self.count_neighbors(
            "UA", self.sequence
        )
        self.ca_count = (
            self.count_neighbors("CA", self.sequence)
            + self.count_neighbors("TG", self.sequence)
            + self.count_neighbors("UG", self.sequence)
        )
        self.gt_count = (
            self.count_neighbors("GT", self.sequence)
            + self.count_neighbors("AC", self.sequence)
            + self.count_neighbors("GU", self.sequence)
        )
        self.ct_count = (
            self.count_neighbors("CT", self.sequence)
            + self.count_neighbors("AG", self.sequence)
            + self.count_neighbors("CU", self.sequence)
        )
        self.ga_count = (
            self.count_neighbors("GA", self.sequence)
            + self.count_neighbors("TC", self.sequence)
            + self.count_neighbors("UC", self.sequence)
        )
        self.cg_count = self.count_neighbors("CG", self.sequence)
        self.gc_count = self.count_neighbors("GC", self.sequence)
        self.gg_count = self.count_neighbors("GG", self.sequence) + self.count_neighbors(
            "CC", self.sequence
        )

        # Calculate IUpac pairs
        self.IUpairVals_min = [0, 0, 0]
        self.IUpairVals_max = [0, 0, 0]
        for j in range(3):  # (var j=0; j<3; j++):
            self.IUpairVals_min[j] = 0
            self.IUpairVals_max[j] = 0

        # *******************************************************
        for i in range(1, len(sequence)):  # first base can not be IUpacbase
            base0 = sequence[i - 1]
            base = sequence[i]
            temp = [0, 0, 0]

            temp = self.calc_IU_pair(base0, base, i, sequence, "min")
            # if (debug) alert("temp "+temp[0]+" "+temp[1]+" "+temp[2]);
            for j in range(3):  # (var j=0; j<3; j++) {
                self.IUpairVals_min[j] += temp[j]

            # if (debug) alert("mim"+self.IUpairVals_min[0]+" "+self.IUpairVals_min[1]+" "+self.IUpairVals_min[2]);
            temp = self.calc_IU_pair(base0, base, i, sequence, "max")
            # if (debug) alert("temp "+temp[0]+" "+temp[1]+" "+temp[2]);
            for j in range(3):
                self.IUpairVals_max[j] += temp[j]
            # if (debug) alert("max"+self.IUpairVals_max[0]+" "+self.IUpairVals_max[1]+" "+self.IUpairVals_max[2]);

    def rt_choice(self, sequence: str) -> str:
        """degenerate nucleotide checker
        `min`: [ATCG]+
        `max`: degenerate nt exists

        Args:
            sequence (str): primer sequence

        Returns:
            str: min or max
        """
        atcg_cnt = [self.count_char(v, sequence) for v in "ATCG"]
        if len(sequence) - sum(atcg_cnt) > 0:
            return "max"
        return "min"

    def count_char(self, base: str, sequence: str) -> int:
        """_summary_
        IUPAC nucleotide code	Base
        A	Adenine
        C	Cytosine
        G	Guanine
        T (or U)	Thymine (or Uracil)
        R	A or G
        Y	C or T
        S	G or C
        W	A or T
        K	G or T
        M	A or C
        B	C or G or T
        D	A or G or T
        H	A or C or T
        V	A or C or G
        N	any base
        . or -	gap

        Args:
            base (str): _description_
            sequence (str): _description_

        Returns:
            int: _description_
        """
        res = sequence.count(base)
        return res

    def cal_adj_melting_temp_short(
        self, choice: str, sequence: str, salt_concentration: float
    ) -> float:
        # pylint: disable=missing-function-docstring
        res = round(
            2 * (len(sequence) - self.e_gc_max)
            + 4 * (self.e_gc_max)
            + 21.6
            + (7.21 * math.log(salt_concentration / 1000)),
            2,
        )
        if choice == "min":
            res = round(
                2 * (len(sequence) - self.e_gc_min)
                + 4 * (self.e_gc_min)
                + 21.6
                + (7.21 * math.log(salt_concentration / 1000)),
                2,  # math.log is a natural logarithm
            )
        return res

    def cal_adj_melting_temp_long(
        self, choice: str, sequence: str, salt_concentration: float
    ) -> float:
        # pylint: disable=missing-function-docstring
        res = round(
            100.5
            + (0.41 * self.gc_val_max)
            - (820 / len(sequence))
            + (7.21 * math.log(salt_concentration / 1000)),
            2,
        )

        if choice == "min":
            res = round(
                100.5
                + (0.41 * self.gc_val_min)
                - (820 / len(sequence))
                + (7.21 * math.log(salt_concentration / 1000)),
                2,
            )
        return res

    def cal_gc(self, sequence: str) -> float:
        """return gc contents depending on the existence of the degenerate nucleotides

        Args:
            sequence (str): primer sequence

        Returns:
            float: return GC contents
        """
        gc_val_min = 0
        gc_val_max = 0
        if len(sequence) > 0:
            gc_val_min = round(100 * self.e_gc_min / len(sequence), 2)
            if self.choice == "min":
                gc_val_max = gc_val_min
            else:
                gc_val_max = round(100 * self.e_gc_max / len(sequence), 2)
        return gc_val_min, gc_val_max

    def cal_adj_melting_temp(self, salt_concentration: float = 50.0) -> float:
        """
        The basic calculation is
        For < 14bp sequences
        Tm= (wA+xT)*2 + (yG+zC)*4 - 16.6*log10(0.050) + 16.6*log10([Na+])
        where w,x,y,z are the number of the bases A,T,G,C in the sequence, respectively.
        The concentration of [Na+] might be micro mole.

        for longer sequences
        Tm= 100.5 + (41 * (yG+zC)/(wA+xT+yG+zC)) - (820/(wA+xT+yG+zC)) + 16.6*log10([Na+])

        이 수식에서 16.6* log10(0.05) 가 -21.6 이고
        log 가 자연로그이다 보니 수식에 맞추어 log10 으로 변화해주기 위해서
        7.21/ln(10.0) 을 해주어서 7.21 이 나오게 되는 것임.

        Args:
            choice (str): when sequence has degenerate nucleotide(s) then set `max` otherwise `min`
            sequence (str): primer sequence itself
            saltConcetration (float): Na+ mM Salt concentration

        Returns:
            float: adjusted melting temperature for a given primer sequence
        """
        if len(self.sequence) < 14:
            res = self.cal_adj_melting_temp_short(self.choice, self.sequence, salt_concentration)
        else:
            res = self.cal_adj_melting_temp_long(self.choice, self.sequence, salt_concentration)
        return res

    def count_neighbors(self, seek_dimer: str, target_seq: str) -> int:
        res = 0
        i = 0
        while i >= 0 and i < len(target_seq):
            offset = target_seq[i:].find(seek_dimer)
            if offset >= 0:
                res += 1
                i = i + 1 + offset
            else:
                i += 1

        return res

    def IsIUpacBase(self, theBase: str) -> bool:
        """if theBase is a degenerate base return True

        Args:
            theBase (str): one base

        Returns:
            bool: True if the base is one of degenerate bases
        """
        if theBase in ["M", "R", "W", "S", "Y", "K", "V", "H", "D", "B", "N"]:
            return True
        return False

    def calc_IU_pair(self, base0, base, i, theSequence, choice):
        IUpacBase = ""
        pair1 = ""
        pair2 = ""
        temp1 = [0, 0, 0]  # new Array(0,0,0)
        temp2 = [0, 0, 0]  # new Array(0,0,0)
        reValue = [0, 0, 0]  # new Array(0,0,0)
        base2 = theSequence[i]
        if self.IsIUpacBase(base0):  # if previous base is IUpacBase, do nothing
            return reValue

        if self.IsIUpacBase(base):
            # if (debug) alert("base0 "+base0+"base"+base+" base2 "+base2);
            if base == "M":
                IUpacBase = "AC"
            elif base == "R":
                IUpacBase = "AG"
            elif base == "W":
                IUpacBase = "AT"
            elif base == "S":
                IUpacBase = "CG"
            elif base == "Y":
                IUpacBase = "CT"
            elif base == "K":
                IUpacBase = "GT"
            elif base == "V":
                IUpacBase = "ACG"
            elif base == "H":
                IUpacBase = "ACT"
            elif base == "D":
                IUpacBase = "AGT"
            elif base == "B":
                IUpacBase = "CGT"
            elif base == "N":
                IUpacBase = "ACGT"

            j = 0
            while IUpacBase[j] != "":
                base = IUpacBase[j]
                # //alert("base choose "+base);
                pair1 = base0 + base
                # //	alert("pair1 "+pair1);
                if pair1 == "AA":
                    temp1[0] = 1.2
                    temp1[1] = 8.0
                    temp1[2] = 21.9
                elif pair1 == "AT":
                    temp1[0] = 0.9
                    temp1[1] = 5.6
                    temp1[2] = 15.2
                elif pair1 == "TA":
                    temp1[0] = 0.9
                    temp1[1] = 6.6
                    temp1[2] = 18.0
                elif pair1 == "CA":
                    temp1[0] = 1.7
                    temp1[1] = 8.2
                    temp1[2] = 21.0
                elif pair1 == "GT":
                    temp1[0] = 1.5
                    temp1[1] = 9.4
                    temp1[2] = 25.5
                elif pair1 == "CT":
                    temp1[0] = 1.5
                    temp1[1] = 6.6
                    temp1[2] = 16.4
                elif pair1 == "GA":
                    temp1[0] = 1.5
                    temp1[1] = 8.8
                    temp1[2] = 23.5
                elif pair1 == "CG":
                    temp1[0] = 2.8
                    temp1[1] = 11.8
                    temp1[2] = 29.0
                elif pair1 == "GC":
                    temp1[0] = 2.3
                    temp1[1] = 10.5
                    temp1[2] = 26.4
                elif pair1 == "GG":
                    temp1[0] = 2.1
                    temp1[1] = 10.9
                    temp1[2] = 28.4

                if base2 == "":
                    for k in range(2):  # (k=0; k<2; k++)
                        temp2[k] = 0.0

                elif not self.IsIUpacBase(base2):
                    pair2 = base + base2
                    # //alert("pair2 "+pair2);
                    if pair2 == "AA":
                        temp2[0] = 1.2
                        temp2[1] = 8.0
                        temp2[2] = 21.9
                    elif pair2 == "AT":
                        temp2[0] = 0.9
                        temp2[1] = 5.6
                        temp2[2] = 15.2
                    elif pair2 == "TA":
                        temp2[0] = 0.9
                        temp2[1] = 6.6
                        temp2[2] = 18.4
                    elif pair2 == "CA":
                        temp2[0] = 1.7
                        temp2[1] = 8.2
                        temp2[2] = 21.0
                    elif pair2 == "GT":
                        temp2[0] = 1.5
                        temp2[1] = 9.4
                        temp2[2] = 25.5
                    elif pair2 == "CT":
                        temp2[0] = 1.5
                        temp2[1] = 6.6
                        temp2[2] = 16.4
                    elif pair2 == "GA":
                        temp2[0] = 1.5
                        temp2[1] = 8.8
                        temp2[2] = 23.5
                    elif pair2 == "CG":
                        temp2[0] = 2.8
                        temp2[1] = 11.8
                        temp2[2] = 29.0
                    elif pair2 == "GC":
                        temp2[0] = 2.3
                        temp2[1] = 10.5
                        temp2[2] = 26.4
                    elif pair2 == "GG":
                        temp2[0] = 2.1
                        temp2[1] = 10.9
                        temp2[2] = 28.4
                elif self.IsIUpacBase(base2):
                    base0 = base
                    base = base2
                    i += 1
                    temp2 = self.CalcIUpair(base0, base, i, theSequence, choice)
                    i = i - 1

                for k in range(3):  # (k=0;k<3;k++) :
                    if j == 0:
                        reValue[k] = temp1[k] + temp2[k]
                    else:
                        if (choice == "max") and (reValue[k] < temp1[k] + temp2[k]):
                            reValue[k] = temp1[k] + temp2[k]
                        elif (choice == "min") and (reValue[k] > temp1[k] + temp2[k]):
                            reValue[k] = temp1[k] + temp2[k]
                j += 1
        return reValue

    def DeltaG(self, ignore_seq_len: bool = False):
        """unit: Kcal/mol

        Returns:
            _type_: deltaG (Kcal/mol)
        """
        if len(self.sequence) > 7 or ignore_seq_len:
            val = -5.0
            # // Helix initiation Free Energy of 5 kcal.
            # // symmetry function: if symmetrical, subtract another 0.4
            val += 1.2 * self.aa_count
            val += 0.9 * self.at_count
            val += 0.9 * self.ta_count
            val += 1.7 * self.ca_count
            val += 1.5 * self.gt_count
            val += 1.5 * self.ct_count
            val += 1.5 * self.ga_count
            val += 2.8 * self.cg_count
            val += 2.3 * self.gc_count
            val += 2.1 * self.gg_count
            if self.choice == "min":
                val += self.IUpairVals_min[0]
            else:
                val += self.IUpairVals_max[0]
            return round(1000 * val, 2) / 1000
        return ""

    def DeltaH(self, ignore_short_len: bool = False):
        """unit: Kcal/mol

        Returns:
            _type_: deltaH (Kcal/mol)
        """
        if len(self.sequence) > 7 or ignore_short_len:
            val = 0.0
            if self.isDeoxy:
                val += 8.0 * self.aa_count
                val += 5.6 * self.at_count
                val += 6.6 * self.ta_count
                val += 8.2 * self.ca_count
                val += 9.4 * self.gt_count
                val += 6.6 * self.ct_count
                val += 8.8 * self.ga_count
                val += 11.8 * self.cg_count
                val += 10.5 * self.gc_count
                val += 10.9 * self.gg_count
            else:
                val += 6.8 * self.aa_count
                val += 9.38 * self.at_count
                val += 7.69 * self.ta_count
                val += 10.44 * self.ca_count
                val += 11.4 * self.gt_count
                val += 10.48 * self.ct_count
                val += 12.44 * self.ga_count
                val += 10.64 * self.cg_count
                val += 14.88 * self.gc_count
                val += 13.39 * self.gg_count
            if self.choice == "min":
                val += self.IUpairVals_min[1]
            else:
                val += self.IUpairVals_max[1]

            return round(1000 * val, 2) / 1000

        return ""

    def DeltaS(self, ignore_short_len: bool = False):
        """unit: cal/(K*mol) where K is Kelvin temperature

        Args:
            choice (_type_): _description_

        Returns:
            _type_: deltaS (cal/(K*mol)) where K is a Kelvin temperature
        """
        if len(self.sequence) > 7 or ignore_short_len:
            val = 0
            if self.isDeoxy:
                val += 21.9 * self.aa_count
                val += 15.2 * self.at_count
                val += 18.4 * self.ta_count
                val += 21.0 * self.ca_count
                val += 25.5 * self.gt_count
                val += 16.4 * self.ct_count
                val += 23.5 * self.ga_count
                val += 29.0 * self.cg_count
                val += 26.4 * self.gc_count
                val += 28.4 * self.gg_count
            else:
                val += 19.0 * self.aa_count
                val += 26.7 * self.at_count
                val += 20.5 * self.ta_count
                val += 26.9 * self.ca_count
                val += 29.5 * self.gt_count
                val += 27.1 * self.ct_count
                val += 32.5 * self.ga_count
                val += 26.7 * self.cg_count
                val += 36.9 * self.gc_count
                val += 32.7 * self.gg_count

            if self.choice == "min":
                val += self.IUpairVals_min[2]
            else:
                val += self.IUpairVals_max[2]
            return round(1000 * val, 2) / 1000
        return ""

    def cal_RlnK(self, primer_concentration: float) -> float:
        K = 10**9 * 1 / primer_concentration  # nanomoles to moles
        R = 1.987
        res = round(R * math.log(K), 3)
        return res

    def NearestNeighborTM(
        self,
        primer_concentration: float = 50,
        salt_concentration: float = 50,
        ignore_seq_len: bool = False,
    ) -> float:
        """For Nearest Neighbor calculation
        deltaH=deltaG + TdeltaS
        delta G =  RT ln ([DNA-primer complex]/([DNA][primer]))
        T = deltaH/(deltaS + R ln([DNA-primer complex]/([DNA][primer])))  + 7.21* ln [Na concentration])
            Assume [DNA-primer complex] = [DNA] so we only need primer concentration
            NOTE: Python math.log() calls a NATURAL LOG, not Log base 10!

        Args:
            primer_concentration (float): nM primer
            salt_concentration (float): mM salt Na+
        Return:
            float : NN Tm
        """
        res = ""
        if len(self.sequence) > 7 or ignore_seq_len:
            K = 1 / (primer_concentration * 1e-9)  # Convert from nanomoles to moles
            R = 1.987  # cal/(mole*K);
            RlnK = R * math.log(K)  # natural log (ln)
            self.RlogK = round(1000 * RlnK, 2) / 1000
            # // Helix initiation Free Energy of 3.4 kcal (Sugimoto et al, 1996)
            # // symmetry function: if symmetrical, subtract another 0.4 - not implemented
            res = 1000 * (
                (self.DeltaH(ignore_short_len=ignore_seq_len) - 3.4)
                / (self.DeltaS(ignore_short_len=ignore_seq_len) + RlnK)
            )
            res += -272.9  # // Kelvin to Centigrade
            res += 7.21 * math.log(salt_concentration / 1000)  # milimoles to moles
            res = round(res, 2)
        else:
            self.RlogK = ""
        return res


if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser(
        prog="cal_melting_temperature.py",
        description="calculate Tm of a given sequence",
    )
    parser.add_argument("sequence", type=str, help="a primer sequence")
    parser.add_argument(
        "--type",
        choices=["adj", "nn"],
        help="{adj|nn} salt(Na+) adjusted Tm or nearest neighbor Tm",
    )
    parser.add_argument("--na", type=float, help="[Na+] mM. default=50(mM)", default=50)
    parser.add_argument(
        "--primer_conc",
        type=float,
        default=50,
        help="primer concentration(nM). default=50(nM). required for nn",
    )

    if len(sys.argv[:1]) < 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(sys.argv[1:])
    runner = MeltingTemp(isdeoxy=True)
    runner.set_sequence(args.sequence)
    if args.type == "adj":
        tmval = runner.cal_adj_melting_temp(salt_concentration=args.na)
    elif args.type == "nn":
        tmval = runner.NearestNeighborTM(
            primer_concentration=args.primer_conc, salt_concentration=args.na
        )
    print(f"{args.sequence} Tm: {tmval}")
