from Bio import SeqIO
from Bio.Seq import Seq

if __name__ == "__main__":
    import sys

    alnfn = sys.argv[1]
    fafn = sys.argv[2]
    contdict = dict()
    with open(fafn, "r") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            name = rec.id
            seq = rec.seq
            desc = rec.description.split(" ")[1:]
            desc = "_".join(desc).replace("HSP:", "")
            desc = desc[: desc.rfind(":")]
            contdict[name.split(":")[0]] = desc
    with open(alnfn, "r") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            name = rec.id.split(":")[0]
            seq = rec.seq
            desc = contdict[name]
            out_rec = SeqIO.SeqRecord(
                seq=Seq(seq), id=f"{rec.id}:{desc}", description=""
            )
            print(out_rec.format("fasta"), end="")
