import sys
from Bio import SeqIO


def stats():
    fastain = sys.argv[1]
    with open(fastain) as handle:
        for rec in SeqIO.FastaIO.FastaIterator(handle):
            print(f"{len(rec.seq)} {rec.description}")


if __name__ == "__main__":
    stats()
