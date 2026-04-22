#!/usr/bin/env python
import sys
from Bio import SeqIO


def print_length(fa_list: list):
    for fa in fa_list:
        with open(fa, "r") as FastaFile:
            for rec in SeqIO.parse(FastaFile, "fasta"):
                name = rec.id
                seq = rec.seq
                seqLen = len(rec)
                print(name, seqLen)


def print_ids(fa_list: list):
    for fa in fa_list:
        with open(fa, "r") as FastaFile:
            for rec in SeqIO.parse(FastaFile, "fasta"):
                print(rec.id)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Process subcommands.")
    subparsers = parser.add_subparsers(title="subcommands", dest="subcommand")

    len_parser = subparsers.add_parser("length", help="print length of each fa id")
    len_parser.add_argument("fa_in", nargs="+", help="input fasta files")

    id_parser = subparsers.add_parser("id", help="print id list of the given fasta file")
    id_parser.add_argument("fa_in", nargs="+", help="input fasta files")

    if len(sys.argv[1:]) < 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    if args.subcommand == "length":
        print_length(args.fa_in)
    elif args.subcommand == "id":
        print_ids(args.fa_in)
