#!/usr/bin/env python
import sys
from Bio import SeqIO
from io import StringIO


def parse_description(desc: str):
    """download_genebank.py 에서 나오는 fasta 파일의 리코드들의 descripttion
    ex::
    CP053068.1:4612617-4614637 HSP:Mycobacterium avium subsp. paratuberculosis strain DSM 44135 chromosome, complete genome:4613617-4613637
    complement 는 reverse complement 이므로 y축 대칭으로 모두 음으로 바꾸고 g_end 를 0 으로 평행이동 시키면 된다
    """
    cont = desc.split(" ")
    g_reg = cont[0].split(":")[1]
    hit_reg = " ".join(cont[1:]).split(":")[-1]
    if "complement(" not in g_reg:
        g_start, g_end = g_reg.split("-")
        hit_start, hit_end = hit_reg.split("-")
        res = (int(hit_start) - int(g_start), int(hit_end) - int(g_start))
    else:
        g_end, g_start = g_reg.replace("complement(", "").split("-")
        hit_end, hit_start = hit_reg.split("-")
        g_end = -int(g_end)
        # g_start = -int(g_start)
        hit_end = -int(hit_end)
        hit_start = -int(hit_start)
        res = (hit_start - g_end, hit_end - g_end)
    # try:

    # except ValueError as ve:
    #     print(hit_reg, ve)
    #     raise ValueError(ve)

    return res


def word_distance(q1, q2):
    cont = [v1 != v2 for v1, v2 in zip(q1, q2)]
    return sum(cont)


if __name__ == "__main__":
    fa = sys.argv[1]
    queryseq = SeqIO.read(StringIO(f">qseq\n{sys.argv[2]}"), "fasta")
    complementaryq = sys.argv[3]
    if eval(complementaryq) == True:
        queryseq = queryseq.reverse_complement()
    fastafile = open(fa, "r")

    for rec in SeqIO.parse(fastafile, "fasta"):
        name = rec.id
        desc = rec.description
        try:
            hit_start, hit_end = parse_description(desc)
        except Exception as e:
            print(desc, e)
            raise e
        seq = rec.seq
        hsp_seq = seq[hit_start - 1 : hit_end]
        print(hsp_seq, queryseq.seq, word_distance(hsp_seq, queryseq.seq), "<- " + desc)

    fastafile.close()
