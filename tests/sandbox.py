# %%
import Bio.Blast.NCBIWWW as blast

res = blast.qblast("blastn", "nt", "gatgctcgcaaccactatcca", short_query=True)
# with open("my_blast.xml","w") as handle:
#     handle.write(res.read())
# res.close()

# %%
### issue #27 의 multiple nt 를 사용하였을 때 결과 리턴 문제 테스트
import Bio.Blast.NCBIWWW as blast

res = blast.qblast("blastn", "nt", "TGCTSGACATCTACCGVAARC", short_query=True)
with open("issue27-blastn.xml", "w") as handle:
    handle.write(res.read())
res.close()
### 위 결과는 잘 나온다. issue #27 확인
# %%
### issue #27 의 multiple nt 를 사용하였을 때 결과 리턴 문제 테스트
### 현재 call_blastn.py 에서는 megablast=True 조건을 주고 있음. 이것이 문제인지 확인
import Bio.Blast.NCBIWWW as blast

res = blast.qblast(
    "blastn", "nt", "TGCTSGACATCTACCGVAARC", short_query=True, megablast=True
)
with open("issue27-blastn.xml", "w") as handle:
    handle.write(res.read())
res.close()
### megablast=True 를 하니 결과가 하나도 나오지 않는다. megablast 옵션때문임.
### megablast option 에 대해서 다시한번 확인해 본다.
### megablast=False 를 하니 결과가 많이 나온다.
### megablast=None 과 False 의 다른점은 없는 것 같다.
# %%
import Bio.Blast.NCBIXML as bxml

blast_records = bxml.parse(open("issue27-blastn.xml", "r"))
for blast_record in blast_records:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            print("****Alignment****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query)
            print(hsp.match)
            print(hsp.sbjct)

# %%

### taxid:1772 taxid:1773 exclude 조건 테스트
import Bio.Blast.NCBIWWW as blast

res = blast.qblast(
    "blastn",
    "nt",
    "gatgctcgcaaccactatcca",
    short_query=True,
    entrez_query="all[filter] NOT txid1772[ORGN] NOT txid1773[ORGN]",
)
with open("blast_exclude_example.xml", "w") as handle:
    handle.write(res.read())

# %%
### taxid:1772 에만 쿼리 테스트
import Bio.Blast.NCBIWWW as blast

res = blast.qblast(
    "blastn",
    "nt",
    "gatgctcgcaaccactatcca",
    short_query=True,
    entrez_query="(txid1772[ORGN])",
    hitlist_size=100,
    megablast=True,
)
with open("blast_include_example.xml", "w") as handle:
    handle.write(res.read())

# %%
import Bio.Blast.NCBIXML as bxml

blast_records = bxml.parse(open("blast_exclude_example.xml", "r"))
for blast_record in blast_records:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            print("****Alignment****")
            print("description:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query)
            print(hsp.match)
            print(hsp.sbjct)
            break
    break


# %%
### Entrez genebank sequence download
from Bio import Entrez

Entrez.email = "io2oi@yahoo.co.kr"  # Always tell NCBI who you are
iids = "CP101406.1,CP097320.1,CP076378.1,CP014860.2,CP092427.1,CP092362.1,CP092365.1,CP023138.1,CP081006.1"
handle = Entrez.efetch(db="nucleotide", id=iids, rettype="fasta", retmode="text")
with open("entrez_genbank.fasta", "w") as outfp:
    outfp.write(handle.read())
handle.close()

# %%
### gb 파일 한개씩 내려받을 때와 여러 개 내려받을 때 시간 차이 계산
import time
from Bio import Entrez
import os


def gb_download_test(gb_ids):
    gb_ids = gb_ids.split(",")
    n_gb_ids = len(gb_ids)
    start = time.time()
    os.makedirs("sandbox_files", exist_ok=True)
    for gbid in gb_ids:
        with open(f"sandbox_files/{gbid}.gb", "w") as handle:
            Entrez.email = "io2oi@yahoo.co.kr"
            enz_handle = Entrez.efetch(
                db="nucleotide", id=gbid, rettype="gb", retmode="text"
            )
            handle.write(enz_handle.read())
            enz_handle.close()
    print(f"개별 QUERY time ({n_gb_ids}): ", time.time() - start)

    gb_ids = ",".join(gb_ids)
    start = time.time()

    with open(f"sandbox_files/gb_ids.gb", "w") as handle:
        Entrez.email = "io2oi@yahoo.co.kr"
        enz_handle = Entrez.efetch(
            db="nucleotide", id=gb_ids, rettype="gb", retmode="text"
        )
        handle.write(enz_handle.read())
        enz_handle.close()
    print(f"동시 QUERY time ({n_gb_ids}): ", time.time() - start)


# %%
iids = "CP092362.1,CP092365.1,CP023138.1"
gb_download_test(iids)
iids = "CP101406.1,CP097320.1,CP076378.1,CP014860.2,CP092427.1"
gb_download_test(iids)
iids = "CP101406.1,CP097320.1,CP076378.1,CP014860.2,CP092427.1,CP092362.1,CP092365.1,CP023138.1,CP081006.1"
gb_download_test(iids)

# %%
### gb id 대신 gi id 도 같은 결과를 주는지 확인
# gi|328456706|gb|CP001662.1|
iids = "328456706,CP001662.1"
gb_download_test(iids)
# %%
iids = "XM_001613789.1"
gb_download_test(iids)
# %%
### BlastCommandline
### 아래 명령어는 그냥 command string 만 만들어 주는 것이다. 사용할 필요가 없다.
from Bio.Blast.Applications import NcbiblastnCommandline

cline = NcbiblastnCommandline(
    query="foo.fasta", db="nt", evalue=0.001, out="foo.xml", outfmt=5
)
print(cline)

# %%
### Entrez 에서 fasta 형식으로 정보를 가지고 올 수 있는 지 확인
from Bio import Entrez
from Bio import SeqIO
from io import StringIO

iid = "LR812045.1"
pred = (
    ">LR812045.1 Danio rerio genome assembly, chromosome: 8"
    "AGTTGACTCTACCCTAGTTGACTCTACCCTAGTTGACTCTACCCCTAGTTGACTCTACCCTAGTTGACTC"
    "TACCCCTAGTTGACTCTACCCTAGTTGACTCTACCCCTAGTTGACTCTACCCCTAGTTGACTCTACCCCT"
    "AGTTGACTCTACCCCTAGTTGACTCTACCCCCTAGTTGACTCTACCCCTAGTTGACTCTACCCCTAGTTG"
    "ACTCTACCCTAGTTGACTCTACCCCTAGTTGACTCTACCCTAGTTGACTCTACCCCTAGTTGACTCTACC"
    "CCGTATAGTTGACTCTACCCCTAGTTGACTCTACCCCTAGTTGACTCTACCCCTAGTTGACTCTACCCCT"
    "AGTTGACTCTACCCCTAGTTGACTCTACCCCTAGTTGACTCTACCCCTAGTTGACTCTACCCCTAGTTGA"
)
Entrez.email = "io2oi@yahoo.co.kr"
obs_res = Entrez.efetch(db="nucleotide", id=iid, rettype="fasta", retmode="text")
obs = SeqIO.read(StringIO(obs_res.read()), "fasta")
print(obs.id)
print(obs.seq[:10])
# %%
