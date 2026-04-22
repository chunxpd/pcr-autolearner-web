#!/bin/bash

# Check if the input file argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <input_seq_file>"
    exit 1
fi

# Check if the input filename ends with ".fq"
if [[ ${1##*.} == "fq" ]]; then
    IN_FQ="$1"  # Use the entire input filename if it's a FASTQ file
else
    IN_FA="$1"  # Use the entire input filename if it's a FASTA file

    # Extract base filename without extension (handle variations)
    base_filename="${IN_FA%.*}"
    IN_FQ="${base_filename}.fq" # Create the output FASTQ filename

    # Convert FASTA to FASTQ with quality scores as '#' (adjust as needed)
    seqtk seq -F '#' "$IN_FA" > "$IN_FQ"
fi

base_fq_name=`basename $IN_FQ`
OUT_BAM=tempdir/${base_fq_name%.fq}.bam
OUT_SAM=tempdir/${base_fq_name%.fq}.sam

REF=references/ref_seqs100.fasta
### bam file output
# bowtie2 -a -p 2 -x $REF -U $IN_FQ | samtools view -bS - > $OUT_BAM
# samtools view -f 4 $OUT_BAM
### sam file output
bowtie2 -a -p 2 -x $REF -U $IN_FQ > $OUT_SAM

echo "OUTPUT file: $OUT_SAM"