import sys
from Bio import SeqIO
from Bio.Seq import MutableSeq


def main():
    # Specify input and output file paths
    # input_file = "input.clustal"
    # output_file = "output.fasta"
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Read the Clustal file and write the sequences in FASTA format
    input_format = "clustal"
    output_format = "fasta"
    gap_character = "-"
    with open(input_file, "r") as input_handle, open(output_file, "w") as output_handle:
        sequences = SeqIO.parse(input_handle, input_format)
        for seq_record in sequences:
            # seq_record.seq = seq_record.seq.tomutable()  # Convert to mutable sequence
            seq_record.seq = MutableSeq(str(seq_record.seq))
            seq_record.seq = seq_record.seq.replace(gap_character, "")  # Replace gaps with '-'
            SeqIO.write(seq_record, output_handle, output_format)


if __name__ == "__main__":
    main()
    # print("Conversion complete.")
