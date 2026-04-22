#!/bin/bash
### https://github.com/io2oi/ctlabs-primerwrapper/issues/183 참고

# Check if the input file argument is provided
if [ $# -eq 0 ]; then
	echo "Usage: $0 <input_file>"
	exit 1
fi

DIRECTORY=references
input_file=$1
IIDS=""
FINAL=ref_seqs.fasta
FINAL100=ref_seqs100.fasta

# Create the DIRECTORY if it doesn't exist
[ ! -d "$DIRECTORY" ] && mkdir -p "$DIRECTORY"

cd $DIRECTORY
# Read each ID from the input file
while IFS= read -r id; do
	# Check if the ID starts with "GC"
	if [[ $id =~ ^GC ]]; then
		# Use NCBI Datasets for GC IDs
		echo "Downloading results for ID: $id (using NCBI Datasets)"
		datasets download genome accession $id --include genome
		unzip ncbi_dataset.zip
		mv ncbi_dataset/data/$id/*.fna ./
		rm -rf ncbi_dataset.zip ncbi_dataset README.md
	else
		IIDS="$IIDS $id"
		# Use esearch and efetch for non-GC IDs
		# echo "Downloading results for ID: $id (using esearch and efetch)"
		# esearch -db nucleotide -query "$id" | efetch -format fasta > "./$id.fasta"
	fi
#   echo "Downloaded results for ID: $id"
done < ../"$input_file"

echo "Downloading results for ID: $IIDS (using esearch and efetch)"
for iid in $IIDS; do
    esearch -db nucleotide -query "$iid" | efetch -format fasta # > "./$id.fasta"
done > "__temp__.fasta"

cat *.fna *.fasta *.fa > $FINAL
cdhit -i $FINAL -o $FINAL100 -c 1.0

files_to_delete=$(find . -type f -name "*" ! -name ref_seqs100.fasta)
rm $files_to_delete

echo "Downloaded results saved in DIRECTORY: $DIRECTORY"
