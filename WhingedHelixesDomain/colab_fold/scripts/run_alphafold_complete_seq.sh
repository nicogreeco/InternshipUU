#!/bin/bash

# Define variables
protein="$1"

fasta_input="/home/nicola/internship/winghel/alignments/${protein}_curated_complete_seq.fasta"
splitted_fasta_dir="/home/nicola/internship/colabfold-1.5.5/input_fastas/${protein}/complete_seq"
output_dir="/home/nicola/internship/colabfold-1.5.5/results/${protein}"

# Create output directory if not exists
mkdir -p "${splitted_fasta_dir}"

# Split fasta into individual files
python -c "
import os
from Bio import SeqIO
fasta_input = '${fasta_input}'
directory_output = '${splitted_fasta_dir}'
for record in SeqIO.parse(fasta_input, 'fasta'):
    output_file = os.path.join(directory_output, record.id + '.fasta')
    with open(output_file, 'w') as output_handle:
        SeqIO.write(record, output_handle, 'fasta')
"

# Count total number of fasta files
total_files=$(ls -1 "${splitted_fasta_dir}"/*.fasta | wc -l)
counter=0

# Run ColabFold on each splitted fasta file
for file in "${splitted_fasta_dir}"/*.fasta; do
    ((counter++))
    echo "Processing ${counter} of ${total_files}"
    base_name=$(basename "$file" .fasta)
    output_dir_file="${output_dir}/${base_name}"
    mkdir -p "${output_dir_file}/results_complete_seq"
    /home/jankees-colabfold-155/colabfold-1.5.5/run_colabfold.sh --templates --amber --use-gpu-relax "$file" "${output_dir_file}/results_complete_seq" |& tee "${output_dir_file}/results_complete_seq/file.log"
done

echo "Completed processing ${total_files} files."
