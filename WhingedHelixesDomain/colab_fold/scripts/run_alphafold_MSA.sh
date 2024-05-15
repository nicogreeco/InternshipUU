#!/bin/bash

# Define variables
protein="$1"
input_dir="/home/nicola/internship/colabfold-1.5.5/results_MSA_complete_seq/${protein}"

# Run ColabFold on each splitted fasta file
for folder in "${input_dir}"/max_*; do
    # Find the .a3m file in the folder
    for file in "${folder}"/*.a3m; do
        echo "Processing $file"
        # Run the prediction command only if the file exists
        if [ -f "$file" ]; then
            /home/jankees-colabfold-155/colabfold-1.5.5/run_colabfold.sh --templates --amber --use-gpu-relax "$file" "${folder}" |& tee "${folder}/file.log"
        fi
    done
done

echo "Completed processing files."
