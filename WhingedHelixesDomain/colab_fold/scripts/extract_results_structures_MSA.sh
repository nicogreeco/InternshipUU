#!/bin/bash

protein="$1"
results_directory="/home/nicola/internship/colabfold-1.5.5/results_MSA_complete_seq/${protein}"
copy_folder="/home/nicola/internship/colabfold-1.5.5/to_download"




# Ensure the copy folder exists
mkdir -p "$copy_folder"/"$1"
rm -rf "$copy_folder"/"$1"/*
mkdir -p "$copy_folder"/"$1"/json

# Iterate over subdirectories in the results directory
for protein_folder in "$results_directory"/max_*; do
    if [ -d "$protein_folder" ]; then  # Check if it's a directory
        folder_name=$(basename "$protein_folder")
        # Copy files containing the specific string in their names and rename them
        for file in "$protein_folder"/*"_relaxed_rank_001_"*; do
            # echo "Found file: $file"
            if [ -f "$file" ]; then  # Check if it is a file
                new_name="${folder_name}_"$2"_1.pdb"  # Modify this line to change the naming rule
                # echo "Copying $file to $copy_folder/$new_name"
                cp "$file" "$copy_folder/"$1"/$new_name"
            fi
        done
        for file in "$protein_folder"/*"_scores_rank_001_"*; do
            # echo "Found file: $file"
            if [ -f "$file" ]; then  # Check if it is a file
                new_name="${folder_name}_"$2"_1.json"  # Modify this line to change the naming rule
                # echo "Copying $file to $copy_folder/$new_name"
                mkdir -p "$copy_folder"/json
                cp "$file" "$copy_folder/"$1"/json/$new_name"
            fi
        done
    fi
done

