#!/bin/bash

# List of proteins to process
proteins=("ask1" "ska1" "ska2" "ska3")

# Loop through each protein and run the analysis script with nohup
for protein in "${proteins[@]}"; do
  echo "Processing $protein"
  nohup ./scripts/run_alphafold_MSA.sh "$protein" > "nohup_output_${protein}_MSA_complete_not_in_my.log"
done

echo "All proteins are processed in the background."


for protein in "${proteins[@]}"; do
  echo "Exporting $protein"
  nohup ./scripts/extract_results_structures_MSA.sh "$protein" "complete_seq_msa"
done


echo "All proteins are copied in download folder."