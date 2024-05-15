#!/bin/bash


# Get the directory where the script is located
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# echo "Working dir: $script_dir"
# Directory containing .hhm files, relative to the script's location
# echo "Please enter the directory containing input .hhm files:"
# read input_dir
input_dir=$1
echo "Input dir: $input_dir"

# Directory where .hhm files will be saved, relative to the script's location
# echo "Please enter the path to the directory where .hhr files will be saved:"
# read output_dir
output_dir=$2
# Check if the output directory exists, create it if not
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

echo "Output dir: $output_dir"

databases_dir="/home/nicola/internship/hhsuite/databases" 
# Count total number of .hhm files
total_files=$(ls -1 "$input_dir"/*.hhm | wc -l)
count=0
echo "Starting processing of $total_files files..."

# Loop through all .hhm files in the input directory
for file in "$input_dir"/*.hhm; do

    # Extract the base name without the extension
    base_name=$(basename "$file" .hhm)

    echo "Searching $base_name ($percent% complete)."

    # Increment counter
    ((count++))

    # Calculate percentage of completion
    percent=$(echo "scale=3; $count/$total_files*100" | bc)

    # Define the output file path with .hhm extension
    output_file="$output_dir/${base_name}.hhr"

    # Check if the output file already exists
    if [ ! -f "$output_file" ]; then
        # Output file does not exist, run hhsearch
        hhsearch -cpu 46 -i "$file" -d $databases_dir/COG_KOG/COG_KOG -d $databases_dir/WHD_proteins_db/WHD_db -o "$output_file" -v 0
    else
        echo "Output file already exists: $base_name"
    fi
done
echo "Search completed."


