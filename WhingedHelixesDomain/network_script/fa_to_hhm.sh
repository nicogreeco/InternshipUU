#!/bin/bash


# Get the directory where the script is located
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Directory containing .fa files, relative to the script's location
input_dir="$script_dir/MSAs/fasta"

# Directory where .hhm files will be saved, relative to the script's location
output_dir="$script_dir/MSAs/hhm"

# Add to the PERL5LIB environment variable the path to hhsuite scripts directory so that addss.pl finds HHPaths
PERL5LIB=$PERL5LIB:/home/nicola/internship/hhsuite/hh-suite/build/scripts

# Path to hhsuite scripts directory
hhscript_dir="/home/nicola/internship/hhsuite/hh-suite/build/scripts"

# Ensure the output directory exists
mkdir -p "$output_dir"

# Loop through all .fa files in the input directory
for fa_file in "$input_dir"/*.fasta; do

    # Extract the base name without the extension
    base_name=$(basename "$fa_file" .fasta)

    # Extract the part of the base name before the first underscore
    short_name=${base_name%%_*}

    echo "Parsing ${short_name}"

    perl "$hhscript_dir/reformat.pl" fas a3m "$fa_file" "$script_dir/MSAs/a3m/${short_name}.a3m" -M 65

    sed -i "1s/^/#$short_name\n/" "$script_dir/MSAs/a3m/${short_name}.a3m"

    hhconsensus -i "$script_dir/MSAs/a3m/${short_name}.a3m" -o "$script_dir/MSAs/a3m/${short_name}.a3m"

    echo "Predicting secondary structure of ${fa_file}"

    perl "$hhscript_dir/addss.pl" "$script_dir/MSAs/a3m/${short_name}.a3m" -a3m

    echo "Making .hhm of ${short_name}"

    # Run hhmake on the input file and output to the defined output file
    hhmake -i "$script_dir/MSAs/a3m/${short_name}.a3m" -o "$output_dir/${short_name}.hhm"
done

echo "Conversion completed."