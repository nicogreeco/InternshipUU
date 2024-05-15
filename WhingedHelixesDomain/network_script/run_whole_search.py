from hhr_parsing_functions import *
from network_creation_functions import *
import subprocess
from collections import defaultdict


for number_of_search in [1,2,3]:

        ## ------- SEARCH 
    input_dir=f'./search_outputs_range/{number_of_search-1}_search/hits_hhm'
    output_dir=f'./search_outputs_range/{number_of_search}_search/not_parsed_hhr'
    if number_of_search != 1:
        subprocess.run(["./hhm_search.sh", input_dir, output_dir])

    ## ------- EXTRACTION OF HITS FROM not_parsed_hhr-- CUTTED PROFILES

    # Directory containing output files from HHsearch.
    hhsearch_output_dir = output_dir
    # hhsearch_output_dir=output_dir
    hits_dict = defaultdict(list)

    # Iterate over each HHsearch output file to extract and compile a list of unique hits in the fisrt search
    for filename in os.listdir(hhsearch_output_dir):
        hhsearch_output_file=os.path.join(hhsearch_output_dir, filename) # Create path to file
        hits_dict=extract_hits_names_and_range(hhsearch_output_file, hits_dict, p_cutoff=40, cols_cutoff=30) # Extract the hits name

    # Removes duplicates from the previous search
    previous_search_queries=[file.split('.')[0].split('_')[0] for file in os.listdir(hhsearch_output_dir)]
    print(f"In total there are {len(hits_dict)} unique hits, before removing the already searched queries")
    hits_dict={key: value for key, value in hits_dict.items() if key not in previous_search_queries} # Removes duplicates from the previous search
    print(f"In total there are {len(hits_dict)} unique hits")

    # Process the ranges
    processed_ranges = process_ranges(hits_dict)
    adjusted_ranges = adjust_ranges(processed_ranges)

    ## ------- EXTRACTION OF .hhm -- CUTTED PROFILES

    # Directory to save the extracted HMM profiles.
    save_directory=f'./search_outputs_range/{number_of_search}_search/hits_hhm'

    # Extracts and saves the HMM profiles for the list of unique hits.
    extract_hmm_profiles_and_cut_range(adjusted_ranges, save_directory, database_path="/home/nicola/internship/hhsuite/databases/COG_KOG/COG_KOG_hhm.ffdata")
