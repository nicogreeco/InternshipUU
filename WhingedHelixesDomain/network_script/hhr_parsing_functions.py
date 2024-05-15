
import os
from tqdm import tqdm
import numpy as np
import math
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import textwrap
import re

# ----- FUNCTIONS -----

## EXTRACTION OF .hhm HITS OF QUERY SEARCH

# Saves a profile's data to a file within a specified directory.
# - profile_name: The name of the profile to save.
# - output_lines: The content of the profile to write to the file.
# - save_directory: The directory where the file will be saved.

def save_profile(profile_name, output_lines, save_directory):
    # Extract the first word of the profile name as the filename
    filename = f"{profile_name}.hhm"
    # Full path to where the file will be saved
    full_path = os.path.join(save_directory, filename)
    with open(full_path, 'w') as output_file:
        output_file.writelines(output_lines)

def cut_and_save_profile(profile_name, save_directory, range_str):

    # Full path to where the hmm file was saved
    full_hhm_path = os.path.join(save_directory, f"{profile_name}.hhm")

    # Full path to where the cutted hmm file will be saved
    output_file_path = os.path.join(save_directory, f"{profile_name}_{range_str}.hhm")

    extract_hmm_range(full_hhm_path, output_file_path, range_str)



# Extracts and saves specified HMM profiles from a database.
# - profile_names: A list of profile names to extract.
# - save_directory: Directory to save the extracted profiles.
# - database_path: Path to the database file containing HMM profiles.

def extract_hmm_profiles(profile_names, save_directory, database_path="/home/nicola/internship/hhsuite/databases/COG_KOG/COG_KOG_hhm.ffdata"):
    
    # Ensure the save directory exists, create if it does not
    os.makedirs(save_directory, exist_ok=True)

    # Markers to identify the start and end of profiles in the database.
    start_marker = "HHsearch"
    end_marker = "//"
    profile_found = False
    output_lines = []
    profile_name = ""

    # First, count the total number of lines for the progress bar (optional, might be skipped for very large files)
    total_lines = sum(1 for _ in open(database_path, 'r'))

    with open(database_path, 'r') as file:
        for line in tqdm(file, total=total_lines, desc="Processing file"):
            if line.startswith(start_marker):                  
                output_lines = [line]  # Start a new profile
            elif line.startswith("NAME ") and line.split()[1] in profile_names:
                profile_name = line.split()[1]
                profile_found = True  # The correct profile has been found
                output_lines.append(line)
            elif line.startswith(end_marker):
                if profile_found:
                    # Finish writing the previous profile and save it
                    output_lines.append(line)
                    save_profile(profile_name, output_lines, save_directory)
                    profile_names.remove(profile_name)
                    # Resets variables for the next profile.
                    output_lines = []
                    profile_found = False
                    profile_name = ""
            elif profile_found:
                output_lines.append(line)

    if len(profile_names) > 0:
        for profile in profile_names:
            print(f"Profile with name {profile} not found.")


def extract_hmm_profiles_and_cut_range(dict_of_profiles, save_directory, database_path="/home/nicola/internship/hhsuite/databases/COG_KOG/COG_KOG_hhm.ffdata"):
    
    # Ensure the save directory exists, create if it does not
    os.makedirs(save_directory, exist_ok=True)

    # Markers to identify the start and end of profiles in the database.
    start_marker = "HHsearch"
    end_marker = "//"
    profile_found = False
    output_lines = []
    profile_name = ""
    profile_names=list(dict_of_profiles.keys())
    # First, count the total number of lines for the progress bar (optional, might be skipped for very large files)
    total_lines = sum(1 for _ in open(database_path, 'r'))

    with open(database_path, 'r') as file:
        for line in tqdm(file, total=total_lines, desc="Processing file"):
            if line.startswith(start_marker):                  
                output_lines = [line]  # Start a new profile
            elif line.startswith("NAME ") and line.split()[1] in profile_names:
                profile_name = line.split()[1]
                profile_found = True  # The correct profile has been found
                output_lines.append(line)
            elif line.startswith(end_marker):
                if profile_found:
                    # Finish writing the previous profile and save it
                    output_lines.append(line)
                    save_profile(profile_name, output_lines, save_directory)
                    for range_str in dict_of_profiles[profile_name]:
                        cut_and_save_profile(profile_name, save_directory, range_str)
                    full_path_to_complete_hhm = os.path.join(save_directory,  f"{profile_name}.hhm")
                    os.remove(full_path_to_complete_hhm)
                    profile_names.remove(profile_name)
                    # Resets variables for the next profile.
                    output_lines = []
                    profile_found = False
                    profile_name = ""
            elif profile_found:
                output_lines.append(line)

    if len(profile_names) > 0:
        for profile in profile_names:
            print(f"Profile with name {profile} not found.")


# Extracts hit names from a specified file.
# - file_path: Path to the file from which to extract hits.
# Returns a list of extracted hit names.

def extract_hits_names(file_path, e_cutoff=math.inf, p_cutoff=0, score_cutoff=0):

    # positions of the variables in the results table
    positions = [
        (4, 11),  # Name
        (35, 40), # Prob
        (41, 48), # E-value
        (49, 56), # P-value
        (57, 63), # Score
        (64, 69), # SS
        (70, 74), # Cols
        (75, 84), # Query HMM start-end
        (85, 100) # Template HMM start-end
    ]
    

    start_table = " No Hit                             Prob E-value P-value"
    hits_names = []
    is_table_started = False

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith(start_table):                  
                is_table_started = True # Marks the start of the Table
                continue # Skip the header line itself
            elif is_table_started:
                if line.startswith("No 1"):
                    break # Stop at the first alignment (table is finished)
                split_line = line.split()

                # Check if the line is empty
                if len(split_line)>1: 
                    # Extract attributes based on positions
                    attributes = [line[start:end].strip() for start, end in positions]
                    if float(attributes[1])>p_cutoff and float(attributes[2])<e_cutoff and float(attributes[4])>score_cutoff:
                        hits_names.append(attributes[0]) # Add the name of the hit to the list
    return hits_names




def extract_hits_names_and_range(file_path, hits_dict, e_cutoff=math.inf, p_cutoff=0, score_cutoff=0, cols_cutoff=0):

    # positions of the variables in the results table
    positions = [
        (4, 11),  # Name
        (35, 40), # Prob
        (41, 48), # E-value
        (49, 56), # P-value
        (57, 63), # Score
        (64, 69), # SS
        (70, 74), # Cols
        (75, 84), # Query HMM start-end
        (85, 100) # Template HMM start-end
    ]


    start_table = " No Hit                             Prob E-value P-value"
    is_table_started = False

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith(start_table):                  
                is_table_started = True # Marks the start of the Table
                continue # Skip the header line itself
            elif is_table_started:
                if line.startswith("No 1"):
                    break # Stop at the first alignment (table is finished)
                split_line = line.split()

                # Check if the line is empty
                if len(split_line)>1: 
                    # Extract attributes based on positions
                    attributes = [line[start:end].strip() for start, end in positions]
                    if float(attributes[1])>p_cutoff and float(attributes[2])<e_cutoff and float(attributes[4])>score_cutoff and int(attributes[6])>cols_cutoff:
                        hits_dict[attributes[0]].append(attributes[8].split('(')[0].strip())
    return hits_dict


def parse_range(range_str):
    """ Parse a range string 'start-end' into a tuple of integers (start, end). """
    start, end = map(int, range_str.split('-'))
    return start, end

def merge_ranges(ranges):
    """ Merge overlapping ranges into the smallest possible number of ranges. """
    # Sort ranges by the start position
    ranges = sorted(ranges, key=lambda x: x[0])
    
    # Initialize the list of merged ranges
    merged_ranges = []
    
    for current_range in ranges:
        if not merged_ranges:
            merged_ranges.append(current_range)
        else:
            last_range = merged_ranges[-1]
            
            # Check for overlap and merge if necessary
            last_start, last_end = last_range
            current_start, current_end = current_range
            if current_start <= last_end:  # There is an overlap
                overlap_length = min(last_end, current_end) - current_start
                # Use the minimum length of the two ranges to calculate the overlap percentage
                total_length = min(current_end - current_start, last_end - last_start)
                if overlap_length / total_length >= 0.6:
                    new_end = max(last_end, current_end)
                    merged_ranges[-1] = (last_start, new_end)  # Update the last range
                else:
                    merged_ranges.append(current_range)
            else:
                merged_ranges.append(current_range)
    
    return merged_ranges


def process_ranges(ranges_dict):
    """ Process each key in the dictionary to merge overlapping ranges. """
    result_dict = {}
    for key, range_list in ranges_dict.items():
        # Parse the ranges into tuples of (start, end)
        parsed_ranges = [parse_range(r) for r in range_list]
        # Merge the ranges
        merged_ranges = merge_ranges(parsed_ranges)
        # Convert ranges back to the original format
        result_dict[key] = [f'{start}-{end}' for start, end in merged_ranges]
    return result_dict


def adjust_ranges(ranges_dict, min_adjustment=10):
    """
    Adjust each range in the dictionary by increasing both the start and end by a given amount.
    The start will not go below zero.

    Args:
    ranges_dict (dict): Dictionary with keys and list of range strings as values.
    adjustment (int): Max amount of the range.

    Returns:
    dict: Updated dictionary with adjusted ranges.
    """
    adjusted_dict = {}
    for key, ranges in ranges_dict.items():
        adjusted_ranges = []
        for range_str in ranges:
            start, end = map(int, range_str.split('-'))
            range_len = end-start
            range_diff = 70 - range_len
            adjustment = range_diff//2
            adjustment = max(adjustment, min_adjustment)
            # Increase the start and end by 'adjustment', ensuring start does not go below zero
            new_start = max(1, start - adjustment)
            new_end = end + adjustment
            adjusted_ranges.append(f'{new_start}-{new_end}')
        adjusted_dict[key] = adjusted_ranges
    return adjusted_dict


## QUERY REUNDANCY - KOG AND COG


# Removes files from the specified directory if their base names (without extensions) match any name in the given list.
# This function is specifically used to clean up the directory by removing COG and KOG profiles that refer to input MSA.
# - directory: The directory from which files will be removed.
# - list_of_COG: A list containing the base names of COG and KOG profiles to be removed.

def remove_COG_files_from_dir_by_names(directory, list_of_COG):
# Remove frpm the folder with .hhm and .hhr the COG and KOG that referes to input MSA
    # List all files in the directory
    for filename in os.listdir(directory):
        # Extract the base name (without extension) of the file
        base_name = os.path.splitext(filename)[0]
        
        # Check if the base name is in the list of names to remove
        if base_name in list_of_COG:
            # Construct the full path to the file
            file_path = os.path.join(directory, filename)
            
            # Remove the file
            os.remove(file_path)
            print(f"Removed file: {file_path}")


# Substitutes a segment of the original string with a given replacement string ('prot'), ensuring the substituted segment
# has a specific length by padding with spaces if necessary. This function is used to substitute the COG and KOG profile name
# of profiles refering to the Query proteins with the protein names.
# - original_string: The line to be modified.
# - prot: The protein name that will replace a segment of the original string.
# Returns the modified string with the specified substitution.

def substitute(original_string, prot): 
    # Define the start and end positions for replacement
    start_pos = 4
    end_pos = 11

    # Calculate the length of the replacement segment
    segment_length = end_pos - start_pos

    # Adjust prot to match the segment length by padding with spaces if necessary
    prot_padded = prot + ' ' * (segment_length - len(prot) + 1)

    # Construct the new string
    new_string = original_string[:start_pos] + prot_padded + original_string[end_pos + 1:]

    return new_string



# Modifies lines in a given file based on a list of OrthoGroups (OG) to be substituted with corresponding Query protein names.
# This function is used to substitute parts of lines where the OG code matches an entry in the given list, using the provided mapping.
# - file_path: Path to the original file containing lines to be modified.
# - save_path: Path where the modified lines will be saved.
# - og_list: A list of KOG and COG OrthoGroups (OG) whose occurrences should be substituted.
# - og_to_prot: A dictionary mapping OGs to their corresponding protein.

def substitute_OG_with_prot(file_path, save_path, og_list, og_to_prot):
    start_table = " No Hit                             Prob E-value P-value"
    is_table_started = False
    modified_lines = []  # Initialize a list to hold all lines, modified or not

    with open(file_path, 'r') as file:
        lines = file.readlines()  # Read all lines into memory

    for line in lines:
        if line.startswith(start_table):
            is_table_started = True  # Marks the start of the Table
            modified_lines.append(line)  # Add the header line to the list
            continue  # Skip the header line itself
        elif is_table_started and line.startswith("No 1"):
            is_table_started=False
        elif is_table_started:
            split_line = line.split()
            if len(split_line) > 1 and split_line[1] in og_list: 
                prot = og_to_prot[split_line[1]]
                line = substitute(line, prot)  # Substitute the part of the line
        modified_lines.append(line)  # Add the (possibly modified) line to the list

    # Write the modified lines to the new file
    with open(save_path, 'w') as new_file:
        new_file.writelines(modified_lines)


# Modifies lines in a given file based on a list of OrthoGroups (OG) to be removed.
# This function is used to remove lines where the OG code matches an entry in the given list.
# - file_path: Path to the original file containing lines to be modified.
# - save_path: Path where the modified lines will be saved.
# - og_list: A list of KOG and COG OrthoGroups (OG) whose occurrences should be removed.

def remove_OG_of_from_list(file_path, save_path, og_list):
    start_table = " No Hit                             Prob E-value P-value"
    is_table_started = False
    modified_lines = []  # Initialize a list to hold all lines, modified or not

    with open(file_path, 'r') as file:
        lines = file.readlines()  # Read all lines into memory

    for line in lines:
        if line.startswith(start_table):
            is_table_started = True  # Marks the start of the Table
            modified_lines.append(line)  # Add the header line to the list
            continue  # Skip the header line itself
        elif is_table_started and line.startswith("No 1"):
            is_table_started=False
        elif is_table_started:
            split_line = line.split()
            if len(split_line) > 1 and split_line[1] in og_list: 
                continue
        modified_lines.append(line)  # Add the (possibly modified) line to the list

    # Create the directory if it does not exist before writing the file
    os.makedirs(os.path.dirname(save_path), exist_ok=True)

    # Write the modified lines to the new file
    with open(save_path, 'w') as new_file:
        new_file.writelines(modified_lines)

# Removes lines from a file that represent self-hits, where the best hit of a query is the query itself.
# This function is used to clean up hhsearch results by excluding self-referential hits.
# - file_path: Path to the file from which self-hits will be removed.

def remove_self_hit(file_path, save_path=None):
    if save_path is None:
        save_path = file_path
    start_table = " No Hit                             Prob E-value P-value"
    is_table_started = False
    modified_lines = []  # Initialize a list to hold all lines, modified or not
    file_prot = os.path.basename(file_path).split('.')[0] # Extract the name of query from filename
    file_prot = file_prot.split('_')[0] # Remove the range if present
    with open(file_path, 'r') as file:
        lines = file.readlines()  # Read all lines into memory

    for line in lines:
        if line.startswith(start_table):                  
            is_table_started = True  # Marks the start of the Table
        elif is_table_started and line.startswith("No 1"):
            is_table_started=False
        elif is_table_started:
            split_line = line.split()
            if len(split_line) > 1 and split_line[1]==file_prot: 
                continue
        modified_lines.append(line)  # Add the (possibly modified) line to the list
    
    with open(save_path, 'w') as file:
        file.writelines(modified_lines)


# Removes lines from a file that do not correspond to profiles hit in the first round of query searches.
# This function is used to filter the hhsearch results, keeping only hits that were also found in an initial search.
# - file_path: Path to the file from which lines will be removed.
# - list_of_hitted_profile_in_query_search: A list of profile names that were hit in the first query search.

def remove_OG_not_in_first_output(file_path, list_of_hitted_profile_in_query_search):
    # Get the names of the profile that wher hitted in the first query search    
    start_table = " No Hit                             Prob E-value P-value"
    is_table_started = False
    modified_lines = []  # Initialize a list to hold all lines, modified or not
    file_prot = os.path.basename(file_path).split('.')[0] # Extract the name of query from filename

    with open(file_path, 'r') as file:
        lines = file.readlines()  # Read all lines into memory
    
    for line in lines:
        if line.startswith(start_table):                  
            is_table_started = True  # Marks the start of the Table
        elif is_table_started and line.startswith("No 1"):
            is_table_started=False
        elif is_table_started:
            split_line = line.split()
            if len(split_line) > 1 and split_line[1]  not in list_of_hitted_profile_in_query_search: 
                continue
        modified_lines.append(line)  # Add the (possibly modified) line to the list
    
    with open(file_path, 'w') as file:
        file.writelines(modified_lines)


import time

# Counts the number of files in he input dorectory every 30 sec

def count_files_in_directory(directory):
    """Counts the number of files in the specified directory."""
    try:
        # List directory contents and filter out directories, counting only files
        num_files = len([f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))])
        return num_files
    except FileNotFoundError:
        print("The specified directory was not found.")
        return None

def monitor_directory(directory, interval=30):
    """Monitors the specified directory, reporting the file count at specified intervals."""
    while True:
        num_files = count_files_in_directory(directory)
        if num_files is not None:
            print(f"{time.ctime()}: {num_files} files present in '{directory}'")
        else:
            # Exit the loop if the directory was not found
            break
        time.sleep(interval)  # Wait for specified interval (in seconds) before next check


def cut_list_of_sequences(lines_file_sequences, start, end):
    # Convert list of lines into a single string formatted as FASTA file content
    fasta_content = "".join(lines_file_sequences)

    # Parse the fasta content using SeqIO
    records = list(SeqIO.parse(StringIO(fasta_content), 'fasta'))

    # Slice each record to the desired range (adjust for zero-based index)
    sliced_records = []
    for record in records:
        # Make sure the range is within the sequence length
        if len(record.seq) >= end:
            sliced_seq = record.seq[start-1:end]  # Convert 1-based index to 0-based
            sliced_record = SeqRecord(sliced_seq, id=record.id, description=record.description)
            sliced_records.append(sliced_record)
        else:
            sliced_seq = record.seq[start-1:]  # Convert 1-based index to 0-based
            sliced_record = SeqRecord(sliced_seq, id=record.id, description=record.description)
            sliced_records.append(sliced_record)
    return sliced_records

def process_NAME_line(NAME_line, range_str):
    NAME_split=NAME_line.split()
    NAME_split[1]=f'{NAME_split[1]}_{range_str}'
    NAME_line=NAME_split[0]+'  '+' '.join(NAME_split[1:])+'\n'
    return NAME_line

def format_seq_records(records, line_length=100):
    """Format sequence records to ensure specific line length in FASTA."""
    formatted_records = []
    for record in records:
        # Wrap the sequence string to the specified line length
        wrapped_sequence = "\n".join(textwrap.wrap(str(record.seq), line_length))
        # Create a new SeqRecord with the wrapped sequence
        new_record = SeqRecord(Seq(wrapped_sequence), id=record.id, description=record.description)
        formatted_records.append(new_record)
    return formatted_records

def write_formatted_fasta(records, output, line_length=100):
    """Write formatted SeqRecord objects to a FASTA file."""
    formatted_records = format_seq_records(records, line_length)
    for record in formatted_records:
        # Write the record header
        output.write(f">{record.id} {record.description}\n")
        # Write the formatted sequence
        output.write(f"{record.seq}\n")

def extract_hmm_range(hmm_file_path, output_file_path, range_str):
    # Parse the range
    start, end = map(int, range_str.split('-'))
    cols = end - start + 1
    with open(hmm_file_path, 'r') as file, open(output_file_path, 'w') as output:
        lines_file = file.readlines() # Read all the lines from the hmm file
        metadata_lines = lines_file[:9]
        metadata_lines[0] = process_NAME_line(metadata_lines[0], range_str)
        output.writelines(metadata_lines) # Write the metadata lines
        lines_file_sequences = [] # Initialize all the lines-storing lists
        lines_file_profile =[]
        are_in_sequences = True # Initialize tokens
        are_in_profile = False # Initialize
        pattern = r'^[A-Z] (\d{1,4}) '
        record_count = 0 # Initialize counter for the sequences (I want just to write the first 4 sequences: SS, ss_conf, Consensus, Main seq)
        null_hmm_lines_count = 0 # Initialize the counter for the null frequencies and hhm header lines in the profile (the first 4 lines have to be always written)
        are_in_range = False
        for line in lines_file[9:]:
            if line.strip() == '#': # Check the presence of the #-line that marks the start of the profile hmm
                are_in_profile = True
                continue
            if are_in_sequences: # Processing the sequences
                if line.startswith('>'): # When we have processed already 4 sequences, write a #-line and exit the sequences parsing
                    if record_count == 4: # Stop processing if already four sequences are written 
                        # lines_file_sequences.append('#\n')
                        are_in_sequences = False
                        continue 
                    record_count += 1
                lines_file_sequences.append(line)
                continue
            if are_in_profile: # Processing the sequences
                if null_hmm_lines_count < 4: # Writing the first four lines of the profile 
                    lines_file_profile.append(line)
                    null_hmm_lines_count += 1
                    continue
                else:
                    match = re.match(pattern, line)
                    if match:
                        number_position = int(match.group(1))
                        are_in_range = start <= number_position <= end
                    if are_in_range:
                        lines_file_profile.append(line)
        if not lines_file_profile[-1]=='//\n':
            lines_file_profile.append('//\n')

        lines_file_sequences_cutted = cut_list_of_sequences(lines_file_sequences, start, end)
        write_formatted_fasta(lines_file_sequences_cutted, output)
        output.write('#\n')
        output.writelines(lines_file_profile)