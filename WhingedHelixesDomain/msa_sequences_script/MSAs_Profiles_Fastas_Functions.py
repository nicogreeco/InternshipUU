import pandas as pd
from Bio import SeqIO
import os
from io import StringIO
import subprocess
import random
import numpy as np



################## ------------------------------------------------------------------------------------------------- ##################
################################################        OLD FUNCTIONS       ###########################################################
################## ------------------------------------------------------------------------------------------------- ##################

# Creates a dictionary mapping species IDs to their respective proteome filenames in a specified directory.
def get_index_of_proteomes(path_to_folder):
    file_index = {file.split('_')[0]: file for file in os.listdir(path_to_folder)}
    return file_index

# Retrieves a BioPython SeqIO.index object for quick sequence access from a proteome file based on a species ID.
# This allows efficient retrieval of sequences from large fasta files by their IDs.
def get_sequence_index_from_id(id, path_to_database='/home/nicola/internship/winghel/EukProt/proteins'):
    file_index = get_index_of_proteomes(path_to_database)
    file = file_index[id]
    path = os.path.join(path_to_database, file)
    seq_index = SeqIO.index(path, 'fasta')
    return seq_index

# Parses the filename of a FASTA file to extract a concise identifier, 
# removing any prefixes or file extensions, and handling special cases like filenames starting with "DNA_".
def parse_file_name(fasta):
    fasta_parsed = fasta.split('.')[0]
    if fasta_parsed.startswith('DNA_'):
        fasta_parsed = fasta_parsed.replace('DNA_', '')
        fasta_parsed = fasta_parsed.split('_')[0]
    return fasta_parsed


################## ------------------------------------------------------------------------------------------------- ##################
#######################################################################################################################################
################## ------------------------------------------------------------------------------------------------- ##################


# Retrieves a BioPython SeqIO.index object for quick sequence access from a proteome file based on a species ID.
# This allows efficient retrieval of sequences from large fasta files by their IDs.
def get_eukarya5_index(path_to_database='/home/max/NOBINFBACKUP/eukarya5_v2/data_set/eukarya_proteomes_reduced99_lt.fa'):
    seq_index = SeqIO.index(path_to_database, 'fasta')
    return seq_index

# Given a query filename, retrieves matching sequence entries from a species' proteome.
# It leverages previously defined functions to efficiently locate and return these sequences.
def get_sequence_from_query(query_file, string=True):
    id = parse_file_name(query_file)
    try:
        matched_sequences = eukarya5_index[id]
    except KeyError:
        print(f"No sequence found for ID: {query_file}")
        return None
    if string:
        return str(matched_sequences.seq)
    else:
        return matched_sequences


def get_species_name(query_file):
    id = parse_file_name(query_file)
    species_id = eukarya_species[eukarya_species['Abbreviation'] == id]['Scientific name']
    return species_id.iloc[0]

def get_protein_id(query_file):
    id = parse_file_name(query_file)
    gene_id = eukarya_metadata[eukarya_metadata['protein_id'] == id]['original_protein_id']
    return gene_id.iloc[0]

def get_gene_id(query_file):
    id = parse_file_name(query_file)
    gene_id = eukarya_metadata[eukarya_metadata['protein_id'] == id]['original_gene_id']
    return gene_id.iloc[0]

def clean_seq(remove_X=True):
    input_str = input("Enter a string: ")
    cleaned_str = input_str.replace("\n", "").replace("-", "").replace(" ", "")
    if 'X' in cleaned_str and remove_X:
        index = cleaned_str.find('X')
        cleaned_str = cleaned_str.replace("X", "")
        print(f'An X was replaced replaced at position {index}')
    print(f'Sequence with length: {len(cleaned_str)}')
    return cleaned_str

def remove_gaps_from_alignment(input_file, output_file):
    # Open the output file in write mode
    with open(output_file, 'w') as output_handle:
        # Iterate over each record in the input FASTA file
        for record in SeqIO.parse(input_file, "fasta"):
            # Remove all gap characters from the sequence
            gapless_sequence = str(record.seq).replace('-', '')
            # Update the record's sequence with the gapless sequence
            record.seq = gapless_sequence
            # Write the updated record to the output file
            SeqIO.write(record, output_handle, "fasta")

def combine_sequences(basenames, directory, output_file):
    """
    Combines sequences from multiple FASTA files into a single output file.
    
    :param basenames: A list of basenames (filename without the .fasta extension) of the files to include.
    :param directory: The directory containing the FASTA files.
    :param output_file: The path to the output FASTA file to create.
    """

    sequences = []  # List to store sequences
    
    # Iterate over each basename and construct the full path to the FASTA file
    for basename in basenames:
        fasta_file = os.path.join(directory, basename + '.fasta')
        
        # Check if the FASTA file exists
        if os.path.isfile(fasta_file):
            # Read the sequence from the file and add it to the list
            for record in SeqIO.parse(fasta_file, "fasta"):
                sequences.append(record)
        else:
            print(f"File {fasta_file} not found.")
    
    # Write all collected sequences to the output file
    with open(output_file, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")
    
    print(f"All sequences have been combined into {output_file}")

def clean_fasta_seq(record, remove_X=True):
    cleaned_str = str(record.seq).replace("-", "").replace(" ", "")
    if 'X' in cleaned_str and remove_X:
        index = cleaned_str.find('X')
        cleaned_str = cleaned_str.replace("X", "")
        print(f'An X was replaced at position {index}')

    if '*' in cleaned_str and remove_X:
        index = cleaned_str.find('*')
        cleaned_str = cleaned_str.replace("*", "")
        print(f'An * was replaced at position {index}')
    return cleaned_str

def get_fasta_sequence(protein, sequence_id):
    fasta_file=f"/home/nicola/internship/winghel/fastas/{protein}/{protein}_conc.fasta"
    found = False
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == sequence_id:
            return(record.seq)
            found = True
            break
    if not found:
        return(f"Sequence ID {sequence_id} not found in the file {fasta_file}.")

def clean_fasta_sequences(fasta_input, is_path=False):
    if is_path:
        # Check if the file exists
        if not os.path.isfile(fasta_input):
            raise FileNotFoundError(f"No file found at the provided path: {fasta_input}")
        
        # Open the file and parse the sequences
        with open(fasta_input, 'r') as file:
            records = list(SeqIO.parse(file, "fasta"))
    else:
        # Parse the sequences from the string
        fasta_io = StringIO(fasta_input)
        records = list(SeqIO.parse(fasta_io, "fasta"))
    
    cleaned_fasta = ""
    for record in records:
        cleaned_sequence = clean_fasta_seq(record)
        record_id_modified = record.id.replace("/", " ")
        cleaned_fasta += f">{record_id_modified}\n{cleaned_sequence}\n"
    
    return cleaned_fasta

def get_genes_ids(fasta_string):
    fasta_io = StringIO(fasta_string)
    records = list(SeqIO.parse(fasta_io, "fasta"))
    gene_ids = ""

    for record in records:
        record_id_modified = record.id.split("/")[0]
        gene_id = get_gene_id(record_id_modified)
        gene_ids += f">{record_id_modified}\n{gene_id}\n"      

    return gene_ids


def get_fasta_sequences(protein, fasta_input, is_path=False):
    if is_path:
        # Check if the file exists
        if not os.path.isfile(fasta_input):
            raise FileNotFoundError(f"No file found at the provided path: {fasta_input}")
        
        # Open the file and parse the sequences
        with open(fasta_input, 'r') as file:
            records = list(SeqIO.parse(file, "fasta"))
    else:
        # Parse the sequences from the string
        fasta_io = StringIO(fasta_input)
        records = list(SeqIO.parse(fasta_io, "fasta"))
    
    fasta_sequences = ""
    for record in records:
        record_id_modified = record.id.split("/")[0]
        fasta_seq = get_fasta_sequence(protein, record_id_modified)
        fasta_sequences += f">{record_id_modified}\n{fasta_seq}\n"
    
    return fasta_sequences
################## ------------------------------------------------------------------------------------------------- ##################
##############################################        Random MSA search       #########################################################
################## ------------------------------------------------------------------------------------------------- ##################


def filter_fasta(fasta_path, ids_to_remove, output_path):
    # Read sequences from the original FASTA file
    sequences = SeqIO.parse(fasta_path, "fasta")
    
    # Filter out sequences whose IDs are in the ids_to_remove list
    filtered_sequences = [seq for seq in sequences if seq.id not in ids_to_remove]
    
    # Write the filtered sequences to a new FASTA file
    SeqIO.write(filtered_sequences, output_path, "fasta")

    
def select_random_sequences(fasta_file, number=4):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    selected = random.sample(sequences, number)
    with open("temp_selected.fasta", "w") as output_file:
        SeqIO.write(selected, output_file, "fasta")
    return "temp_selected.fasta"

def create_msa(input_fasta):
    msa_output = "temp_msa.fasta"
    subprocess.run(["mafft", "--globalpair", "--maxiterate", "1000", input_fasta], stdout=open(msa_output, "w"))
    return msa_output

def build_hmm(msa_fasta):
    hmm_output = "temp_profile.hmm"
    subprocess.run(["hmmbuild", hmm_output, msa_fasta])
    return hmm_output

def run_hmmsearch(hmm_profile, target_fasta):
    search_output = "temp_hmmsearch_results.txt"
    # Direct standard output to the file instead of using --tblout
    subprocess.run(["hmmsearch", hmm_profile, target_fasta], stdout=open(search_output, "w"))
    return search_output

def get_sequence_ids(fasta_file):
    sequence_ids = [record.id.split('/')[0] for record in SeqIO.parse(fasta_file, "fasta")]
    return sequence_ids

def get_hmm_search_hits(hmm_file_path):
    hit_id_list=[]
    hit_score_list=[]
    hit_Eval_list=[]
    table_started=False
    with open(hmm_file_path, 'r') as f:
        for line in f:
            if line.startswith('    ------- ------ -----    '):
                table_started=True
                continue
            elif line.startswith('\n'): 
                table_started=False
            elif table_started:
                hit_id_list.append(line[60:72])
                hit_score_list.append(float(line[13:18]))
                hit_Eval_list.append(float(line[4:11]))

    df = pd.DataFrame({'hit_id': hit_id_list, 'Eval': hit_Eval_list, 'score': hit_score_list})
    return df

def integrate_hmm_scores(counts_of_hits_df, hmm_search_hits_df):
    # Initialize a single-row DataFrame with NaNs
    new_row = pd.DataFrame(index=[0], columns=counts_of_hits_df.columns)
    new_row.iloc[0] = np.nan  # Fill the row with NaN
    
    # Loop through each hit in the hmm_search_hits_df DataFrame
    for index, row in hmm_search_hits_df.iterrows():
        hit_id = row['hit_id'].strip()  # Ensure there's no trailing whitespace
        hit_score = row['score']
        
        # If the hit_id matches a column in counts_of_hits_df, update the new_row with the score
        if hit_id in new_row.columns:
            new_row.at[0, hit_id] = hit_score
        else:
            # If the column does not exist, add it to the DataFrame
            new_row[hit_id] = pd.Series(hit_score, index=[0])
    
    # Concatenate this new_row to counts_of_hits_df
    updated_counts_of_hits_df = pd.concat([counts_of_hits_df, new_row], ignore_index=True)
    
    return updated_counts_of_hits_df


def integrate_hmm_scores_to_csv(counts_of_hits_csv_path, hmm_search_hits_df):
    if os.path.exists(counts_of_hits_csv_path):
        counts_of_hits_df = pd.read_csv(counts_of_hits_csv_path)
    else:
        unique_hits = hmm_search_hits_df['hit_id'].unique()
        counts_of_hits_df = pd.DataFrame(columns=unique_hits)
    
    # Use DataFrame to handle the merge and avoid explicit loops
    new_row = hmm_search_hits_df.set_index('hit_id')['score'].to_frame().T
    counts_of_hits_df = pd.concat([counts_of_hits_df, new_row], ignore_index=True, sort=False).fillna(0)
    
    counts_of_hits_df.to_csv(counts_of_hits_csv_path, index=False)

def add_hmm_search_results_to_csv(csv_file_path, precision, recall):
    # Read the CSV file as if it were a TSV
    df = pd.read_csv(csv_file_path, sep='\t')

    # The new data to add (as a dictionary where keys are column names and values are the data)
    new_data = {'Precision': [precision], 'Recall': [recall]}

    # Create a DataFrame from the new data
    new_df = pd.DataFrame(new_data)

    # Append the new data to the DataFrame using pd.concat
    df = pd.concat([df, new_df], ignore_index=True)

    # Save the updated DataFrame back to a CSV file, maintaining the TSV format
    df.to_csv(csv_file_path, sep='\t', index=False)

    
def save_hmm_queries_to_csv(query_csv_path, msa_fasta):
    query_ids = get_sequence_ids(msa_fasta)
    new_row_df = pd.DataFrame([query_ids], columns=[f'Query{num}' for num in range(1,11)])

    # Check if the file exists
    if os.path.exists(query_csv_path):
        # Append without reading the entire file
        new_row_df.to_csv(query_csv_path, mode='a', header=False, index=False)
    else:
        # Create new file with header
        new_row_df.to_csv(query_csv_path, index=False)
################## ------------------------------------------------------------------------------------------------- ##################
################################################        Objects       ###########################################################
################## ------------------------------------------------------------------------------------------------- ##################


eukarya5_index = get_eukarya5_index()
eukarya_species = pd.read_csv('/home/nicola/internship/winghel/eukarya_species.tsv', sep='\t', usecols=[1,2,4,5])
eukarya_metadata = pd.read_csv('/home/nicola/internship/winghel/eukarya_proteomes_metadata.tsv', sep='\t', usecols=[0,1,2,8,9,10])

