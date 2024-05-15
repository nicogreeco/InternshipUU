import pandas as pd
from Bio import SeqIO
import pandas as pd
import random
import subprocess
import os
import numpy as np


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
    # Check if the file exists and is not empty
    if os.path.exists(counts_of_hits_csv_path):
        counts_of_hits_df = pd.read_csv(counts_of_hits_csv_path)
    else:
        unique_hits = hmm_search_hits_df['hit_id'].unique()
        counts_of_hits_df = pd.DataFrame(columns=unique_hits)
#        counts_of_hits_df.loc[0] = pd.NA  # Add an initial row filled with NA values

    # This should not raise an error since counts_of_hits_df now definitely exists
    new_row = pd.DataFrame(0, index=[0], columns=counts_of_hits_df.columns)  # Initialize with 0 for simplicity
    for _, row in hmm_search_hits_df.iterrows():
        hit_id = row['hit_id'].strip()
        hit_score = row['score']
        if hit_id in new_row.columns:
            new_row[hit_id] = hit_score
        else:
            # If the hit_id does not exist in the DataFrame, add it as a new column
            new_row[hit_id] = hit_score

    # Concatenate this new row to counts_of_hits_df
    updated_counts_of_hits_df = pd.concat([counts_of_hits_df, new_row], ignore_index=True)

    # Save the updated DataFrame back to the CSV file, overwriting it
    updated_counts_of_hits_df.to_csv(counts_of_hits_csv_path, index=False)


# Main workflow function
def main():
    protein='ska3'
    queries_fasta = f"/home/nicola/internship/winghel/alignments/{protein}_all_seq_only_WHD.fasta"
    target_fasta = f"/home/nicola/internship/winghel/fastas/{protein}/{protein}_conc.fasta"
    csv_file_path = f"/home/nicola/internship/winghel/random_search_banchmark/hmm_search_result_{protein}_notquery.csv"
    
    # Step 1: Select 5 random sequences and create a temp FASTA
    selected_fasta = select_random_sequences(queries_fasta)

    # Step 2: Create MSA from selected sequences
    msa_fasta = create_msa(selected_fasta)

    # Step 3: Build HMM profile from MSA
    hmm_profile = build_hmm(msa_fasta)

    # Step 4: Run hmmsearch against target dataset
    hmm_search_results = run_hmmsearch(hmm_profile, target_fasta)

    # Step 5: Add scores to CSV
    hmm_hits_df = get_hmm_search_hits(hmm_search_results)
    integrate_hmm_scores_to_csv(csv_file_path, hmm_hits_df)

    # Cleanup temporary files
    os.remove("temp_selected.fasta")
    os.remove("temp_msa.fasta")
    os.remove("temp_profile.hmm")
    os.remove("temp_hmmsearch_results.txt")

if __name__ == "__main__":
    main()


