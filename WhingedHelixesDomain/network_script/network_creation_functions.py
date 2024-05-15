import os
import csv
import pandas as pd

## FUNCTIONS

def remove_brakets(range_str):
    range_= range_str.split('(')[0]
    return range_


def calculate_overlap(range1, range2):
    """
    Calculate the overlap between two ranges.
    
    Args:
    range1 (tuple): First range as a tuple (start, end).
    range2 (tuple): Second range as a tuple (start, end).
    
    Returns:
    int: The overlap between the two ranges, 0 if there is no overlap.
    """
    # Extract start and end points from ranges
    start1, end1 = range1
    start2, end2 = range2
    
    # Calculate the maximum of the start points and the minimum of the end points
    max_start = max(start1, start2)
    min_end = min(end1, end2)
    
    # Calculate overlap
    if min_end > max_start:
        # There is an overlap
        overlap = min_end - max_start
    else:
        # No overlap
        overlap = 0
    return overlap



# Extracts an edge from a line of hhsearch output, including specified edge features.
def get_edge_from_line(query_profile_name, line, WHD_range=False, edge_feature=""):
    # Define the positions of each column in the line
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

    # Extract attributes based on positions
    attributes = [line[start:end].strip() for start, end in positions]

    attributes[8]=str(remove_brakets(attributes[8])).strip()

    if WHD_range:
        if attributes[0] not in ['ask1', 'dam1', 'ska1', 'ska2', 'ska3', 'spc19', 'spc34']:
            attributes[0]=f'{attributes[0]}_{attributes[8]}'
        

    
    # Combine the query profile name with the extracted attributes, separated by tabs
    edge = query_profile_name + "\t" + "\t".join(attributes) + "\t" + edge_feature
    
    return edge

# Extracts edges from a file, selecting the best and second-best hits and any number of additional hits specified.
def get_edges_from_file_n(file_path, filename, hits_list, additional_hits=0):
    with open(file_path, 'r') as file:
        query_profile_name = filename.split('.')[0]  # Extract the name of the profile from filename
        for current_line_number, line in enumerate(file, 1):  # Iterate on lines
            if current_line_number == 10:  # Get the line with the best hit
                edge = get_edge_from_line(query_profile_name, line, edge_feature='BH')
                hits_list.append(edge)
            elif current_line_number == 11:  # Get the line with the second best hit
                edge = get_edge_from_line(query_profile_name, line, edge_feature='SBH')
                hits_list.append(edge)
            elif additional_hits > 0 and current_line_number > 11:
                # Only process additional hits if additional_hits > 0
                if line.strip() == "":
                    break  # Exit loop early if an empty line is encountered
                edge = get_edge_from_line(query_profile_name, line, edge_feature='H')
                hits_list.append(edge)
                # If we've reached the limit of additional hits to process, break out of the loop
                if current_line_number >= 11 + additional_hits:
                    break

    return hits_list

# Extracts edges from a file with specified E-value and alignment length cutoffs.
def get_edges_from_file_cut_off(file_path, filename, hits_list, E_cut_off=0.001, min_alig_len=10):
    with open(file_path, 'r') as file:
        query_profile_name = filename.split('.')[0]  # Extract the name of the profile from filename
        for current_line_number, line in enumerate(file, 1):  # Iterate on lines
            if current_line_number >= 10 and line.strip() == "":
                break  # Exit loop early if an empty line is encountered
            if current_line_number == 10:  # Get the line with the best hit
                edge = get_edge_from_line(query_profile_name, line, edge_feature='BH')
                # print('Line N: '+str(current_line_number)+'; Edge: '+edge)
                hits_list.append(edge)
            elif current_line_number == 11:  # Get the line with the second best hit
                edge = get_edge_from_line(query_profile_name, line, edge_feature='SBH')
                # print('Line N: '+str(current_line_number)+'; Edge: '+edge)
                hits_list.append(edge)
            elif current_line_number > 11:
                E_val=float(line[41:48].strip())
                alig_len=float(line[70:74].strip())
                if E_val<=E_cut_off and alig_len>=min_alig_len:
                    edge = get_edge_from_line(query_profile_name, line, edge_feature='H')
                    # print('Line N: '+str(current_line_number)+'; Edge: '+edge)
                    hits_list.append(edge)
                # If we've reached the limit of additional hits to process, break out of the loop
        
    return hits_list


def get_edges_from_file_cut_off_and_range(file_path, filename, hits_list, E_cut_off=35, min_alig_len=10, Prob_cut_off = 20):
    with open(file_path, 'r') as file:
        query_profile_name = filename.split('.')[0]  # Extract the name of the profile from filename
        for current_line_number, line in enumerate(file, 1):  # Iterate on lines
            if current_line_number >= 10 and line.strip() == "":
                break  # Exit loop early if an empty line is encountered
            if current_line_number == 10:  # Get the line with the best hit
                edge = get_edge_from_line(query_profile_name, line, WHD_range=True, edge_feature='BH')
                # print('Line N: '+str(current_line_number)+'; Edge: '+edge)
                hits_list.append(edge)
            elif current_line_number == 11:  # Get the line with the second best hit
                edge = get_edge_from_line(query_profile_name, line, WHD_range=True, edge_feature='SBH')
                # print('Line N: '+str(current_line_number)+'; Edge: '+edge)
                hits_list.append(edge)
            elif current_line_number > 11:
                E_val=float(line[41:48].strip())
                Prob=float(line[35:40].strip())
                alig_len=float(line[70:74].strip())
                if E_val<=E_cut_off and alig_len>=min_alig_len and Prob>=Prob_cut_off:
                    edge = get_edge_from_line(query_profile_name, line, WHD_range=True, edge_feature='H')
                    # print('Line N: '+str(current_line_number)+'; Edge: '+edge)
                    hits_list.append(edge)
                # If we've reached the limit of additional hits to process, break out of the loop
        
    return hits_list

# Filters edges to remove duplicates, keeping only the most significant edge (lowest E-value) for each query-hit pair.
def filter_edges(input_file_path, output_file_path):
    unique_edges = {}  # Dictionary to hold the best edge for each query-hit pair

    # Read the input file
    with open(input_file_path, 'r') as infile:
        reader = csv.DictReader(infile, delimiter='\t')
        for row in reader:
            query = row['Query']
            hit = row['Hit']
            e_value = float(row['E-value'].strip())
            # Create a unique key for each query-hit pair
            key = (query, hit)

            # Check if the query-hit pair is already in the dictionary and compare E-values correctly
            if key not in unique_edges or float(unique_edges[key]['E-value']) > e_value:
                # Store the row with E-value converted back to string for consistency
                row['E-value'] = str(e_value)
                unique_edges[key] = row  # Update with the new best edge

    # Write the filtered dataset to the output file
    with open(output_file_path, 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames, delimiter='\t')
        writer.writeheader()
        for edge in unique_edges.values():
            writer.writerow(edge)

# Function to add an 'IsQuery' column to the network, indicating whether the 'Query' is a query protein.
def add_IsQuery_column(path):
    network = pd.read_csv(path, sep='\t')
    network['IsQuery']= ~network['Query'].str.contains('COG') & ~network['Query'].str.contains('KOG')
    network.to_csv(path, sep='\t', index=False)

def add_IsThird_column(csv_path, query_list):
    # Load the CSV file
    network = pd.read_csv(csv_path, sep='\t')
    network['IsThird'] = network['Query'].apply(lambda x: x in query_list)
    network.to_csv(csv_path, sep='\t', index=False)


### WHD-Only functions

def has_significant_overlap(range1, range2):
    start1, end1 = map(int, range1.split('-'))
    start2, end2 = map(int, range2.split('-'))
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    if overlap_start <= overlap_end:  # Overlap exists
        overlap_length = overlap_end - overlap_start
        range1_length = end1 - start1
        if range1_length > 0:
            return overlap_length / range1_length > 0.7  # Adjust threshold as needed
    return False

def update_dataframe_ranges(df, protein_ranges):
    def update_column(value):
        if '_' not in value:
            return value  # No range present, return original
        
        # Extract protein name and range
        protein_name, range_str = value.split('_')
        if protein_name in protein_ranges:
            for dict_range in protein_ranges[protein_name]:
                if has_significant_overlap(range_str, dict_range):
                    # print(f'{protein_name}_{range_str} --> {protein_name}_{dict_range}')
                    return f"{protein_name}_{dict_range}"  # Return the dictionary range
            print(f'{protein_name}_{range_str} has no significant overlap found with {protein_ranges[protein_name]} indicate row to be removed')    
            return None  # No significant overlap found, indicate row to be removed
        else:
            # print(f'{protein_name} not in dictionary, indicate row to be removed')    
            return None  # Protein not in dictionary, indicate row to be removed

    # print('Hit:')
    df['Hit'] = df['Hit'].apply(update_column)
    
    # Drop rows where either Query or Hit have been set to None
    df.dropna(subset=['Hit'], inplace=True)
    return df