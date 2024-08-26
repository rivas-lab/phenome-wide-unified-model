# metareg_prep.py

# PREREQUISITES: data access to the following datasets
# 1. Constraint probabilities data located at: 
# /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm/results/HMM_rgc_0.9_over20_chr2_predictions_rgc_wes.tsv.gz
# 2. Alpha Missense scores located in:
# /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/data/AlphaMissense_hg38.tsv.gz
# OUT: a tsv file ready to be analyzed by unified_reg.py


import pandas as pd
import argparse
import json

# preprocess the TSV file:

# split 'locus' column into 'CHR' and 'POS' to match constraint probabilities columns
def reformat_locus(data):
    data[['chr', 'pos']] = data['locus'].str.split(':', expand=True)
    data['pos'] = data['pos'].astype(int)
    return data


# split 'alleles' column into 'REF' and 'ALT' to match constraint probabilities columns
def reformat_alleles(data):
    # extracts REF and ALT from the alleles column
    def extract_alleles(alleles_str):
        # Convert string representation of list to actual list
        alleles_list = json.loads(alleles_str)
        if len(alleles_list) != 2:
            raise ValueError(f"Unexpected alleles format: {alleles_str}")
        return pd.Series(alleles_list, index=['ref', 'alt'])
    
    # Apply the extraction function to the alleles column
    data[['ref', 'alt']] = data['alleles'].apply(extract_alleles)
    return data


# Merge to create new columns in TSV file

# adds the prob_0 column from constraint probabilities data (merge on CHR and POS (locus))
def add_constraint_prob_data(data, constraint_probabilities):
    constraint_prob_0 = constraint_probabilities[['chr', 'pos', 'prob_0']]
    return pd.merge(data, constraint_prob_0, on=['chr', 'pos'], how='left')


# adds the am_pathogenicity column from alpha missense dataset (merge on CHR, POS, REF, ALT)
def add_alpha_missense_data(data, alpha_missense):
    alpha_missense.rename(columns={'#CHROM': 'chr', 'POS': 'pos', 'REF': 'ref', 'ALT': 'alt'}, inplace=True)
    alpha_missense_filtered = alpha_missense[['chr', 'pos', 'ref', 'alt', 'am_pathogenicity']]
    return pd.merge(data, alpha_missense_filtered, on=['chr', 'pos', 'ref', 'alt'], how='left')


# adds a pLoF indicator with 'LC' and 'pLoF' mapping to 1, otherwise 0
def add_plof_ind_data(data):
    # Define the condition for setting the indicator to 1
    plof_indicators = ['pLoF', 'LC']
    condition = data['annotation'].isin(plof_indicators)

    # Create the new indicator column based on the condition
    data['pLoF_indicator'] = condition.astype(int)
    return data

# adds a missense indicator
def add_missense_ind_data(data):
    # Define the condition for setting the indicator to 1
    missense_indicators = ['missense']
    condition = data['annotation'].isin(missense_indicators)

    # Create the new indicator column based on the condition
    data['missense_indicator'] = condition.astype(int)
    return data

def main(dataset_path, output_path):
    # Define the file paths
    constraint_probabilities_path = '/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm/results/HMM_rgc_0.9_over20_chr2_predictions_rgc_wes.tsv.gz'
    alpha_missense_path = '/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/data/AlphaMissense_hg38.tsv.gz'
    
    # Load the datasets for constraint probabilities, alpha missense
    constraint_probabilities = pd.read_csv(constraint_probabilities_path, compression='gzip', sep='\t', comment='#')
    print("Constraint Probabilities Data:")
    print(constraint_probabilities.iloc[0:5])

    alpha_missense = pd.read_csv(alpha_missense_path, compression='gzip', sep='\t', skiprows=3)
    print("\nAlpha Missense Data:")
    print(alpha_missense.iloc[0:5])

    # we will call the astrazeneca file for this phenotype 'data'
    data = pd.read_csv(dataset_path, sep='\t', comment='#')
    print("\nOriginal Dataset:")
    print(data.iloc[0:5])

    # preprocess data
    data = reformat_locus(data)
    data = reformat_alleles(data)

    # add columns: prob_0, am_pathogenicity, pLoF_indicator
    data = add_constraint_prob_data(data, constraint_probabilities)
    data = add_alpha_missense_data(data, alpha_missense)
    data = add_plof_ind_data(data)
    data = add_missense_ind_data(data)

    # Write the final DataFrame to a TSV file
    data.to_csv(output_path, sep='\t', index=False)
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Join Constraint Probabilities Data with astrazeneca by-phenotype-data.")
    parser.add_argument('dataset_path', type=str, help='Path to the astrazeneca dataset file')
    parser.add_argument('output_path', type=str, help='Path to save the processed output file')
    
    args = parser.parse_args()
    main(args.dataset_path, args.output_path)
