# PREREQUISITES: access to tsv file for desired phenotype which has been pre-processed by metareg_prep.py
# IN: gene name, phenotype name, path to phenotype pre-processed tsv file
# OUT: saves the results of a genotype-phenotype pair to a csv file with name <geno>_<pheno>.tsv

# note that 'input_file' param actually determines the phenotype being analyzed for
# and the 'phenotype' parameter simply sets the phenotype name for the output path

import pandas as pd
import numpy as np
import sys
import os

def main(gene, phenotype, input_file):
    # Read the TSV file
    df = pd.read_csv(input_file, sep='\t')

    # Filter the data for just the selected gene
    filtered_df = df[df['gene'] == gene]
    
    if filtered_df.empty:
        print(f"No data found for gene: {gene} and phenotype: {phenotype}")

    # Select required columns and categories
    output_df = filtered_df[['gene', 'pos', 'markerID', 'BETA', 'SE', 'prob_0', 'am_pathogenicity', 'pLoF_indicator', 'missense_indicator']]
    
    # Ensure probabilities are not exactly 1 or 0
    epsilon = 1e-2
    output_df['prob_0'] = np.clip(output_df['prob_0'], epsilon, 1 - epsilon)
    output_df['am_pathogenicity'] = np.clip(output_df['am_pathogenicity'], epsilon, 1 - epsilon)

    # Apply log transformations
    output_df[['log_constraint', 'log_pathogenicity']] = -np.log1p(-(output_df[['prob_0', 'am_pathogenicity']]))
    
    # Impute missing values with column means
    output_df[['log_constraint', 'log_pathogenicity']] = output_df[['log_constraint', 'log_pathogenicity']].fillna(output_df[['log_constraint', 'log_pathogenicity']].mean())
    output_df[['pLoF_indicator', 'missense_indicator']] = output_df[['pLoF_indicator', 'missense_indicator']].fillna(0)

    # Remove rows with missing effect_size or var_effect_size
    output_df = output_df.dropna(subset=['BETA', 'SE'])
    output_df = output_df[output_df['SE'] != 0]

    # Export to CSV
    directory_path = "/scratch/groups/mrivas/larissal/unified_model/geno_pheno/"
    new_filename = f"{gene}_{phenotype}.tsv"
    output_path = os.path.join(directory_path, new_filename)

    output_df.to_csv(output_path, sep='\t', index=False)
    print(f"Data successfully exported to {output_path}")
    print(output_df.head())

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <gene> <phenotype> <input_file>")
        sys.exit(1)
    
    gene = sys.argv[1]
    phenotype = sys.argv[2]
    input_file = sys.argv[3]
    
    main(gene, phenotype, input_file)
