# unified_reg_MAF.05.py

# PREREQUISITES: you need to pass in the filepath to a tsv file which has been pre-processed by metareg_prep.py
# OUT: saves the results of a meta regression analysis (p values and coefficients) to a compressed TSV file
# be sure to pass in a .tsv.gz file for the output file, and not just a .tsv

import pandas as pd
import numpy as np
import os
from tqdm import tqdm
from statsmodels.regression.linear_model import WLS
from statsmodels.tools.tools import add_constant
import logging


def main(results_path, output_path):
    # Read the data from the file
    input_df = pd.read_csv(results_path, sep='\t')

    # Initialize lists to store results
    meta_model_results = []

    # frequency threshold
    input_df = input_df[input_df['AF'] <= 0.05]

    # Ensure probabilities are not exactly 1 or 0
    epsilon = 1e-2
    input_df['prob_0'] = np.clip(input_df['prob_0'], epsilon, 1 - epsilon)
    input_df['am_pathogenicity'] = np.clip(input_df['am_pathogenicity'], epsilon, 1 - epsilon)

    # Apply log transformations
    input_df[['log_constraint', 'log_pathogenicity']] = -np.log1p(-(input_df[['prob_0', 'am_pathogenicity']]))

    # Impute missing values with column means
    input_df[['log_constraint', 'log_pathogenicity']] = input_df[['log_constraint', 'log_pathogenicity']].fillna(input_df[['log_constraint', 'log_pathogenicity']].mean())
    input_df[['pLoF_indicator', 'missense_indicator']] = input_df[['pLoF_indicator', 'missense_indicator']].fillna(0)

    # Remove rows with missing effect_size or var_effect_size
    input_df = input_df.dropna(subset=['BETA', 'SE'])
    input_df = input_df[input_df['SE'] != 0]

    # Group data by gene
    grouped_gene_data = input_df.groupby(['gene'])

    # Loop over each gene group and build a meta-regression model
    for gene, gene_data in tqdm(grouped_gene_data, desc="Processing genes", unit="gene"):
        
        if gene_data[['log_constraint', 'log_pathogenicity', 'pLoF_indicator', 'missense_indicator', 'BETA', 'SE']].isnull().any().any():
            continue

        # Meta-regression model for the gene
        X = add_constant(gene_data[['log_constraint', 'log_pathogenicity', 'pLoF_indicator', 'missense_indicator']])
        y = gene_data['BETA']
        weights = 1 / gene_data['SE']**2

        try:
            model = WLS(y, X, weights=weights, missing='drop').fit()
            
            # Set up logging basic config
            pheno_code = os.path.basename(results_path).replace(".genebass.tsv.gz", "")
            log_dir = "/scratch/groups/mrivas/larissaredo/unified_model"
            log_file = os.path.join(log_dir, f"{pheno_code}.log")
            logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

            # Perform diagnostics
            logging.info(f"Gene {gene[0]}")
            logging.info(model.summary().as_text())

            # Append results to the list: coefficients and p-values
            meta_model_results.append({
                'gene': gene[0],
                # model p-value
                'p_model': model.f_pvalue,
                # log constraint
                'coef_log_constraint': model.params['log_constraint'],
                'p_log_constraint': model.pvalues['log_constraint'],
                # log pathogenicity
                'coef_log_pathogenicity': model.params['log_pathogenicity'],
                'p_log_pathogenicity': model.pvalues['log_pathogenicity'],
                # pLoF
                'coef_pLoF_indicator': model.params['pLoF_indicator'],
                'p_pLoF_indicator': model.pvalues['pLoF_indicator'],
                # missense
                'coef_missense_indicator': model.params['missense_indicator'],
                'p_missense_indicator': model.pvalues['missense_indicator'],
                # constant
                'coef_constant': model.params['const'],
                'p_constant': model.pvalues['const']
            })
        except Exception as e:
            logging.error(f"An error occurred with gene {gene}: {e}")

    # Create a DataFrame from the results
    unified_model_df = pd.DataFrame(meta_model_results)

    # Save the results to a compressed CSV file
    unified_model_df.to_csv(output_path, index=False, compression='gzip', sep='\t')


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python unified_reg.py <results_path> <output_path>")
        sys.exit(1)
    
    results_path = sys.argv[1]
    output_path = sys.argv[2]
    main(results_path, output_path)
