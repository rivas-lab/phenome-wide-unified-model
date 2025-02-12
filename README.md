# Phenome-wide unified model

## Table of Contents
1. [Download genebass files](https://github.com/rivas-lab/phenome-wide-unified-model/blob/main/README.md#1-download-genebass-files)
2. [Process phenotype descriptions](https://github.com/rivas-lab/phenome-wide-unified-model/blob/main/README.md#2-process-phenotype-descriptions)
3. [Prepare Files for Meta-Regression](https://github.com/rivas-lab/phenome-wide-unified-model/blob/main/README.md#3-prepare-files-for-meta-regression) (```metareg_prep.py```)
4. [Identify Continuous Phenotypes](https://github.com/rivas-lab/phenome-wide-unified-model/blob/main/README.md#4-identify-continuous-phenotypes)
5. [Job Scripting for ```metareg_prep.py```]()
6. Unified Meta-Regression Model


# 1. Download genebass files

## Prerequisites: 
- Install hail and pandas for python3
- Allow hail to read from Google Cloud Storage. Hail provides documentation for doing this [here](https://hail.is/docs/0.2/cloud/google_cloud.html).

```python3
# reads genebass files from GCS and converts them to TSV files.

import hail as hl
import pandas as pd

hl.init()

gene_burden_results = hl.read_matrix_table('gs://ukbb-exome-public/500k/results/results.mt')
single_variant_results = hl.read_matrix_table('gs://ukbb-exome-public/500k/results/variant_results.mt')
phenotype_metadata = hl.read_table('gs://ukbb-exome-public/500k/results/pheno_results.ht')

scratch_path = '/scratch/groups/mrivas/larissal' # replace with path to where you want to save these files

gene_burden_table = gene_burden_results.rows()
gene_burden_table.export(scratch_path + 'gene_burden_results.tsv')

single_variant_results_table = single_variant_results.rows()
single_variant_results_table.export(scratch_path + 'single_variant_results.tsv')

phenotype_metadata_df = phenotype_metadata.to_pandas()
phenotype_metadata_df.to_csv(scratch_path + 'pheno_results.tsv', sep='\t', index=False)
```

# 2. Process phenotype descriptions

## Prerequisites:
- Install hail.
- Have access to the phenotype file as a TSV.

```python3
# process_phenotypes.py


import sys
import hail as hl
import os

def process_phenotypes(start_index, end_index):
    hl.init()

    # Read the table (TSV file)
    path = '/scratch/groups/mrivas/larissal/pheno_results.tsv'
    table = hl.import_table(path, impute=True)
    
    # Collect unique phenotype codes
    pheno_codes_set = table.aggregate(hl.agg.collect_as_set(table.phenocode))
    
    # Convert the set to a list and sort it
    pheno_codes_list = sorted(list(pheno_codes_set))

    # Process the subset of phenotypes using the list
    for pheno_code in pheno_codes_list[start_index:end_index]:
        file_path = f"/scratch/groups/mrivas/larissal/genebassout/{pheno_code}.genebass.tsv.gz"
        if os.path.exists(file_path):
            print(f"{pheno_code} EXISTS, SKIPPING")
            continue
        # Filter the table based on phenotype code
        subset_table = table.filter(table.phenocode == pheno_code)
        # Export the filtered table
        subset_table.export(file_path, header=True)

if __name__ == "__main__":
    start_index = int(sys.argv[1])
    end_index = int(sys.argv[2])
    process_phenotypes(start_index, end_index)
```

# 3. Prepare Files for Meta-Regression (```metareg_prep.py```)
OUTPUT: a TSV file ready to be analyzed by unified_reg.py
## Prerequisites:
Access to the following datasets
1. Constraint probabilities - data located at ``` /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm/results/HMM_rgc_0.9_over20_chr2_predictions_rgc_wes.tsv.gz ```
2. Alpha Missense scores located in:
``` /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/data/AlphaMissense_hg38.tsv.gz ```
3. Downloaded genebass files for phenotypes

```python3
# metareg_prep.py


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
```

# 4. Identify Continuous Phenotypes
Helpful when we want to analyze only continuous phenotypes.

## Prerequisities:
1. Downloaded genebass file ```pheno_results.tsv```

```python3
# get_continuous_phenos.py
# creates a text file containing every continuous phenotype code

import pandas as pd

def get_continuous_phenos(path):
    df = pd.read_csv(path, sep='\t')

    # Filter rows where 'trait_type' is 'continuous'
    continuous_df = df[df['trait_type'] == 'continuous']
    continuous_phenocodes = set(continuous_df['phenocode'])
    
    # Export the set of continuous phenotypes to a text file
    with open('/scratch/groups/mrivas/larissaredo/continuous_phenos.txt', 'w') as f:
        for pheno in continuous_phenocodes:
            f.write(f"{pheno}\n")

if __name__ == "__main__":
    # edit path to match location of phenoresults.tsv
    my_path = "/scratch/groups/mrivas/larissaredo/pheno_results.tsv"
    continuous_phenos = get_continuous_phenos(my_path)
```

# 5. Job Scripting for ```metareg_prep.py```
Push the file pre-processing step across all continuous phenotypes.

## Prerequisities:
1. Working version of ```metareg_prep.py```
2. A file ```continous_phenos.txt``` containing a list of every continuous phenotype (use ```get_continuous_phenos.py```)
3. A job script (here it is named ```metareg_prep.sh```) for submitting each phenotype to be processed.

```bash
#!/bin/bash

# Load the necessary module
module load python/3.9

# Directory containing the input files
DIRECTORY="/scratch/groups/mrivas/genebassout"

# Directory to save the output files
OUTPUT_DIR="/scratch/groups/mrivas/larissaredo/prepped_files"

# Path to the file containing continuous phenotypes
CONTINUOUS_PHENOS_FILE="/scratch/groups/mrivas/larissaredo/continuous_phenos.txt"

# Ensure the continuous phenotype file exists
if [ ! -f "$CONTINUOUS_PHENOS_FILE" ]; then
  echo "Error: Continuous phenotypes file not found at $CONTINUOUS_PHENOS_FILE."
  exit 1
fi

# Loop over each file in the directory with .tsv extension
for file in "$DIRECTORY"/*.tsv.gz; do
  # Extract the phenotype code from the file name
  pheno=$(basename "$file" .tsv)
  # removing the suffix .genebass because i want files to match phenocode name
  pheno=$(echo "$pheno" | cut -d'.' -f1)
  
  # Check if the phenotype code is in the file of continuous phenotypes
  if grep -Fxq "$pheno" "$CONTINUOUS_PHENOS_FILE"; then
    output_file="$OUTPUT_DIR/$(basename "$file")"
    sbatch metareg_prep.sh "$file" "$output_file"
    echo "Phenotype $pheno is continuous. Submitting job..."
  else
    echo "Skipping: $(basename "$file"). Reason: Not a continuous phenotype"
  fi 
done
```


# II. Unified Meta Regression Model
