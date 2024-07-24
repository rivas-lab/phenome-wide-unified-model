# Phenome-wide unified model

# Download genebass files

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

# Process phenotype descriptions

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
