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
