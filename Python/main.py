import pandas as pd
import os

src_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(src_dir, '..', 'data/')
global_dir = os.path.join(src_dir, '..')

############################# Data Preprocessing ###############################

metadata = pd.read_csv(
    f'{data_dir}metadata/biome_classification.csv').filter(
        items=['samples', 'habitat', 'ecosystem'])

annotation = pd.read_csv(f'{data_dir}taxon_abundances/kraken_custom_phyla.csv')

merged_data = pd.merge(annotation, metadata, on='samples', how='inner')

ecosystems = merged_data['ecosystem'].unique()

for ecosystem in ecosystems:
    subset = merged_data[merged_data['ecosystem'] == ecosystem]
    habitats = subset['habitat'].unique()
    
    for habitat in habitats:
        sub_subset = subset[subset['habitat'] == habitat]
        filename = f"{data_dir}community_subsets/{ecosystem}.{habitat}.csv"
        if sub_subset.shape[0] < 12:
            print(f"Skipping {ecosystem}.{habitat} because it has fewer than 12 samples.")
        else:
            sub_subset.to_csv(filename, index=False)