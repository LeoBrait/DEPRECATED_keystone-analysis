import pandas as pd

metadata = pd.read_csv('data/metadata/biome_classification.csv').filter(
    items=['samples', 'habitat', 'ecosystem'])
annotation = pd.read_csv('data/taxon_abundances/kraken_custom_phyla.csv')

merged_data = pd.merge(annotation, metadata, on='samples', how='inner')

ecosystems = merged_data['ecosystem'].unique()

for ecosystem in ecosystems:
    subset = merged_data[merged_data['ecosystem'] == ecosystem]
    habitats = subset['habitat'].unique()
    for habitat in habitats:
        sub_subset = subset[subset['habitat'] == habitat]
        filename = f"data/community_subsets/{ecosystem}.{habitat}.csv"
        if sub_subset.shape[0] < 12:
            print(f"Skipping {filename} because it has fewer than 12 samples.")
        else:
            sub_subset.to_csv(filename, index=False)