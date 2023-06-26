import pandas as pd
import os
from datetime import datetime
from numpy.linalg import norm
import numpy as np

# Paths
this_dir = os.path.dirname(os.path.abspath(__file__))
global_dir = os.path.join(this_dir, '../..')
data_dir = os.path.join(global_dir, 'data/')


pre_process_time = datetime.now()



################################################################
# Read the data
kraken_custom_phyla = pd.read_csv(f'{data_dir}taxon_abundances/kraken_custom_phyla.csv')

# Define the desired sum value
desired_sum = 10000000000000000000  # 1 billion

# Calculate the row sums of the data
data_rowsum = np.sum(kraken_custom_phyla.iloc[:, 1:], axis=1)

# Calculate the scaling factor to achieve the desired sum
scaling_factor = desired_sum / np.sum(data_rowsum)

# Multiply each relative abundance by the scaling factor
kraken_custom_phyla.iloc[:, 1:] *= scaling_factor

# Calculate the row sums of the transformed data
data_rowsum_transformed = np.sum(kraken_custom_phyla.iloc[:, 1:], axis=1)

# Normalize the transformed data by dividing each element by its corresponding row sum
kraken_custom_phyla.iloc[:, 1:] = kraken_custom_phyla.iloc[:, 1:].div(data_rowsum_transformed, axis=0)

# Verify the sum of each column
column_sums = np.sum(kraken_custom_phyla.iloc[:, 1:], axis=0)
print(column_sums)







merged_data = pd.read_csv(
    f'{data_dir}metadata/biome_classification.csv').filter(
        items=['samples', 'habitat', 'ecosystem']).merge(
            kraken_custom_phyla,
            on='samples',
            how='inner')

#summarize the number of rows per habitat
habitat_counts = merged_data.groupby(
    ['ecosystem','habitat']).count()['samples']
os.makedirs(f'{data_dir}summaries/', exist_ok=True)
habitat_counts.to_csv(f'{data_dir}summaries/habitat_counts.csv')

# Subset data by ecosystem and habitat
ecosystems = merged_data['ecosystem'].unique()

for ecosystem in ecosystems:
    subset = merged_data[merged_data['ecosystem'] == ecosystem]
    habitats = subset['habitat'].unique()

    for habitat in habitats:
        sub_subset = subset[subset['habitat'] == habitat]
        filename = f"{data_dir}community_subsets/{ecosystem}.{habitat}.tsv"
  
        #TODO: remove this when correlation function is fixed to accept tsv files
        filename_tsv = f"{data_dir}community_subsets_raw/{ecosystem}.{habitat}.csv"

        if sub_subset.shape[0] < 3:
            print(f"Skipping {ecosystem}.{habitat} because it has"
                  "fewer than 3 samples.")
        else:
            column_sums = sub_subset.sum()
            zero_sum_columns = column_sums[column_sums == 0].index
            sub_subset = sub_subset.drop(zero_sum_columns, axis=1)
            sub_subset = sub_subset.drop(['habitat', 'ecosystem'], axis=1)

            #TODO: remove this when correlation function is fixed to accept
            #transposed tsv files
            os.makedirs(f'{data_dir}community_subsets_raw/', exist_ok=True)
            sub_subset.to_csv(filename_tsv, sep=',', index=False)

            #prepare data for sparcc (transpose and remove metadata columns)
            os.makedirs(f'{data_dir}community_subsets/', exist_ok=True)
            sub_subset = sub_subset.transpose().reset_index()
            sub_subset.columns = sub_subset.iloc[0]
            sub_subset = sub_subset.drop([0]).reset_index(drop=True).rename(
            columns={sub_subset.columns[0]: "#OTU ID"})
            sub_subset.to_csv(filename, sep='\t', index=False)
        continue
print(f"Data preprocessing took {datetime.now() - pre_process_time}")