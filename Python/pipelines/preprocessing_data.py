import pandas as pd
import os
from datetime import datetime
from numpy.linalg import norm
import sys
import json

#heritage from the main
analysis_frame=str(sys.argv[1])
multi_const = int(sys.argv[2])
minimum_samples = int(sys.argv[3])
annotated_table = str(sys.argv[4])
metadata_table = str(sys.argv[5])


# Paths
this_dir = os.path.dirname(os.path.abspath(__file__))
global_dir = os.path.join(this_dir, '../..')
data_dir = os.path.join(global_dir, 'data/')

# create dir
os.makedirs(f'{data_dir}/{analysis_frame}/'
              'community_subsets/', exist_ok=True)
os.makedirs(f'{data_dir}/{analysis_frame}/'
              'community_subsets_raw/', exist_ok=True)


pre_process_time = datetime.now()

#TODO: DECIDE THE FACTOR OF MULTIPLICATION TO AVOID FLOATING POINT ERRORS
#multiply each relative abundance for 1000
kraken_custom_phyla = pd.read_csv(
    f'{data_dir}/taxon_abundances/{annotated_table}')
kraken_custom_phyla.iloc[:,1:] = kraken_custom_phyla.iloc[:,1:].apply(
    lambda x: x * multi_const)

merged_data = pd.read_csv(
    f'{data_dir}/metadata/{metadata_table}').filter(
        items=['samples', 'habitat', 'ecosystem']).merge(
            kraken_custom_phyla,
            on='samples',
            how='inner')

#summarize the number of rows per habitat
habitat_counts = merged_data.groupby(
    ['ecosystem','habitat']).count()['samples']
os.makedirs(f'{data_dir}/{analysis_frame}/summaries/', exist_ok=True)
habitat_counts.to_csv(f'{data_dir}{analysis_frame}/'
                        f'summaries/habitat_counts.csv')

# Subset data by ecosystem and habitat
ecosystems = merged_data['ecosystem'].unique()

for ecosystem in ecosystems:
    subset = merged_data[merged_data['ecosystem'] == ecosystem]
    habitats = subset['habitat'].unique()

    for habitat in habitats:
        sub_subset = subset[subset['habitat'] == habitat]
        filename = (f"{data_dir}/{analysis_frame}/community_subsets/"
                        f"{ecosystem}.{habitat}.tsv")
  
        #TODO: remove this when correlation function accepts tsv files
        filename_tsv = (f"{data_dir}/{analysis_frame}/community_subsets_raw/"
                            f"{ecosystem}.{habitat}.csv")

        if sub_subset.shape[0] < minimum_samples:
            print(f"Skipping {ecosystem}.{habitat} because it has"
                  f"fewer than {minimum_samples} samples.")
        else:
            column_sums = sub_subset.sum()
            zero_sum_columns = column_sums[column_sums == 0].index
            sub_subset = sub_subset.drop(zero_sum_columns, axis=1)
            sub_subset = sub_subset.drop(['habitat', 'ecosystem'], axis=1)

            #remove undesired taxa
            if habitat == "estuarine_seawater":
              sub_subset = sub_subset.drop(['Myxococcota'], axis=1)
            if habitat == "salt_lake":
              sub_subset = sub_subset.drop(
                [    'Candidatus Azambacteria',             'Thermomicrobiota',
                   'Candidatus Liptonbacteria',      'Candidatus Tagabacteria',
                    'Candidatus Terrybacteria',    'Candidatus Spechtbacteria',
                   'Candidatus Microgenomates',      'Candidatus Wallbacteria',
                    'Candidatus Coatesbacteria', 'Candidatus Blackallbacteria',
                    'Candidatus Calescamantes',                'Chrysiogenota',
                    'Coprothermobacterota',        'Candidatus Diapherotrites',
                    'Candidatus Aenigmarchaeota', 'Candidatus Thermoplasmatota',
                    'Candidatus Thorarchaeota'], 
                axis=1)

            #TODO: remove this when correlation function is fixed to accept
            #transposed tsv files
            os.makedirs(f'{data_dir}/{analysis_frame}/community_subsets_raw/',
                        exist_ok=True)
            sub_subset.to_csv(filename_tsv, sep=',', index=False)

            #prepare data for sparcc (transpose and remove metadata columns)
            os.makedirs(f'{data_dir}/{analysis_frame}/community_subsets/',
                        exist_ok=True)
            sub_subset = sub_subset.transpose().reset_index()
            sub_subset.columns = sub_subset.iloc[0]
            sub_subset = sub_subset.drop([0]).reset_index(drop=True).rename(
            columns={sub_subset.columns[0]: "#OTU ID"})
            sub_subset.to_csv(filename, sep='\t', index=False)
        continue
print(f"Data preprocessing took {datetime.now() - pre_process_time}")
