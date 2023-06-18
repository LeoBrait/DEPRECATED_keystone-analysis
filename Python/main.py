import pandas as pd
import subprocess
import os
import sys

src_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(src_dir, '..', 'data/')
global_dir = os.path.join(src_dir, '..')

############################# Data Preprocessing ###############################

metadata = pd.read_csv(f'{data_dir}metadata/biome_classification.csv').filter(
        items=['samples', 'habitat', 'ecosystem'])
annotation = pd.read_csv(f'{data_dir}taxon_abundances/kraken_custom_phyla.csv')

merged_data = pd.merge(annotation, metadata, on='samples', how='inner')

ecosystems = merged_data['ecosystem'].unique()

for ecosystem in ecosystems:
    subset = merged_data[merged_data['ecosystem'] == ecosystem]
    habitats = subset['habitat'].unique()
    
    for habitat in habitats:
        sub_subset = subset[subset['habitat'] == habitat]
        filename = f"{data_dir}community_subsets/{ecosystem}.{habitat}.tsv"
        if sub_subset.shape[0] < 12:
            print(f"Skipping {ecosystem}.{habitat} because it has"
                  "fewer than 12 samples.")
        else:
            #drop columns with sum equal to 0
            column_sums = sub_subset.sum()
            zero_sum_columns = column_sums[column_sums == 0].index
            sub_subset = sub_subset.drop(zero_sum_columns, axis=1)

            sub_subset = sub_subset.drop(['habitat', 'ecosystem'], axis=1)
            sub_subset.to_csv(filename, sep='\t', index=False)

############################# Fastspar #########################################

# # fastspar_dir = f'{data_dir}/fastspar_networks/'
# community_subsets = os.listdir(f'{data_dir}community_subsets/*.csv')


# for subset in community_subsets:
#     net_output = f'{fastspar_dir}/{subset}/sparcc_data/'
#     os.mkdir(net_output)

#     community_matrix = pd.read_csv(f'{data_dir}community_subsets/{subset}')
#     community_matrix.index.name = "#OTU ID"
#     community_matrix.columns.name = "#OTU ID"
#     community_matrix.T.to_csv(f'{net_output}/{subset}.tsv', sep='\t')

#     #transposing the data
#     x = pd.read_csv(f, index_col=0)

#     # fastspar implementation requirements (BIOME tsv format)
#     x.index.name = "#OTU ID"
#     x.columns.name = "#OTU ID"

#     x.T.to_csv(o,sep='\t')


#     community_matrix = pd.read_csv(f'{data_dir}community_subsets/{subset}')
#     community_matrix = community_matrix.set_index('samples').T
    

#     fastspar_command = (
        # f""
        "fastspar "
        "-c ${subset} "
        "-r ${net_output}/cor.tsv "
        "-a ${net_output}/cov.tsv "
        "-t 4 "
        "-s 1 "
        "-i 3000 "
        "-x 50 "
        "-e 0.1 "
        # "-y")

# bash_command = '''''
#!/bin/bashls
# [[ -z $1 ]] && echo 'Missing input file!' && exit 1

# nthreads=4 # WARNING: It does contribute to randomness!
# # instead of a single seed, I'm using a random variable for all used fastspar steps

# iter=3000
# xter=50
# corThr=0.1
# seed=1


# input=${1##*/}


# base=output/transposed/$input/sparcc/sparcc_data
# anot=anot.tsv

# # creating final folder
# rm -rf $base && mkdir -p $base

# # creating input file for sparcc
# ls
# echo -e 'Converting annotation file to fastspar input format\n'
# python3 csvtool_pandas.py community_matrix/${input} ${base}/${anot}

# # computing the original data
# fastspar -c ${base}/$anot -r ${base}/cor.tsv -a ${base}/cov.tsv -t $nthreads -s $seed -i $iter -x $xter -e $corThr -y
# '''