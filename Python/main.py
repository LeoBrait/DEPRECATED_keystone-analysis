#!/usr/bin/env python3
try:
    import pandas as pd
    import subprocess
    import os
    import glob
    import sys
    import numpy
    #import gephitools
    import scipy
    import matplotlib
    import mpl_toolkits
    import networkx
    import PIL
except ImportError as e:
    print(e)
    sys.exit()
from datetime import datetime

# Paths
src_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(src_dir, '..', 'data/')
global_dir = os.path.join(src_dir, '..')

# Custom functions
sys.path.append(f'{src_dir}/src')
from xonsh_py import Logger, lsgrep, existOldData, mkdir_p, cpr, rmr, sexec
from correlation import main as correlation_main
from identify_keystones import identify_keystones
sys.stderr = sys.stdout = Logger()

############################# Data Preprocessing ###############################

metadata = pd.read_csv(f'{data_dir}metadata/biome_classification.csv').filter(
        items=['samples', 'habitat', 'ecosystem'])
annotation = pd.read_csv(f'{data_dir}taxon_abundances/kraken_custom_phyla.csv')

merged_data = pd.merge(annotation, metadata, on='samples', how='inner')

# Subset data by ecosystem and habitat
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
            #drop columns with sum equal to 0 (irrelevant taxa)
            column_sums = sub_subset.sum()
            zero_sum_columns = column_sums[column_sums == 0].index
            sub_subset = sub_subset.drop(zero_sum_columns, axis=1)

            #prepare data for sparcc (transpose and remove metadata columns)
            sub_subset = sub_subset.drop(['habitat', 'ecosystem'], axis=1)
            sub_subset = sub_subset.transpose().reset_index()
            sub_subset.columns = sub_subset.iloc[0]
            sub_subset = sub_subset.drop([0]).reset_index(drop=True).rename(
                columns={sub_subset.columns[0]: "#OTU ID"})

            sub_subset.to_csv(filename, sep='\t', index=False)

#delete karst porous table
#TODO: remove this when karst porous is fixed
if os.path.exists(f'{data_dir}community_subsets/groundwater.karst-porous.tsv'):
    os.remove(f'{data_dir}community_subsets/groundwater.karst-porous.tsv')

############################# Fastspar #########################################

startTime = datetime.now()
#create output generel directorie
networks_dir = f'{data_dir}/fastspar_networks/'
os.makedirs(f'{networks_dir}', exist_ok=True)

#list jobs
community_subsets = glob.glob(f'{data_dir}community_subsets/*.tsv')


for subset_path in community_subsets:
    
    subset_name = (
        subset_path.split('/')[-1].split('.')[0] + 
        '.' + 
        subset_path.split('/')[-1].split('.')[1])
    
    #prepare specific output directories
    os.makedirs(f'{networks_dir}{subset_name}', exist_ok=True)
    os.makedirs(f'{networks_dir}{subset_name}/sparcc', exist_ok=True)
    net_output = f'{networks_dir}/{subset_name}/sparcc/sparcc_data/'
    os.makedirs(net_output, exist_ok=True)

    #create sparcc shell command
    fastspar_command = (
        f"fastspar "
        "-c {subset_path} "
        "-r {net_output}/cor.tsv "
        "-a {net_output}/cov.tsv "
        "-t 4 "
        "-s 1 "
        "-i 3000 "
        "-x 50 "
        "-e 0.1 "
        "-y > "
        "{net_output}/log.txt")
    fastspar_command = fastspar_command.replace("{subset_path}", subset_path)
    fastspar_command = fastspar_command.replace("{net_output}", net_output)

    #run fastspar
    if existOldData(f'{net_output}/cor.tsv'):
        print('There is an already computed SparCC matrix for %s' %subset_name)
    else:
        subprocess.run(fastspar_command, shell=True, capture_output=True)
        print(f'fastspar of {subset_name} finished.')
        print(
            'Progress: %s of %s' 
            % (
                community_subsets.index(subset_path) + 1,
                len(community_subsets)))
        print('\nElapsed time: ', (datetime.now()-startTime))
 

################ Preprocessing for Keystones Identification ####################

for subset_path in community_subsets:
    
    #TODO: correct this for coherence with the rest of the code
    subset_name = (
        subset_path.split('/')[-1].split('.')[0] +
        '.' +
        subset_path.split('/')[-1].split('.')[1])

    meta = "none"
    
    # checking if there is a computed sparcc matrix
    print("running keystones preprocessing for:", subset_name)
    if existOldData(f'{data_dir}/fastspar_networks/'+subset_name+'/sparcc/sparcc_data/cor.tsv'):
        print('There is an already computed SparCC matrix for %s. Using it for the analysis.\n'%subset_name)
        spcc_corr_mat = f'{data_dir}/fastspar_networks/'+subset_name+'/sparcc/sparcc_data/cor.tsv'
    else:
        print(f'NO SPARCC MATRIX FOUND! JUMPING {subset_name}')
        continue

    # Create output directory
    mkdir_p([f'{data_dir}fastspar_networks/'+subset_name])


    coSparCC = correlation_main(subset_path, meta, subset_name, spcc_corr_mat)

    # running code-integration steps
    # this step is responsible for running the CNM and LIASP algorithms
    # and identifying the most important nodes in the network (the keystones)
    # identify_keystones("output/transposed/"+fname+"/sparcc", coSparCC, fname)

    # copy the results to the corresponding environment directory
    cpr(lsgrep('out',['']),f'{data_dir}fastspar_networks/'+subset_name+'/sparcc/')
    # delete the temporary 'out' directory
    rmr(['out'])

# record end time   
print('\nTotal execution time: ', (datetime.now()-startTime))



################## Keystones identification ####################################
startTime = datetime.now() # record start time

print("\n\nStarting analysis through environments.\n\n")
for level in ['phyla']:

    # create output directory
    # mkdir -p @('output/transposed_all_environments/'+level)
    mkdir_p(['output/transposed_all_environments/'])

    # run the identification of keystones
    print("Concat keystones for %s." % level)
    sexec('./src/concat_keystones.py '+level)

    # Generate the  keystones heatmap (total effect)
    # print("Heatmap keystones for %s." % level)
    # sexec('./src/heatmap_keystones.py '+level+' 1')

    # Generate the  keystones heatmap (indirect effect)
    # print("Heatmap keystones (indirect) for %s." % level)
    # sexec('./src/heatmap_keystones.py '+level+' 2')

print('\nPost-results execution time: ', (datetime.now()-startTime))