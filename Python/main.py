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

src_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(src_dir, '..', 'data/')
global_dir = os.path.join(src_dir, '..')

sys.path.append(f'{src_dir}/src')
from xonsh_py import Logger, lsgrep, existOldData, mkdir_p, cpr, rmr, sexec
sys.stderr = sys.stdout = Logger() # check src/xonsh_py.py for details
from datetime import datetime
from correlation import main as correlation_main
from identify_keystones import identify_keystones

startTime = datetime.now()

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

            #prepare data for sparcc
            sub_subset = sub_subset.drop(['habitat', 'ecosystem'], axis=1)
            sub_subset = sub_subset.transpose().reset_index()
            sub_subset.columns = sub_subset.iloc[0]
            sub_subset = sub_subset.drop([0]).reset_index(drop=True).rename(
                columns={sub_subset.columns[0]: "#OTU ID"})

            sub_subset.to_csv(filename, sep='\t', index=False)

############################# Fastspar #########################################

networks_dir = f'{data_dir}/fastspar_networks/'
community_subsets = glob.glob(f'{data_dir}community_subsets/*.tsv')


for subset_path in community_subsets:
    
    #get subset name
    subset_name = (
        subset_path.split('/')[-1].split('.')[0] + 
        '.' + 
        subset_path.split('/')[-1].split('.')[1])

    os.makedirs(f'{networks_dir}{subset_name}', exist_ok=True)
    net_output = f'{networks_dir}/{subset_name}/sparcc_data/'
    os.makedirs(net_output, exist_ok=True)


    fastspar_command = (
        f"fastspar "
        "-c {subset_path} "
        "-r {net_output}/cor.tsv "
        "-a {net_output}/cov.tsv "
        "-t 4 "
        "-s 1 "
        "-i 30 " #TODO: change this to 3000
        "-x 50 "
        "-e 0.1 "
        "-y > "
        "{net_output}/log.txt")
    fastspar_command = fastspar_command.replace("{subset_path}", subset_path)
    fastspar_command = fastspar_command.replace("{net_output}", net_output)
    result = subprocess.run(fastspar_command, shell=True, capture_output=True)
    print(result.stdout)

########################## Keystone Identification #############################
for f in community_subsets:
    
    #TODO: correct this for coherence with the rest of the code
    fname = (
        f.split('/')[-1].split('.')[0] +
        '.' +
        f.split('/')[-1].split('.')[1])

    

    # checking if the actual file has metadata
    # TODO: solve this artifact
    meta = lsgrep('metadata',[fname])
    meta = meta[0] if meta else "none"

    print('Community: %s\nMetadata: %s\n' % (f, meta))
    print("fname: ", fname)
    # checking if there is a computed sparcc matrix
    if existOldData('output/transposed/'+fname+'/sparcc/sparcc_data/cor.tsv'):
        print('There is an already computed SparCC matrix for %s. Using it for the analysis.\n'%fname)
        spcc_corr_mat = 'output/transposed/'+fname+'/sparcc/sparcc_data/cor.tsv'
    else:
        print(f'NO SPARCC MATRIX FOUND! JUMPING {fname}')
        continue

    # Create output directory
    mkdir_p(['output/transposed/'+fname])

    # run the analysis
    # f : community matrix
    # meta : metadata
    # fname : environment name
    # spcc_corr_mat : sparcc matrix
    coSparCC = correlation_main(f, meta, fname, spcc_corr_mat)

    # running code-integration steps
    # this step is responsible for running the CNM and LIASP algorithms
    # and identifying the most important nodes in the network (the keystones)
    # identify_keystones("output/transposed/"+fname+"/sparcc", coSparCC, fname)

    # copy the results to the corresponding environment directory
    cpr(lsgrep('out',['']),'output/transposed/'+fname+'/sparcc/')
    # delete the temporary 'out' directory
    rmr(['out'])

# record end time   
print('\nTotal execution time: ', (datetime.now()-startTime))

# running all analysis script
startTime = datetime.now() # record start time
#run the script that will generate the results run_keystone_analysis.py


print("if error, try to Manually exclude karst from the analysis ")
sexec('./run_keystone_analysis.py') # Keystones analysis
#TODO: Maybe this should be replaced with -> exec(open('run_keystone_analysis.py').read())
# This should solve the problem for Windows users
#TODO: Manually excluded karst from the analysis to perform withou errors

print('\nPost-results execution time: ', (datetime.now()-startTime)) # record end time
