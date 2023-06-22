#!/usr/bin/env python3
import sys
try:
    import pandas as pd
    import subprocess
    import os
    import glob
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
from numpy.linalg import norm

# Paths
src_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(src_dir, '../..', 'data/')
global_dir = os.path.join(src_dir, '../..')

# Custom functions
sys.path.append(f'{src_dir}/src')
from xonsh_py import Logger, lsgrep, existOldData, mkdir_p, cpr, rmr, sexec, cat
from correlation import main as correlation_main
from identify_keystones import identify_keystones
sys.stderr = sys.stdout = Logger()
total_time = datetime.now()

############################# Data Preprocessing ###############################

merged_data = pd.read_csv(
    f'{data_dir}metadata/biome_classification.csv').filter(
        items=['samples', 'habitat', 'ecosystem']).merge(
            pd.read_csv(f'{data_dir}taxon_abundances/kraken_custom_phyla.csv'),
            on='samples',
            how='inner')

#summarize the number of rows per habitat
habitat_counts = merged_data.groupby(['ecosystem','habitat']).count()['samples']
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

############################# Fastspar #########################################



#*************************** Bootstraping **************************************
startTime = datetime.now()

bootstrap_dir = f'{data_dir}/performance_iteractions_fastspar/'
os.makedirs(f'{bootstrap_dir}', exist_ok=True)

bootstrap_subsets = [
    f'{data_dir}community_subsets/'
        'animal_host-associated.aqueous_humour.tsv', #N=8
    f'{data_dir}community_subsets/'
        'animal_host-associated.animal_feces.tsv',   #N=675
    f'{data_dir}community_subsets/'
        'saline_water.coastal_seawater.tsv',         #N=286
    f'{data_dir}community_subsets/'
        'saline_water.hypersaline_water.tsv',        #N=16
    f'{data_dir}community_subsets/'
        'soil.savanna_soil.tsv',                     #N=21
    f'{data_dir}community_subsets/'
        'soil.tundra_soil.tsv',                      #N=3
    f'{data_dir}community_subsets/'
        'groundwater.porous_contaminated.tsv',       #N=48
    f'{data_dir}community_subsets/'
        'groundwater.mine.tsv']                      #N=3

processes = []
for subset_path in bootstrap_subsets:

    subset_name = (
        subset_path.split('/')[-1].split('.')[0] +
        '.' +
        subset_path.split('/')[-1].split('.')[1])
    
    #prepare specific output directories
    os.makedirs(f'{bootstrap_dir}/{subset_name}', exist_ok=True)

    iteractions = [300,  400,  500,  1000, 
                   1500, 3000, 3500, 4000, 
                   5000, 6000, 7000, 8000, 9000]
    
    for iteraction in iteractions:
        iter_name = str(iteraction)
        iteraction_dir = (f'{bootstrap_dir}/{subset_name}/{iter_name}')
        
        
        if(os.path.exists(f'{iteraction_dir}')):
            print(f'{subset_name} and {iteraction} already exists')
        else:
            os.makedirs(f'{iteraction_dir}', exist_ok=True)
    

            for seed in range(0, 20):
                    out_cor = f'{iteraction_dir}/cor_'f'{seed}''.cor'
                    out_cov = f'{iteraction_dir}/cov_'f'{seed}''.cov'
                    remove = iteraction / 2.5
                    remove = int(remove)
                
                    network_time = datetime.now()
                    completed_process = subprocess.run(['fastspar',
                        '-c', f'{subset_path}',
                        '-r', f'{out_cor}',
                        '-a', f'{out_cov}',
                        '-t', '5',
                        '-s', f'{seed}',
                        '-i', f'{iteraction}',
                        '-x', f'{remove}',
                        '-e', '0.1',
                        '-y'],
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE, text=True)
                    
                    log_filepath = f'{iteraction_dir}/log_'f'{seed}''.txt'
                
                    output = completed_process.stdout
                    with open(log_filepath, 'w') as log_file:
                        log_file.write(output)
                
                    network_time = datetime.now() - network_time
                    time_filepath = f'{iteraction_dir}/time_'f'{seed}''.txt'
                    with open(time_filepath, 'w') as time_file:
                        time_file.write(str(network_time))
                
                    print(f"Finished {subset_name} with "
                      f"{iteraction} iteractions and seed {seed}")
                    print("elapsed time in bootstrap: ",
                      datetime.now() - startTime)
         

#*************************** Similarity of matrices ****************************          
'''
def cosine_similarity(matrix1, matrix2):
    flattened1 = matrix1.flatten()
    flattened2 = matrix2.flatten()
    dot_product = np.dot(flattened1, flattened2)
    norm_product = norm(flattened1) * norm(flattened2)
    similarity = dot_product / norm_product
    return similarity

def calculate_similarity(matrices):
    num_matrices = len(matrices)
    similarity_matrix = numpy.zeros((num_matrices, num_matrices))

    for i in range(num_matrices):
        for j in range(i, num_matrices):  # Only calculate upper triangular matrix
            similarity = cosine_similarity(matrices[i], matrices[j])
            similarity_matrix[i, j] = similarity
            similarity_matrix[j, i] = similarity  # Assign to both symmetric positions

    return similarity_matrix


similarity_matrix = calculate_similarity(matrices)
print("Similarity matrix:")
print(similarity_matrix)


networks_dir = f'{data_dir}/fastspar_networks/'
os.makedirs(f'{networks_dir}', exist_ok=True)

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
        "-i 5000 "
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
        print('nElapsed time: ', (datetime.now()-startTime))
 

################ Preprocessing for Keystones Identification ####################
startTime = datetime.now()

#TODO: remove this when correlation function is fixed to accept tsv files
community_subsets = glob.glob(f'{data_dir}community_subsets_raw/*.csv')

for subset_path in community_subsets:

    #TODO: correct this for coherence with the rest of the code
    subset_name = (
        subset_path.split('/')[-1].split('.')[0] +
        '.' +
        subset_path.split('/')[-1].split('.')[1])

    meta = "none"
    
    # checking if there is a computed sparcc matrix
    if existOldData(f'{networks_dir}'
                    +subset_name
                    +'/sparcc/sparcc_data/cor.tsv'):
        
        print('Using %s for the analysis correlation.\n'%subset_name)
        spcc_corr_mat = (f'{networks_dir}'
                        +subset_name
                        +'/sparcc/sparcc_data/cor.tsv')
    else:
        print(f'NO SPARCC MATRIX FOUND! JUMPING {subset_name}')
        continue

    if existOldData(f'{networks_dir}'+subset_name+'/sparcc/keystones.csv'):        
        print(f'Keystones of {subset_name} already exist')
        #TODO: reactivete the line below if you want to delete 
        #the keystones tables
        #os.remove(f'{networks_dir}'+subset_name+'/sparcc/keystones.csv')
        
    else:    
        coSparCC = correlation_main(
            subset_path,
            meta,
            subset_name,
            spcc_corr_mat)
        print(f'correlattion of {subset_name} finished.')
        print('Progress: %s of %s' 
                % (
                community_subsets.index(subset_path) + 1,
                len(community_subsets)))
        print('\nElapsed time: ', (datetime.now()-startTime))

        #TODO: for some reason the code iteracts using the out folder as a
        # temporary folder, and the cpr function copies it to the output folder
        cpr(
            lsgrep('out',['']),
            f'{data_dir}fastspar_networks/'+subset_name+'/sparcc/')
        rmr(['out'])
        continue



################## Keystones identification ####################################
startTime = datetime.now()

print("\n\nStarting analysis through environments.\n\n")
os.makedirs(f'{data_dir}/final_keystones_table/', exist_ok=True)

files = glob.glob(f'{data_dir}fastspar_networks/*')
df = pd.DataFrame()

# reading all files
for i in files:
    try:
        peak = cat(i+'/sparcc/raw_data/cnm_highest_peak.txt').strip()
    except:
        print('Could not read data of '+i)
        continue

    # appending new table indexed by the Taxon name
    if os.path.exists(i+'/sparcc/figures/0p%s/keystones.csv'%peak):
        df = pd.concat(
            [df,
             pd.read_csv(
                 i+'/sparcc/figures/0p%s/keystones.csv'%(peak),
                 index_col=0)],
            sort=False)
    else:
        print('There is no data for %s\n'%i)

df.to_csv(f'{data_dir}/final_keystones_table/keystones.csv')
    
    #TODO: heatmaps are not working due to the lack of the grephi module
    # Generate the  keystones heatmap (total effect)
    # print("Heatmap keystones for %s." % level)
    # sexec('./src/heatmap_keystones.py '+level+' 1')

    # # Generate the  keystones heatmap (indirect effect)
    # print("Heatmap keystones (indirect) for %s." % level)
    # sexec('./src/heatmap_keystones.py '+level+' 2')
print('total time: ', (datetime.now()-startTime))
'''