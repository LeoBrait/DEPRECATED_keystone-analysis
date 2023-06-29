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

analysis_frame=str(sys.argv[1])


# Paths
src_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(src_dir, '../..', 'data/')
global_dir = os.path.join(src_dir, '../..')

# Custom functions
sys.path.append(f'{src_dir}/../src')
from xonsh_py import Logger, lsgrep, existOldData, mkdir_p, cpr, rmr, sexec, cat
from correlation import main as correlation_main
from identify_keystones import identify_keystones
#sys.stderr = sys.stdout = Logger()
total_time = datetime.now()



################ Preprocessing for Keystones Identification ####################
startTime = datetime.now()
networks_dir=(f"{data_dir}{analysis_frame}/'fastspar_correlations/")


#TODO: remove this when correlation function is fixed to accept tsv files
community_subsets = glob.glob(f'{data_dir}{analysis_frame}/'
                                  'community_subsets_raw/*.csv')

for subset_path in community_subsets:

    #TODO: correct this for coherence with the rest of the code
    subset_name = (
        subset_path.split('/')[-1].split('.')[0] +
        '.' +
        subset_path.split('/')[-1].split('.')[1])

    meta = "none"

    # checking if there is a computed sparcc matrix
    if existOldData(f'{networks_dir}/{subset_name}/cor_{subset_name}'):
        
        print('Using %s for the analysis correlation.\n'%subset_name)
        spcc_corr_mat = (f'{networks_dir}'
                        +subset_name
                        +f'/cor_{subset_name}')
    else:
        print(f'NO SPARCC MATRIX FOUND! JUMPING {subset_name}')
        continue



    if existOldData(
        f'{data_dir}/{analysis_frame}/'
          f'parcial_keystones_results/{subset_name}/keystones.csv'):

        print(f'Keystones of {subset_name} already exist')

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
        os.makedirs(f'{data_dir}/{analysis_frame}/'
                      f'final_keystones_table/', exist_ok=True)
        cpr(
            lsgrep('out',['']),
            f'{data_dir}/{analysis_frame}/'
              f'parcial_keystones_results/'+subset_name)
        rmr(['out'])
        continue



################## Keystones identification ####################################
startTime = datetime.now()

print("\n\nStarting analysis through environments.\n\n")
os.makedirs(f'{data_dir}/{analysis_frame}/'
              f'final_keystones_table/', exist_ok=True)

files = glob.glob(f'{data_dir}/{analysis_frame}/'
                    'parcial_keystones_results/*')
df = pd.DataFrame()

# reading all files
for i in files:
    try:
        peak = cat(i+'/raw_data/cnm_highest_peak.txt').strip()
    except:
        print('Could not read data of '+i)
        continue

    # appending new table indexed by the Taxon name
    if os.path.exists(i+'/figures/0p%s/keystones.csv'%peak):
        df = pd.concat(
            [df,
             pd.read_csv(
                 i+'/figures/0p%s/keystones.csv'%(peak),
                 index_col=0)],
            sort=False)
    else:
        print('There is no data for %s\n'%i)

df.to_csv(f'{data_dir}/{analysis_frame}/'
            f'final_keystones_table/keystones.csv')
    
#     #TODO: heatmaps are not working due to the lack of the grephi module
#     # Generate the  keystones heatmap (total effect)
#     # print("Heatmap keystones for %s." % level)
#     # sexec('./src/heatmap_keystones.py '+level+' 1')

#     # # Generate the  keystones heatmap (indirect effect)
#     # print("Heatmap keystones (indirect) for %s." % level)
#     # sexec('./src/heatmap_keystones.py '+level+' 2')
# print('total time: ', (datetime.now()-startTime))