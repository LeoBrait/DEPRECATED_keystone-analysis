#!/usr/bin/env python3

import sys
sys.path.append('src')

from xonsh_py import Logger, lsgrep, existOldData, mkdir_p, cpr, rmr, sexec
sys.stderr = sys.stdout = Logger() # check src/xonsh_py.py for details
from datetime import datetime

from correlation import main as correlation_main
from identify_keystones import identify_keystones

# checking for dependencies
try:
    import numpy
    import pandas
    import gephitools
    import scipy
    import matplotlib
    import mpl_toolkits
    import networkx
    import PIL
except ImportError as e:
    print(e)
    sys.exit()

# record start time
startTime = datetime.now()

# loops through all the environments' community matrices
# the community matrices are stored in the 'community_matrix' directory 
# the file name should contain 'phyla' in accordance with the level of organization
# being considered.
# the metadata is stored in the 'metadata' directory and should have an identical
# folder structure.
for f in lsgrep('community_matrix',['phyla']):

    # get the environment name
    #ATTENTION: this is likely to break if the file name is not in the expected format
    fname = f[:-4].split('/')[1].split('_')[0]

    # checking if the actual file has metadata
    meta = lsgrep('metadata',[fname])
    meta = meta[0] if meta else "none"

    print('Community: %s\nMetadata: %s\n' % (f, meta))

    # checking if there is a computed sparcc matrix
    if existOldData('output/environments/'+fname+'/sparcc_data/cor.tsv'):
        print('There is an already computed SparCC matrix for %s. Using it for the analysis.\n'%fname)
        spcc_corr_mat = 'output/environments/'+fname+'/sparcc_data/cor.tsv'
    elif existOldData('input/environments/'+fname+'/sparcc_data/cor.tsv.gz'):
        #BUG
        pass
    else:
        print(f'NO SPARCC MATRIX FOUND! JUMPING {fname}')
        continue

    # Create output directory
    mkdir_p(['output/environments/'+fname])

    # run the analysis
    # f : community matrix
    # meta : metadata
    # fname : environment name
    # spcc_corr_mat : sparcc matrix
    coSparCC = correlation_main(f, meta, fname, spcc_corr_mat)

    # running code-integration steps
    # this step is responsible for running the CNM and LIASP algorithms
    # and identifying the most important nodes in the network (the keystones)
    identify_keystones("sparcc", coSparCC, fname)

    # copy the results to the corresponding environment directory
    cpr(lsgrep('out',['']),'output/environments/'+fname+'/')
    # delete the temporary 'out' directory
    rmr(['out'])

# record end time   
print('\nTotal execution time: ', (datetime.now()-startTime))

# running all analysis script
startTime = datetime.now() # record start time
sexec('./run_keystone_analysis.py') # Keystones analysis
print('\nPost-results execution time: ', (datetime.now()-startTime)) # record end time
