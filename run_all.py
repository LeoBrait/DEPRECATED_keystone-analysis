#!/usr/bin/python3

import os
import sys
sys.path.append('src')
from xonsh_py import Logger, lsgrep, existOldData, mkdir_p, cpr, rmr, sexec
sys.stderr = sys.stdout = Logger() # check src/xonsh_py.py for details
from datetime import datetime

from correlation import main as correlation_main

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
for f in lsgrep('community_matrix',['phyla']):

    # get the environment name
    fname = f[:-4].split('/')[1].split('_')[0]

    # checking if the actual file has metadata
    meta = lsgrep('metadata',[fname])
    meta = meta[0] if meta else "none"

    print('Community: %s\nMetadata: %s\n' % (f, meta))

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
    correlation_main(f, meta, fname, spcc_corr_mat)

    # copy the results to the corresponding environment directory
    cpr(lsgrep('out',['']),'output/transposed/'+fname+'/sparcc/')
    # delete the temporary 'out' directory
    rmr(['out'])

# record end time   
print('\nTotal execution time: ', (datetime.now()-startTime))

# running all analysis script
startTime = datetime.now() # record start time
sexec('./run_keystone_analysis.py')
print('\nPost-results execution time: ', (datetime.now()-startTime)) # record end time
