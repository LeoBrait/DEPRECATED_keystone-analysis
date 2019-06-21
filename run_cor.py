#!/usr/bin/python3

import os
import sys
sys.path.append('src')
from xonsh_py import *
sys.stderr = sys.stdout = Logger() # check src/xonsh_py.py for details
from datetime import datetime

from correlation import main

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

startTime = datetime.now()

for f in lsgrep('community_matrix',['phyla']):
    if 'animal.manure' not in f:
        continue

    fname = f[:-4].split('/')[1].split('_')[0]

    # checking if the actual file has metadata
    meta = lsgrep('metadata',[fname])
    meta = meta[0] if meta else "none"

    print('Community: %s\nMetadata: %s\n' % (f, meta))
    # checking if it has an already computed sparcc matrix
    if existOldData('output/transposed/'+fname+'/sparcc/sparcc_data/cor.tsv'):
        print('There is an already computed SparCC matrix for %s. Using it for the analysis.\n'%fname)
        spcc_corr_mat = 'output/transposed/'+fname+'/sparcc/sparcc_data/cor.tsv'
    else:
        print('NO SPARCC MATRIX FOUND! JUMPING')
        continue

    mkdir_p(['output/transposed/'+fname])

    main(f, meta, fname, spcc_corr_mat)

    cpr(lsgrep('out',['']),'output/transposed/'+fname+'/sparcc/')
    rmr(['out'])

print('\nTotal execution time: ', (datetime.now()-startTime))

# running all analysis script
startTime = datetime.now()
sexec('./run_all_analysis.py')
print('\nPost-results execution time: ', (datetime.now()-startTime))
