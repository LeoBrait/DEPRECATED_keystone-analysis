#!/usr/bin/python3

import os
import sys
sys.path.append('src')
import pandas as pd

from xonsh_py import cat, lsgrep

level = sys.argv[1]

files = lsgrep('output/transposed',[level])

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
        df = df.append(pd.read_csv(i+'/sparcc/figures/0p%s/keystones.csv'%(peak),index_col=0),sort=False)
    else:
        print('There is no data for %s\n'%i)
df.to_csv('output/transposed_all_environments/%s/keystones.csv'%(level))
