#!/usr/bin/env python3

import os
import pandas as pd
import sys
import os

sys.path.append("src")
from xonsh_py import cat, lsgrep

level = sys.argv[1]

files = lsgrep('output/environments',[level])

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
        df = pd.concat([df, pd.read_csv(i+'/figures/0p%s/keystones.csv'%(peak),index_col=0)],sort=False)
    else:
        print('There is no data for %s\n'%i)
df.to_csv('output/all_environments/%s/keystones.csv'%(level))
