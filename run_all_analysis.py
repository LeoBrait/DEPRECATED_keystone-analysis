#!/usr/bin/python3

import sys
sys.path.append('src')
from xonsh_py import *

# this script aims to execute all analysis that needs more than one Ecossistem/habitat to work

print("\n\nStarting analysis through environments.\n\n")
# steps that uses all the generated data (usually used for analysis among all environments and taxonomic levels)
for level in ['phyla']:

    # mkdir -p @('output/transposed_all_environments/'+level)
    mkdir_p(['output/transposed_all_environments/'+level])

    print("Concat keystones for %s." % level)
    sexec('./src_all/concat_keystones.py '+level)

    print("Heatmap keystones for %s." % level)
    sexec('./src_all/heatmap_keystones.py '+level+' 1')

    print("Heatmap keystones (indirect) for %s." % level)
    sexec('./src_all/heatmap_keystones.py '+level+' 2')
