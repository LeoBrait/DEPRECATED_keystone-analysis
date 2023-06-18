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

entrada = "oi"
subprocess.run(['echo', f'{entrada}'])
