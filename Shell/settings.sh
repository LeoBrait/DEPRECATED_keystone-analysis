#!/bin/bash


# Package manager path
package_manager="miniconda3"
source ~/$package_manager/etc/profile.d/conda.sh

# Performance analysis
iterations=(300 400 500 600 700 800 900 1000
1000 2000 3000 4000 5000 6000 7000 8000 9000 10000)

seeds=(1 2 3 4 5 6 7 8 9 10
11 12 13 14 15 16 17 18 19 20
21 22 23 24 25 26 27 28 29 30
31 32 33 34 35 36 37 38 39 40
41 42 43 44 45 46 47 48 49 50)

synthetic_communities=250

#fastspar settings
definitive_iter=10000
remove=15


# Computational resources
parallel=40
