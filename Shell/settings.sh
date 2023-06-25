#!/bin/bash


# Package manager path
package_manager="miniconda3"
source ~/$package_manager/etc/profile.d/conda.sh

# Performance analysis
iterations=( 300   400)
seeds=(1  2)
synthetic_communities=100

#fastspar settings
remove=15
definitive_iter=300

# Computational resources
parallel=40
