#!/bin/bash


# Package manager path
package_manager="miniconda3"
source ~/$package_manager/etc/profile.d/conda.sh

# Performance analysis
iterations=( 300   400)
seeds=(1  2)
definitive_iter=4000
synthetic_communities=100

# Computational resources
parallel=40