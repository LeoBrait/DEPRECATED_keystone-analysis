#!/bin/bash
#Author: Bright Mage

echo "
################################################################################
#                                     Setup                                    #
################################################################################
Start time: $(date "+%Y-%m-%d %H:%M:%S")
"
mkdir -p logs
mkdir -p Shell/jobs
source Shell/settings.sh

echo "
starting analysis with the following parameters:
"

echo "package manager
  ${package_manager}"

echo "iterations "
for iteration in "${iterations[@]}"; do
    echo "  $iteration"
done

echo "iterations test subsets "
for subset in "${iterations_test_subsets[@]}"; do
    echo "  $subset"
done

echo "seeds
  ${#seeds[@]}"


echo "synthetic communities 
  $synthetic_communities
definitive iteration
  $definitive_iter
multiplier constant
  $multi_const
remove correlates
  $remove
parallel processes
  $parallel"

echo "
####################################################################### FASTSPAR
"

source Shell/pipelines/calculating_fastspar.sh
echo "End Of Fastspar pipeline: $(date "+%Y-%m-%d %H:%M:%S")"

echo "
######################################################## Keystone Identification
"

conda activate pyshell_biome_keystones
python3 Python/pipelines/identifying_keystones.py