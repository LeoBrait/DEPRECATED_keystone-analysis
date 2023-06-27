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
remove correlates
  $remove
parallel processes
  $parallel"

echo "
################################################################################
#                               FASTSPAR PIPELINE                              #
################################################################################
Start time: $(date "+%Y-%m-%d %H:%M:%S")
"

source Shell/pipelines/calculating_fastspar.sh

#end time
echo "End Of Fastspar pipeline: $(date "+%Y-%m-%d %H:%M:%S")"