#!/bin/bash
#Author: Bright Mage

echo "
################################################################################
#                                     Setup                                    #
################################################################################
Start time: $(date "+%Y-%m-%d %H:%M:%S")
"

source Shell/settings.sh
mkdir -p "logs"
mkdir -p Shell/jobs
mkdir -p data/$frame_analysis


echo "
$frame_analysis frame with the following parameters:"


echo "
//Inputs//
"
echo "annotated table
  $annotated_table"
echo "metadata_table
  $metadata"

echo "
//Environment//
"
echo "package manager
  ${package_manager}"

echo "
//Preprocess data//
"
echo "multiplier constant
  $multiplicative_const"
echo "number of samples
  $minimum_samples"

echo "
//Performance analysis//
"
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

echo "
//Final analysis//
"
echo "synthetic communities 
  $synthetic_communities"
echo "definitive iteration
  $definitive_iter"
echo "remove correlates
  $remove"

echo "
//Computational resources//
"
echo "parallel processes
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
python3 Python/pipelines/identifying_keystones.py $frame_analysis