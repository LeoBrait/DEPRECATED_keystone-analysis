#!/bin/bash
#Author: Bright Mage

echo "
################################################################################
#                                     Setup                                    #
################################################################################
Start time: $(date "+%Y-%m-%d %H:%M:%S")
"

source Shell/settings.sh
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
echo "samples minimum threshold
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
################################################################## POSPROCESSING
"

conda activate R_biome_keystones
Rscript R/pipelines/plotting_iterations_performance.R $frame_analysis