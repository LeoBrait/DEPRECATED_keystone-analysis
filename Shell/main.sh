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
mkdir -p data/$analysis_frame


echo "
$analysis_frame frame with the following parameters:"


echo "
//Inputs//
"
echo "annotated table absolute
  $annotated_table_absolute"
echo "annotated table relative
  $annotated_table_relative"
echo "metadata_table
  $metadata_table"

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
    echo "  $iteration "
done

echo "iterations test subsets "
for subset in "${iterations_test_subsets[@]}"; do
    echo "  $subset "
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

source Shell/pipelines/calculating_fastspar.sh $analysis_frame
echo "End Of Fastspar pipeline: $(date "+%Y-%m-%d %H:%M:%S") "

echo "
######################################################## Keystone Identification
"

conda activate pyshell_biome_keystones
python3 Python/pipelines/identifying_keystones.py $analysis_frame

echo "
################################################################## Taxa Grouping
"
conda deactivate
Rscript R/pipelines/verifying_taxa_grouping.R \
  $analysis_frame \
  $annotated_table_relative \
  $metadata_table \
  $parallel \
  $minimum_samples