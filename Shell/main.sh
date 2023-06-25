#!/bin/bash
#Author: Bright Mage

# Environment Settings
package_manager="miniconda3"
source ~/$package_manager/etc/profile.d/conda.sh

# Computational resources
parallel=40
definitive_iter=4000

echo "
################################################################################
#                              Preprocessing Data                              #
################################################################################
Start time: $(date "+%Y-%m-%d %H:%M:%S")"


conda activate pyshell_biome_keystones
python3 Python/pipelines/data_preprocessing.py


echo "
################################################################################
#                     performancing fastspar iterations                        #
################################################################################
Start time: $(date "+%Y-%m-%d %H:%M:%S")"


# Environment ******************************************************************
communities_path=data/community_subsets
iterations_test_subsets=( 
    "${communities_path}/animal_host-associated.aqueous_humour.tsv" #N=8
    "${communities_path}/animal_host-associated.animal_feces.tsv"   #N=675
    "${communities_path}/saline_water.coastal_seawater.tsv"         #N=286
    "${communities_path}/saline_water.hypersaline_water.tsv"        #N=16
    "${communities_path}/soil.savanna_soil.tsv"                     #N=21
    "${communities_path}/soil.tundra_soil.tsv"                      #N=3
    "${communities_path}/groundwater.porous_contaminated.tsv"       #N=48
    "${communities_path}/groundwater.mine.tsv")                     #N=3

iterations=(   300   400   500  1000
              1500  3000  3500  4000)

seeds=( 1  2  3  4  5  6  7  8  9 10
       11 12 13 14 15 16 17 18 19 20
       21 22 23 24 25 26 27 28 29 30
       31 32 33 34 35 36 37 38 39 40
       41 42 43 44 45 46 47 48 49 50)

performance_dir=data/performance_fastspar_iterations
mkdir -p "${performance_dir}"

# Job List *********************************************************************
for subset_path in "${iterations_test_subsets[@]}"
do
    #parse the file name without the .tsv
    filename=$(basename -- "$subset_path")
    filename="${filename%.*}"
    habitat_dir="${performance_dir}/${filename}"

    for iteration in "${iterations[@]}"
    do
        iteration_dir="${habitat_dir}/${iteration}"

        for seed in "${seeds[@]}"
        do
            #Environment for each iteration and seed
            out_cor="${iteration_dir}/cor_${seed}.cor"
            out_cov="${iteration_dir}/cov_${seed}.cov"
            log="${iteration_dir}/log_${seed}.txt"
            time_var="${iteration_dir}/time_${seed}.txt"
            remove=15
            
            #check the pre-existence of each output
            if [ -f "${out_cor}" ]; then
                echo "echo The Sparcc for: " \
                        "${filename} with ${iteration} iterations," \
                            "seed ${seed} was already done" > general_log.txt
            else                
                #create the jobs file in jobs folder
                echo "echo The Sparcc for:" \
                        " ${filename} with ${iteration} iterations," \
                            "seed ${seed} is running..."
                echo "mkdir -p ${habitat_dir}"
                echo "mkdir -p ${iteration_dir}"
                echo "fastspar "\
                    "-c ${subset_path} "\
                    "-r ${out_cor} "\
                    "-a ${out_cov} "\
                    "-t 2 "\
                    "-s $seed "\
                    "-i $iteration "\
                    "-x $remove "\
                    "-e 0.1 "\
                    "-y > ${log}"
                echo "echo ${filename} with ${iteration} iterations," \
                        " and seed ${seed} done!"
                echo
            fi
        done
    done
done > Shell/jobs/performance_iterations.txt

# Run the jobs 
xargs -P $parallel -I {} bash -c "{}" < Shell/jobs/performance_iterations.txt

# Maxrix similarities
# conda activate R_biome_keystones
# Rscript R/pipelines/measuring_matrix_similarities.R

echo "
################################################################################
#                          Fastspar for all communities                        #
################################################################################
Start time: $(date "+%Y-%m-%d %H:%M:%S")"


conda activate pyshell_biome_keystones
communities_path=data/community_subsets
tablenames=($(\ls ${communities_path}))
fastspar_dir=data/fastspar_correlations
mkdir -p "${fastspar_dir}"

# Job List
for tablename in "${tablenames[@]}";
do
    #parse the file name without the .tsv
    filename=$(basename -- "$tablename")
    filename="${filename%.*}"

    #create the jobs file in jobs folder
    echo "echo The Sparcc for: ${filename} is running..."
    echo "mkdir -p ${fastspar_dir}/${filename}"
    echo "fastspar "\
            "-c ${communities_path}/${tablename} "\
            "-r ${fastspar_dir}/${filename}/${filename}.cor "\
            "-a ${fastspar_dir}/${filename}/${filename}.cov "\
            "-t 2 "\
            "-s 1 "\
            "-i $definitive_iter "\
            "-x 15 "\
            "-e 0.1 "\
            "-y > ${fastspar_dir}/${filename}/${filename}.log"
    echo "echo ${filename} done!"
    echo
done > Shell/jobs/fastspar_all_communities.txt


echo "
################################################################################
#                              Fastspar P-values                               #
################################################################################
Start time: $(date "+%Y-%m-%d %H:%M:%S")"

# Generate Synthetic habitats
if [ ! -d "data/synthetic_habitats" ]; then
    source Shell/pipelines/creating_fake_habitats.sh
    else
        echo "Seems that the synthetic habitats are already created.
                jumping to the next step..."
fi

# Run fastspar with synthetic communities
source Shell/pipelines/fastspar_synthetics.sh

# # Test real correlations vs. fake ********************************************

# Environment
# communities_path=data/community_subsets
# tablenames_real=($(\ls ${communities_path}))
# p_values_dir=data/fastspar_pvalues
# mkdir -p "${p_values_dir}"

# # Job List
# for tablename_real in "${tablenames_real[@]}";
# do
#     #parse the file name without the .tsv
#     filename=$(basename -- "$tablename_real")
#     filename="${filename%.*}"

#     #create the jobs file in jobs folder
#     echo "echo The p-values for: ${filename} is running..."
#     echo "python3 Python/pipelines/fastspar_pvalues.py "\
#             "-c ${communities_path}/${tablename_real} "\
#             "-r ${p_values_dir}/${filename}.cor "\
#             "-n 1000 "\
#             "-p ${p_values_dir}/${filename}/ "\
#             "-o ${p_values_dir}/${filename}.tsv "\
#             "-t 2"



#fastspar_pvalues -c ${base}/$anot -r ${path}/${sseed}/cor.tsv -n 1000 -p ${path}/${sseed}/cor/ -o ${path}/${sseed}/pval.tsv -t $nthreads
echo "End time: $(date "+%Y-%m-%d %H:%M:%S")"