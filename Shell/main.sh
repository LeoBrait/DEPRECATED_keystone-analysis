#!/bin/bash
#Author: Bright Mage

# Environment Settings
package_manager="miniconda3"
source ~/$package_manager/etc/profile.d/conda.sh

# Computational resources
parallel=40

######################### Data pre-process #####################################


conda activate pyshell_biome_keystones
python3 Python/pipelines/data_preprocessing.py

######################### performance fastspar iterations ######################

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
              1500  3000  3500  4000
              5000  6000  7000  8000
              9000  12000 14000 16000)

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


# Run the jobs *****************************************************************

xargs -P $parallel -I {} bash -c "{}" < Shell/jobs/performance_iterations.txt

# Maxrix similarities **********************************************************

# conda activate R_biome_keystones
# Rscript R/pipelines/measuring_matrix_similarities.R


######################### Fastspar P-values ####################################

# Environment ******************************************************************
# conda activate pyshell_biome_keystones
# source ~/Shell/pipelines/fastspar_pvalues.sh

confirm() {
    read -N 1 REPLY
    echo
    if [[ "$REPLY" = "y" || "$REPLY" = "Y" ]]; then
        "$@"
    else
        echo "Cancelled by user"
        exit 1
    fi
}
echo "Do you want to create fake habitats?"
confirm source ~/Shell/pipelines/creating_fake_habitats.sh


## ask to run the line bellow
# echo "Do you want to run a line?"
# read -p "Enter the line number: " line
# source ~/Shell/pipelines/fastspar_pvalues.sh

# echo "Do you want to run a line?"
# read -p "Enter the line number: " line
# #do fastspar on fake habitats
# for seed in {1..2}; do
#     sseed=$(printf "%03d" $seed)
#     for i in {0..999}; do
#         if [[ ! -f ${path}/${sseed}/cor/${i}.tsv ]]; then
#             fastspar -c ${path}/${sseed}/anot/a_${i}.tsv -r ${path}/${sseed}/cor/${i}.tsv -a ${path}/${sseed}/cov/${i}.tsv -t $nthreads -s 1 -i $iter -x $xter -e $corThr -y
#         fi
#     done
# done

# # Test real correlations vs. fake *********************************************
# for seed in {1..2}; do
#     sseed=$(printf "%03d" $seed)
#     fastspar_pvalues -c ${base}/$anot -r ${path}/${sseed}/cor.tsv -n 1000 -p ${path}/${sseed}/cor/ -o ${path}/${sseed}/pval.tsv -t $nthreads
# done