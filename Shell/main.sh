#!/bin/bash


#### Data pre-process
python3 Python/pipelines/data_preprocessing.py

################ performance fastspar iterations ##############################

# Environment ******************************************************************
communities_path=data/community_subsets
booststrap_subsets=( 
    "${communities_path}/animal_host-associated.aqueous_humour.tsv" #N=8
    "${communities_path}/animal_host-associated.animal_feces.tsv"   #N=675
    "${communities_path}/saline_water.coastal_seawater.tsv"         #N=286
    "${communities_path}/saline_water.hypersaline_water.tsv"        #N=16
    "${communities_path}/soil.savanna_soil.tsv"                     #N=21
    "${communities_path}/soil.tundra_soil.tsv"                      #N=3
    "${communities_path}/groundwater.porous_contaminated.tsv"       #N=48
    "${communities_path}/groundwater.mine.tsv")                     #N=3

iterations=( 300  400  500 1000
             1500 3000 3500 4000
             5000 6000 7000 8000 9000)

seeds=( 1  2  3  4  5  6  7  8  9 10
       11 12 13 14 15 16 17 18 19 20)

performance_dir=data/performance_fastspar_iterations
mkdir -p "${performance_dir}"

# Job List *********************************************************************
for subset_path in "${booststrap_subsets[@]}"
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
                echo "echo The job for: " \
                        "${filename} with ${iteration} iteration," \
                            "seed ${seed} was already done" > general_log.txt
            else                
                #create the jobs file in jobs folder
                echo "echo The job for:" \
                " ${filename} with ${iteration}, seed ${seed} is running..."
                echo "mkdir -p ${habitat_dir}"
                echo "mkdir -p ${iteration_dir}"
                echo "fastspar " \
                    "-c ${subset_path} "\
                    "-r ${out_cor} " \
                    "-a ${out_cov} " \
                    "-t 2 " \
                    "-s $seed " \
                    "-i $iteration " \
                    "-x $remove " \
                    "-e 0.1 " \
                    "-y > ${log} "
                echo "echo ${filename} with ${iteration}, seed ${seed} done!"
                echo
            fi
        done
    done
done > Shell/jobs/performance_iterations.txt


# Run the jobs *****************************************************************

xargs -P 20 -I {} bash -c "{}" < Shell/jobs/performance_iterations.txt

# Maxrix similarities **********************************************************

######################### Fastspar P-values ####################################

# # Generate fake habitats
# iter=3000
# path=${base}/boot_iter_${iter}

# for seed in {1..2}; do
#     sseed=$(printf "%03d" $seed)
#     mkdir -p ${path}/${sseed}/anot ${path}/${sseed}/cor ${path}/${sseed}/cov
#     fastspar_bootstrap -c ${base}/$anot -n 1000 -p ${path}/${sseed}/anot/a -t $nthreads -s $seed
# done

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