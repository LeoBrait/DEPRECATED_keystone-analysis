#!/bin/bash


#### Data pre-process
python3 Python/pipelines/data_preprocessing.py

################ performance fastspar iteractions ##############################

# Environment **************************
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

iteractions=( 300  400  500 1000
             1500 3000 3500 4000
             5000 6000 7000 8000 9000)

seeds=( 1  2  3  4  5  6  7  8  9 10
       11 12 13 14 15 16 17 18 19 20)

performance_dir=data/performance_fastspar_iteractions
mkdir -p "${performance_dir}"

# Job List ***************************
job_creation_time=$(date +%Y-%m-%d_%H-%M-%S)
for subset_path in "${booststrap_subsets[@]}"
do
    #parse the file name without the .tsv
    filename=$(basename -- "$subset_path")
    filename="${filename%.*}"
    habitat_dir="${performance_dir}/${filename}"

    for iteraction in "${iteractions[@]}"
    do
        iteraction_dir="${habitat_dir}/${iteraction}"

        for seed in "${seeds[@]}"
        do
            out_cor="${iteraction_dir}/cor_${seed}.cor"
            out_cov="${iteraction_dir}/cov_${seed}.cov"
            log="${iteraction_dir}/log_${seed}.txt"
            time_var="${iteraction_dir}/time_${seed}.txt"
            if [ -f "${out_cor}" ]; then
                echo "echo The job for: " \
                    "${filename} with ${iteraction} an ${seed} was already done"
            else

                #Environment for each iteraction and seed
                remove=$(echo "scale=2; $iteraction / 10" | bc | cut -d '.' -f 1)
            
                #create the jobs
                echo "echo The job for:" \
                " ${filename} with ${iteraction} an ${seed} is running"
                echo "mkdir -p ${habitat_dir}"
                echo "mkdir -p ${iteraction_dir}"
                echo "fastspar " \
                    "-c ${subset_path} "\
                    "-r ${out_cor} " \
                    "-a ${out_cov} " \
                    "-t 2 " \
                    "-s ${seed} " \
                    "-i ${iteraction} " \
                    "-x ${remove} " \
                    "-e 0.1 " \
                    "-y > ${log} "
                echo
            fi
        done
    done
done > Shell/jobs/performance_iteractions.txt


# Run the jobs 5 process
xargs -P 5 -I {} bash -c "{}" < Shell/jobs/performance_iteractions.txt