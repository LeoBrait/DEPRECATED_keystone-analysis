#!/bin/bash
#Author: Bright Mage

#redundant line for modularization
source Shell/settings.sh


echo "
################################################################################
#                              Preprocessing Data                              #
################################################################################
Start time: $(date "+%Y-%m-%d %H:%M:%S")
"

conda activate pyshell_biome_keystones
python3 Python/pipelines/preprocessing_data.py \
    $frame_analysis \
    $multiplicative_const \
    $minimum_samples \
    $annotated_table \
    $metadata_table \


echo "
################################################################################
#                     performancing fastspar iterations                        #
################################################################################
Start time: $(date "+%Y-%m-%d %H:%M:%S")
"

performance_dir="data/${frame_analysis}/performance_fastspar_iterations"
mkdir -p "${performance_dir}"

# Job List
for subset_path in "${iterations_test_subsets[@]}"
do
    #parse names
    table_name=$(basename -- "$subset_path")
    habitat_name="${table_name%.*}"

    for iteration in "${iterations[@]}"
    do
        iteration_dir="${performance_dir}/${habitat_name}/${iteration}"

        for seed in "${seeds[@]}"
        do
           

            #check the pre-existence of each output
            if [ -f "${iteration_dir}/cor_${seed}" ];
            then
                echo "echo The Sparcc for: " \
                        "${habitat_name} with ${iteration} iterations," \
                            "seed ${seed} was already"\
                                        "done" >> logs/jump_log.txt
            else                
                mkdir -p "${performance_dir}/${habitat_name}"
                mkdir -p "${iteration_dir}"

                #create the jobs file in jobs folder
                echo "echo The Sparcc for:" \
                         "${habitat_name} with ${iteration} iterations," \
                            "seed ${seed} is running... | fold -w 80" 
                echo "fastspar "\
                    "-c ${subset_path} "\
                    "-r ${iteration_dir}/cor_${seed}"\
                    "-a ${iteration_dir}/cov_${seed}"\
                    "-t 2 "\
                    "-s $seed "\
                    "-i $iteration "\
                    "-x $remove "\
                    "-e 0.1 "\
                    "-y > ${iteration_dir}/log_${seed}.txt"
                echo "echo ${habitat_name} with ${iteration} iterations"\
                        "and seed ${seed} done! | fold -w 80"
                echo
            fi
        done
    done
done > Shell/jobs/performance_iterations.txt

# Run the jobs 
xargs -P $parallel -I {} bash -c "{}" < Shell/jobs/performance_iterations.txt


echo "
################################################################################
#                          Fastspar for all communities                        #
################################################################################
Start time: $(date "+%Y-%m-%d %H:%M:%S")"


communities_path="data/${frame_analysis}/community_subsets"
tablenames=($(\ls ${communities_path}))
fastspar_dir="data/${frame_analysis}/fastspar_correlations"
mkdir -p "${fastspar_dir}"

# Job List
for tablename in "${tablenames[@]}";
do
    #parse the file name without the .tsv
    habitat_name="${tablename%.*}"
    mkdir -p "${fastspar_dir}/${habitat_name}"

    if [ -f "${fastspar_dir}/${habitat_name}/cor_${habitat_name}" ];
    then
        echo "echo The Sparcc for: ${habitat_name}"\
        "was already done" >> logs/jump_log.txt
    else
    #create the jobs file in jobs folder
        echo "echo The Sparcc for: ${habitat_name} is running..."
        echo "fastspar "\
            "-c ${communities_path}/${tablename} "\
            "-r ${fastspar_dir}/${habitat_name}/cor_${habitat_name}"\
            "-a ${fastspar_dir}/${habitat_name}/cov_${habitat_name}"\
            "-t 2 "\
            "-s 1 "\
            "-i $definitive_iter "\
            "-x $remove "\
            "-e 0.1 "\
            "-y > ${fastspar_dir}/${habitat_name}/log_${habitat_name}.txt"
        echo "echo ${habitat_name} done!"
        echo
    fi
done > Shell/jobs/fastspar_all_communities.txt

xargs -P $parallel -I {} bash -c "{}" < Shell/jobs/fastspar_all_communities.txt

echo "
################################################################################
#                          Generate Synthetic habitats                         #
################################################################################
Start time: $(date "+%Y-%m-%d %H:%M:%S")"

conda activate pyshell_biome_keystones

#Generate Synthetic habitats
if [ ! -d "data/${frame_analysis}/synthetic_habitats" ]; then

    # Environment
    communities_path="data/${frame_analysis}/community_subsets"
    tablenames_real=($(\ls ${communities_path}))
    general_synthetics_dir="data/${frame_analysis}/synthetic_habitats"
    mkdir -p "${general_synthetics_dir}"
    mkdir -p "logs/synthetic_habitats_logs"

    # Job List
    for real_tablename in "${tablenames_real[@]}"; 
    do
        #parse the file name without the .tsv
        habitat="${real_tablename%.*}"
        mkdir -p "${general_synthetics_dir}/${habitat}"

        #create the jobs file in jobs folder
        echo "echo fake_habitat for: ${habitat} is running..."
        echo "fastspar_bootstrap" \
            "-c ${communities_path}/${real_tablename}" \
            "-n $synthetic_communities" \
            "-p ${general_synthetics_dir}/${habitat}/synt_" \
            "-t 2 " \
            "-s 1 > logs/synthetic_habitats_logs/log_${habitat}.txt"
        echo "echo fake ${habitat} done!"
        echo
    done > Shell/jobs/synthetic_habitats.txt

# Run the jobs
xargs -P $parallel -I {} bash -c "{}" < Shell/jobs/synthetic_habitats.txt

else
    echo "Seems that the synthetic habitats are already created.
                jumping to the next step..."
fi


echo "
################################################################################
#                        Fastspar for synthetic habitats                       #
################################################################################
Start time: $(date "+%Y-%m-%d %H:%M:%S")"

# input
general_synthetics_dir="data/${frame_analysis}/synthetic_habitats"
synt_habitats_dirs=($(ls ${general_synthetics_dir}))

# output
mkdir -p "data/${frame_analysis}/synthetic_fastspar"

for synt_habitat_dir in "${synt_habitats_dirs[@]}";
do
    tables=($(ls ${general_synthetics_dir}/${synt_habitat_dir}))

    for table in "${tables[@]}";
    do
        #parse the file name without the .tsv
        filename=$(basename -- "$table")
        filename="${filename%.*}"
        habitat=$(basename -- "${synt_habitat_dir}")

        #get table number
        table_number=$(echo "$filename" | awk -F'_' '{print $NF}')

        #Environment
        synt_fastspar_dir="data/${frame_analysis}/synthetic_fastspar/${habitat}"

    if [ -f "${synt_fastspar_dir}/cor_${habitat}_${table_number}" ]; then
        echo "echo The Sparcc for: " \
                "${habitat} " \
                   "table ${table_number} was already done" >> logs/jump_log.txt
    else
        mkdir -p ${synt_fastspar_dir}

        #create the jobs file in jobs folder
         echo "echo The Sparcc for:" \
                " synthetic ${habitat}," \
                    "table ${table_number} is running..."
        echo "fastspar "\
                "-c ${general_synthetics_dir}/${synt_habitat_dir}/${table}"\
                "-r ${synt_fastspar_dir}/cor_${habitat}_${table_number}"\
                "-a ${synt_fastspar_dir}/cov_${habitat}_${table_number}"\
                "-t 2 "\
                "-s 1 "\
                "-i $definitive_iter "\
                "-x $remove "\
                "-e 0.1 "\
                "-y > ${synt_fastspar_dir}/log_${habitat}_${table_number}.txt"
        echo "echo table ${table_number} of ${habitat}" \
                         " done! | fold -w 80"
        echo       
    fi
    done
done > Shell/jobs/fake_fastspar.txt

# Run the jobs
xargs -P $parallel -I {} bash -c "{}" < Shell/jobs/fake_fastspar.txt


echo "
################################################################################
#                          P-values: real vs synthetic                         #
################################################################################
Start time: $(date "+%Y-%m-%d %H:%M:%S")"

#Environment
p_values_dir="data/${frame_analysis}/fastspar_pvalues"
mkdir -p "${p_values_dir}"

#real habitats
communities_path="data/${frame_analysis}/community_subsets"
tablenames_real=($(\ls ${communities_path}))

#synthetic habitats
general_synthetics_dir="data/${frame_analysis}/synthetic_habitats"
synt_fastspar_dir="data/${frame_analysis}/synthetic_fastspar"

#real correlations
fastspar_dir="data/${frame_analysis}/fastspar_correlations"
habitat_dirs=($(ls ${fastspar_dir}))

for tablename_real in "${tablenames_real[@]}";
do  
    habitat_name="${tablename_real%.*}"

    if [ -f "${p_values_dir}/${habitat_name}/real.tsv" ]; then
        echo "echo The p-values for: " \
                "${habitat_name} " \
                   "was already done" >> logs/jump_log.txt
    else
        #create the jobs file in jobs folder
        echo "echo The p-values for:" \
                " ${habitat_name} is running..."
        echo "fastspar_pvalues "\
                "-c ${communities_path}/${tablename_real}"\
                "-r ${fastspar_dir}/${habitat_name}/cor_${habitat_name}"\
                "-p ${synt_fastspar_dir}/${habitat_name}/cor_"\
                "-n $synthetic_communities"\
                "-o ${p_values_dir}/${habitat_name}/real.tsv"\
                "-t 2"
        echo "echo ${habitat_name} done!"
        echo
    fi
done > Shell/jobs/fastspar_pvalues.txt

# Run the jobs
xargs -P $parallel -I {} bash -c "{}" < Shell/jobs/fastspar_pvalues.txt

echo "
################################################################################
#                       Generating performance tables                          #
################################################################################
Start time: $(date "+%Y-%m-%d %H:%M:%S")"


conda activate R_biome_keystones
Rscript R/pipelines/measuring_matrix_similarities.R $frame_analysis