#!/bin/bash

# Environment ******************************************************************
communities_path=data/community_subsets
tablenames_real=($(\ls ${communities_path}))
general_synthetics_dir=data/synthetic_habitats
mkdir -p "${general_synthetics_dir}"

# Job creation *****************************************************************
for single_tablename_real in "${tablenames_real[@]}"; 
do
    #parse the file name without the .tsv
    habitat=$(basename -- "${single_tablename_real}")
    habitat="${habitat%.*}"

    synt_habitat_dir="${general_synthetics_dir}/${habitat}"

    #Environment for each iteration and seed
    log="${synt_habitat_dir}/log.txt"
                     
    #create the jobs file in jobs folder
    echo "echo fake_habitat for: ${habitat} is running..."
    echo "mkdir -p ${synt_habitat_dir}"
    echo "fastspar_bootstrap" \
            "-c ${communities_path}/${single_tablename_real}" \
            "-n 250" \
            "-p ${synt_habitat_dir}/fake_${habitat}" \
            "-t 2 " \
            "-s 1 > ${log}"
    echo "echo fake ${habitat} done!"
    echo
done > Shell/jobs/synthetic_habitats.txt

#Run the jobs *****************************************************************

xargs -P $parallel -I {} bash -c "{}" < Shell/jobs/synthetic_habitats.txt