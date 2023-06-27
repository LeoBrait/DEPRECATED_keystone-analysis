communities_path=data/community_subsets
tablenames=($(\ls ${communities_path}))
fastspar_dir=data/fastspar_correlations
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