

# Test real correlations vs. fake *********************************************
fake_communities_path=data/fake_habitats
fake_habitats=($(ls ${fake_communities_path}))
seeds=(1 2)
mkdir -p data/fake_fastspar

for fake_habitat in "${fake_habitats[@]}";
do
    #parse the file name without the .tsv
    filename=$(basename -- "${fake_habitat}")
    filename="${filename%.*}"
    fake_fastspar_dir="data/fake_fastspar/${fake_habitat}"


    for seed in "${seeds[@]}"; 
    do
        #Environment for each iteration and seed
        log="${fake_fastspar_dir}/log_${seed}.txt"
        fake_cor="${fake_fastspar_dir}/cor_${seed}"
        fake_cov="${fake_fastspar_dir}/cov_${seed}"
        remove=15
        if [ -f "${fake_cor}" ]; then
                echo "echo The Sparcc for: " \
                        "${filename} with  iterations, " \
                            "seed ${seed} was already done" > general_log.txt
            else
                for table in {0..500}
                do             
                #create the jobs file in jobs folder
                echo "echo The Sparcc for:" \
                        " ${filename}," \
                            "seed ${seed} is running..."
                echo "mkdir -p ${fake_fastspar_dir}"
                echo "fastspar "\
                    "-c ${fake_communities_path}/${filename}__$table"\
                    "-r ${fake_cor}_table "\
                    "-a ${fake_cov}_table "\
                    "-t 2 "\
                    "-s $seed "\
                    "-i 16000 "\
                    "-x $remove "\
                    "-e 0.1 "\
                    "-y > ${log}"
                echo "echo table ${table} of ${filename}" \
                        " and seed ${seed} done!"
                echo
                done         
            fi
    done
done > Shell/jobs/fake_fastspar.txt

# Run the jobs *****************************************************************
xargs -P $parallel -I {} bash -c "{}" < Shell/jobs/fake_fastspar.txt