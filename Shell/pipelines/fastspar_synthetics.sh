
# Environment ******************************************************************
#!/bin/bash
general_synthetics_dir=data/synthetic_habitats
synt_habitats_dirs=($(ls ${general_synthetics_dir}))
mkdir -p data/synthetic_fastspar

for synt_habitat_dir in "${synt_habitats_dirs[@]}";
do
    #parse the file name without the .tsv
    habitat=$(basename -- "${synt_habitat_dir}")
    habitat="${habitat%.*}"
    mkdir -p "data/synthetic_fastspar/${habitat}"
    
    tables=($(ls ${general_synthetics_dir}/${synt_habitat_dir}))
    for table in "${tables[@]}";
    do
        #Environment
        echo $table
        synt_fastspar_dir="data/synthetic_fastspar/${habitat}"
        remove=15
        mkdir -p "${synt_fastspar_dir}"
        synt_cor="${synt_fastspar_dir}/cor_${table}.cor"
        synt_cov="${synt_fastspar_dir}/cov_${table}.cov"
        log="${synt_fastspar_dir}/log_${table}.txt"

    #if [ -f "${synt_cor}" ]; then
        # echo "echo The Sparcc for: " \
        #         "${habitat} with  iterations, " \
        #             "table ${table} was already done" > general_log.txt
    #else         
        #create the jobs file in jobs folder
        # echo "echo The Sparcc for:" \
        #         " ${habitat}," \
        #             "table ${table} is running..."
        # echo "fastspar "\
        #             "-c ${general_synthetics_dir}/${synt_habitat_dir}/"\
        #             "-r ${fake_cor}_table "\
        #             "-a ${fake_cov}_table "\
        #             "-t 2 "\
        #             "-s $seed "\
        #             "-i 16000 "\
        #             "-x $remove "\
        #             "-e 0.1 "\
        #             "-y > ${log}"
        #         echo "echo table ${table} of ${filename}" \
        #                 " and seed ${seed} done!"
        #         echo       
        #fi
    done
done #> Shell/jobs/fake_fastspar.txt

# Run the jobs *****************************************************************
#xargs -P $parallel -I {} bash -c "{}" < Shell/jobs/fake_fastspar.txt