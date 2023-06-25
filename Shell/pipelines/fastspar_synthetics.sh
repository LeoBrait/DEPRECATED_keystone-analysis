
# Environment ******************************************************************
#!/bin/bash

# input
general_synthetics_dir=data/synthetic_habitats
synt_habitats_dirs=($(ls ${general_synthetics_dir}))

# output
mkdir -p data/synthetic_fastspar

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
        synt_fastspar_dir="data/synthetic_fastspar/${habitat}"
        remove=15
        synt_cor="${synt_fastspar_dir}/${habitat}_${table_number}.cor"
        synt_cov="${synt_fastspar_dir}/${habitat}_${table_number}.cov"
        log="${synt_fastspar_dir}/log_${habitat}"

    if [ -f "${synt_cor}" ]; then
        echo "echo The Sparcc for: " \
                "${habitat} " \
                   "table ${table_number} was already done" > general_log.txt
    else         
        #create the jobs file in jobs folder
         echo "echo The Sparcc for:" \
                " synthetic ${habitat}," \
                    "table ${table_number} is running..."
        echo "mkdir -p ${synt_fastspar_dir}"
        echo "fastspar "\
                     "-c ${general_synthetics_dir}/${synt_habitat_dir}/${table}"\
                     "-r ${synt_cor}"\
                     "-a ${synt_cov}"\
                     "-t 2 "\
                     "-s 1 "\
                     "-i 4000 "\
                     "-x $remove "\
                     "-e 0.1 "\
                     "-y > ${log}_${table_number}.txt"
                 echo "echo table ${table_number} of ${filename}" \
                         " done!"
                 echo       
        fi
    done
done > Shell/jobs/fake_fastspar.txt



# Run the jobs *****************************************************************
xargs -P $parallel -I {} bash -c "{}" < Shell/jobs/fake_fastspar.txt