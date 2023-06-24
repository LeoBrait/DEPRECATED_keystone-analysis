communities_path=data/community_subsets

general_fake_dir=data/fake_habitats
mkdir -p "${general_fake_dir}"

#list of achives inside the folder
subsets=($(\ls ${communities_path}))

for subset_path in "${subsets[@]}"; 
do
    #parse the file name without the .tsv
    filename=$(basename -- "${subset_path}")
    filename="${filename%.*}"
    fake_habitat_dir="${general_fake_dir}/${filename}"

        #Environment for each iteration and seed
        log="${fake_habitat_dir}/log.txt"
                     
        #create the jobs file in jobs folder
         echo "echo fake_habitat for:" \
                " ${filename} seed: ${seed} is running..."
         echo "mkdir -p ${fake_habitat_dir}"
         echo "fastspar_bootstrap" \
                "-c ${communities_path}/${subset_path} " \
                "-n 500 " \
                "-p ${fake_habitat_dir}/fake_${filename}_" \
                "-t 2 " \
                "-s 1 > ${log}"
        echo "echo fake ${filename} done!"
        echo
done > Shell/jobs/fake_habitats.txt

#Run the jobs *****************************************************************

xargs -P $parallel -I {} bash -c "{}" < Shell/jobs/fake_habitats.txt