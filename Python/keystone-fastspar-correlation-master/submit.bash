#!/bin/bash

# number of fastspar executions per submission file
# each done by each input matrix in community_matrix
nsubs=80

#for level in class family phyla order genus; do
for level in phyla;
	do echo "${level}"

    # getting reverse size-ordered list of input files
	list=( ls -Sr community_matrix/*csv)

	n=0
	for w in $(seq 0 $nsubs ${#list[@]});
	do

		sublist=${list[@]:$w:$nsubs}

		j=$(printf "%03d" ${n})

        # change the submission file header below as you wish
        # this is a PBS header used by "qsub", you can edit it
        # for SLURM requirements for instance
	echo "#!/bin/bash
#PBS -N SparCC_${j}
#PBS -V
#PBS -l nodes=1:ppn=48
#PBS -q exec_8_600

cd $(pwd)

" > jobs/subs_${level}_${j}.sh

		l=$w
		for m in ${sublist[@]}; do
            # comment the line below if you want sequencial execution
			# echo "./run_fastspar.sh ${m} &> outs/${level}_$(printf "%03d" $l).out &" >> jobs/subs_${level}_${j}.sh

            # comment the line below if you want "parallel" execution
			echo "./run_fastspar.sh ${m} &> outs/${level}_$(printf "%03d" $l).out" >> jobs/subs_${level}_${j}.sh
			let "l++"
		done

		echo "wait" >> jobs/subs_${level}_${j}.sh

		let "n++"
	done

done
