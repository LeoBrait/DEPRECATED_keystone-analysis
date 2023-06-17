#!/bin/bash
[[ -z $1 ]] && echo 'Missing input file!' && exit 1

nthreads=5

xter=10
corThr=0.1

input=${1##*/}

base=output/transposed/$input/sparcc/sparcc_data
anot=anot.tsv

# creating final folder
rm -rf $base && mkdir -p $base

# creating input file for sparcc
echo -e 'Converting annotation file to fastspar input format\n'
python3 csvtool_pandas.py community_matrix/${input} ${base}/${anot}

# real cor data
for iter in 125 250 500 750 1500 3000 4500; do
    path=${base}/iter_${iter}
    mkdir -p $path
    for seed in {0..20}; do
        sseed=$(printf "%03d" $seed)
        if [[ ! -f ${path}/cor_${sseed}.tsv ]]; then
            fastspar -c ${base}/$anot -r ${path}/cor_${sseed}.tsv -a ${path}/cov_${sseed}.tsv -t $nthreads -s $seed -i $iter -x $xter -e $corThr -y
        fi
    done
done

# bootstrapping
iter=3000
path=${base}/boot_iter_${iter}

for seed in {1..2}; do
    sseed=$(printf "%03d" $seed)
    mkdir -p ${path}/${sseed}/anot ${path}/${sseed}/cor ${path}/${sseed}/cov
    fastspar_bootstrap -c ${base}/$anot -n 1000 -p ${path}/${sseed}/anot/a -t $nthreads -s $seed
done

# not random seeded fake correlation
for seed in {1..2}; do
    sseed=$(printf "%03d" $seed)
    for i in {0..999}; do
        if [[ ! -f ${path}/${sseed}/cor/${i}.tsv ]]; then
            fastspar -c ${path}/${sseed}/anot/a_${i}.tsv -r ${path}/${sseed}/cor/${i}.tsv -a ${path}/${sseed}/cov/${i}.tsv -t $nthreads -s 1 -i $iter -x $xter -e $corThr -y
        fi
    done
done

# getting p-value
for seed in {1..2}; do
    sseed=$(printf "%03d" $seed)
    # creating a single real computed correlation
    fastspar -c ${base}/$anot -r ${path}/${sseed}/cor.tsv -a ${path}/${sseed}/cov.tsv -t $nthreads -s 1 -i $iter -x $xter -e $corThr -y
    fastspar_pvalues -c ${base}/$anot -r ${path}/${sseed}/cor.tsv -n 1000 -p ${path}/${sseed}/cor/ -o ${path}/${sseed}/pval.tsv -t $nthreads
done