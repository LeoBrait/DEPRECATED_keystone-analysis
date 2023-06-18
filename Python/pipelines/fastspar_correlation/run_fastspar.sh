#!/bin/bashls
[[ -z $1 ]] && echo 'Missing input file!' && exit 1

nthreads=4 # WARNING: It does contribute to randomness!
# instead of a single seed, I'm using a random variable for all used fastspar steps

iter=3000
xter=50
corThr=0.1
seed=1


input=${1##*/}


base=output/transposed/$input/sparcc/sparcc_data
anot=anot.tsv

# creating final folder
rm -rf $base && mkdir -p $base

# creating input file for sparcc
ls
echo -e 'Converting annotation file to fastspar input format\n'
python3 csvtool_pandas.py community_matrix/${input} ${base}/${anot}

# computing the original data
fastspar -c ${base}/$anot -r ${base}/cor.tsv -a ${base}/cov.tsv -t $nthreads -s $seed -i $iter -x $xter -e $corThr -y
