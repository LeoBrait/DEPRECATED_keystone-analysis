# Please refer to the [WIKI](https://github.com/MeirellesLab/keystone-analysis/wiki) for details

# Keystone Analysis

This repository holds all codes to reproduce the major part of the Keystones determination procedures described in the Methods section of the paper.


## Citation

***

## Description

The entire analysis aims to determine which are the Keystones nodes based on multiple metrics for a network constructed from a correlation matrix. It compares many metrics used on the paper and rank the Keystones through a continuous range of Keystoneness. It also employs the LIASP metric, proposed [here](~Paper de Roberto.~).


## Overview

![](https://github.com/MeirellesLab/keystone-webcontent/blob/master/major-fluxogram.svg)

__Check the sidebar links for detailed information about this flowchart.__

## Requirements

1. Data preparation

The program inputs are basically two:
    - A occurence matrix, which the first column reffer to the samples(i.e, communities) and the other columns reffer to the OTUs. The values are the abundances. The ocurrence matrix can be found [here](data/taxon_abundances). Extract the given files just right there or use your own data.
    - A metadata matrix, which the first column reffer to the samples(i.e, communities) and the other columns reffer to the metadata. In this pipeline we use habitat and ecosystem as metadata. The metadata matrix can be found [here](data/metadata). Extract the given files just right there or use your own data.

2. Software dependencies

This program relies on a Python-3 interpreter and the anaconda package manager. Also, it was tested and developed on Ubuntu-18-04(WSL) with the Anaconda3 environment manager but probably can be run on other Linux distributions and MacOS.
All python packages required by the program are listed in `requirements.txt`. In order to create and activate environment with all dependencies, you should install [anaconda 3](https://www.anaconda.com/) on your linux machine and run:

```bash
bash install_dependencies.sh
conda activate keystones
```
## Usage

After installing all requirements and setting your files in the right directories, you can run the analysis just by doing:

```sh
Python3 Python/main.py
```

This script will index all input matrices you gave, execute fastspar with bootstraping, check for the matrix of edges. If the matrix of edges is found, the script proceeds to execute [`correlation.py`]((A)-Keystone-Analysis-Program) for each input matrix.


## Contact

You may contact the corresponding author Pedro for questions related to the paper itself through the email pedrommeirelles@gmail.com. 

For questions related to this repository and its auxiliary ones, you may also contact the developers Bertolino (jgabbc@hotmail.com), Flávia Mayumi (flaviamayumi.rh@gmail.com), .

## Contributors

Original develepers of ths program:
- Bertolino - @bertolinocastro 
- Flávia Ruziska - @flaviamayumi
- Rafael Menezes - @r-menezes

Later developers:
- Leonardo Brait
- Felipe Alexandre

## License

~Not defined yet~

