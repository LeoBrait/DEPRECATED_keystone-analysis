# Please refer to the [WIKI](https://github.com/MeirellesLab/keystone-analysis/wiki) for details

# Keystone Analysis

This repository holds all codes to reproduce the major part of the Keystones determination procedures described in the Methods section of the paper.

This repository generates only the published results and a few more needed to make them by other scripts.

## Citation

***

## Description

The entire analysis aims to determine which are the Keystones nodes based on multiple metrics for a network constructed from a correlation matrix. It compares many metrics used on the paper and rank the Keystones through a continuous range of Keystoneness. It also employs the LIASP metric, proposed [here](~Paper de Roberto.~).

### Keystone Analysis Repositories
- [Main code](https://github.com/MeirellesLab/keystone-analysis)
- [Input Community Matrices](https://github.com/MeirellesLab/keystone-analysis-community-input)
- [Input SparCC Correlation Matrices](https://github.com/MeirellesLab/keystone-fastspar-correlation)
- [Development repo](https://bitbucket.org/bertolinocastro/the_model/)

## Overview

![](https://github.com/MeirellesLab/keystone-webcontent/blob/master/major-fluxogram.svg)

__Check the sidebar links for detailed information about this flowchart.__

## Requirements

This program relies on an updated Python-3 interpreter and a Fortran compiler. It was tested and developed on Ubuntu-18-04(WSL) with the Anaconda3 environment manager.
All python packages required by the program are listed in `requirements.txt`. In order to create and activate environment with all dependencies, you should run:


```bash
bash install_dependencies.sh
conda activate keystones
```

## Data Preparation

The first step is to prepare the data to be analyzed. Execute the main.R file to generate the data to be analyzed. From the metadata and tax occurence files, it will generate the files that will be used in the analysis. The files will be generated in the folder `data/`. The files generated are:

You should have your input matrices copied into `community_matrix`. Their full description is deeply discussed [here]((1)-Community-Matrices). You may use the ones from the paper that you can find [here](https://github.com/MeirellesLab/keystone-analysis-community-input). Just download the repository and copy the files.

You may also have metadata matrices about the nodes for your input matrix. Copy them inside `metadata` folder. They are well described [here]((3)-Metadata-Matrices).

Finally, you should copy the correlation matrices computed by your desired criteria to be the edges of the network. You can create these edges by whatever mean you want. For the paper, we created them as a SparCC correlation matrix. You should read [this repository](https://github.com/MeirellesLab/keystone-fastspar-correlation) to get the same data as we did. These files are well described [here]((2)-SparCC-Correlation-Matrices).

## Usage

After installing all requirements and setting your files in the right directories, you can run the analysis just by doing:

```sh
Python3 Python/main.py
```

This script will index all input matrices you copied, check if they have metadata based on their filename and check for the matrix of edges (SparCC correlation matrix in our case). If the matrix of edges is found, the script proceeds to execute [`correlation.py`]((A)-Keystone-Analysis-Program) for each input matrix. After all, the script will call [`run_keystone_analysis.py` script](https://github.com/MeirellesLab/keystone-analysis/blob/master/run_keystone_analysis.py), which runs the scripts that need the outputs from all input matrices.

I suggest you to read the [`run_all.py`](https://github.com/MeirellesLab/keystone-analysis/blob/master/run_all.py) script and modify it based on your need.

## Contact

You may contact the corresponding author Pedro for questions related to the paper itself through the email pedrommeirelles@gmail.com. 

For questions related to this repository and its auxiliary ones, you may also contact the developers Bertolino (jgabbc@hotmail.com) and Flávia Mayumi (flaviamayumi.rh@gmail.com).

## Contributors

This program was developed by:
- Bertolino - @bertolinocastro
- Flávia Ruziska - @flaviamayumi
- Rafael Menezes - @r-menezes
- Leonardo Brait
- Felipe Alexandre

## License

~Not defined yet~

