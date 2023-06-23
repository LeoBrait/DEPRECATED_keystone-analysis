# Please refer to the [WIKI](https://github.com/MeirellesLab/keystone-analysis/wiki) for details

# Keystone Analysis

This repository holds all codes to reproduce the major part of the Keystones determination procedures described in the Methods section of the paper.


## Citation

***

## Description

The entire analysis aims to determine which are the Keystones nodes based on multiple metrics for a network constructed from a correlation matrix. It compares many metrics used on the paper and rank the Keystones through a continuous range of Keystoneness. It also employs the LIASP metric, proposed [here](~Paper de Roberto.~).

## Requirements

1. Data preparation

The program inputs are basically two:
- A occurence matrix, which the first column reffer to the samples(i.e, communities) and the other columns reffer to the OTUs. The values are the abundances. The ocurrence matrix can be found [here](data/taxon_abundances). Extract the given files just right there or use your own data.  
- A metadata matrix, which the first column reffer to the samples(i.e, communities) and the other columns reffer to the metadata. In this pipeline we use habitat and ecosystem as metadata. The metadata matrix can be found [here](data/metadata). Extract the given files just right there or use your own data.

Once in the keystone-analysis directory you can run

```sh
unzip data/taxon_abundances/kraken_relative_matrix_biome_db.zip
```

2. Software dependencies

This program relies on a Python-3 interpreter and the anaconda (or miniconda) package manager. Also, it was tested and developed on Ubuntu-18-04(WSL) with the Anaconda3 environment manager but probably can be run on other Linux distributions and MacOS.  
All python packages required by the program are listed in `requirements.txt`. In order to create and activate environment with all dependencies, you should install [anaconda 3](https://www.anaconda.com/) or miniconda on your linux machine and run:

```sh
bash Shell/install_dependencies.sh
```
* This is necessay just for the first time you run the program.
* **Troubleshooting:** If you have problems to find your conda installation and profiles, you can try the following:

```sh
which conda
```



## Usage

After installing all requirements and setting your files in the right directories, you should check for the environment settings in the file [main.sh](Shell/main.sh). In this program we use environment switchs between R and Python environments, so please make sure your package manager and profiles path are correctly assigned. After that, you can run the analysis just by doing:

```sh
bash Shell/main.sh
```

This script will index all input matrices you gave, execute fastspar with bootstraping, check for the matrix of edges. If the matrix of edges is found, the script proceeds to execute [`correlation.py`]((A)-Keystone-Analysis-Program) for each input matrix.

## Output

# Description of Variables in the keystones.csv file

| Variable            | Description                                                                                                                                                                                                                               |
| ------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Taxon               | Tthe taxon identification                                                                                                                                                                                                                 |
| Ecosystem           | Ecosystem associated with Habitat                                                                                                                                                                                                         |
| Habitat             | Habitat at which samples were taken                                                                                                                                                                                                       |
| Abundance           | Mean relative number of reads of taxon in samples from Habitat, in percentage. Column should add to in each Habitat.                                                                                                                      |
| AbundanceAbsolute   | Mean relative number of reads of taxon in samples from Habitat.                                                                                                                                                                           |
| Marker              | Marker symbol. Used for diagnostic plots.                                                                                                                                                                                                 |
| Median              | Median number of reads for taxon in samples from Habitat                                                                                                                                                                                  |
| Prevalence          | Proportion of samples from Habitat in which taxon was detected                                                                                                                                                                            |
| Sem                 | Standard error of the number of reads of taxon in samples from Habitat                                                                                                                                                                    |
| Std                 | Standard deviation of the number of reads of taxon in samples from Habitat                                                                                                                                                                |
| BC                  | Betweeness Centrality. Defined as the number of shortest paths between all pairs of nodes in the Habitat network that pass through the node representing the Taxon.                                                                       |
| D                   | Degree. Defined as the number of edges connected to the node representing the Taxon in the Habitat network.                                                                                                                               |
| DxCC                | Product of Degree and Closeness Centrality. Closeness centrality is defined as the inverse of the sum of the shortest path distances from the node representing the Taxon to all other nodes in the Habitat network.                      |
| LIASP               | Largest Impact on Average Shortest Path. Defined in the main text.                                                                                                                                                                        |
| LIASPdir            | Direct component of LIASP. Defined in the main text.                                                                                                                                                                                      |
| LIASPindir          | Indirect component of LIASP. Defined in the main text.                                                                                                                                                                                    |
| effOriginal         | Efficiency of the network representing the Habitat before the Taxon is removed. Defined in the main text.                                                                                                                                 |
| effRemoved          | Efficiency of the network representing the Habitat after the Taxon is removed. Defined in the main text.                                                                                                                                  |
| totalEffChange      | Total change in the efficiency the network representing the Habitat when the Taxon is removed. Defined in the main text.                                                                                                                  |
| dirEffChange        | Direct component of the change in the efficiency the network representing the Habitat when the Taxon is removed. Defined in the main text.                                                                                                |
| indirEffChange      | Indirect component of the change in the efficiency the network representing the Habitat when the Taxon is removed. Defined in the main text.                                                                                              |
| VARIABLE_isKeystone | Boolean (1/0) variable indicating if the value of VARIABLE for this taxon is above the treshold given by MEDIAN(VARIABLE) + 2*STANDARD_DEVIATIONS(VARIABLE) [1] or not [0], when considering the distribution of VARIABLE for the Habitat |
| VARIABLE_scaled     | The scaled value of the variable, with the largest occurrence in the Habitat assigned to 1 and remaining occurrences scaled as VALUE/MAX(VARIABLE)                                                                                        |
| VARIABLE_rank       | The rank of the Taxon when all Taxa in the Habitat are ordered by increasing values of the variable, with the smallest value observed in the Habitat assigned to 1                                                                        |
| VARIABLE_zScore     | The z-score of the variable, scaled as (VALUE - MEAN(VARIABLE))/STANDARD_DEVIATION(VARIABLE)                                                                                                                                              |
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

