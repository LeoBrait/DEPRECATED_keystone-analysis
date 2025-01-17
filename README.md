# Please refer to the [WIKI](https://github.com/MeirellesLab/keystone-analysis/wiki) for details

# Keystone Analysis

This repository holds all codes to reproduce the major part of the Keystones determination procedures described in the Methods section of the paper.
##

## Citation

***

## Description

The entire analysis aims to determine which are the Keystones nodes based on multiple metrics for a network constructed from a correlation matrix. It compares many metrics used on the paper and rank the Keystones through a continuous range of Keystoneness. It also employs the LIASP metric, proposed [here](~Paper de Roberto.~).

## Requirements

1. Data preparation for keystone analysis

The program inputs are basically two:
- A occurence matrix, which the first column reffers to the samples(i.e, communities) and the other columns reffer to the OTUs. The values are the abundances. The ocurrence matrix can be found [here](data/taxon_abundances). Extract the given files just right there or use your own data.  
- A metadata matrix, which the first column reffer to the samples(i.e, communities) and the other columns reffer to the metadata. In this pipeline we use habitat and ecosystem as metadata. The metadata matrix can be found [here](data/metadata). Extract the given files just right there or use your own data.

To extract the referenced files, once in the keystone-analysis directory you can run:

```sh
unzip data/taxon_abundances/kraken_biomedb_absolute.zip -d data/taxon_abundances/
unzip data/taxon_abundances/kraken_biomedb_relative.zip -d data/taxon_abundances/
```

2. Software dependencies

This program relies on a Python-3 interpreter and the anaconda (or miniconda) package manager. Also, it was tested and developed on Ubuntu-18-04(WSL) with the Anaconda3 environment manager but probably can be run on other Linux distributions and MacOS.  
All python packages required by the program are listed in [requirements file](requirements/requirements_pyshell.txt). In order to create and activate environment with all dependencies, you first need to install [anaconda 3](https://www.anaconda.com/) or miniconda on your linux machine, then run:

```sh
bash Shell/installation/install_dependencies.sh
```
* This is necessay just for the first time you run the program.
* **Troubleshooting:** If you have problems to find your conda installation and profiles, you can try the following:

```sh
which conda
```

3. (Optional) Data preparation for posprocessing

If you want to proceed to the postprocess and generate the figures of the paper, you will need to extract the other archives:

```sh
unzip data/taxon_abundances/kraken_stdb_absolute.zip -d data/taxon_abundances/
unzip data/taxon_abundances/kraken_stdb_relative.zip -d data/taxon_abundances/
```


## Usage

1. Main program

After installing all requirements and setting your files in the right directories, you should check for the environment settings in the file [settings.sh](Shell/settings.sh). In this program we use environment switchs between R and Python environments, so please make sure your package manager and profiles path are correctly assigned. After that, you can run the analysis just by doing:

```sh
nohup bash Shell/main.sh
```

- Please note that the main program calls an [python script](Python/pipelines/preprocessing_data.py) to prepare the data for the analysis. The script have a chunk "#remove undesired taxa". This chunk is setted according to the feedbacks from the "Fastspar for all communities" in the [fastspar's pipeline](Shell/pipelines/calculating_fastspar.sh). Plese, analyse the nohup file and ajust it according to your data.

2. Posprocess Material

(not ready yet)
Some other features are very data sensitive, and could break the main analyses. This is the case of most of our data vizualisation/posprocessing resources. To maintain the program safe and easy, we decided to keep them in modularized scripts. For customized analyses, you will need to adapt the scripts to your data. But if you only need to reprode the results of the [paper](), you can run the following protocol:

First, you must install "Arial" font in your system to run the scripts.

```sh
sudo apt install ttf-mscorefonts-installer
sudo fc-cache -f
```

This next step depends on a file named "radiations.csv", the file serves to analyse the microbial groups by their radiation(CPR, DPANN or Bonafide).
The file must be in  "data/<frame_analysis>/radiations.csv". If you only want to reproduce the results of the paper, you can just run:

```sh
bash Shell/pipelines/posprocessing.sh 
```

## Outputs

This pipeline is results sensitive. It means that most of its processes can be paused then continued at any time. The pipeline also accepts a custom name for the analysis frame, that can be setted in settings script, which will allow analysis of different datasets without overwriting previous results.
The following directories will be generated in the data folder after running the code:

- **data/<analyses_frame>/summaries/**  
    Contains the number of samples for each habitat and ecosystem.

- **data/<analyses_frame>/performance_fastspar_iterations/** 
    Contains sparcc correlatiom matrices, covariance matrices and logs with different seeds and iteractions.

- **data/<analyses_frame>/fastspar_correlations/**
    Contains the final sparcc correlation matrices, covariance matrices and logs.

- **data/<analyses_frame>/synthetic_habitats/**  
    Contains tables for synthetic habitats generated by the fastspar_bootstrap function.

- **data/<analyses_frame>/synthetic_fastspar**
    The Sparcc's cor, cov and log files for each synthetic habitat.

If you want to regenerate any sttep of your results, you can simply delete one or all of the directories above and run the code again.


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

## Troubleshooting
Some R pckges can be problematic to install and the Anaconda environment manager sometimes does not perceive it. If you have any problems with the Anaconda envrironment solving, try to run the R-libraries troubleshooting trhough:

```sh
bash nohup Shell/pipelines/r_troubleshooting.sh
```

This will try to install all packages needed in a dummy folder outside the conda enviroment and produces a report of errors: "nohup.out". You should read the report and try to solve the issues of your machine manually, installing the missing modules into your machine. This could need several runs. After this procedure, you will need to reinstall the R conda environment.

## Contact

You may contact the corresponding author Pedro Meirelles for questions related to the paper itself through the email pedrommeirelles@gmail.com. 

For questions related to this repository and its auxiliary ones, you may also contact the developers:
Felipe Alexandre (felipe.as.barbosa99@gmail.com)
Leonardo Brait (leonardobrait@gmail.com )
Rafael Menezes (menezes.santos.rafael@gmail.com )
Gabriel Bertolino (jgabbc@hotmail.com)
Flávia Mayumi (flaviamayumi.rh@gmail.com), .

## Contributors

- Bertolino - @bertolinocastro 
- Flávia Ruziska - @flaviamayumi
- Rafael Menezes - @r-menezes
- Leonardo Brait - 
- Felipe Alexandre - @Felipe-Alexandre

## License

~Not defined yet~

