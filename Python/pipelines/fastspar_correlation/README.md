# Scripts for SparCC correlation used by [Keystones Analysis](https://github.com/MeirellesLab/keystone-analysis)

Here you find the scripts executed to generate the SparCC correlation matrices used for the [Keystones Analysis](https://github.com/MeirellesLab/keystone-analysis) published on the [paper]() ~not yet published~.

## Citation

***

### Keystone Analysis Repositories
- [Main code](https://github.com/MeirellesLab/keystone-analysis)
- [Input Community Matrices](https://github.com/MeirellesLab/keystone-analysis-community-input)
- [Input SparCC Correlation Matrices](https://github.com/MeirellesLab/keystone-fastspar-correlation)
- [Development repo](https://bitbucket.org/bertolinocastro/the_model/)

### Requisites

- You must download the [input matrices](https://github.com/MeirellesLab/keystone-analysis-community-input) and copy them to the [`community_matrix/` folder](community_matrix/). Feel free to use your own matrices. Just make sure they fit the FastSpar input requirements.
- You must download the [FastSpar](https://github.com/scwatts/fastspar) binaries and copy them to the root of this repository. [This commit version](https://github.com/scwatts/fastspar/tree/d0c83b324645722f1ea57689aca305b2b4b204c1) was used on the paper.

### Execution

1. Copy all your inputs to the `community_matrix` folder
2. Make the scripts executable: `chmod +x submit.bash run_fastspar.sh`
3. Execute `./submit.bash`
4. Enter in `jobs` folder and execute each submission script. You can just type `for i in *; do qsub ${i}; done`
    1. You can check the execution status in the respective `logfile` inside `outs` folder
5. After finished, copy the `output` folder to the Keystone Analysis [root directory](https://github.com/MeirellesLab/keystone-analysis)

### Scripts description

#### `submit.bash`

This script creates a submission script with a PBS header to be submitted and executed by a cluster, although it can be run in any environment with `Bash` cli/interpreter.

This script is suitable for the data used in the Paper (found in [this repo](https://github.com/MeirellesLab/keystone-analysis-community-input)).

#### `run_fastspar.sh`

This script reshapes the matrices found in `community_matrix` to the FastSpar requirements and then execute the FastSpar computation.
It also maintains the same folder structure used by the [main repo](https://github.com/MeirellesLab/keystone-analysis) so that you just need to copy the `output` folder.

###### You can add `Boostrapping` to the `run_fastspar.sh` script as described in the [FastSpar repository](https://github.com/scwatts/fastspar), but be aware that it was not done in the Paper and that the Keystones Analysis is very sensitive to the  `Boostrapping` randomness.

#### `csvtool_pandas.py`

This script just reshapes the input matrix to the FastSpar requirements.

## Contact

You may contact the corresponding author Pedro for questions related to the paper itself through the email pedrommeirelles@gmail.com. 

For questions related to this repository and its auxiliary ones, you may also contact the developers Bertolino (jgabbc@hotmail.com) and Flávia Mayumi (flaviamayumi.rh@gmail.com).

## Contributors

This program was developed by:
- Bertolino - @bertolinocastro

## License

~Not defined yet~
