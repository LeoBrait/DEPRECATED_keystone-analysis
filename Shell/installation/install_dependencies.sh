#!/bin/bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda create --name pyshell_biome_keystones --file requirements/requirements_pyshell.txt


