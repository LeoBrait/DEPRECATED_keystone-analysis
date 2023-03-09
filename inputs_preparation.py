#!/usr/bin/python3.7

#Split input tables

# If you have a matrix with all annotated samples results and want to perform keystones analysis
# run this code to generate the proper inputs to keystones analysis scripts

#This script will split your matrix in a new matrix for each habitat

import pandas as pd
import os

path = os. getcwd()

os.chdir(path)

# Read data
df = pd.read_csv('phyla.general.relative.matrix_2022_08_19.csv')

# Split table by habitats in column "Level 3"
df_split = [g for _, g in df.groupby('Level_3')]

#Saving new dataframes
for i, split in enumerate(df_split):
    df_split[i].to_csv(f"{df_split[i]['Level_2'].iloc[0]}." + f"{df_split[i]['Level_3'].iloc[0]}.csv", index=False)
