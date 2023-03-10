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
df = pd.read_csv('phyla-relative-ab-keystones-validation.csv')

# Split table by habitats in column "Level 3"
df_split = [g for _, g in df.groupby('habitat')]

#Saving new dataframes
for i, split in enumerate(df_split):
    df_split[i].to_csv(f"{df_split[i]['Ecosystem'].iloc[0]}." + f"{df_split[i]['habitat'].iloc[0]}.csv", index=False)
