#!/usr/bin/env python3

# setting src to the local path in order to import the other files
from sys import path
path.append('src')

from datetime import datetime
import math
import numpy as np
import pandas as pd
import gephitools as gp
from scipy import stats as sp
import os
from xonsh_py import mkdir_p, rmr, cd

def existOldData(path=None):
    """
    This function checks if the path exists
    """
    if not isinstance(path,str):
        raise Exception("Wrong path data type!")
    return os.path.exists(path)

def debug(string):
    """
    Print current time for debugging purposes
    """
    print("[%s] %s"%(datetime.now().strftime('%H:%M:%S'), string))

# a function that converts raw counts to relative ones (%)
# WARNING: only one axis sums to 100%
def relative(comm_mat):
    mat = np.zeros(shape=(comm_mat.shape))
    for i in range(comm_mat.shape[1]):
        s = sum(comm_mat[:,i])
        for j in range(comm_mat.shape[0]):
            mat[j,i] = comm_mat[j,i]/s*100.
    return mat

def main(inFile, inMeta, host, spcc_backlog):
    """
    This function is responsible for running the correlation analysis

    The step by step of the analysis is:
    1. Read the communit_matrix file
    2. Remove empty rows and columns from the community matrix
    3. Compute the co-occurrence matrices
    4. Compute the correlation matrix
    5. Compute the SparCC matrix
    6. Compute the CNM and LIASP algorithms
    7. Compute the keystone nodes
    8. Compute the Venn diagram
    9. Compute the Gephi files
    10. Compute the networkx files
    11. Compute the NGA files

    Parameters
    ----------
    inFile : string
        Path to the calculated community matrix file. The file must be a csv.
        The community matrix is a file where the rows are the samples and the columns are the taxa.

    inMeta : string
        Path to the metadata file.

    host : string
        Name of the ecosystem and habitat.

    spcc_backlog : string
        Path to the SparCC correlation matrix.
        The correlation matrix is a file where the rows and columns are the taxa.
        The values are the correlations between the taxa.

    rareAbund : float
        Minimum abundance of a node to be considered in the analysis.

    Returns
    -------
    None
    """

    print()
    debug("Reading input...")

    # input metadata file
    hasMeta = True if not "none" in inMeta else False

    # loading the community_matrix file
    # The community matrix is a file where the rows are the samples and the columns are the nodes
    # we transpose the matrix in order to have the nodes as rows and the samples as columns
    comm_matrix = pd.read_csv(inFile, skipinitialspace=True, index_col=0).T.sort_index()

    # checking for empty rows and columns before start the cnm_analysis
    deleted_indexes=[]
    for i in range(comm_matrix.shape[0]): # iterating over rows
        if sum(comm_matrix.iloc[i]) <= 0.:
            debug("WARNING: An empty row was detected! Deleting %s row..."%comm_matrix.index[i])
            deleted_indexes.append(comm_matrix.index[i])

    # deleting empty rows and resetting the list of empty indexes
    comm_matrix.drop(index=deleted_indexes,inplace=True)
    deleted_indexes=[]
    
    for i in range(comm_matrix.shape[1]): # iterating over columns
        if sum(comm_matrix.iloc[:,i]) <= 0.:
            debug("WARNING: An empty column was detected! Deleting %s column..."%comm_matrix.keys()[i])
            deleted_indexes.append(comm_matrix.keys()[i])
    comm_matrix.drop(columns=deleted_indexes,inplace=True)

    # getting labels from the file
    labelSamp = list(comm_matrix.index)

    # getting matrix' shape
    nsamp = comm_matrix.shape[0]
    nspc =  comm_matrix.shape[1]

    debug("This dataset has %d Nodes with length %d"%(nsamp,nspc))

    # renaming repeated columns
    for i in range(len(labelSamp)-1):
        k=0
        for j in range(i+1,len(labelSamp)):
            if labelSamp[i] == labelSamp[j]:
                labelSamp[j] += '_'+str(k+1)
                k += 1
                debug("There was an exchange: %s and %s were duplicated at positions %d and %d!"%(vec[i],vec[j],i,j))

    # creating a copy of the anotation data for parallelization purpose
    anotShared = comm_matrix.values # original data (without modifications)

    # ABUNDANCE and DESCRIPTIVE STATISTICS
    # abundance of counts for each node
    abund    = [np.mean(ser) for ser in anotShared]
    # log of abundance of counts for each node
    abundLog = [math.log10(ser) for ser in abund]
    # converting raw counts to relative ones (%)
    relAnot = relative(anotShared)
    # mean of relative abundance
    abundRel = [np.mean(ser) for ser in relAnot] # Taxon mean abundance
    # standard deviation of the mean for relative abundance
    stdAbund = [np.std(ser,ddof=1) for ser in relAnot]
    # standard error of the mean for relative abundance
    semAbund = [sp.sem(ser,axis=None) for ser in relAnot]
    # median of relative abundance
    abundMdn = [np.median(ser) for ser in relAnot]
    # computing prevalence for each node
    prevalence = [sum(map(lambda x: x>0., ser))/len(ser) for ser in anotShared]

    debug("Computing co-ocurrences types...")

    # CO-OCCURRENCE MATRICES
    # declaring auxiliary a, b, and c matrices for further use
    # the matrices are used to calculate the frequency of co-occurrences
    # of each node with the rest of the nodes
    # a if A and B co-occurred
    # b if A xor B occurred (A or B occurred, but not both)
    # c if both A and B did not occur
    a = np.zeros((nsamp, nsamp))
    b = np.zeros((nsamp, nsamp))
    c = np.zeros((nsamp, nsamp))

    # loop over the nodes
    for i in range(nsamp-1):
        A = anotShared[i] # number of reads for node i (row) in each sample (column)
        boolA = A > 0.
        for j in range(i+1,nsamp):
            B = anotShared[j] # number of reads for node j (row) in each sample (column)
            boolB = B > 0.
            # counts the number of samples where A and B following the pattern presented above
            a[i][j] = np.sum(boolA & boolB)
            b[i][j] = np.sum(boolA ^ boolB)
            c[i][j] = np.sum(~boolA & ~boolB)

    debug("Getting metadata file...")

    # generating the metadata file
    # TODO: the filtering method needs to filter the rows from the df based on the
    # labelSamp and also needs to add the remaining indexes that exist in the labelSamp list
    if hasMeta:
        meta = pd.read_csv(inMeta, skipinitialspace=True)
        meta.set_index(meta.keys()[0], inplace=True)
        # removing rows that didn't match the ones received in the community file
        meta = meta.filter(items=labelSamp, axis='index')
        # checking if all the indexes were deleted by the previous method
        if meta.shape[0] < 1:
            meta = pd.DataFrame(data=np.zeros(shape=(nsamp,2)), index=labelSamp, columns=['abundance', 'abundance (log10)'])
        # sorting the metadata dataframe in order to match the same nodes as in community file
        meta.sort_index(inplace=True)
    else:
        meta = pd.DataFrame(data=np.zeros(shape=(nsamp,2)), index=labelSamp, columns=['abundance', 'abundance (log10)'])

    # adding the abundance and descriptive statistics columns to the metadata
    meta['abundance']         = abund
    meta['abundance (log10)'] = abundLog
    meta['abundance (%)']     = abundRel
    meta['std'] = stdAbund
    meta['sem'] = semAbund
    meta['med'] = abundMdn
    meta['prevalence'] = prevalence

    # creating folders to output files
    rmr(['out']) # removing the folder if it already exists
    mkdir_p(['out']) # recreating the folder for the output files
    cd('out') # changing the current directory to the output folder

    # creating the output folders
    # $[mkdir -p "raw_data" "gephi_data" "cnm_data" "nga_data" "figures" "sparcc_data" "liasp_data" "nx_data" "matrices"]
    mkdir_p(["raw_data", "gephi_data", "cnm_data", "nga_data", "figures", "sparcc_data", "liasp_data", "nx_data", "matrices"])

    debug("Getting already computed SparCC matrix...")

    # reading the data with all rows and columns sorted
    coSpccPar = pd.read_csv('../'+spcc_backlog,sep='\t',index_col=0).sort_index().sort_index(axis=1).values
    # coSpccPar.drop(index=deleted_indexes,inplace=True)
    # coSpccPar.drop(columns=deleted_indexes,inplace=True)

    debug("Converting calculated data to DataFrame...")

    # inserting calculated correlations inside a DataFrame
    coSparCC = pd.DataFrame(data=coSpccPar, index=labelSamp, columns=labelSamp)

    debug("Saving non-filtered correlation data...")

    # printing non filtered correlation matrix to file
    coSparCC.to_csv("raw_data/correlation_matrix.csv")

    # printing non filtered correlation matrix to file as gephi format
    gp.printCorr("gephi_data/correlation_matrix.csv", np.array(coSpccPar), a, b, c)

    # printing nodes' labels using gephi format
    gp.printNodes("gephi_data/nodes.csv", meta)

    debug("Executing Critical Network Method (CNM) and Largest Influence on Average Shortest Path (LIASP) algorithms...")

    cd('..')
    debug("Exiting...")

    return coSparCC

if __name__ == '__main__':
    assert len(sys.argv) > 3, "Missing arguments!\nUsage: ./main <community> <metadata> <host> <spcc_corrpath>\n"
    _ = main(sys.argv[1],sys.argv[2],sys.argv[3])
