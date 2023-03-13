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
import integration as itg
from sys import exit
import os
from xonsh_py import *

def existOldData(path=None):
    if not isinstance(path,str):
        raise Exception("Wrong path data type!")
    return os.path.exists(path)

def debug(string):
    print("[%s] %s"%(datetime.now().strftime('%H:%M:%S'), string))

def main(inFile, inMeta, host, spcc_backlog, rareAbund=0.1):
    print()
    debug("Reading input...")

    # input metadata file
    hasMeta = True if not "none" in inMeta else False

    # getting the data from file
    anot = pd.read_csv(inFile, skipinitialspace=True, index_col=0).T.sort_index()

    # checking for empty rows and columns before start the cnm_analysis
    delete=[]
    for i in range(anot.shape[0]): # iterating over rows
        if sum(anot.iloc[i]) <= 0.:
            debug("WARNING: An empty row was detected! Deleting %s row..."%anot.index[i])
            delete.append(anot.index[i])

    # deleting empty rows and resetting the list of empty indexes
    anot.drop(index=delete,inplace=True)
    delete=[]
    
    for i in range(anot.shape[1]): # iterating over columns
        if sum(anot.iloc[:,i]) <= 0.:
            debug("WARNING: An empty column was detected! Deleting %s column..."%anot.keys()[i])
            delete.append(anot.keys()[i])
    anot.drop(columns=delete,inplace=True)

    # getting labels from the file
    labelSamp = list(anot.index)
    labelSpc  = list(anot.keys())

    # getting matrix' shape
    nsamp = anot.shape[0]
    nspc =  anot.shape[1]

    debug("This dataset has %d Nodes with length %d"%(nsamp,nspc))

    # renaming repeated columns
    for i in range(len(labelSamp)-1):
        k=0
        for j in range(i+1,len(labelSamp)):
            if labelSamp[i] == labelSamp[j]:
                labelSamp[j] += '_'+str(k+1)
                k += 1
                debug("There was an exchange: %s and %s were duplicated at positions %d and %d!"%(vec[i],vec[j],i,j))

    # a function that converts raw counts to relative ones (%)
    # WARNING: It only converts TRANSPOSED raw matrices to relative!!!
    # *******  That's because only one axis sums to 100%
    def relative(anota):
        mat = np.zeros(shape=(anota.shape))
        for i in range(anota.shape[1]):
            s = sum(anota[:,i])
            for j in range(anota.shape[0]):
                mat[j,i] = anota[j,i]/s*100.
        return mat

    # creating a copy of the anotation data for parallelization purpose
    anotShared = anot.values # original data (without modifications)

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

    # if wanted, removing the nodes where the abundance is too low.
    # here we define as low values lower than 0.1 accordingly to this reference:
    # doi:10.1038/nrmicro3400
    # for instance, in the case of holobionts, I used the 0.1 since the quantities of phyla
    # for <0.1 and <0.01, respectively, were 119 and 66 within a total of 169 phyla.
    # TODO: Insert treatment here if the program is running in SparCC mode
    # ****  for compute the rare vs. abund regarding to the abundRel object.
    if False:
        i = 0; toDel = []
        abd = abundRel
        for x in list(map(lambda x: x < rareAbund, abd)):
            # deleting values where the abundance was different from the expected
            if x: toDel.append(i)
            i+=1

        debug("Removing "+str(len(toDel))+" nodes which have unexpected abundance")
        anotShared = np.delete(anotShared, toDel, axis=0)

        # correcting the boundaries for the rest of the program
        nsamp       = anotShared.shape[0]
        labelSamp   = np.delete(labelSamp, toDel, axis=0)
        abund       = list(np.delete(abund, toDel, axis=0))
        abundLog    = list(np.delete(abundLog, toDel, axis=0))
        prevalence  = list(np.delete(prevalence, toDel, axis=0))

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
    # TODO: URGENT!
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
    # $[mkdir -p "sparcc/raw_data" "sparcc/gephi_data" "sparcc/cnm_data" "sparcc/nga_data" "sparcc/figures" "sparcc/sparcc_data" "sparcc/liasp_data" "sparcc/nx_data" "sparcc/matrices"]
    mkdir_p(["sparcc/raw_data", "sparcc/gephi_data", "sparcc/cnm_data", "sparcc/nga_data", "sparcc/figures", "sparcc/sparcc_data", "sparcc/liasp_data", "sparcc/nx_data", "sparcc/matrices"])

    debug("Getting already computed SparCC matrix...")
    # reading the data with all rows and columns sorted
    coSpccPar = pd.read_csv('../'+spcc_backlog,sep='\t',index_col=0).sort_index().sort_index(axis=1).values
    coSpccPar.drop(index=delete,inplace=True)
    coSpccPar.drop(columns=delete,inplace=True)

    debug("Converting calculated data to DataFrame...")

    # inserting calculated correlations inside a DataFrame
    coSparCC = pd.DataFrame(data=coSpccPar, index=labelSamp, columns=labelSamp)

    debug("Saving non-filtered correlation data...")

    # printing non filtered correlation matrix to file
    coSparCC.to_csv("sparcc/raw_data/correlation_matrix.csv")

    # printing non filtered correlation matrix to file as gephi format
    gp.printCorr("sparcc/gephi_data/correlation_matrix.csv", np.array(coSpccPar), a, b, c)

    # printing nodes' labels using gephi format
    gp.printNodes("sparcc/gephi_data/nodes.csv", meta)

    debug("Executing Critical Network Method (CNM) and Largest Influence on Average Shortest Path (LIASP) algorithms...")

    # running code-integration steps
    # this step is responsible for running the CNM and LIASP algorithms
    # and identifying the most important nodes in the network (the keystones)
    itg.run("sparcc", coSparCC, host)

    cd('..')
    debug("Exiting...")

if __name__ == '__main__':
    assert len(sys.argv) > 3, "Missing arguments!\nUsage: ./main <community> <metadata> <host> <spcc_corrpath>\n"
    main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
