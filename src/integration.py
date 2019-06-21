from xonsh_py import *

"""

File used to hold the functions to integrate CNM and Newman-Girvan Analyses.

"""

import numpy as np
import pandas as pd
from keystone_indexes_functions import *
from metadata_topython_functions import *
from table_keystones import *
import matplotlib.pyplot as plt
from detect_peaks import detect_peaks

# this function aims to execute the following steps:
# - create the apropriate arg file for cnm code
# - call the executable for the proper correlation matrix
# - generate visual outputs (3rd party codes) and select the desired peaks
# - run the nga algorithm for the data obtained
def run(path=None,df=None,host=None):
    if not isinstance(path, str):
        raise Exception("Wrong path data type!")
    if not isinstance(df, pd.DataFrame):
        raise Exception("Wrong df data type!")
    if not isinstance(host , str):
        raise Exception("Wrong host data type!")

    # writing the matrix as tsv with the desired format
    df.to_csv(path+"/cnm_data/correlation_matrix.dat", sep='\t',index=False,header=False)

    # copying executable temporarily to the folder
    cp(['../cnm/a.out'], path+'/cnm_data/')
    cp(['../liasp/a.out'], path+'/edm_data/')
    oldDir = pwd()

    argBase =   "%d %d\n"\
                "%3.2f 0.01 %d\n"\
                "correlation_matrix.dat\n"\
                "p_envi_%3.2f.dat\n"\
                "n_envi_%3.2f.dat\n"\
                "a_envi_%3.2f.dat\n"\
                "-1 -1 -1\n"

    n = df.shape[0]

    # creating right args file
    echo(path+"/cnm_data/redecrit1mc13.dat", argBase % (n,n,0.,100,1.,1.,1.))

    # running CNM for the whole space
    cd(path+'/cnm_data')
    sexec('./a.out')
    cd(oldDir)

    # choosing peaks
    peaksL, high = deltaFuncAndPeaksSelector(path)
    echo(path+'/raw_data/cnm_highest_peak.txt', high)

    # doing only the highest peak in order to improve performance as we are not using else peaks
    peaksL = [high]

    cd(path+'/cnm_data')
    # running CNM for each peak
    for peak in peaksL:
        peak = peak/100.
        # creating right args file
        echo("redecrit1mc13.dat", argBase%(n,n,peak-.01,1,peak,peak,peak))
        sexec('./a.out')
    rm(['a.out'])
    cd(oldDir)

    # starting to produce euclidian distance method outputs
    # argfile for edm
    argBase = "%d 1000\n"\
            "a_envi_%3.2f.dat\n"\
            "eudist_%3.2f.dat\n"\
            "-1 -1 -1\n"

    # $[cd @(path+'/edm_data')]
    cd(path+'/edm_data')
    for peak in peaksL:
        cp(['../cnm_data/a_envi_%3.2f.dat'%(peak/100.)], '.')
        echo("liasp2.dat", argBase%(n,peak/100.,peak/100.))
        sexec('./a.out')
        rm(['a_envi_%3.2f.dat'%(peak/100.)])
    rm(['a.out'])
    cd(oldDir)

    # executing the keystone finder
    for peak in peaksL:
        mkdir_p([path+'/figures/0p'+str(peak)])
        find_keystones_envi(path,peak,host)

# ----------------------
# Pasted  code from graphics_cnm_analysis.py
# ----------------------

"""
Created on Wed Mar 14 13:44:45 2018

@author: flavia
"""


def deltaFuncAndPeaksSelector(way=None):
    if not isinstance(way,str):
        raise Exception("Wrong way data type!")

    # # Part 1:
    # # Read data of distance between netorks in function of threshold
    geral = way+'/cnm_data/p_envi_1.00.dat'
    f = np.loadtxt(geral, unpack=True)
    x = f[0]
    y = f[7]
    tam = len(x)

    cutoff=sorted(y)[-int(len(y)*.05)]

    argpeaks = detect_peaks(y,mpd=1,mph=cutoff)

    # getting the highest peak separately or 0 if there's no peak
    if len(argpeaks) < 1:
        argpeaks = [0]
    high = np.where(y == max(y[argpeaks]))[0][0]

    return argpeaks, high

# --------------------
# Pasted from find_keystones_all.py
# --------------------

"""
Spyder Editor

"""

def find_keystones_envi(base,peak,host):
    if not isinstance(base,str):
        raise Exception("Wrong base data type!")

    rare = [0.1,1,5]

    #In this line I'm assuming that the ecological features has the name of the directory - class + _ecol_features_class.csv, but we have to see if this pattern has been followed and we can change the directorie if necessarie
    N, Marker, Taxa, Abundance, AbundanceAbsolute, Prevalence, Sem, Std, Med = metadata(base+"/gephi_data/nodes.csv")

    df = pd.read_csv(base+'/raw_data/correlation_matrix.csv', index_col=0)
    correlation_matrix_pos = df.values
    np.fill_diagonal(correlation_matrix_pos,0) # removing auto-edges that may exist

    #Construct network, positve network and negative
    thr = peak/100.
    correlation_matrix_pos[correlation_matrix_pos<thr]=0

    # Saving filtered matrices for this peak (it'll be used on other programs)
    df[:] = correlation_matrix_pos
    df.to_csv(base+'/matrices/filtered_correlation_matrix_%3.2f_pos.csv'%thr)

    #correlation_matrix_neg = np.absolute(correlation_matrix_neg)
    G_pos = nx.from_numpy_matrix(correlation_matrix_pos)

    # saving adjacency matrix from graphs
    adj = (correlation_matrix_pos != 0.)*1
    np.savetxt(base+'/matrices/adjacency_matrix_%3.2f_pos.txt'%thr,adj,fmt='%d')

    #create dictionaire with the taxa and relabel nodes
    dict_taxa = dict(enumerate(Taxa))
    G_pos = nx.relabel_nodes(G_pos, dict_taxa)

    #create dictionaire with Abundance, AbundanceAbsolute, Marker
    dict_abundance =  dict(zip(Taxa, Abundance))
    dict_absolute = dict(zip(Taxa, AbundanceAbsolute))
    dict_marker = dict(zip(Taxa,Marker))
    dict_prevalence = dict(zip(Taxa,Prevalence))
    dict_sem = dict(zip(Taxa,Sem))
    dict_std = dict(zip(Taxa,Std))
    dict_med = dict(zip(Taxa,Med))

    #set Abundance, AbundanceAbsolute, Marker, Prevalence, Std and Sem as attributtes to the nodes
    nx.set_node_attributes(G_pos,values=dict_abundance,name='abundance')
    nx.set_node_attributes(G_pos,values=dict_absolute,name='abundance_absolute')
    nx.set_node_attributes(G_pos,values=dict_marker,name='marker')
    nx.set_node_attributes(G_pos,values=dict_prevalence,name='prevalence')
    nx.set_node_attributes(G_pos,values=dict_sem,name='sem')
    nx.set_node_attributes(G_pos,values=dict_std,name='std')
    nx.set_node_attributes(G_pos,values=dict_med,name='med')

    # some information to further use
    subgraphs = find_subgraphs(G_pos,Taxa)

    # #Calculate d, bc, cc and d*cc
    dd_pos = G_pos.degree()
    d_pos = {i:j for (i,j) in dd_pos}
    bc_pos = nx.betweenness_centrality(G_pos)
    cc_pos = nx.closeness_centrality(G_pos)
    dxcc_pos=prod_dict(d_pos,cc_pos,Taxa)

    maxN = sum([len(i) for i in subgraphs])
    d_div_pos = {i:d_pos[i]/maxN for i in d_pos} # a new step that tries to 'normalise' the degree
    dxcc_div_pos=prod_dict(d_div_pos,cc_pos,Taxa)

    #Create names based on directory
    prefix_pos = "positive_" + base + '_0p%d'%peak

    #Euclidean distance divided by higher subgraph diameter and by node mean shortestpath
    func = lambda x:350000*x
    euclidians_distance_divided_and_msp(base,base+'/edm_data/eudist_%3.2f.dat'%thr,G_pos,rare,prefix_pos,Taxa,peak)

    #Metrics Literature
    metrics = ['BC','D','DxCC','Ddiv','DxCCdiv']
    funcs = [lambda x: x*4000+10,lambda x: x*2.5,lambda x: x*2.5,lambda x: x*5, lambda x: x*5]
    list_dict_metrics_pos = [bc_pos,d_pos,dxcc_pos,d_div_pos,dxcc_div_pos]
    for i,metric in enumerate(metrics):
        func = funcs[i]
        construct_rank_dict_set_to_G(G_pos,Taxa,list_dict_metrics_pos[i],metric)
        save_central_values(base,G_pos,metric,prefix_pos,peak)

    nx.write_gpickle(G_pos,base+"/nx_data/"+prefix_pos+"_G.gpickle")

    # running the table_keystones procedure
    table_analysis(base, peak, G_pos,host)
