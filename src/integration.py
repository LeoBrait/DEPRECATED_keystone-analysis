"""

File used to hold the functions to integrate Critical Network Method (CNM) and Newman-Girvan Analyses (NGA).

"""

from xonsh_py import cp, echo, sexec, mkdir_p, rm, cd, pwd
import numpy as np
import pandas as pd
from keystone_indexes_functions import construct_rank_dict_set_to_G, euclidians_distance_divided_and_msp, find_subgraphs, prod_dict, save_central_values
from metadata_topython_functions import metadata
from table_keystones import table_analysis
from detect_peaks import detect_peaks
import networkx as nx

# this function aims to execute the following steps:
# - create the apropriate arg file for CNM code
# - call the executable for the proper correlation matrix
# - generate visual outputs (3rd party codes) and select the desired peaks
# - run the NGA algorithm for the data obtained
def run(path=None, df=None, host=None):
    if not isinstance(path, str):
        raise Exception("`path` must be a string. It stores the path to the containing folder.")
    if not isinstance(df, pd.DataFrame):
        raise Exception("`df` must be a pandas DataFrame containing the correlation matrix") 
    if not isinstance(host , str):
        raise Exception("`host` must be a string containing the name of the environment")

    # writing the matrix as tsv with the desired format
    df.to_csv(path+"/cnm_data/correlation_matrix.dat", sep='\t',index=False,header=False)

    # copying executable temporarily to the folder
    cp(['../cnm/cnm.out'], path+'/cnm_data/')
    cp(['../liasp/liasp.out'], path+'/liasp_data/')
    oldDir = pwd() # saving the path to current directory

    # ------------------------------ CNM --------------------------------------------
    # The Critical Network Method (CNM) is a method to identify a critical subnetwork.
    # The algorithm removes the least important links in a network until the network is disconnected.
    # By comparing the dissimilarity between consecutive networks, the algorithm identifies the critical
    # treshold, which is the point where the dissimilarity between the networks is the highest.
    # The critical network is the subnetwork in which all links below the critical threshold are removed.

    # The following code:
    # 1) generate auxiliary file with inputs for CNM code
    # 2) runs the CNM program (a compiled Fortran code)
    # 3) identifies the peaks in the sequence of dissimilarity values
    # 4) selects the highest peak as the critical threshold
    # 5) generates the critical network

    # 
    # generate auxiliary file with inputs for CNM code
    argBaseCNM =   "%d %d\n"\
                "%3.2f 0.01 %d\n"\
                "correlation_matrix.dat\n"\
                "p_envi_%3.2f.dat\n"\
                "n_envi_%3.2f.dat\n"\
                "a_envi_%3.2f.dat\n"\
                "-1 -1 -1\n" # "-1 -1 -1" indicates the end of the input file and is required by the CNM code

    # number of nodes
    n = df.shape[0]

    # creating right args file
    echo(path+"/cnm_data/redecrit1mc13.dat", argBaseCNM % (n,n,0.,100,1.,1.,1.))

    # running CNM for the whole space
    cd(path+'/cnm_data')
    # running CNM
    sexec('./cnm.out')
    # go to the original directory
    cd(oldDir)

    # choosing peaks
    peaks_list, high = deltaFuncAndPeaksSelector(path)

    # going a step before the peak in order to construct the networks before they break apart
    peaks_list[:] = [x-1 for x in peaks_list]
    high -= 1

    # Save the highest peak in a file
    echo(path+'/raw_data/cnm_highest_peak.txt', high)

    # doing only the highest peak in order to improve performance as we are not using else peaks
    peaks_list = [high]

    cd(path+'/cnm_data')
    # running CNM for each peak
    for peak in peaks_list:
        peak = peak/100.
        # generate auxiliary file with inputs for CNM code
        echo("redecrit1mc13.dat", argBaseCNM%(n,n,peak-.01,1,peak,peak,peak))
        # execute CNM program
        sexec('./cnm.out')
    # remove the temporary executable
    rm(['cnm.out'])
    # go back to the original directory
    cd(oldDir)

    # ------------------------------ LIASP ------------------------------------------- 
    # The Largest Influence on Average Shortest Path (LIASP) is a metric to identify
    # the most influential nodes in a network. The measure is based on the efficiency
    # of the network, which is the inverse of the average shortest path length.
    # We use LIASP to identify the most influential nodes in the critical network.
    # The algorithm identifies the direct and indirect influence of each node in the network.
    # The direct is the decrease in efficiency of the network due to the 0 efficiency of
    # the connections between the removed node and its neighbors after the node is removed.
    # The indirect is the decrease in efficiency of the network due to the increase in
    # the shortest path length among all the nodes in the network.

    # The following code:
    # 1) generate auxiliary file with inputs for LIASP code
    # 2) runs the LIASP program (a compiled Fortran code)

    # argfile for liasp
    argBaseLIASP = "%d 1000\n"\
            "a_envi_%3.2f.dat\n"\
            "eudist_%3.2f.dat\n"\
            "-1 -1 -1\n" # "-1 -1 -1" indicates the end of the input file and is required by the LIASP code

    # entering the LIASP folder
    cd(path+'/liasp_data')
    # running LIASP code for each peak
    for peak in peaks_list:
        # copy the correlation matrix to the folder
        cp(['../cnm_data/a_envi_%3.2f.dat'%(peak/100.)], '.')
        # generate auxiliary file with inputs for LIASP code
        echo("liasp2.dat", argBaseLIASP%(n,peak/100.,peak/100.))
        # execute LIASP program
        sexec('./liasp.out')
        # clean temporary argument file
        rm(['a_envi_%3.2f.dat'%(peak/100.)])
    # remove the temporary executable
    rm(['liasp.out'])
    # return to the original directory
    cd(oldDir)

    # ---------------------------- KEYSTONES ---------------------------------------- 
    # Keystones are the most influential nodes in a network. They are identified by
    # using the LIASP metric. We use a script to identify the keystones in the critical
    # network. The script reads the output of the LIASP code and selects the nodes with
    # the highest LIASP metric, identifying as keystones the nodes with LIASP above the
    # mean + 2 standard deviations of the distribution of LIASP values.

    # The following code:
    # 1) runs the find_keystones.py script to identify the keystones in the critical network

    # executing the keystone finder
    for peak in peaks_list:
        # creating the folder for the peak
        mkdir_p([path+'/figures/0p'+str(peak)])
        # calling the function to find the keystones
        find_keystones_envi(path,peak,host)

# ----------------------
# Pasted  code from graphics_cnm_analysis.py
# ----------------------

"""
Created on Wed Mar 14 13:44:45 2018

@author: flavia
"""


def deltaFuncAndPeaksSelector(path):
    if not isinstance(path,str):
        raise Exception("Wrong path data type!")

    # # Part 1:
    # # Read data of distance between netorks in function of threshold
    geral = path+'/cnm_data/p_envi_1.00.dat'
    dat_f = np.loadtxt(geral, unpack=True)
    # treshold = dat_f[0] # valor do threshold
    dissim = dat_f[7] # dissimilarity between network i and network i+1

    # set minimum peak height so that identified peaks are among the 5% highest peaks
    minimum_peak_height = sorted(dissim)[-int(len(dissim)*.05)]

    # detect peaks in the dissimilarity curve
    # mpd = 1 means that the minimum distance between peaks is 1 (no consecutive peaks)
    argpeaks = detect_peaks(dissim, mpd=1, mph=minimum_peak_height )

    # getting the highest peak separately or 0 if there's no peak
    if len(argpeaks) < 1:
        argpeaks = [0]
    high = np.where(dissim == max(dissim[argpeaks]))[0][0]

    return argpeaks, high

# --------------------
# Pasted from find_keystones_all.py
# --------------------

"""
Spyder Editor

"""

def find_keystones_envi(base,peak,host):
    if not isinstance(base,str):
        raise Exception("`base` path must be a string")

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

    # Calculate the degree (dd_pos), betweeness centrality (bc_pos), cc and d*cc
    dd_pos = G_pos.degree()
    d_pos = {i:j for (i,j) in dd_pos}
    bc_pos = nx.betweenness_centrality(G_pos)
    cc_pos = nx.closeness_centrality(G_pos)
    dxcc_pos=prod_dict(d_pos,cc_pos,Taxa)

    maxN = sum([len(i) for i in subgraphs])
    d_div_pos = {i:d_pos[i]/maxN for i in d_pos} # a new step that tries to 'normalise' the degree
    dxcc_div_pos = prod_dict(d_div_pos,cc_pos,Taxa)

    #Create names based on directory
    prefix_pos = "positive_" + base + '_0p%d'%peak

    #Euclidean distance divided by higher subgraph diameter and by node mean shortestpath
    func = lambda x:350000*x
    euclidians_distance_divided_and_msp(base,base+'/liasp_data/eudist_%3.2f.dat'%thr,G_pos,rare,prefix_pos,Taxa,peak)

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
