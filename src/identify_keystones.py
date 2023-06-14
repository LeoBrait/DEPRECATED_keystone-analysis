"""

File used to hold the functions to integrate Critical Network Method (CNM) and Newman-Girvan Analyses (NGA).

"""

from xonsh_py import echo, mkdir_p
import numpy as np
import pandas as pd
from keystone_indexes_functions import construct_rank_dict_set_to_G, LIASP_dissimilarity, find_subgraphs, prod_dict, save_central_values
from metadata_topython_functions import metadata
from table_keystones import table_analysis
from detect_peaks import detect_peaks
import networkx as nx

# this function aims to execute the following steps:
# - create the apropriate arg file for CNM code
# - call the executable for the proper correlation matrix
# - generate visual outputs (3rd party codes) and select the desired peaks
# - run the NGA algorithm for the data obtained
def identify_keystones(path=None, df=None, host=None):
    if not isinstance(path, str):
        raise Exception("`path` must be a string. It stores the path to the containing folder.")
    if not isinstance(df, pd.DataFrame):
        raise Exception("`df` must be a pandas DataFrame containing the correlation matrix") 
    if not isinstance(host , str):
        raise Exception("`host` must be a string containing the name of the environment")

    # choosing peaks
    peaks_list, high = deltaFuncAndPeaksSelector(df)

    # going a step before the peak in order to construct the networks before they break apart
    peaks_list[:] = [x-1 for x in peaks_list]
    high -= 1

    # Save the highest peak in a file
    echo(path+'/raw_data/cnm_highest_peak.txt', high)

    # doing only the highest peak in order to improve performance as we are not using other peaks
    peaks_list = [high]

    # executing the keystone finder
    for peak in peaks_list:
        # creating the folder for the peak
        mkdir_p([path+'/figures/0p'+str(peak)])
        # calling the function to find the keystones
        find_keystones_envi(path,peak,host)

def drop_weights(G):
    """
    Drop the weights from a networkx weighted graph.
    
    Parameters
    ----------
    G : networkx graph
        Graph of the ecosystem. Modified in place.

    Returns
    -------
    None

    Source
    ------
    https://stackoverflow.com/a/74540827/13845224
    """
    for node, edges in nx.to_dict_of_dicts(G).items():
        for edge, attrs in edges.items():
            attrs.pop('weight', None)

def treshold_efficiency_dissimilarities(G, tresholds):
    """
    The algorithm removes the least important links in a network in a sequential way following tresholds values defined in `tresholds`.
    
    Returns a list with the dissimilarities between networks for each treshold value.

    Parameters
    ----------

    G : networkx graph
        Graph of the ecosystem.

    tresholds : list
        List of tresholds to be used in the algorithm.

    Returns
    -------
    dissimilarities : list
        List of dissimilarities between networks for each treshold value.
    """

    # creating the list of efficiencies
    efficiencies = []

    # store the efficiency of the original network
    efficiencies.append(nx.global_efficiency(G))

    # creating a copy of the graph
    Gc = G.copy()

    # loop over the tresholds
    for t in tresholds:

        # removing the edges with weight below the treshold
        Gc.remove_edges_from([(u,v) for (u,v,d) in Gc.edges(data=True) if d['weight'] < t])

        # calculating the efficiency of the network
        eff = nx.global_efficiency(Gc)

        # storing the efficiency
        efficiencies.append(eff)

    # calculating the dissimilarities between networks
    dissimilarities = -np.diff(np.array(efficiencies))

    # print(efficiencies)
    # input()
    # print(dissimilarities)

    return dissimilarities


def deltaFuncAndPeaksSelector(df):
    """
    This function aims to select the peaks in the dissimilarity curve.

    Parameters
    ----------

    df : pandas DataFrame
        Correlation matrix of the ecosystem.

    Returns
    -------
    argpeaks : list
        List of the peaks in the dissimilarity curve.

    high : int
        Highest peak in the dissimilarity curve.
    """

    if not isinstance(df,pd.DataFrame):
        raise Exception("Wrong path data type!")

    # converting the dataframe to a networkx graph
    G = nx.from_pandas_adjacency(df)
    
    # defining the tresholds
    tresholds = np.linspace(0,1.,101)[:-1]

    # calculating the dissimilarities
    # measure the time the function takes to run
    # start = time.time()
    dissim = treshold_efficiency_dissimilarities(G, tresholds)
    # end = time.time()
    # print("Time to calculate the dissimilarities: %f"%(end-start))

    # set minimum peak height so that identified peaks are among the 5% highest peaks
    minimum_peak_height = sorted(dissim)[-int(len(dissim)*.05)]

    # detect peaks in the dissimilarity curve
    # mpd = 1 means that the minimum distance between peaks is 1 (no consecutive peaks)
    # start = time.time()
    argpeaks = detect_peaks(dissim, mpd=1, mph=minimum_peak_height )
    # end = time.time()
    # print("Time to detect peaks: %f"%(end-start))

    # getting the highest peak separately or 0 if there's no peak
    if len(argpeaks) < 1:
        argpeaks = [0]
    high = np.where(dissim == max(dissim[argpeaks]))[0][0]

    return argpeaks, high

def find_keystones_envi(base, peak, host):
    """
    This function aims to find the keystones in the ecosystem.

    Parameters
    ----------

    base : string
        Path to the base directory of the ecosystem.

    peak : string
        Peak of the Critical Network to be considered.

    host : string
        Name of the ecosystem and habitat.

    """
    if not isinstance(base,str):
        raise Exception("`base` path must be a string")

    #In this line I'm assuming that the ecological features has the name of the directory - class + _ecol_features_class.csv, but we have to see if this pattern has been followed and we can change the directorie if necessarie
    N, Marker, Taxa, Abundance, AbundanceAbsolute, Prevalence, Sem, Std, Med = metadata(base+"/gephi_data/nodes.csv")

    # Read the correlation matrix from raw data
    df = pd.read_csv(base+'/raw_data/correlation_matrix.csv', index_col=0)
    correlation_matrix_pos = df.values
    np.fill_diagonal(correlation_matrix_pos,0) # removing auto-edges that may exist

    # Construct the Critical Network from the correlation matrix and the threshold
    threshold = peak/100.
    correlation_matrix_pos[correlation_matrix_pos<threshold] = 0

    # Saving filtered matrices for this peak (it'll be used on other programs)
    df[:] = correlation_matrix_pos
    df.to_csv(base+'/matrices/filtered_correlation_matrix_%3.2f_pos.csv'%threshold)

    # Transform the filtered correlation matrix into a Networkx graph
    G_pos = nx.from_numpy_matrix(correlation_matrix_pos)

    # saving adjacency matrix for the graphs
    adj = (correlation_matrix_pos != 0.)*1
    np.savetxt(base+'/matrices/adjacency_matrix_%3.2f_pos.txt'%threshold,adj,fmt='%d')

    # create dictionaire with the taxa and relabel nodes
    dict_taxa = dict(enumerate(Taxa))
    G_pos = nx.relabel_nodes(G_pos, dict_taxa)

    # create dictionaires for Abundance, AbundanceAbsolute, Marker, Prevalence, Sem, Std, Med
    dict_abundance =  dict(zip(Taxa, Abundance))
    dict_absolute = dict(zip(Taxa, AbundanceAbsolute))
    dict_marker = dict(zip(Taxa,Marker))
    dict_prevalence = dict(zip(Taxa,Prevalence))
    dict_sem = dict(zip(Taxa,Sem))
    dict_std = dict(zip(Taxa,Std))
    dict_med = dict(zip(Taxa,Med))

    # set Abundance, AbundanceAbsolute, Marker, Prevalence, Std and Sem as attributtes of the nodes
    nx.set_node_attributes(G_pos,values=dict_abundance,name='Abundance')
    nx.set_node_attributes(G_pos,values=dict_absolute,name='AbundanceAbsolute')
    nx.set_node_attributes(G_pos,values=dict_marker,name='Marker')
    nx.set_node_attributes(G_pos,values=dict_prevalence,name='Prevalence')
    nx.set_node_attributes(G_pos,values=dict_sem,name='Sem')
    nx.set_node_attributes(G_pos,values=dict_std,name='Std')
    nx.set_node_attributes(G_pos,values=dict_med,name='Median')

    # Create dictionaries for the degree (dd_pos), betweeness centrality (bc_pos),
    # closeness centrality (cc_pos), and the product of degree and closeness centrality
    # of the nodes in the critical network
    # The dictionaries have the structure: {node: value}
    dd_pos = G_pos.degree()
    d_pos = {i:j for (i,j) in dd_pos}
    bc_pos = nx.betweenness_centrality(G_pos)
    cc_pos = nx.closeness_centrality(G_pos)
    dxcc_pos = prod_dict(d_pos, cc_pos, Taxa)

    # Create dictionaries for "normalized" degree (ndd_pos), betweeness centrality (nbc_pos),
    # maxN = sum([len(i) for i in subgraphs])
    # d_div_pos = {i:d_pos[i]/maxN for i in d_pos} # a new step that tries to 'normalise' the degree
    # dxcc_div_pos = prod_dict(d_div_pos,cc_pos,Taxa)

    # Create names based on directory
    prefix_pos = "positive_" + base + '_0p%d'%peak

    # calculate and store LIASP dissimilarity
    LIASP_dissimilarity(base, G_pos, prefix_pos, Taxa, peak)
    # base, G, prefix, Taxa, peak

    #Metrics Literature
    metrics = ['BC','D','DxCC']#,'Ddiv','DxCCdiv']
    list_dict_metrics_pos = [bc_pos,d_pos,dxcc_pos]#,d_div_pos,dxcc_div_pos]
    for i,metric in enumerate(metrics):
        construct_rank_dict_set_to_G(G_pos,Taxa,list_dict_metrics_pos[i],metric)
        save_central_values(base,G_pos,metric,prefix_pos,peak)

    nx.write_gpickle(G_pos,base+"/nx_data/"+prefix_pos+"_G.gpickle")

    # running the table_keystones procedure
    table_analysis(base, peak, G_pos, host)
