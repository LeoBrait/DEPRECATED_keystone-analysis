# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 09:07:32 2018
Modified on Sat Nov 19 09:30:00 2022

@author: Flavia Mayumi
modifications: Rafael Menezes
"""
import numpy as np
import networkx as nx
import pandas as pd
from scipy import stats
from collections import defaultdict
from functools import lru_cache

def find_subgraphs(GG,Taxa):
    '''
    GG is a networkx graph
    Taxa is a list of taxa

    the function returns a list of the subgraphs of GG, sorted by size (from largest to smallest)
    '''
    G = GG.copy()
    dis = list(nx.isolates(G))
    # node_list = list(Taxa)
    # for el in dis:
    #     node_list.remove(el)
    G.remove_nodes_from(dis)
    return sorted(nx.connected_components(G), key=len, reverse=True)

def max_diameter_subgraphs(GG,Taxa):
    resp = find_subgraphs(GG,Taxa)
    diameters = []
    for el in resp:
        g = GG.subgraph(el)
        sp = nx.diameter(g)
        diameters.append(sp)
    return max(np.array(diameters))

def get_values(G,par):
    dict_par = G.nodes(par)
    return np.array([j for (i,j) in dict_par])

def get_node_at(G,atr,value):
    node = [x for x,y in G.nodes(data=True) if y[atr]==value]
    return node[0]

def rank_x(x):
    x =np.array(x)
    a = np.array([[y,i] for i, y in enumerate(x.argsort()[::-1])])
    b = a[a[:,0].argsort(),1]
    b = b+1
    return b

def construct_rank_dict_set_to_G(G, Taxa, metric_dict, metric):
    '''
    G is a networkx graph
    Taxa is a list of taxa (nodes of G)
    metric_dict is a dictionary of values of the metric with the taxa as keys
    metric is the name of the metric

    the function:
    1) constructs a dictionary with the ranks of the values of the metric
    2) sets the dictionary as an attribute of the nodes of G
    '''
    # construct the dictionary with the ranks of the values of the metric
    metric_values = np.array([metric_dict[t] for t in Taxa])
    rank_k = np.argsort(metric_values) # returns an array with the indices of the sorted values
    rank_k = np.argsort(rank_k) + 1 # returns an array with the ranks of the values of the metric
    dict_k = dict(zip(Taxa, metric_values))
    dict_rank_k = dict(zip(Taxa, rank_k))

    # set the dictionary as an attribute of the nodes of G
    nx.set_node_attributes(G,values=dict_k,name=metric)
    nx.set_node_attributes(G,values=dict_rank_k,name=metric+"_rank")

def _construct_rank_dict_set_to_G(G,Taxa,dict_metric,metric):
    nx.set_node_attributes(G,values=dict_metric,name=metric)
    k = []
    for taxon in Taxa:
        k.append(dict_metric[taxon])
    rank_k = rank_x(k)
    dict_rank_k = dict(zip(Taxa, rank_k))
    nx.set_node_attributes(G,values=dict_rank_k,name=metric+"_rank")
    return k

def filter_disconnected_nodes(G):
    g = G.copy()
    toDel = [node for node,deg in g.degree if deg == 0]
    g.remove_nodes_from(toDel)
    return g

def drop_weights(G):
    '''Drop the weights from a networkx weighted graph.'''
    # https://stackoverflow.com/a/74540827/13845224
    for node, edges in nx.to_dict_of_dicts(G).items():
        for edge, attrs in edges.items():
            attrs.pop('weight', None)

def calc_equivalent_increase_in_path_length(eff1, eff2):
    '''
    given the efficiencies of two graphs, the function returns the equivalent increase in path length.

    The equivalent increase in path length is the increase in path of the network 1 that would be necessary to obtain the same efficiency of network 2.

    the function returns the equivalent increase in path length of the efficiency matrices
    '''

    return eff1/eff2 - 1

# @lru_cache(maxsize=256)
def generate_efficiency_matrix(G):
    '''
    G is a networkx graph

    the function returns the efficiency matrix of G
    '''
    # get the distance matrix
    dist = nx.floyd_warshall_numpy(G)

    # transform the distance matrix into an efficiency matrix
    np.fill_diagonal(dist, 1)
    eff = 1/dist
    np.fill_diagonal(eff, 0)

    return eff

def LIASP(G, node, node_idx = None):
    '''
    G is a networkx graph
    node is the index of a node of G (integer)
    node_label is the label of the node (string)

    the function returns the contribution to dissimilarity from direct interactions of the removal of node
    '''
    # get the index of the node
    if node_idx is None:
        node_idx = list(G.nodes()).index(node)
    
    import time
    # start = time.time()

    # remove weights from the graph
    G = G.copy()
    drop_weights(G)

    # number of nodes in the graph
    n = len(G.nodes())

    # get the efficiency matrix of G
    eff_mat_original = generate_efficiency_matrix(G)

    # remove all edges connected to node from G
    G_ = G.copy()
    edges_incident_to_node = list(G_.edges(node))
    G_.remove_edges_from(edges_incident_to_node)

    # get the efficiency matrix of G_
    eff_mat_removed = generate_efficiency_matrix(G_)
    
    # calculate the efficiencies of original and modified graphs
    eff_original = nx.global_efficiency(G)
    eff_removed = nx.global_efficiency(G_)

    # calculate the difference between the efficiency matrices
    delta_eff = eff_mat_original - eff_mat_removed
    try:
        assert np.all(delta_eff >= 0), "Error in LIASP Routine: For some pair of nodes, the efficiency of their interaction in the original graph is smaller than the efficiency of their interaction in the graph without node" + str(node)
        assert np.all(eff_mat_original >= 0), "Error in LIASP Routine: For some pair of nodes, the efficiency of their interaction in the original graph is smaller than 0 in node" + str(node)
        assert np.all(eff_mat_removed >= 0), "Error in LIASP Routine: For some pair of nodes, the efficiency of their interaction in the graph without node is smaller than 0 in node" + str(node)
    except AssertionError as e:
        print(e)
        print("Efficiency matrix of original graph:")
        print(eff_mat_original)
        print("Efficiency matrix of graph without node:")
        print(eff_mat_removed)
        print("Delta efficiency matrix:")
        print(delta_eff)
        raise

    # calculate the total dissimilarity and the contribution to dissimilarity from direct interactions of the removal of node
    LIASP = eff_original - eff_removed#np.sum(delta_eff)/(n*(n-1))
    # assert np.allclose(LIASP, eff1 - eff2), "Error in LIASP Routine: LIASP is not equal to the difference in global efficiency of the graph without node" + str(node)
    LIASP_dir = 2*np.sum(delta_eff[node_idx, :])/(n*(n-1))

    # create a copy of delta_eff, and remove lines and columns corresponding to node
    delta_eff_ = delta_eff.copy()
    delta_eff_ = np.delete(delta_eff_, node_idx, 0)
    delta_eff_ = np.delete(delta_eff_, node_idx, 1)

    # calculate the contribution to dissimilarity from indirect interactions of the removal of node
    LIASP_indir = np.sum(delta_eff_)/(n*(n-1))

    # Sanity check
    try:
        # assert LIASP >= 0, "Error in LIASP Routine: LIASP is negative for node " + str(node)
        # print(f"LIASP: {LIASP}, LIASP_dir: {LIASP_dir}, LIASP_indir: {LIASP_indir}")
        # input()
        assert all(np.array([LIASP, LIASP_dir, LIASP_indir]) >= 0), "Error in LIASP Routine: LIASP, LIASP_dir, LIASP_indir are not all positive for node " + str(node)
        assert np.allclose(LIASP - LIASP_dir, LIASP_indir) , "Error in LIASP Routine: direct and indirect components of LIASP do not sum to the total value for graph without node " + str(node)
        assert all(np.array([LIASP, LIASP_dir, LIASP_indir]) < 1.), "Error in LIASP Routine: LIASP, LIASP_dir, LIASP_indir are not all smaller then 1 for node " + str(node)
    except AssertionError as e:
        print(e)
        print("LIASP: ", LIASP)
        print("LIASP_dir: ", LIASP_dir)
        print("LIASP_indir: ", LIASP_indir)
        print("Efficiency matrix of original graph:")
        print(eff_mat_original)
        print("Efficiency matrix of graph without node:")
        print(eff_mat_removed)
        print("Delta efficiency matrix:")
        print(delta_eff)
        raise

    path_increase = LIASP/eff_removed
    path_increase_dir = LIASP_dir/eff_removed
    path_increase_indir = LIASP_indir/eff_removed

    # end = time.time()
    # print("Time to calculate the LIASP indexes: %f"%(end-start))

    return {'LIASP': LIASP,
            'LIASPdir': LIASP_dir,
            'LIASPindir': LIASP_indir,
            'effOriginal': eff_original,
            'effRemoved': eff_removed,
            'pathIncrease': path_increase,
            'pathIncreaseDir': path_increase_dir,
            'pathIncreaseIndir': path_increase_indir,
            'node': node}


def LIASP_dissimilarity(base, G, prefix, Taxa, peak):

    dictionaries = []
    for node in G.nodes():

        # calculate LIASP for node
        LIASP_dict = LIASP(G, node)

        # add EDpDM and iEDpDM to the dictionary for backwards compatibility
        LIASP_dict['EDpDM'] = LIASP_dict['LIASP']
        LIASP_dict['iEDpDM'] = LIASP_dict['LIASPindir']

        # store the dictionary
        dictionaries.append(LIASP_dict)

    # convert a list of dictionaries to a dictionary of lists
    # https://stackoverflow.com/a/11450683/13845224
    metrics = defaultdict(list)
    for node_dict in dictionaries:
        for key, val in node_dict.items():
            metrics[key].append(val)

    # add the metrics to the graph
    for metric in metrics:
        # if metric is 'node' go to the next iteration
        if metric == 'node':
            continue

        # construct the rank of the metric and add it to the graph
        construct_rank_dict_set_to_G(G,Taxa,dict(zip(metrics['node'], metrics[metric])),metric)

        # save information
        save_central_values(base, G, metric, prefix, peak)

def _LIASP_dissimilarity(base,file_eu,G,rare,prefix,Taxa,peak):
    # read LIASP data from the output of the Fortran code (liasp2.for)
    data = np.loadtxt(file_eu).T
    diss_total = data[-4] # total dissimilarity
    diss_dir = data[-3] # contribution to dissimilarity from direct interactions
    eff_orig = data[-2] # efficiency of the original network
    eff_mod = data[-1] # efficiency of the modified network

    # remove the first element of the arrays (0.)
    diss_total = diss_total[1:]
    diss_dir = diss_dir[1:]
    eff_orig = eff_orig[1:]
    eff_mod = eff_mod[1:]

    # define the indirect contribution as the difference between the full and direct contribution
    diss_ind = diss_total - diss_dir

    # calculate the equivalent increase in path length
    # if the efficiency of the modified network is 0, the equivalent increase in path length is set to 0
    
    path_inc_total = np.zeros(len(diss_total))
    path_inc_dir = np.zeros(len(diss_dir))
    path_inc_ind = np.zeros(len(diss_ind))

    for idx, efficiency in enumerate(eff_mod):
        if efficiency > 0:
            path_inc_total[idx] = diss_total[idx]/efficiency
            path_inc_dir[idx] = diss_dir[idx]/efficiency
            path_inc_ind[idx] = diss_ind[idx]/efficiency

    # Define the metrics that are going to be added to the graph
    metrics = {"LIASP" : diss_total,
            "EDpDM" : diss_total, # for back compatibility
            "LIASP_dir" : diss_dir,
            "LIASP_indir" : diss_ind,
            "iEDpDM" : diss_ind, # for back compatibility
            "path_increase" : path_inc_total,
            "path_increase_dir" : path_inc_dir,
            "path_increase_indir" : path_inc_ind}

    # add the metrics to the graph
    for metric in metrics:
        construct_rank_dict_set_to_G(G,Taxa,dict(zip(Taxa,metrics[metric])),metric)
        save_central_values(base, G, metric, prefix, peak)

def save_central_values(base,G,metric,prefix,peak):
    # G_ = filter_disconnected_nodes(G)
    nnk = get_values(G,metric)
    k = np.array(nnk)
    number_keystones(base,k,metric,prefix,peak)

def prod_dict(x1,x2,Taxa):
    '''
    x1 and x2 are dictionaries

    the function returns a dictionary with the product of the values on dictionaries x1 and x2
    '''
    prod = []
    for taxon in Taxa:
        prod.append(x1[taxon]*x2[taxon])
    return dict(zip(Taxa, prod))

def number_keystones(base,k,metric,prefix,peak):
    X = k.copy()
    X.sort()
    mx = np.mean(X)
    stdx = np.std(X, ddof=1)
    medx = np.median(X)
    cutoff1 = mx + stdx
    cutoff2 = mx + 2*stdx
    cutoff3 = medx + stdx
    cutoff4 = medx + 2*stdx
    Nkey1 = np.sum(X>cutoff1)
    Nkey2 = np.sum(X>cutoff2)
    Nkey3 = np.sum(X>cutoff3)
    Nkey4 = np.sum(X>cutoff4)
    f = open(base+'/figures/0p%d/'%peak+metric+'_number_keystones_'+prefix + ".csv",'w')
    f.write("mean,{:.7e}\n".format(mx))
    f.write("median,{:.7e}\n".format(medx))
    f.write("std,{:.7e}\n".format(stdx))
    f.write("number of keystones (mean + std),{}\n".format(Nkey1))
    f.write("number of keystones (mean + 2std),{}\n".format(Nkey2))
    f.write("number of keystones (median + std),{}\n".format(Nkey3))
    f.write("number of keystones (median + 2std),{}\n".format(Nkey4))
    f.close()
    return mx, medx, stdx