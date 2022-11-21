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

def construct_dic_rank_dict_set_to_G(G,Taxa,k,metric):
    rank_k = rank_x(k)
    dict_k = dict(zip(Taxa, k))
    dict_rank_k = dict(zip(Taxa, rank_k))
    nx.set_node_attributes(G,values=dict_k,name=metric)
    nx.set_node_attributes(G,values=dict_rank_k,name="rank_"+metric)

def construct_rank_dict_set_to_G(G,Taxa,dict_metric,metric):
    nx.set_node_attributes(G,values=dict_metric,name=metric)
    k = []
    for taxon in Taxa:
        k.append(dict_metric[taxon])
    rank_k = rank_x(k)
    dict_rank_k = dict(zip(Taxa, rank_k))
    nx.set_node_attributes(G,values=dict_rank_k,name="rank_"+metric)
    return k

def filter_disconnected_nodes(G):
    g = G.copy()
    toDel = [node for node,deg in g.degree if deg == 0]
    g.remove_nodes_from(toDel)
    return g


def LIASP_dissimilarity(base,file_eu,G,rare,prefix,Taxa,peak):
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
    metrics = {"LIASP_total" : diss_total,
            "EDpDM" : diss_total, # for back compatibility
            "LIASP_dir" : diss_dir,
            "LIASP_ind" : diss_ind,
            "iEDpDM" : diss_ind, # for back compatibility
            "path_inc_total" : path_inc_total,
            "path_inc_dir" : path_inc_dir,
            "path_inc_indir" : path_inc_ind}

    # add the metrics to the graph
    for metric in metrics:
        construct_rank_dict_set_to_G(G,Taxa,dict(zip(Taxa,metrics[metric])),metric)
        save_central_values(base, G, metric, prefix, peak)

def save_central_values(base,G,metric,prefix,peak):
    G_ = filter_disconnected_nodes(G)
    nnk = get_values(G_,metric)
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

# A function to save the LIASP information in the graph G
# information saved: ['LIASP-total', 'LIASP-direct', 'LIASP-indirect', 'pathIncrease-total', 'pathIncrease-direct', 'pathIncrease-indirect', 'keystone']
def save_keystone_info_in_G(G, Taxa, LIASP, pathIncrease, keystone):
    # Create dictionaries
    dict_LIASP_total = dict(zip(Taxa, LIASP[:,0]))
    dict_LIASP_direct = dict(zip(Taxa, LIASP[:,1]))
    dict_LIASP_indirect = dict(zip(Taxa, LIASP[:,2]))
    dict_pathIncrease_total = dict(zip(Taxa, pathIncrease[:,0]))
    dict_pathIncrease_direct = dict(zip(Taxa, pathIncrease[:,1]))
    dict_pathIncrease_indirect = dict(zip(Taxa, pathIncrease[:,2]))
    dict_keystone = dict(zip(Taxa, keystone))

    # Add attributes to the graph
    nx.set_node_attributes(G,values=dict_LIASP_total,name="LIASP-total")
    nx.set_node_attributes(G,values=dict_LIASP_direct,name="LIASP-direct")
    nx.set_node_attributes(G,values=dict_LIASP_indirect,name="LIASP-indirect")
    nx.set_node_attributes(G,values=dict_pathIncrease_total,name="pathIncrease-total")
    nx.set_node_attributes(G,values=dict_pathIncrease_direct,name="pathIncrease-direct")
    nx.set_node_attributes(G,values=dict_pathIncrease_indirect,name="pathIncrease-indirect")
    nx.set_node_attributes(G,values=dict_keystone,name="keystone")

def save_habitat_keystone_table(path, G, Taxa, habitat):
    # TODO: improve documentation
# I think the best thing would be to have a {habitat}_LIASP.csv file with each taxon as a row and columns being: [LIASP-total, LIASP-direct, LIASP-indirect, pathIncrease-total, pathIncrease-direct, pathIncrease-indirect, keystone]
# keystone is a binary variable that is 1 if the taxon is a keystone and 0 otherwise    
# I will need to add the information about the keystone status to the table
    # initialize the table
    df = pd.DataFrame(index=Taxa, columns=['LIASP-total', 'LIASP-direct', 'LIASP-indirect', 'pathIncrease-total', 'pathIncrease-direct', 'pathIncrease-indirect', 'keystone'])

    # fill the table
    for taxon in Taxa:
        df.loc[taxon, 'LIASP-total'] = G.nodes[taxon]['LIASPefficiency']
        df.loc[taxon, 'LIASP-direct'] = G.nodes[taxon]['LIASPdirect']
        df.loc[taxon, 'LIASP-indirect'] = G.nodes[taxon]['LIASPindirect']
        df.loc[taxon, 'pathIncrease-total'] = G.nodes[taxon]['pathIncrease']
        df.loc[taxon, 'pathIncrease-direct'] = G.nodes[taxon]['pathIncreaseDirect']
        df.loc[taxon, 'pathIncrease-indirect'] = G.nodes[taxon]['pathIncreaseIndirect']
        df.loc[taxon, 'keystone'] = G.nodes[taxon]['keystone']

    # save the table
    df.to_csv(path + f'{habitat}_LIASP.csv')

    return df