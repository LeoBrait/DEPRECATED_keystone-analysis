# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 09:07:32 2018

@author: flavia
"""
import numpy as np
import networkx as nx
from scipy import stats

def find_subgraphs(GG,Taxa):
    G = GG.copy()
    dis = list(nx.isolates(G))
    node_list = list(Taxa)
    for el in dis:
        node_list.remove(el)
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

def euclidians_distance_divided_and_msp(base,file_eu,G,rare,prefix,Taxa,peak):
    diam = max_diameter_subgraphs(G,Taxa)
    data = np.loadtxt(file_eu).T
    dist = data[-2] # node full contribution
    dist = dist[1:]/diam
    dist2 = data[-1] # node indirect contribution
    dist2 = dist2[1:]/diam

    for i,taxon in enumerate(Taxa):
        length = nx.single_source_shortest_path_length(G,taxon)
        scale = np.mean(np.array(list(length.values())))
        dist[i] = float(dist[i])/scale if float(dist[i]) > 0. else float(dist[i])
        dist2[i] = float(dist2[i])/scale if float(dist2[i]) > 0. else float(dist2[i])

    # full contribution
    metric = "EDpDM"
    construct_dic_rank_dict_set_to_G(G,Taxa,dist,metric)
    save_central_values(base,G,metric,prefix,peak)

    # indirect contribution
    metric = "iEDpDM"
    construct_dic_rank_dict_set_to_G(G,Taxa,dist2,metric)
    save_central_values(base,G,metric,prefix,peak)

def save_central_values(base,G,metric,prefix,peak):
    G_ = filter_disconnected_nodes(G)
    nnk = get_values(G_,metric)
    k = np.array(nnk)
    number_keystones(base,k,metric,prefix,peak)

def prod_dict(x1,x2,Taxa):
    prod = []
    for taxon in Taxa:
        prod.append(x1[taxon]*x2[taxon])
    return dict(zip(Taxa, prod))

def number_keystones(base,k,metric,prefix,peak):
    X = k.copy()
    X.sort()
    mx = np.mean(X)
    stdx = np.std(X,ddof=1)
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
