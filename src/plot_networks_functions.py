#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 14:27:16 2018

@author: flavia
"""

import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import community
from networkx.algorithms.community import greedy_modularity_communities
import pandas as pd

def find_modules_nga(GG,Taxa):
    G = GG.copy()
    dis = list(nx.isolates(G))
    node_list = list(Taxa)
    for el in dis:
        node_list.remove(el)
    G.remove_nodes_from(dis)
    communities_generator = community.girvan_newman(G)
    top_level_communities = next(communities_generator)
    next_level_communities = next(communities_generator)
    groups = sorted(map(sorted, next_level_communities))
    return groups

def find_modules_modularity(GG,Taxa):
    G = GG.copy()
    dis = list(nx.isolates(G))
    node_list = list(Taxa)
    for el in dis:
        node_list.remove(el)
    communities_generator = greedy_modularity_communities(G)
    groups = list(greedy_modularity_communities(G))
    return groups

def getKeystonenessCutts(base,peak,mode):
    df = pd.read_csv(base+'/figures/0p%d/%s_number_keystones_positive_%s_0p%d.csv'%(peak,mode,base,peak), index_col=0, header=None)
    mea = df.loc['mean']
    med = df.loc['median']
    std = df.loc['std']

    return mea,med,std

def plot_network_modules_sizes_positive(base,GG,subgraphs,Taxa,prefix,metric,function,peak):
    mea,med,std = getKeystonenessCutts(base,peak,metric)

    G = GG.copy()
    node_list = list(Taxa)
    plt.figure()
    dis = list(nx.isolates(G))
    for el in dis:
        node_list.remove(el)
    G.remove_nodes_from(dis)
    pos=nx.spring_layout(G,k=0.7)
    edges = G.edges()
    nodes = G.nodes(data=True)
    weights = [2.5*G[u][v]['weight'] for u,v in edges]
    nx.draw_networkx_edges(G, pos,width=weights, edgelist=edges,edge_color='k',alpha=0.2)
    for i,el in enumerate(subgraphs):
        list_el = list(el)
        for noh in list_el:
            color = plt.cm.RdPu(nodes[noh]['abundance'])
            nds = function(nodes[noh][metric])
            # plotting different shape when the node is keystone
            noh_shape = 'v' if nodes[noh][metric] > float(med)+float(std) else 'o'
            nx.draw_networkx_nodes(G, pos,  nodelist=[noh],node_size = nds,node_color=[color],alpha = .8,node_shape=noh_shape)

    rank_metric = G.nodes("rank_"+metric)
    dict_rank_metric = {i:j for (i,j) in rank_metric}
    nx.draw_networkx_labels(G, pos, labels = dict_rank_metric,font_size=6,font_color='k')
    plt.axis('off')
    plt.savefig(base+'/figures/0p%d/'%peak+metric+"_"+prefix+'_network.pdf',dpi=150)
