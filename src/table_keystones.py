# -*- coding: utf-8 -*-
"""
Spyder Editor

"""
import numpy as np
import networkx as nx
import pandas as pd
from scipy import stats

rare = [0.1,1,5]

def calculateZScoreMax(G,metric,mea,med,std):
    taxaG = list(G.nodes)
    listMetric = np.array([G.nodes[x][metric] for x in taxaG])
    maxMetric = max(listMetric)
    zscoreMetric = stats.zscore(listMetric)
    dict_zscoreMetric = dict(zip(taxaG, zscoreMetric))
    dict_maxMetric = dict(zip(taxaG, listMetric/maxMetric))
    iskeeeys = list(map(lambda x: int(x > mea+2*std), listMetric))
    dict_isKeystone = dict(zip(taxaG, iskeeeys))
    nx.set_node_attributes(G,values=dict_zscoreMetric,name="zscore"+metric)
    nx.set_node_attributes(G,values=dict_maxMetric,name=metric+"scaled")
    nx.set_node_attributes(G,values=dict_isKeystone,name=metric+"isKey")

def getKeystonenessCutts(base,peak,mode):
    df = pd.read_csv(base+'/figures/0p%d/%s_number_keystones_positive_%s_0p%d.csv'%(peak,mode,base,peak), index_col=0, header=None)
    mea = df.loc['mean']
    med = df.loc['median']
    std = df.loc['std']

    return mea,med,std

# creating a csv table to save the desired metric as keystoneness and a Venn diagram for caparision
def table_analysis(base, peak, G, host):
    tmp = host.split('-')
    envi0 = ' '.join('-'.join(tmp[:-1]).split('.'))
    envi1 = ' '.join(tmp[-1].split('.')[:-1])

    metrics = ['EDpDM','iEDpDM','D','BC','DxCC','Ddiv','DxCCdiv']

    # for mode in ['ED','D','BC']:
    for mode in metrics:
        # getting the mean, median and std of the environment metric
        mea, med, std = getKeystonenessCutts(base,peak,mode)
        calculateZScoreMax(G, mode, mea, med, std) # setting zscore and scaled metric as node attr

    dataN = {'EDpDM':[],'iEDpDM':[],'D':[],'BC':[],'DxCC':[],'Ddiv':[],'DxCCdiv':[]}

    #keystone or non keystone = se é keystone por media + sigma, se é keystone por mediana + sigma
    df = pd.DataFrame(columns=['Ecosystem', 'Habitat', 'Taxon', 'Abundance', 'Std', 'Sem', 'Prevalence', 'Median']
    + [x+y for x in metrics for y in ['', '_Z-score', '_Scaled', '_isKeystone', '_rank']])

    for j,i in enumerate(G):
        node = G.nodes[i]
        abd = node['abundance']
        ast = node['std']
        sem = node['sem']
        pre = node['prevalence']
        amd = node['med']

        df.loc[j,:] = (envi0, envi1, i, abd, ast, sem, pre, amd) +\
            tuple(node[y] for mode in metrics for y in [mode,'zscore'+mode,mode+'scaled',mode+'isKey','rank_'+mode])

        # adding to the list of dicts if it's keystone in any metric
        for mode in metrics:
            if node[mode+'isKey']:
                dataN[mode].append(i)

    for mode in metrics:
        df.sort_values(by=[mode, 'Abundance'], ascending=False, inplace=True)
        df[mode+'_rank'] = df[mode].rank(ascending=False,method='min').astype(int)

    df.sort_values(by=['Taxon'], ascending=True, inplace=True)
    df.set_index('Taxon', inplace=True)
    df.to_csv(base+'/figures/0p%d/keystones.csv'%peak)
