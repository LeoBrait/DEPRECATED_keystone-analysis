# -*- coding: utf-8 -*-
"""
Spyder Editor

"""
import numpy as np
import networkx as nx
import pandas as pd
from scipy import stats

rare = [0.1,1,5]

def identify_keystones_and_add_metric_props(G, metric):
    # get list of taxa (node labels) and list of metric values
    taxaG = list(G.nodes())
    listMetric = np.array([G.nodes[x][metric] for x in taxaG])

    # If the metric is not numerical, skip it
    if not np.issubdtype(listMetric.dtype, np.number):
        return

    # get max value of metric
    maxMetric = max(listMetric)
    maxMetric = maxMetric if maxMetric > 0 else 1 # avoid division by zero
    dict_maxMetric = dict(zip(taxaG, listMetric/maxMetric))

    #  calculate mean, median and std of metric
    medianMetric = np.median(listMetric)
    meanMetric = np.mean(listMetric)
    stdMetric = np.std(listMetric)

    # calculate zscore of metric
    zscoreMetric = stats.zscore(listMetric)
    dict_zscoreMetric = dict(zip(taxaG, zscoreMetric))
    
    # calculate if it's keystone a taxa is keystone for some metric
    # if its metric value is higher than the median + 2*std
    iskeystone = list(map(lambda x: int(x > medianMetric + 2*stdMetric), listMetric))
    dict_isKeystone = dict(zip(taxaG, iskeystone))

    # set node attributes
    nx.set_node_attributes(G,values=dict_zscoreMetric,name=metric+"_zScore")
    nx.set_node_attributes(G,values=dict_maxMetric,name=metric+"_scaled")
    nx.set_node_attributes(G,values=dict_isKeystone,name=metric+"_isKeystone")

def _calculateZScoreMax(G,metric,mea,med,std):
    taxaG = list(G.nodes)
    listMetric = np.array([G.nodes[x][metric] for x in taxaG])
    maxMetric = max(listMetric)
    maxMetric = maxMetric if maxMetric > 0 else 1
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
def _table_analysis(base, peak, G, host):
    tmp = host.split('-')
    envi0 = ' '.join('-'.join(tmp[:-1]).split('.'))
    envi1 = ' '.join(tmp[-1].split('.')[:-1])

    metrics = ['LIASP', 'LIASP_dir', 'LIASP_indir', 'path_increase', 'EDpDM','iEDpDM','D','BC','DxCC','Ddiv','DxCCdiv']

    # for mode in ['ED','D','BC']:
    for mode in metrics:
        # getting the mean, median and std of the environment metric
        mea, med, std = getKeystonenessCutts(base,peak,mode)
        identify_keystones_and_add_metric_props(G, mode, mea, med, std) # setting zscore and scaled metric as node attr

    dataN = {'LIASP':[], 'LIASP_dir':[], 'LIASP_indir':[], 'path_increase':[], 'EDpDM':[],'iEDpDM':[],'D':[],'BC':[],'DxCC':[],'Ddiv':[],'DxCCdiv':[]}

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

# creating a csv table to save the desired metric as keystoneness and a Venn diagram for caparision
def table_analysis(base, peak, G, host):
    # get ecosystem and habitat names from "host" string
    tmp = host.split('-')
    ecosystem = ' '.join('-'.join(tmp[:-1]).split('.'))
    habitat = ' '.join(tmp[-1].split('.')[:-1])
    
    # set ecosystem and habitat names as node attributes
    nx.set_node_attributes(G,values=ecosystem, name="Ecosystem")
    nx.set_node_attributes(G,values=habitat, name="Habitat")

    # add Taxon name (node label) as node attribute
    taxa = list(G.nodes(data=False))
    nx.set_node_attributes(G, values=dict(zip(taxa, taxa)), name="Taxon")

    # read the list of node attributes by accessing an arbitrary node
    arbitrary_node = next(iter(G.nodes(data=True)))
    metrics = list(arbitrary_node[1].keys())

    # setting zscore and scaled metric as node attr
    for mode in metrics:
        identify_keystones_and_add_metric_props(G, mode)

    # All relevant metrics are stored as node attributes
    # We can now create a dataframe from the graph
    df = pd.DataFrame.from_dict(dict(G.nodes(data=True)), orient='index')

    # sort the dataframe
    df = dataframe_aesthetics(df)

    # save dataframe to csv
    df.to_csv(base+'/figures/0p%d/keystones.csv'%peak)
    df.to_csv(base+'/keystones.csv')

def save_keystones_table(base, peak, G, host):
    # construct a pandas dataframe with the information on the nodes
    df = pd.DataFrame(dict(G.nodes(data=True)))

    # sort the dataframe items by the node label
    df.sort_values(by=['Taxon'], ascending=True, inplace=True)

    # sort the dataframe
    df = dataframe_aesthetics(df)

    # save dataframe to csv
    df.to_csv(base+'/figures/0p%d/keystones.csv'%peak)

def dataframe_aesthetics(df):
    # sort the dataframe columns
    cols = df.columns.tolist()
    # remove Habitat and Ecosystem from the list of columns
    cols.remove('Habitat')
    cols.remove('Ecosystem')
    # Sort the columns alphabetically
    cols.sort()
    # add Habitat and Ecosystem to the beginning of the list of columns
    cols = ['Ecosystem', 'Habitat'] + cols
    # set the dataframe columns to the sorted list of columns
    df = df[cols]
    # set the dataframe index to the node label
    df.set_index('Taxon', inplace=True)

    return df