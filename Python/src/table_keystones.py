import numpy as np
import networkx as nx
import pandas as pd
from scipy import stats

# rare = [0.1,1,5]

# creating a csv table to save the desired metric as keystoneness and a Venn diagram for caparision
def table_analysis(base, peak, G, host):
    """
    Create a csv table with the keystone information that is stored on the nodes of the graph G.

    Parameters
    ----------
    base : string
        Path to the base directory of the ecosystem.

    peak : string
        Peak of the Critical Network to be considered.

    G : networkx graph
        Graph of the ecosystem. Modified in place.

    host : string
        Name of the ecosystem and habitat.
    """
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
    for metric in metrics:
        identify_keystones_and_add_metric_props(G, metric)

    # All relevant metrics are stored as node attributes
    # We can now create a dataframe from the graph
    df = pd.DataFrame.from_dict(dict(G.nodes(data=True)), orient='index')

    # sort the dataframe
    df = dataframe_aesthetics(df)

    # save dataframe to csv
    df.to_csv(base+'/figures/0p%d/keystones.csv'%peak)
    df.to_csv(base+'/keystones.csv')

def identify_keystones_and_add_metric_props(G, metric):
    """
    Identify keystone taxa and add metric properties to the graph.

    Parameters
    ----------
    G : networkx graph
        Graph of the ecosystem. Modified in place.

    metric : string
        Name of the metric to be added to the graph

    Returns
    -------
    None.
    """
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

    # calculate mean, median and std of metric
    medianMetric = np.median(listMetric)
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


def dataframe_aesthetics(df):
    """
    Sort the dataframe columns and set the dataframe index to the node label.

    Parameters
    ----------
    df : pandas dataframe
        Dataframe with the keystone information.

    Returns
    -------
    df : pandas dataframe
        Updated dataframe with the keystone information.
    """
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