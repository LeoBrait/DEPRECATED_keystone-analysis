# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 09:07:32 2018

@author: flavia
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
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

def transform(m):
    if (m=='^'):
        m = 'D'
    elif (m=='>'):
        m = 's'
    else:
        m = 'o'
    return m

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

def plot_rare(axs,G,metric,rare):
    nabundance = get_values(G,"abundance")
    nk = np.max(get_values(G,metric))
    min_abundance = np.min(nabundance)
    max_abundance = np.max(nabundance)
    for ax in axs:
        ax.set_xlim([0.8*min_abundance, 1.2*max_abundance])
        ax.plot([rare,rare],[0,nk*1.1],'-.k',linewidth=0.5)
        ax.set_ylim([0,1.1*nk])
    return min_abundance, max_abundance

def write_index(base,G,tam,metric,rank_metric,prefix,peak):
    nodes = G.nodes()
    f = open(base+"/figures/0p%d/"%peak+metric+"_"+prefix+"_reordered.csv",'w')
    f.write("Rank,Taxon,"+metric+",abundance,Marker\n")
    for el in range(1,tam+1,1):
        taxon = get_node_at(G,rank_metric,el)
        k = nodes[taxon][metric]
        abundance = nodes[taxon]["abundance"]
        mark = nodes[taxon]["marker"]
        f.write("{},{},{:.7e},{:.7e},{}\n".format(el,taxon,k,abundance,mark))

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

def index(base,G,rare,prefix,metric,Taxa,peak):
    rank_metric = "rank_"+metric
    tam = len(G)
    fig, ax = plt.subplots()
    nodes = G.nodes()
    bound = int(tam/4)
    for taxon in Taxa:
        k_i = nodes[taxon][metric]
        abundance_i = nodes[taxon]["abundance"]
        rank_ki = nodes[taxon][rank_metric]
        if (k_i > 0):
            ax.scatter(abundance_i,k_i,edgecolors='w',marker=nodes[taxon]["marker"],s=100,color='#DCDCDC')
            ax.annotate(nodes[taxon][rank_metric], (abundance_i, k_i),color='k',horizontalalignment='center', 			verticalalignment='center',fontweight='black',fontsize=6)
    ax.set_xlabel('Abundance (\%)')
    ax.set_ylabel(metric)
    ax.set_xscale('log')
    nabundance = get_values(G,"abundance")
    nnk = get_values(G,metric) # ADDED, but prob temporarily
    nk = np.max(nnk)
    min_abundance = np.min(nabundance)
    max_abundance = np.max(nabundance)
    plt.xlim([0.8*min_abundance, 1.2*max_abundance])
    plt.plot([rare,rare],[0,nk*1.1],'-.k',linewidth=0.5)
    ax.set_ylim([0,1.1*nk])
    fig.tight_layout()
    fig.savefig(base+"/figures/0p%d/"%peak+metric+"_"+prefix+".pdf",dpi=150)
    plt.close()
    write_index(base,G,tam,metric,rank_metric,prefix,peak)
    # TODO: CORRIGIR vvvvvvvvv
    cdf_indexes(base,nnk,metric,prefix,peak)

def index_split(base,G,rare,prefix,metric,Taxa,peak):
    rank_metric = "rank_"+metric
    tam = len(G)
    fig,axs = plt.subplots(2,1)
    ax = axs[0]
    ax1= axs[1]
    nodes = G.nodes()
    bound = int(tam/4)
    for taxon in Taxa:
        k_i = nodes[taxon][metric]
        abundance_i = nodes[taxon]["abundance"]
        rank_k_i = nodes[taxon][rank_metric]
        marker_i = nodes[taxon]["marker"]
        if (k_i > 0):
            if (marker_i == 'D' or marker_i == 's' or marker_i == 'o'):
                ax.scatter(abundance_i,k_i,edgecolors='w',marker=marker_i,s=100,color='#DCDCDC')
                ax.annotate(nodes[taxon][rank_metric], (abundance_i, k_i),color='k',horizontalalignment='center', 			verticalalignment='center',fontweight='black',fontsize=6)
            else:
                marker_i = transform(marker_i)
                ax1.scatter(abundance_i,k_i,edgecolors='w',marker=marker_i,s=100,color='#DCDCDC')
                ax1.annotate(nodes[taxon][rank_metric], (abundance_i, k_i),color='k',horizontalalignment='center', 			verticalalignment='center',fontweight='black',fontsize=6)
    ax1.set_xlabel('Abundance (\%)')
    ax.set_ylabel(metric)
    ax1.set_ylabel('Candidatus ' + metric)
    ax.set_xscale('log')
    ax1.set_xscale('log')
    min_abundance,max_abundance = plot_rare([ax,ax1],G,metric,rare)
    ax1.plot(-10,1,'o',color = 'gray',label = r'$f < 50\%$')
    ax1.plot(-10,2,'s',color = 'gray',label = r'$50\% < f < 75\%$')
    ax1.plot(-101,3,'D',color = 'gray',label = r'$75\% < f$')
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25),
          frameon=False, ncol=5)
    fig.tight_layout()
    fig.savefig(base+"/figures/0p%d/"%peak+metric+"_"+prefix+"_split.pdf",dpi=150)
    plt.close()

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
    index(base,G,rare,prefix,metric,Taxa,peak)
    index_split(base,G,rare,prefix,metric,Taxa,peak)
    histograms(base,dist,metric,prefix,peak)

    # indirect contribution
    metric = "iEDpDM"
    construct_dic_rank_dict_set_to_G(G,Taxa,dist2,metric)
    index(base,G,rare,prefix,metric,Taxa,peak)
    index_split(base,G,rare,prefix,metric,Taxa,peak)
    histograms(base,dist2,metric,prefix,peak)

    return dist

def histograms(base,k,metric,prefix,peak):
    k = np.array(k)
    plt.figure()
    plt.title(metric)
    args = np.nonzero(k)
    nbins = max(np.sqrt(len(args)),10)
    plt.hist(k[args],bins=nbins)
    plt.tight_layout()
    plt.savefig(base+'/figures/0p%d/'%peak+metric+"_histogram_"+prefix+'.pdf',dpi=150)
    plt.close()

def find_neig_key(base,G,Taxa,metric,prefix,peak):
    nodes = G.nodes()
    tam = len(G)
    bound = int(tam/4)
    nei = open(base+"/figures/0p%d/"%peak+metric+"_"+prefix+"_neig.csv",'w')
    nei.write("rank in "+metric+",taxon,neig\n")
    for i in range(1,bound+1,1):
        taxon = get_node_at(G,"rank_"+metric,i)
        nei.write("{},{},".format(i,taxon))
        par_neig = G.edges(taxon)
        neig = [j for (i,j) in par_neig]
        for el in neig:
            nei.write("{},".format(el))
        nei.write("\n")

def prod_dict(x1,x2,Taxa):
    prod = []
    for taxon in Taxa:
        prod.append(x1[taxon]*x2[taxon])
    return dict(zip(Taxa, prod))

def number_keystones(base,k,metric,prefix,peak):
    X = k.copy()
    X.sort()
    X = X[np.nonzero(X)]
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

def pdf_cdf(data):
    X = data.copy()
    X.sort()
    X = X[np.nonzero(X)]
    kernel = stats.gaussian_kde(X)
    X2 = np.linspace(min(X),max(X),1000)
    Y2 = kernel(X2)
    Y = np.arange(1.,len(X)+1.)  / len(X)
    n = np.flatnonzero(X[1:] == X[:-1])
    X = np.delete(X, n)
    Y = np.delete(Y, n)
    return X, Y,X2,Y2

def cdf_indexes(base,k,metric,prefix,peak):
    k = np.array(k)
    mx, medx, stdx = number_keystones(base,k,metric,prefix,peak)
    plt.figure()
    plt.title(metric)
    x,y,x2,y2 = pdf_cdf(k)
    plt.plot(x,y)
    plt.xlabel(metric)
    plt.ylabel("Cumulative"+metric)
    plt.tight_layout()
    plt.savefig(base+'/figures/0p%d/'%peak+metric+"_cdf_"+prefix+'.pdf',dpi=150)
    plt.close()
    dymin = min(y2)
    dymax = max(y2)
    plt.figure()
    plt.xlabel(metric)
    plt.ylabel("PDF of " +metric)
    plt.plot([mx,mx],[dymin,dymax],'r',label = 'mean')
    plt.plot([mx+stdx,mx+stdx],[dymin,dymax],'r--',label='mean + std')
    plt.plot([mx+2*stdx,mx+2*stdx],[dymin,dymax],'r-.',label = 'mean +  2std')
    plt.plot([medx,medx],[dymin,dymax],'k',label = 'median')
    plt.plot([medx+stdx,medx+stdx],[dymin,dymax],'k--',label = 'median + std')
    plt.plot([medx+2*stdx,medx+2*stdx],[dymin,dymax],'k-.',label = 'median + 2std')
    plt.ylim([dymin,dymax])
    plt.plot(x2,y2,'orange')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(base+'/figures/0p%d/'%peak+metric+"_pdf_"+prefix+'.pdf',dpi=150)
    plt.close()
