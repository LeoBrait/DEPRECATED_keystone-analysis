from xonsh_py import *
import concurrent.futures as cf
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Arial'
import threading as thread

from sys import exit

def main():
    sp = '/sparcc/sparcc_data'
    for env in lsgrep('output/transposed',['']):
        for iter in lsgrep(env+sp,['boot']):
            if existOldData(iter+'/001/cor.tsv') and existOldData(iter+'/002/cor.tsv') and \
                existOldData(iter+'/001/pval.tsv') and existOldData(iter+'/002/pval.tsv'):

                hard_work(iter)

            else:
                print('There is no data sufficient for '+iter)
                continue

def hard_work(iter):
    print('Starting: '+iter)

    cor1 = pd.read_csv(iter+'/001/cor.tsv', sep='\t', index_col=0)
    cor2 = pd.read_csv(iter+'/002/cor.tsv', sep='\t', index_col=0)

    pvl1 = pd.read_csv(iter+'/001/pval.tsv', sep='\t', index_col=0)
    pvl2 = pd.read_csv(iter+'/002/pval.tsv', sep='\t', index_col=0)

    cor1[pvl1 > .01] = 0. # filtering cells where pval smaller than 1%
    cor2[pvl2 > .01] = 0. # filtering cells where pval smaller than 1%

    cor1.to_csv(iter.replace('transposed','boot_divergence')+'/filtered_001.tsv',sep='\t')
    cor2.to_csv(iter.replace('transposed','boot_divergence')+'/filtered_002.tsv',sep='\t')

    print('Finishing: '+iter)

main()