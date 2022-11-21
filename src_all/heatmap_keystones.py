#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

fontsiz=3.75
fontszleg=4
plt.rcParams['font.family'] = ['Arial']

experimentalTaxa = ['Chromatium okenii', 'Chromatiaceae', 'Desulfosporosinus', 'Peptococcaceae']

# colors = [
# '#a6cee3', # nonIdentified
# '#1f78b4', # A
# '#b2df8a', # B
# '#33a02c', # AB
# '#fb9a99', # C
# '#e31a1c', # AC
# '#fdbf6f', # BC
# '#ff7f00', # ABC
# '#b15928'] # inexistent

colors = [
[1.,1.,1.], # nonIdentified
[0.12156863, 0.47058824, 0.70588235], # A
[0.69803922, 0.8745098 , 0.54117647], # B
[0.2       , 0.62745098, 0.17254902], # AB
[0.98431373, 0.60392157, 0.6       ], # C
[0.89019608, 0.10196078, 0.10980392], # AC
[0.99215686, 0.74901961, 0.43529412], # BC
[1.        , 0.49803922, 0.        ], # ABC
[0.85, 0.85, 0.85]] # inexisten

multi  = np.array([1,2,4])

def mapColor(val):
    return colors[val]

level = sys.argv[1]

if sys.argv[2] == '1':
    metrics = ['EDpDM', 'BC', 'Ddiv'] # for full contribution
    indirect = False
elif sys.argv[2]== '2':
    metrics = ['iEDpDM', 'BC', 'Ddiv'] # for indirect contribution
    indirect = True
else:
    print('Wrong liasp option as arg2')
    sys.exit()

df = pd.read_csv('output/transposed_all_environments/%s/keystones.csv'%level)

# dropping to improve performance
df.drop(columns=set(df.keys())-set(['Ecosystem','Habitat','Taxon','Abundance']+[i+'_isKeystone' for i in metrics]),inplace=True)

# making one column index
df['Index'] = df['Ecosystem']+'@'+df['Habitat']
df.set_index('Index', inplace=True)

allEnvi = np.array(sorted(set(df.index)))
allTaxn = np.array(sorted(set(df['Taxon'])))

candidatus = np.array(sorted([x for x in allTaxn if 'Candidatus' in x or 'candidate' in x]))
bonafide   = np.array(sorted(set(allTaxn) - set(candidatus)))

groupedTaxa = [bonafide, candidatus]

def getHabitatGroup(allEnvi):
    dic = []; ret = []
    for i in allEnvi:
        eco,hab = i.split('@')
        if not eco in dic:
            dic.append(eco)
            ret.append([i])
        else:
            ret[dic.index(eco)].append(i)

    return ret

groupedEnvi = getHabitatGroup(allEnvi)

X = {}

# getting data as dict
for i in allEnvi:
    data = df.loc[i].copy().set_index('Taxon')
    X[i] = {j:{k:int(data.loc[j,k+'_isKeystone']) for k in metrics} for j in data.index}

# getting all medians abundance as a dict
A = {}
df.set_index('Taxon', inplace=True)
df.drop(columns=set(df.columns)-set(['Taxon','Abundance']),inplace=True)
for i in allTaxn:
    A[i] = np.mean(df.loc[i].values)

# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------

def ranking(val):
    if val in [1,2,4]:
        return 1
    elif val in [3,5,6]:
        return 2
    elif val == 7:
        return 3
    elif val == -1:
        return -0.0000001 # trying to put the grays at the end without penalty to the ranking procedure
    else: # when val is 0
        return 0

def make(allEnvi, allTaxn, X, M):
    for ii,i in enumerate(allEnvi):
        for jj,j in enumerate(allTaxn):
            if j in X[i]:
                val = np.array([X[i][j][k] for k in metrics])
                M[ii,jj] = sum(val*multi)
            else:
                M[ii,jj] = -1 # inexistent

# --------

def paint(Z,M):
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            Z[i,j] = mapColor(int(M[i,j]))

# ---------

def renameTaxa(vec,level):
    dic = {'Unclassified %s'%level:'U.%c'%level[0], 'Candidatus':'Ca.', 'candidate':'cn.'}
    for i in range(len(vec)):
        for j in dic:
            vec[i] = vec[i].replace(j,dic[j])

def renameEnvi(vec,level):
    ret = vec[0].split('@')[0]
    for i in range(len(vec)):
        vec[i] = vec[i].replace('@',' ')
        # vec[i] = vec[i].split('@')[1]
    return ret

# returning the colors for x-Axis labels
def highlightLiterature(vec,lit):
    ret = []
    for i in vec:
        o = len(ret)
        for j in lit:
            if i in j or j in i:
                ret.append('r')
                break
        if o >= len(ret):
            ret.append('k')
    return ret

# ---------
def calculateAccuracy(mat,out):
    for i in range(len(out)):
        out[i] = sum((mat == i).flatten())

# 0 none, ed 1 , bc 2, d 4
# 0 -> allone, 1 -> bc, 2 -> D
# concordance with other metrics
def calculateConcordance(mat,out1,out2):
    for i in range(mat.shape[1]):
        p = np.sum(mat[:,i] == 1) # how many times it is only ED key
        p1 = np.sum((mat[:,i] == 1) | (mat[:,i] == 0)) # how many times it's not BC and not D

        q = np.sum((mat[:,i] == 3) | (mat[:,i] == 7)) # how many times it's ED and BC
        q1 = np.sum((mat[:,i] == 2) | (mat[:,i] == 6) | (mat[:,i] == 3) | (mat[:,i] == 7)) # how many times it's BC

        r = np.sum((mat[:,i] == 5) | (mat[:,i] == 7)) # how many times it's ED and D
        r1 = np.sum((mat[:,i] == 4) | (mat[:,i] == 6) | (mat[:,i] == 5) | (mat[:,i] == 7)) # how many times it's D

        out1[0,i] += p; out2[0,i] += p1
        out1[1,i] += q; out2[1,i] += q1
        out1[2,i] += r; out2[2,i] += r1

def calculateOccurrence(mat,outC,outR):
    for i in range(mat.shape[1]):
        outC[0,i] += np.sum((mat[:,i] == 1) | (mat[:,i] == 3) | (mat[:,i] == 5) | (mat[:,i] == 7)) # counting ED
        outC[1,i] += np.sum((mat[:,i] == 2) | (mat[:,i] == 3) | (mat[:,i] == 6) | (mat[:,i] == 7)) # counting BC
        outC[2,i] += np.sum((mat[:,i] == 4) | (mat[:,i] == 5) | (mat[:,i] == 6) | (mat[:,i] == 7)) # counting D

    for i in range(mat.shape[0]):
        outR[0,i] += np.sum((mat[i,:] == 1) | (mat[i,:] == 3) | (mat[i,:] == 5) | (mat[i,:] == 7)) # counting ED
        outR[1,i] += np.sum((mat[i,:] == 2) | (mat[i,:] == 3) | (mat[i,:] == 6) | (mat[i,:] == 7)) # counting BC
        outR[2,i] += np.sum((mat[i,:] == 4) | (mat[i,:] == 5) | (mat[i,:] == 6) | (mat[i,:] == 7)) # counting D
# ---------

stats = np.zeros(8)

maxn = len(allTaxn)
maxm = len(allEnvi)

MAT = np.zeros((maxm,maxn))

statsTax = np.zeros((2,3,maxn))
statsOcuTax = np.zeros((3,maxn))
statsOcuEco = np.zeros((3,maxm))

def makeData():
    cumM = 0
    for i, myEnvi_ in enumerate(groupedEnvi):
        myEnvi = np.array(myEnvi_)
        # getting matrix shape
        m = len(myEnvi)

        cumN = 0
        for j, myTaxn_ in enumerate(groupedTaxa):
            myTaxn = np.array(myTaxn_)

            # getting matrix shape
            n = len(myTaxn)

            # unicell matrix to hold the i,j value
            M = np.zeros(shape=(m,n))

            # computing all cells positions and axes ranks
            make(myEnvi, myTaxn, X, M)

            # inserting actual data into grand matrix MAT
            MAT[cumM:cumM+m,cumN:cumN+n] = M

            cumN += n
        cumM += m

makeData()

calculateAccuracy(MAT,stats)

# ----- we've got the total data ------

# ranking vertically (regarding with the groups)
def rankEcosys():
    cumM = 0
    l = 0
    for myEnvi in groupedEnvi:
        m = len(myEnvi)

        rankEnv = np.zeros(m)

        for i,k in enumerate(range(cumM,cumM+m)):
            rankEnv[i] = sum([ranking(j) for j in MAT[k,:]])

        indEnv = np.argsort(rankEnv)[::-1]

        # reordering step
        myEnvi = np.array(myEnvi)[indEnv]
        groupedEnvi[l] = myEnvi

        indEnv += cumM # pointing the indexes to the correct location at MAT

        MAT[cumM:cumM+m,:] = MAT[indEnv,:]

        cumM += m
        l += 1

# ranking horizontally
def rankTaxa():
    ret = np.zeros(len(allTaxn))
    cumN=0
    l=0
    for myTaxn in groupedTaxa:
        n = len(myTaxn)

        rankTax = np.zeros(n)

        for i,k in enumerate(range(cumN,cumN+n)):
            rankTax[i] = sum([ranking(j) for j in MAT[:,k]])

        indTax = np.argsort(rankTax)[::-1]

        # reordering step
        myTaxn = np.array(myTaxn)[indTax]
        groupedTaxa[l] = myTaxn
        ret[cumN:cumN+n] = [A[x] for x in groupedTaxa[l]]

        indTax += cumN # pointing the indexes to the correct location at MAT

        MAT[:,cumN:cumN+n] = MAT[:,indTax]

        cumN += n
        l += 1
    return ret

rankEcosys()
abund = rankTaxa()

# ----- we've got the total data reordered

calculateConcordance(MAT,statsTax[0],statsTax[1])
calculateOccurrence(MAT,statsOcuTax,statsOcuEco)
# reordering the statsOcuEco to fit the desired pattern in plot
l=[];c=0
for i in groupedEnvi:
    l.append(list(range(c,len(i)+c)))
    c+=len(i)
statsOcuEco[:,:] = statsOcuEco[:,np.concatenate(l[::-1])]

width = maxn*0.075
height = maxm*0.075
fig = plt.figure(0,figsize=(width,height))

propy = .60
propx = .60

# The dimensions [left, bottom, width, height] of the new axes. All quantities are in fractions of figure width and height.

cumM = 0
cumuY=.0
for i, myEnvi_ in enumerate(groupedEnvi):
    cumN = 0
    cumuX=.15
    for j, myTaxn_ in enumerate(groupedTaxa):
        myEnvi = np.array(myEnvi_)
        myTaxn = np.array(myTaxn_)

        # getting matrix shape
        m,n = len(myEnvi), len(myTaxn)

        ax = fig.add_axes([cumuX, cumuY, propx*n/maxn, propy* m/maxm])

        # unicell matrix to hold the i,j value
        M = MAT[cumM:cumM+m,cumN:cumN+n]

        # converting data to numpy matrix (we'll fill with rgb values)
        Z = np.zeros(shape=(m,n,3))
        paint(Z,M)
        im = ax.imshow(Z)

        # We want to show all ticks...
        if j == 0:
            eco = renameEnvi(myEnvi,level)
            ax.set_yticks(np.arange(len(myEnvi)))
            ax.set_yticklabels(myEnvi,fontsize=fontsiz)
        else:
            ax.set_yticks([])
            ax.set_yticklabels([])

        if i == 0:
            renameTaxa(myTaxn,level)
            labelColors = highlightLiterature(myTaxn,experimentalTaxa)
            ax.set_xticks(np.arange(len(myTaxn)))
            labels = ax.set_xticklabels(myTaxn,fontsize=fontsiz)
            for u,uu in zip(labels,labelColors):
                u.set_color(uu)
        else:
            ax.set_xticks([])
            ax.set_xticklabels([])

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=90,  va='center', ha="right", rotation_mode="anchor")
        plt.setp(ax.get_yticklabels(), rotation=180, va='center', ha="left", rotation_mode="anchor")

        # smoothing spines and ticks color
        plt.setp(ax.spines.values(),color='#939393',linewidth=.25)
        plt.setp([ax.get_xticklines(),ax.get_yticklines()],color='#939393')

        cumN += n
        cumuX += propx*n/maxn + .00

    cumM += m
    cumuY += propy*m/maxm + .00

fig1 = plt.figure(1)
s = sum(stats)
metrics = ['LIASP','Betweenness Centrality','Degree'] # renaming labels
labels = ['Non keystone (#%d %.2f%%)'%(stats[0],100.*stats[0]/s),
    metrics[0] + ' (#%d %.2f%%)'%(stats[1],100.*stats[1]/s),
    metrics[1] + ' (#%d %.2f%%)'%(stats[2],100.*stats[2]/s),
    metrics[0]+' & '+metrics[1] + ' (#%d %.2f%%)'%(stats[3],100.*stats[3]/s),
    metrics[2] + ' (#%d %.2f%%)'%(stats[4],100.*stats[4]/s),
    metrics[0]+' & '+metrics[2] + ' (#%d %.2f%%)'%(stats[5],100.*stats[5]/s),
    metrics[1]+' & '+metrics[2] + ' (#%d %.2f%%)'%(stats[6],100.*stats[6]/s),
    metrics[0]+' & '+metrics[1]+' & '+metrics[2] + ' (#%d %.2f%%)'%(stats[7],100.*stats[7]/s),
    'Absent'+ ' (#%d)'%(maxn*maxm-s)]

patches = [mpatches.Patch(facecolor=colors[i], label=labels[i], edgecolor='k', linewidth=.1) for i in [1,2,4,3,5,6,-1,0,7]]
leg = fig1.legend(handles=patches, frameon=False, fontsize=fontszleg, ncol=3)

indirects = '_indirect' if indirect else ''
fig1.savefig('output/transposed_all_environments/%s/legend%s.png'%(level,indirects),dpi=1000,bbox_inches='tight', pad_inches=.01)
fig1.savefig('output/transposed_all_environments/%s/legend%s.svg'%(level,indirects),dpi=300,bbox_inches='tight', pad_inches=.01)
fig1.savefig('output/transposed_all_environments/%s/legend%s.pdf'%(level,indirects),dpi=300,bbox_inches='tight', pad_inches=.01)

# --------------------------------

# plotting the Taxa abundance above the graph# for the next plot, getting the Abundances as a matrix
cumuX=.15 # repointing to left border
ax = fig.add_axes([cumuX, cumuY+.01, propx, propy*1./maxm])
im = ax.imshow([np.log10(abund)],cmap=plt.cm.Purples)

ax.set_yticks(np.arange(1))
ax.set_yticklabels(['Log10 of Mean Abundance through habitats'], fontsize=fontsiz,
    rotation=180, va='center', ha="left", rotation_mode="anchor")
ax.set_xticks([])
ax.set_xticklabels([])
# smoothing spines and ticks color
plt.setp(ax.spines.values(),color='#939393',linewidth=.25)
plt.setp([ax.get_xticklines(),ax.get_yticklines()],color='#939393')
cumuY += propy*1./maxm + .01

# --------------------------------

# plotting the marginal distributions (above the "mosaic")
cumuX=.15 # repointing to left border
ax = fig.add_axes([cumuX, cumuY+.01, propx, propy*8./maxm])
bins = np.arange(0, maxn)
lasti = 0
ticks_pos = []
for i,clr in zip(statsOcuTax[::-1],np.array(colors)[[1,2,4]][::-1]):
    ax.bar(bins,i,bottom=lasti,align='edge',linewidth=0,width=1,color=clr,alpha=0.6)
    ticks_pos.append(lasti+float(max(i))/2)
    lasti += max(i)
ax.set_yticks(ticks_pos)
ax.set_yticklabels(['Sum of '+x for x in metrics[::-1]], # setting inversed order labels as They're rotated
    fontsize=fontsiz, rotation=180, va='center', ha='left', rotation_mode='anchor')
# ax.set_ylabel('Sum of keystones', fontsize=fontsiz, rotation=180, va='center', ha='left', rotation_mode='anchor')
ax.set_xticks([])
ax.set_xlim([0,maxn])
plt.setp(ax.spines.values(),color='#939393',linewidth=0)
plt.setp([ax.get_xticklines(),ax.get_yticklines()],color='#939393')

# plotting the marginal distributions (aside the "mosaic")
cumuX=.15 + propx + .005 # Repositioning the left
cumuY=.0
ax = fig.add_axes([cumuX, cumuY, propx*8./maxn, propy]) # plotting a 'rotated' axis
bins = np.arange(0, maxm)[::-1]
lasti = 0
ticks_pos = []
for i,clr in zip(statsOcuEco[::-1],np.array(colors)[[1,2,4]][::-1]):
    ax.barh(bins,i,left=lasti,height=1,edgecolor=(.7,.7,.7,.1),linewidth=.05,facecolor=clr,alpha=0.6,align='edge')
    ticks_pos.append(lasti+float(max(i))/2)
    lasti += max(i)
ax.set_yticks([])
ax.set_ylim([0,maxm])
ax.set_xticks(ticks_pos)
ax.set_xticklabels(['Sum of '+x for x in metrics[::-1]], # setting inversed order labels to match the bar colors
    fontsize=fontsiz, rotation=90, va='center', ha='right', rotation_mode='anchor')
# ax.set_xlabel('# of keystones', fontsize=fontsiz, rotation=90)
plt.setp(ax.spines.values(),color='#939393',linewidth=0)
plt.setp([ax.get_xticklines(),ax.get_yticklines()],color='#939393')

# --------------------------------

fig.savefig('output/transposed_all_environments/%s/liasp%s_validation.png'%(level,indirects),dpi=1000,bbox_inches='tight', pad_inches=.01)
fig.savefig('output/transposed_all_environments/%s/liasp%s_validation.svg'%(level,indirects),dpi=300,bbox_inches='tight', pad_inches=.01)
fig.savefig('output/transposed_all_environments/%s/liasp%s_validation.pdf'%(level,indirects),dpi=300,bbox_inches='tight', pad_inches=.01)
