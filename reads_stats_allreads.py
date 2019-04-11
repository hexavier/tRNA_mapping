#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 14:34:06 2019

@author: xhernandez
"""
import sys
import os
import pandas as pd
import re
import numpy as np
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns

def transformdata(data,transf):
    data = data.iloc[0:66,:]
    if transf=="log":
        outdata = data.apply(np.log)
        # Remove inf values
        outdata.replace([np.inf, -np.inf], np.nan,inplace=True)
    elif transf=="arcsinh":
        outdata = data.apply(np.arcsinh)
    elif transf=="rel":
        data = data.iloc[0:61,:]
        # Compute relative data
        outdata = pd.DataFrame(columns=data.columns,index=data.index[0:61])
        aa = list(set([s[0:3] for s in outdata.index]))
        for n in aa:
            idx = [n in s for s in [l[0:3] for l in data.index]]
            total = data.loc[idx,:].sum()
            outdata.loc[data.index[idx],:] = data.loc[data.index[idx],:]/total
            iszero = (total==0)
            if any(iszero):
                outdata.loc[data.index[idx],iszero] = 1.0/sum(idx)
        outdata.iloc[:,:] = np.float64(outdata)
    else:
        outdata=data
        
    return outdata

#%% Load datasets

# Cluster mapping
text_file = open(sys.argv[1], "r")
#text_file = open("../../../Data/Genomes/H.sapiens/hg38.tRNAscan_clusterInfo.fa", "r")
clusters = text_file.read().split('\n')
text_file.close()
seqs = clusters[1::2]
clusters = clusters[0:-1:2]

# tRNAs
text_file = open(sys.argv[2], "r")
#text_file = open("../../../Data/Genomes/H.sapiens/tRNAs.txt", "r")
trnas = text_file.read().split('\n')
text_file.close()
trnas.remove("")

# Load RPM data
trna = pd.read_csv(sys.argv[3],index_col=0)
#trna = pd.read_csv("results_nomod.csv",index_col=0)

#Map cluster id to trnas
cluster_seq = dict(zip([re.findall("cluster[0-9]+",s)[0] for s in clusters],seqs))
mapping = {re.findall("cluster[0-9]+",s)[0]: re.findall("(?<=-)i?[A-Z]{1}[a-zC]+[ACTGN]{3}(?=\()",s)[0] for s in clusters}

# Load datasets
ids = [ name for name in os.listdir(".") if os.path.isdir(name) ]
ids_decap = [s[:34] for s in ids]


#%% Map samples
metadata = pd.read_csv(sys.argv[4]).loc[:,["ID","Source"]]
#metadata = pd.read_csv("../../../Data/NGS/TCGA_pilot/metadata.txt").loc[:,["ID","Source"]]
metadata.drop_duplicates(inplace=True)
metadata.set_index("ID",inplace=True)
metadata = metadata.loc[ids_decap,:]
id_names = list(set([s1+"-"+s2 if sum([s2 in n for n in ids])==1 else s1+"-"+s2+"1" for s1,s2 in zip(metadata.loc[:,"Source"],ids)]))
id_names.extend([s[:-1]+"2" for s in id_names if s[-7:]=="_mirna1"])

# Initialize dataframe
rlen = pd.DataFrame(columns=id_names, index=list(range(1,(max([len(s) for s in seqs])+1))))
cov = pd.DataFrame(columns=id_names, index=trnas)
depth = pd.DataFrame(columns=id_names, index=trnas)

#%% Extract data from datasets
clustlist = list(mapping.keys())
ids=np.array(ids)
rangecov_temp = pd.DataFrame(columns=id_names, index=clustlist)
depthcov_temp = pd.DataFrame(columns=id_names, index=clustlist)
depthperbase_temp = dict.fromkeys(id_names)
listdepth = []
for n,d in enumerate(id_names):
    fileid = ids[np.bool_([s in d for s in ids])][0]
    # Coverage
    cov_filename = os.path.join(fileid,str("%s.depth.allreads.txt" % fileid))
    coverage = pd.read_table(cov_filename,header=None)
    # Initialize structure of perbase computation
    norm_clust_depth = pd.DataFrame(columns=clustlist,index=np.arange(0, 1.05, 0.05))
    for c in [s for s in list(set(coverage.iloc[:,0]))]:
        cov_c = coverage.loc[(coverage.iloc[:,0]==c),:]
        rangecov_temp.loc[c,d]= cov_c.shape[0]/len(cluster_seq[c])
        depthcov_temp.loc[c,d]= cov_c.sum(0)[2]
        normrange = cov_c.iloc[:,1]/len(cluster_seq[c])
        norm_clust_depth.loc[:,c]=[np.sum(cov_c.loc[np.logical_and(normrange<=i,normrange>i-0.05),2]) for i in norm_clust_depth.index]
        
    # Merge clusters into trnas
    norm_trna_depth = pd.DataFrame(columns=trnas,index=np.arange(0, 1.05, 0.05))
    for t in trnas:
        trna_clust = [c for c in mapping.keys() if mapping[c]==t]
        if depthcov_temp.loc[trna_clust,d].sum(0)>0:
            cov.loc[t,d]=(rangecov_temp.loc[trna_clust,d]*(depthcov_temp.loc[trna_clust,d]/
                   depthcov_temp.loc[trna_clust,d].sum(0))).sum(0)
            norm_trna_depth.loc[:,t] = norm_clust_depth.loc[:,trna_clust].sum(1)
    depthperbase_temp[d] = norm_trna_depth
    listdepth.append(norm_trna_depth)
    
    # Read lengths
    rlen_filename = os.path.join(fileid,str("%s.readlengths.allreads.txt" % fileid))
    readlength = pd.read_table(rlen_filename,header=None)
    rlen.loc[readlength.iloc[:,0].values,d] = readlength.loc[:,1].values/readlength.loc[:,1].sum(0)

depth = pd.concat(listdepth, keys=depthperbase_temp)

#%%######################### Plotting #################################
pdf = matplotlib.backends.backend_pdf.PdfPages("reads_stats.allreads.pdf")
# Colors
rlen_source = [re.findall("[A-Za-z0-9_+]+(?=-)",s)[0] for s in rlen.columns]
cov_source = [re.findall("[A-Za-z0-9_+]+(?=-)",s)[0] for s in cov.columns]
labels = list(set(rlen_source))
cmap = plt.cm.get_cmap('jet') # set color map
colors = cmap(np.arange(len(labels))/(len(labels)-1))
handles = [mlines.Line2D([], [], marker="o", alpha=0.4, color=c) for c in colors]
#%% Plot read length
fig = plt.figure(figsize = (14,8))
ax = plt.subplot()
for label, color in zip(labels,colors):
    idx = np.bool_([label==s for s in rlen_source])
    ax.plot(rlen.index,rlen.loc[:,idx],c=color,alpha=0.4)

ax.set_xlabel('Read length')
ax.set_ylabel('Fraction of reads')
fig.legend(handles,labels,loc="center right")
pdf.savefig(fig, bbox_inches="tight")
#%% Plot coverage
fig = plt.figure(figsize = (14,8))
ax = plt.subplot()
dotdf = pd.DataFrame(columns=["Fraction of covered length","source","tRNA"])
for label in labels:
    idx = np.bool_([label==s for s in cov_source])
    source_df = cov.loc[:,idx]
    for col in source_df.columns:
        dotdf = dotdf.append(pd.DataFrame({"Fraction of covered length":source_df.loc[:,col],"source":[label]*source_df.shape[0],
                                   "tRNA":source_df.index}))
sns.stripplot(x='tRNA', y='Fraction of covered length', hue="source", data=dotdf, 
              palette=dict(zip(labels,colors)), jitter=True, edgecolor='none', alpha=0.40)
ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
ax.legend_.remove()
fig.legend(handles,labels,loc="center right")
pdf.savefig(fig, bbox_inches="tight")
#%% Plot coverage vs RPMs
rpm = transformdata(trna,"arcsinh")
fig = plt.figure(figsize = (14,8))
ax = plt.subplot()
for label, color in zip(labels,colors):
    covidx = np.bool_([label==s for s in cov_source])
    rpmidx = cov.columns[covidx]
    for txt in cov.index:
        ax.scatter(rpm.loc[txt,rpmidx],cov.loc[txt,covidx],alpha=0)
        for i in range(sum(covidx)):
            coord = [rpm.loc[txt,rpmidx].values[i], cov.loc[txt,covidx].values[i]]
            if not any(np.isnan(coord)):
                ax.annotate(txt,coord, color=color, ha='center', va='center', fontsize=9)
ax.set_xlabel('arcsinh(RPM)')
ax.set_ylabel('Fraction of covered length')
fig.legend(handles,labels,loc="center right")
pdf.savefig(fig, bbox_inches="tight")
#%% Plot depth across tRNA length
for t in trnas:
    fig = plt.figure(figsize = (14,8))
    ax = plt.subplot()
    for label, color in zip(labels,colors):
        idx = rlen.columns[np.bool_([label==s for s in rlen_source])]
        ax.plot(np.arange(0, 1.05, 0.05),depth.loc[idx,t].unstack(level=0),c=color,alpha=0.4)
    
    ax.set_xlabel('Relative position in tRNA sequence')
    ax.set_ylabel('Number of reads')
    ax.set_title(t)
    fig.legend(handles,labels,loc="center right")
    pdf.savefig(fig, bbox_inches="tight")

#%% Write output
cov.to_csv("coverage.allreads.csv")
rlen.to_csv("readlength.allreads..csv")
pdf.close()
