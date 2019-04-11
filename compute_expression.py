#!/usr/bin/env python

import sys
import pandas as pd
import re

#%% Upload expression data, ASE, trnas and cluster mapping
expression = pd.read_table(sys.argv[1],header=None, index_col=0)
#expression = pd.read_table("TCGA-ZU-A8S4-11A-11R-A41D-13_mirna.expression.txt",header=None, index_col=0)
expression.drop("*",inplace=True)
expression = expression[expression.loc[:,2] != 0]

ase = pd.read_table(sys.argv[2],index_col="contig")
#ase = pd.read_table("TCGA-ZU-A8S4-11A-11R-A41D-13_mirna.ASE.csv",index_col="contig")

text_file = open(sys.argv[3], "r")
#text_file = open("Data/Genomes/H.sapiens/hg38.tRNAscan_clusterInfo.fa", "r")
clusters = text_file.read().split('\n')
text_file.close()
seqs = clusters[1::2]
clusters = clusters[0:-1:2]

text_file = open(sys.argv[4], "r")
#text_file = open("Data/Genomes/H.sapiens/tRNAs.txt", "r")
trnas = text_file.read().split('\n')
text_file.close()

#%% Compute trna expression
#Map cluster id to trnas
cluster_seq = dict(zip([re.findall("cluster[0-9]+",s)[0] for s in clusters],seqs))
mapping = {re.findall("cluster[0-9]+",s)[0]: re.findall("(?<=-)i?[A-Z]{1}[a-zC]+[ACTGN]{3}(?=\()",s)[0] for s in clusters}

# Initialize structures
rpm = pd.DataFrame(0,index = trnas, columns=['rpm','rpm_nomod'])
tot_reads = sum(expression.loc[:,2])

# Analyze expression of clusters
for c in expression.index:
    rpm.loc[mapping[c],"rpm_nomod"] += expression.loc[c,2]*1000000/tot_reads
    if c not in ase.index:
        rpm.loc[mapping[c],"rpm"] += expression.loc[c,2]*1000000/tot_reads
    elif mapping[c]=="UndetNNN":
        rpm.loc[mapping[c],"rpm"] += expression.loc[c,2]*1000000/tot_reads
    else:
        # Localize anticodon
        acod = re.findall(mapping[c][-3:],cluster_seq[c],re.I)
        pos = re.search(mapping[c][-3:],cluster_seq[c],re.I).start()+1
        trim = 0; end5 = 0; end3 = 1
        while len(acod)!=1:
            acod = re.findall(mapping[c][-3:],cluster_seq[c][(trim+end5):-(trim+end3)],re.I)
            pos = re.search(mapping[c][-3:],cluster_seq[c][(trim+end5):-(trim+end3)],re.I).start()+trim+1
            if end3==0:
                end3 += 1
            elif end5==0:
                end5 += 1
            else:
                trim += 1; end5 = 0; end3 = 0
        
        # Get modifications in anticodon
        ase_c=ase.loc[[c],:]
        ase_acod = ase_c.loc[[p in range(pos,pos+3) for p in ase_c.position],]
        
        if ase_acod.empty:
            # Modification does not affect anticodon loop
            rpm.loc[mapping[c],"rpm"] += expression.loc[c,2]*1000000/tot_reads
        else:
            rpms = expression.loc[c,2]*1000000/tot_reads
            
            # Calculate expression of alternative
            alt_rpms = rpms*(sum(ase_acod.loc[:,"altCount"])/
                             sum(ase_acod.loc[:,"rawDepth"]))
            new_acod = list(acod[0])
            for p in ase_acod.position:
                # Determine new trna
                acod_pos = p-pos # will be 0,1 or 2 depending of the position
                new_acod[acod_pos] = ase_acod.loc[ase_acod.position==p,"altAllele"][0]
            new_acod = ''.join(new_acod).upper()
            new_trna = [s for s in trnas if new_acod in s]
            if new_trna:
                new_trna = new_trna[0]
            else:
                new_trna = "UndetNNN"
                
            # Add expression to new trna
            rpm.loc[new_trna,"rpm"] += alt_rpms
                
            # Subtract alternative expression from reference
            rpm.loc[mapping[c],"rpm"] += (rpms - alt_rpms)

rpm.iloc[66,:] = tot_reads; rpm.index = list(rpm.index[:66])+["READS"]

#%% Write results
rpm.to_csv(sys.argv[5])
