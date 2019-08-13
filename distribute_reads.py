#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 11:30:38 2019

@author: xhernandez
"""

import sys
import re


#%% Load clusters
# Cluster info
text_file = open(sys.argv[2], "r")
#text_file = open("../../../../Data/Genomes/H.sapiens/hg38.tRNAscan_clusterInfo.fa", "r")
clusters = text_file.read().split('\n')
text_file.close()
seqs = clusters[1::2]
clusters = clusters[0:-1:2]
mapping = {re.findall("cluster[0-9]+",s)[0]: re.findall("(?<=-)i?[A-Z]{1}[a-zC]+[ACTGN]{3}(?=\()",s)[0] for s in clusters}

#%% Loop through sam file

out_file = open(sys.argv[3], "w")
multimappers = open(sys.argv[4], "w")

ambmap = []
r = str()
lines = []
# Sorted SAM file
#with open("LA1_31413_CAGATC_dist.nSorted.sam", "r") as in_file:
with open(sys.argv[1], "r") as in_file:
    for line in in_file:
        # Print headers
        if line[0]=="@":
            out_file.write(line)
        # Print unique reads
        elif "NH:i:1\t" in line:
            out_file.write(line)
        # Print reads that map to mutiple clusters but sharing the same isoacceptor
        else:
            l = line.split('\t')
            if l[0] == r:
                ambmap.append(mapping[l[2]])
                r = l[0]
                outline = line
                lines.append(line)
            else:
                if len(set(ambmap))==1:
                    out_file.write(outline)
                elif len(set(ambmap))>1:
                    for item in lines:
                        multimappers.write(item)
                ambmap = [mapping[l[2]]]
                r = l[0]
                lines = [line]
            
        
out_file.close()
multimappers.close()
