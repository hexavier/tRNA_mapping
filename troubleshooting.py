#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 08:41:10 2019

@author: xhernandez
"""

import os
import re

# List datasets
ids = os.listdir(".")

# Distinguish downloaded from fastq.gz
bams = [n[0] for n in [re.findall("[a-z0-9]+-[a-z0-9]+-[a-z0-9]+-[a-z0-9]+-[a-z0-9]+",s) for s in ids] if n]
fastq = [n[0] for n in [re.findall("TCGA-[A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]+_mirna[AB]*",s) for s in ids] if n]

# Find missing fastq
missing = []
while bams:
    b = bams[0]
    files = os.listdir(b)
    tcga = [re.findall("TCGA-[A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]+_mirna[AB]*",s)[0] for s in files if s[-3:]=="bam"][0]
    infastq = [tcga in s for s in fastq]
    idxfastq = [n for n in range(len(infastq)) if infastq[n]]
    if not idxfastq:
        missing.append(tcga)
        bams.remove(b)
    else:
        bams.remove(b)
        fastq.remove(fastq[idxfastq[0]])

missing= list(set(missing))

#%% Write output
with open('troubleshoot_files.txt', 'w') as f:
    for item in missing:
        f.write("%s\n" % item)