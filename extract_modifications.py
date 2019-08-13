# -*- coding: utf-8 -*-

import pandas as pd
import os
import sys
import numpy as np

# Load datasets
ids = [ name for name in os.listdir(".") if os.path.isdir(name) ]
ids_decap = [s[:-1] if s[-6:] in ["mirnaA","mirnaB"] else s for s in ids]

# Map samples
metadata = pd.read_csv(sys.argv[1]).loc[:,["ID","Source"]]
#metadata = pd.read_csv("../../../Data/NGS/TCGA_pilot/metadata.txt").loc[:,["ID","Source"]]
metadata.drop_duplicates(inplace=True)
metadata.set_index("ID",inplace=True)
metadata = metadata.loc[ids_decap,:]
id_names = list(set([s1+"-"+s2 if sum([s2 in n for n in ids])==1 else s1+"-"+s2+"1" for s1,s2 in zip(metadata.loc[:,"Source"],ids)]))
id_names.extend([s[:-1]+"2" for s in id_names if s[-7:]=="_mirna1"])

# Initialize dataframe
merged = pd.DataFrame(columns=['contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount',
       'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth',
       'otherBases', 'improperPairs', 'sample', 'modif_id'])
ids=np.array(ids)

# Extract data from datasets
for n,d in enumerate(id_names):
    fileid = ids[np.bool_([s in d for s in ids])][0]
    filename = os.path.join(fileid,str("%s.ASE.csv" % fileid))
    temp = pd.read_csv(filename,sep="\t")
    uniqueids = temp.apply(lambda x: str("%s-%i-%sto%s" % (x.loc["contig"],x.loc["position"],x.loc["refAllele"],x.loc["altAllele"])),axis=1)
    temp.loc[:,'sample'] = d
    temp.loc[:,'modif_id'] = uniqueids
    merged = merged.append(temp, ignore_index=True)

# Read output
merged.to_csv("modifications.csv")
