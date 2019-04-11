# -*- coding: utf-8 -*-

import pandas as pd
import os
import sys
import numpy as np

# Load datasets
ids = [ name for name in os.listdir(".") if os.path.isdir(name) ]
ids_decap = [s[:34] for s in ids]

# Map samples
metadata = pd.read_csv(sys.argv[1]).loc[:,["ID","Source"]]
#metadata = pd.read_table("../../../Data/NGS/TCGA_pilot/metadata.txt").loc[:,["File Name","Project"]]
metadata.drop_duplicates(inplace=True)
metadata.set_index("ID",inplace=True)
metadata = metadata.loc[ids_decap,:]
id_names = list(set([s1+"-"+s2 if sum([s2 in n for n in ids])==1 else s1+"-"+s2+"1" for s1,s2 in zip(metadata.loc[:,"Source"],ids)]))
id_names.extend([s[:-1]+"2" for s in id_names if s[-7:]=="_mirna1"])

# Initialize dataframe
mod = pd.DataFrame(columns=id_names)
nomod = pd.DataFrame(columns=id_names)
ids=np.array(ids)

# Extract data from datasets
for n,d in enumerate(id_names):
    fileid = ids[np.bool_([s in d for s in ids])][0]
    filename = os.path.join(fileid,str("%s.RPM.csv" % fileid))
    temp = pd.read_csv(filename,index_col=0)
    mod.loc[:,d] = temp.loc[:,"rpm"]
    nomod.loc[:,d] = temp.loc[:,"rpm_nomod"]

mod.set_index(temp.index)
nomod.set_index(temp.index)

# Read output
nomod.to_csv("results_nomod.csv")
mod.to_csv("results_mod.csv")
