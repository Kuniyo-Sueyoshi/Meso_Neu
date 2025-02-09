#import pandas as pd
import re, gzip, math, glob
import numpy as np
import pandas as pd 
from fn_make_SegChr import make_SegChr

# Dataframe with index of Segmented hg38 genome 
InFile = "./chr_len.tsv"
OutFile = "./chr_segmented.tsv"

with open(InFile, "rt") as fr:
    next(fr)
    ChrSeg_ls = [] # Initialize
    for line in fr:
        items = line.strip().split("\t")
        CHR = int(items[0])
        ChrLen = int(items[1])
        ChrSegs, ChrSegLens = make_SegChr(CHR = CHR, start = 1, end = ChrLen)
        ChrSeg_ls += ChrSegs
    d_seg_all = pd.DataFrame(index=ChrSeg_ls, columns=[])
    d_seg_all.index.name = "Segmented_Chr"

d_seg_all.to_csv(OutFile, mode='w', index = True)
    
