#import pandas as pd
import re, gzip, math, glob, os
import numpy as np
import pandas as pd 
from fn_make_SegChr import make_SegChr

# Make cnv_facet-TCN dataframe segmented by window
dir="/data/share/WGS/Meso/result/hg38/cnv_facets/"
#InFiles = glob.glob(dir+"*M10_T.vcf.gz")
InFiles = glob.glob(dir+"M03_2_T.vcf.gz")
#InFiles = glob.glob(dir+"M03_cell_*.vcf.gz")
print(InFiles)
for InFile in InFiles:
    sample = re.sub(r".*(M.*)\.vcf\.gz", r"\1", InFile)
    print(sample)
    OutFile = "./facetTCN_SegmentedChr_" + sample + ".csv"
    #if os.path.isfile(OutFile):
    #    continue
    tmpFile = "./tmp/tmp_facetTCN_SegmentedChr_" + sample + ".csv"
    with gzip.open(InFile, "rt") as fr:
        n = 0
        for line in fr:
            if re.compile("^##ploidy=").search(line):
                ploidy=line.strip().replace("##ploidy=", "") # Ploidy of the tumor bulk
                ploidy=round(float(ploidy)) # baseline ploidy = 2 or 3
            if re.compile("^#").search(line): # comment lines
                continue
            else:
                items = line.strip().split("\t")
                FILTER=items[6]
                INFOs=items[7].split(";")
                TCN=int(INFOs[11].replace("TCN_EM=", "")) # Total copy number
                mTCN = int(TCN-ploidy) if TCN==ploidy else int((TCN - ploidy)/abs(TCN - ploidy)) # mTCN=-1, 0(=Gross ploidy), 1 
                CHR=items[0]
                if FILTER=="PASS" and TCN != ploidy and CHR != "chrX":
                    CHR=int(items[0].replace("chr", ""))
                    POS=int(items[1])
                    END=int(INFOs[2].replace("END=", ""))
                    ChrSegs, ChrSegLens =  make_SegChr(CHR = CHR, start = POS, end = END)
                    d_seg_tcn = pd.DataFrame(
                        data={
                            "Segmented_Chr": np.array(ChrSegs),
                            "TCN": np.array([TCN for i in range(len(ChrSegs))]),
                            "mTCN": np.array([mTCN for i in range(len(ChrSegs))]),
                            "Segmented_Chr_Len": np.array(ChrSegLens)
                        }
                    )
                    if n == 0:
                        d_seg_tcn.to_csv(tmpFile, mode='w', index = False)
                    else:
                        d_seg_tcn.to_csv(tmpFile, mode='a', index = False, header=False)
                    n = n + 1
    # Summarize the duplicated Segmeted_Chr
    d = pd.read_csv(tmpFile)
    d["TCN"] = d["TCN"]*d["Segmented_Chr_Len"]
    d["mTCN"] = d["mTCN"]*d["Segmented_Chr_Len"]
    dg = d.groupby("Segmented_Chr")[["TCN","mTCN"]].sum().round({"TCN":3,"mTCN":3})
    #dg = d.groupby("Segmented_Chr")["mTCN"].sum().to_frame("mTCN").round({"mTCN":3})
    print(dg.head)    #
    d_seg_all = pd.read_csv("chr_segmented.tsv").set_index("Segmented_Chr")  # Empty temlate of Segmented chromosome
    dg_all = pd.merge(d_seg_all, dg, how='outer', left_index=True, right_index=True, sort=False)
    dg_all.fillna({"TCN":ploidy, "mTCN":0}).to_csv(OutFile, mode='w', index = True)
    
