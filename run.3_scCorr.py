#!usr/bin/env python
import pandas as pd 

ddir="/data/sueyoshi/mpm/scRNAseq/inferCNV/Ref.M12-15/"
dict_sm_wgs = {
    "M01":"M01_N",
    "M02":"M02_T",
    "M03_2":"M03_2_T",
    "M04":"M04_T",
    "M05":"M05_T",
    "M06":"M06_T",
    "M07":"M07_T",
    "M08":"M08_T",
    "M09":"M09_T",
    "M10":"M10_T",
    "M11":"M11_T"
}

for sm in dict_sm_wgs.keys():
    print(sm)
    sm_wgs = dict_sm_wgs.get(sm)

    InFile_iCNV = "iCNV.HMM_SegmetedChr_" + sm + ".csv" # cell_group-based iCNV
    InFile_mTCN = "facetTCN_SegmentedChr_" + sm_wgs + ".csv" # CNV
    InFile_ann = ddir + sm + "/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings" # cell - cell_group tbl 
    OutFile = "sc-iCNV.Corr_" + sm + ".csv"
    # Read files
    d_iCNV = pd.read_csv(InFile_iCNV).set_index("Segmented_Chr")
    d_mTCN = pd.read_csv(InFile_mTCN).set_index("Segmented_Chr")
    cell_ann = pd.read_csv(InFile_ann, sep = "\t").set_index("cell_group_name")
    # Correlation of "single-cell iTCN" and "WGS-mTCN"
    corr_gp = d_iCNV.corrwith(d_mTCN.mTCN)
    d_corr_gp = pd.DataFrame(corr_gp).set_axis(['Corr'], axis=1)
    # "cell_group_name" to "cell" barcode
    d_corr_cell = pd.merge(cell_ann, d_corr_gp, right_index=True, left_index=True)
    d_corr_cell.to_csv(OutFile,  mode='w', index = False)

    # single-cell x iCNV matrix
    OutFile2 = "iCNV.HMM_SegmetedChr_sc_" + sm + ".csv" # single-cell-based iCNV"
    d_sc_iCNV_T = pd.merge(cell_ann, d_iCNV.T, right_index=True, left_index=True)
    d_sc_iCNV_T.to_csv(OutFile2,  mode='w', index = False)
