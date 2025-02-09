### Arrange CNA patterns of different resourses; genome (WGS) and mRNAs (inferCNV, i6) for comparative analysis

1. mTCN (modifide total copy number)
Fisrt, incremental classification of somatic CNAs (TCN in VCF files) and infer CNV profiles (i6 HMM prediction) were simplifide to obtain mTCN. mTCN ∈{ -1, 0, 1} as -1, any loss/deletion events; 1, any gain/amplification event; 0, TCN is equivalent to estimated ploidy.

2. Binning
Next, fragment length-weighted average of mTCNs were calculated per bins of genomic regions by 1M-base-intervals.
