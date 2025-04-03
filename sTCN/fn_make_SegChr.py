import math

def make_SegChr(CHR, start, end, window=1000000, PseudoChrNum=1000000000): # Return list of segmented chrmoseme location
    cont_start = (CHR*PseudoChrNum + start)/window
    cont_end = (CHR*PseudoChrNum + end)/window
    seg_start = math.floor(cont_start)
    seg_end = math.ceil(cont_end)
    frac_left = math.ceil(cont_start) - cont_start # Fraction of left end (ex. cont_start = 1011.3; frac_left = 0.7)
    frac_right = cont_end - math.floor(cont_end) # Fraction of right end (ex. cont_end = 1014.6; frac_right = 0.6)
    Segmented_Chr = [i for i in range(seg_start, seg_end)]
    if int(cont_start) == int(cont_end):
        Segmented_Chr_Len = [cont_end - cont_start]
    else:
        Segmented_Chr_Len = [frac_left] + [1 for i in range(seg_start+1, seg_end-1)] + [frac_right]
    return Segmented_Chr, Segmented_Chr_Len
