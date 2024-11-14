import os
import sys
sys.path.append("..")
import pickle
import time
import copy
import warnings
import pysam
import pandas as pd
import numpy as np
from dataclasses import dataclass
from intervaltree import IntervalTree, Interval
from clusterSig import readSig, clusterSig
from scipy.stats import ttest_ind
from sklearn.cluster import DBSCAN
from collections import Counter
####################################################
def seqDupInfo(contig, start, end):
    dup_bed = "hg38.inv1.bed.gz"
    tabx_repeat = pysam.TabixFile(dup_bed)

    repeat_starts, repeat_ends, repeat_lens = [], [], []
    for lx in tabx_repeat.fetch(reference=contig, start=start, end=end):
        line = [int(i) if i >0 else i for i,v in enumerate(lx.strip().split('\t'))]
        repeat_starts.extend([line[1],line[2]])
        repeat_ends.extend([line[1]+line[3], line[2]+line[3]])
        repeat_lens.append(line[3])

    return repeat_starts, repeat_ends, repeat_lens

def seqInvInfo(contig, start, end):
    inv_bed = "hg38.inv2.bed.gz"
    tabx_inv_repeat = pysam.TabixFile(inv_bed)

    repeat_starts, repeat_ends, repeat_lens = [], [], []
    for lx in tabx_inv_repeat.fetch(reference=contig, start=start, end=end):
        line = [int(i) if i >0 else i for i,v in enumerate(lx.strip().split('\t'))]
        print(line)
        if line[3] == 1 and line[4] == 2:
            repeat_starts.extend([line[1],line[2]-line[5]])
            repeat_ends.extend([line[1]+line[5], line[2]])
            repeat_lens.append(line[5])
        else:
            repeat_starts.extend([line[1]-line[5], line[2]])
            repeat_ends.extend([line[1], line[2]+line[5]])
            repeat_lens.append(line[5])

    return repeat_starts, repeat_ends, repeat_lens
####################################################
'''
print("a")
#bed_gz = "chrY.11000001.12000000.plk"
win = 1000000
inv_sigs_file = "tmp/tmp/SIM/SIM.INV.sigs"
inv_sub_sigs_file = "tmp/tmp/SIM/SIM.sub.INV.sigs"

#chr1    26639336        26648754        inversion       None    0
contig = "chr1"
start = 26638143
end =  26657722
#chr1    26638143        26657722        INV     1
#chr1    26646421        26657722        INV     1
#chr1    26647390        26650167        INV     1
#chr1    26638149        26650140        INV     7
#chr1    26641632        26646423        INV     3
#chr1    26646437        26650135        INV     6

#hg38.inv2.bed.gz
inv_bed = "hg38.inv2.bed.gz"
tbx = pysam.TabixFile(inv_bed)

clusters = pickle.load(open("tmp/tmp/SIM/chr1.26000001.27000000.plk","rb"))


cluster_D = clusters['D']
cluster_I = clusters["I"]
cluster_L = clusters["L"]
cluster_R = clusters["R"]


repeat_starts, repeat_ends, repeat_lens = seqDupInfo(contig, start, end)
repeat_inv_starts, repeat_inv_ends, repeat_inv_lens = seqInvInfo(contig, start, end)

print("repeat_starts:", repeat_starts)
print("repeat_ends:", repeat_ends)
print("#####################################################################")
print("repeat_inv_starts:", repeat_inv_starts)
print("repeat_inv_ends:", repeat_inv_ends)
'''
