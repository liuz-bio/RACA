import pysam
from intervaltree import IntervalTree, Interval
from sklearn.cluster import KMeans
import numpy as np
from collections import Counter

def ref_sigs(contig, start, end):
    refBed = '/home/liuz/work_space/project/liuz/plotSV/plotSV/callSV/tmp/tmp/HG38/HG38.bed.gz'
    handle = pysam.pysam.TabixFile(refBed)
    sigTs = {'1':'L', '2':'D', '3':'I', '4':'R'}
    refBed_sigs = {'1': IntervalTree(), '2': IntervalTree(), '3': IntervalTree(), '4': IntervalTree()}
    reads = handle.fetch(contig, start, end)
    for read in reads:
        read = read.strip().split('\t')
        for i in read[10:]:
            i = i.strip().split(',')
            s = int(i[0])
            e = int(i[1])
            refBed_sigs[sigTs[i[2]]].addi(s, e)
    return refBed_sigs

def kmeans(svlens):
    kms = KMeans(n_clusters=2)
    kms.fit(svlens.reshape(-1, 1))
    mks_cl = Counter(kms.labels_)
    if mks_cl[0] > mks_cl[1]:
        return svlens[kms.labels_==0]
    elif mks_cl[0] < mks_cl[1]:
        return svlens[kms.labels_==1]
    else:
        return svlens


# def ref_sigs(contig, start, end):
#     #le /home/liuz/work_space/project/liuz/plotSV/plotSV/callSV/tmp/tmp/HG38/HG38.bed.gz|cut -f 1,11-1000000|grep -vP '\t$'|sed "s/\t/|/g"|while read id;do(i=${id%%|*};ii=${id#*|}; echo $ii|sed "s/^/$i\t/g"|sed "s/|/\n$i\t/g"|sed "s/,/\t/g");done|bgzip -@ 10 >ref.HG38.bed.gz
#     refBed = '/home/liuz/work_space/project/liuz/plotSV/plotSV/callSV/ref.HG38.bed.gz'
#     handle = pysam.pysam.TabixFile(refBed)
#     refBed_sigs = {'1': IntervalTree(), '2': IntervalTree(), '3': IntervalTree(), '4': IntervalTree()}
#     reads = handle.fetch(contig, start, end)
#     for read in reads:
#         i = read.strip().split('\t')
#         s = int(i[0])
#         e = int(i[1])
#         refBed_sigs[i[2]].addi(s, e)
#     return refBed_sigs

