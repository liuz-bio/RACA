import gzip
import json
import numpy as np
import pandas as pd
import copy
import time
import os
from sklearn.cluster import DBSCAN
from collections import Counter
#L:0 R:1
def appendRecord(chrom, query_start, query_end, is_reverse, line9, recordHash, allDict, queryDict):
    for lx in line9:
        lx = lx.strip(",").split(',')
        if lx[2] =='1':
            if chrom in allDict:
                try:
                    allDict[chrom]['G'].append(int(lx[1]))
                    allDict[chrom]['Q'].append(int(query_start))
                    allDict[chrom]['RL'].append(0)#L
                    queryDict[query_start] = [query_end, recordHash, is_reverse]
                except KeyError:
                    allDict[chrom]['G'] = [int(lx[1])]
                    allDict[chrom]['Q'] = [int(query_start)]
                    allDict[chrom]['RL'] = [0]#L
                    queryDict[query_start] = [query_end, recordHash, is_reverse]
            else:
                allDict[chrom] = {'G':[int(lx[1])]}
                allDict[chrom]['Q'] = [int(query_start)]
                allDict[chrom]['RL']  = [0]#L
                queryDict[query_start] = [query_end, recordHash, is_reverse]
        elif lx[2] == '4':
            if chrom in allDict:
                try:
                    allDict[chrom]['G'].append(int(lx[0]))
                    allDict[chrom]['Q'].append(int(query_end))
                    allDict[chrom]['RL'].append(1)#R
                    queryDict[query_end] = [query_start, recordHash, is_reverse]
                except KeyError:
                    allDict[chrom]['G'] = [int(lx[0])]
                    allDict[chrom]['Q'] = [int(query_end)]
                    allDict[chrom]['RL'] = [1]#R
                    queryDict[query_end] = [query_start, recordHash, is_reverse]
            else:
                allDict[chrom] = {'G':[int(lx[0])]}
                allDict[chrom]['Q'] = [int(query_end)]
                allDict[chrom]['RL'] = [1]#R
                queryDict[query_end] = [query_start, recordHash, is_reverse]
    return allDict, queryDict


#allDict: {'chr1': {'genome': [2033372], 'query': [9046], 'RL': ['L']}, 'chr4': {'genome': [71592083], 'query': [58287], 'RL': ['L']}}
#allDict: {'chr1': {'genome': [2033372], 'query': [1167], 'RL': ['R']}, 'chr4': {'genome': [71592081], 'query': [3446], 'RL': ['R']}}
#DEL 0, INS 1, INV 2, DUP 3, BND 4
def dealBND(allDict, queryDict, records, sub_records, query_length, reverseLists, outBnds):
    breakpointDeviation = 100
    reverseDeviation = 100

    try:
        reverseL = reverseLists[0] / reverseLists[1]
    except ZeroDivisionError:
        reverseL = 0

    chroms = list(allDict.keys())
    if len(chroms) < 1:
        return []

    for i, CHR1 in enumerate(chroms):
        for CHR2 in chroms[i:]:
            query_CHR1_size = len(allDict[CHR1]['Q'])
            query_CHR2_size = len(allDict[CHR2]['Q'])
            for i1 in range(query_CHR1_size):
                for i2 in range(i1+1, query_CHR2_size):
                    #CHR1==CHR2: chrX chrX 0 1 153066025 153084040 7780 25773
                    #print("CHR1==CHR2:",CHR1, CHR2, allDict[CHR1]['RL'][i1], allDict[CHR2]['RL'][i2], allDict[CHR1]['G'][i1], allDict[CHR2]['G'][i2], allDict[CHR1]['Q'][i1], allDict[CHR2]['Q'][i2])
                    CHR1_RL_i1 =  allDict[CHR1]['RL'][i1]
                    CHR2_RL_i2 = allDict[CHR2]['RL'][i2]
                    CHR1_G_i1 = allDict[CHR1]['G'][i1]
                    CHR2_G_i2 = allDict[CHR2]['G'][i2]
                    CHR1_Q_i1 = allDict[CHR1]['Q'][i1]
                    CHR2_Q_i2 = allDict[CHR2]['Q'][i2]

                    if CHR1 == CHR2:
                        if CHR1_G_i1 == CHR2_G_i2:
                            continue

                        if CHR1_RL_i1 == 0:#L
                            if CHR2_RL_i2 == 0:#L L
                                if abs(CHR1_G_i1 - CHR2_G_i2) > breakpointDeviation:
                                    a = sorted([CHR1_G_i1, CHR2_G_i2])
                                    if abs(CHR1_Q_i1 + CHR2_Q_i2 - query_length) <= reverseDeviation:
                                        records[2].setdefault(CHR1, []).append([a[0], a[1], 2]) #INV
                                    else:
                                        sub_records[2].setdefault(CHR1, []).append([a[0], a[1], 2]) #INV
                                        #print("INV test:", a[0], a[1], 2, CHR1_Q_i1,  CHR2_Q_i2, query_length)
                            else: #L R
                                if CHR1_G_i1 > CHR2_G_i2:
                                    if abs(CHR1_G_i1 - CHR2_G_i2) > breakpointDeviation:
                                        a = sorted([CHR1_G_i1, CHR2_G_i2])
                                        if abs(CHR1_Q_i1 + CHR2_Q_i2 - query_length) <= reverseDeviation:
                                            records[0].setdefault(CHR1, []).append([a[0], a[1], 0]) #DEL
                                        else:
                                            sub_records[0].setdefault(CHR1, []).append([a[0], a[1], 2]) #DEL
                                    elif abs(CHR1_Q_i1 - CHR2_Q_i2) > breakpointDeviation:
                                        a = sorted([CHR1_G_i1, CHR2_G_i2])
                                        records[1].setdefault(CHR1, []).append([a[0], a[1], 1]) #INS
                                else:#CHR1==CHR2: chrX chrX 0 1 153066025 153084040 7780 25773
                                    #print('R>LDUP:', CHR1, CHR2, CHR1_G_i1, CHR2_G_i2, CHR1_Q_i1, CHR2_Q_i2, CHR1_RL_i1, CHR2_RL_i2, "query_length:", query_length)
                                    if abs(CHR1_Q_i1 - CHR2_Q_i2) <= breakpointDeviation:
                                        if abs(CHR1_G_i1 - CHR2_G_i2) > breakpointDeviation:
                                            a = sorted([CHR1_G_i1, CHR2_G_i2])
                                            records[3].setdefault(CHR1, []).append([a[0], a[1], 3]) #DUP
                                    elif abs(queryDict[CHR1_Q_i1][0] - query_length) <= 50 and abs(queryDict[CHR2_Q_i2][0]) <= 50:
                                        a = sorted([CHR1_G_i1, CHR2_G_i2])
                                        records[3].setdefault(CHR1, []).append([a[0], a[1], 3])  # DUP
                                    elif reverseL != 0 and queryDict[CHR1_Q_i1][1]==queryDict[CHR2_Q_i2][1]:
                                        a = sorted([CHR1_G_i1, CHR2_G_i2])
                                        records[2].setdefault(CHR1, []).append([a[0], a[1], 2])
                        elif CHR2_RL_i2 == 0:# R L
                            if CHR1_G_i1 > CHR2_G_i2:
                                #R >  L
                                #print('R>LDUP:', CHR1, CHR2, CHR1_G_i1, CHR2_G_i2, CHR1_Q_i1, CHR2_Q_i2, CHR1_RL_i1, CHR2_RL_i2, "query_length:", query_length)
                                if abs(CHR1_Q_i1 - CHR2_Q_i2) <= breakpointDeviation:
                                    if abs(CHR1_G_i1 - CHR2_G_i2) > breakpointDeviation:
                                        a = sorted([CHR1_G_i1, CHR2_G_i2])
                                        records[3].setdefault(CHR1, []).append([a[0], a[1], 3]) #DUP
                                elif abs(queryDict[CHR2_Q_i2][0] - query_length) <= 50 and abs(queryDict[CHR1_Q_i1][0]) <= 50:
                                    a = sorted([CHR1_G_i1, CHR2_G_i2])
                                    records[3].setdefault(CHR1, []).append([a[0], a[1], 3]) #DUP
                                elif reverseL != 0 and queryDict[CHR1_Q_i1][1]==queryDict[CHR2_Q_i2][1]:
                                    a = sorted([CHR1_G_i1, CHR2_G_i2])
                                    records[2].setdefault(CHR1, []).append([a[0], a[1], 2])
                            else: #R < L
                                if abs(CHR1_G_i1 - CHR2_G_i2) > breakpointDeviation:
                                    a = sorted([CHR1_G_i1, CHR2_G_i2])
                                    if abs(CHR1_Q_i1 + CHR2_Q_i2 - query_length) <= reverseDeviation:
                                        records[0].setdefault(CHR1, []).append([a[0], a[1], 0])
                                    else:
                                        sub_records[2].setdefault(CHR1, []).append([a[0], a[1], 2]) #INV
                                elif abs(CHR1_Q_i1 - CHR2_Q_i2) > breakpointDeviation:
                                    a = sorted([CHR1_G_i1, CHR2_G_i2])
                                    records[1].setdefault(CHR1, []).append([a[0], a[1], 1])
                        else:#CHR1==CHR2: chrX chrX 0 1 153066025 153084040 7780 25773
                            if abs(CHR1_G_i1 - CHR2_G_i2) > breakpointDeviation:
                                a = sorted([CHR1_G_i1, CHR2_G_i2])
                                if abs(CHR1_Q_i1 + CHR2_Q_i2 - query_length) <= reverseDeviation:
                                    records[2].setdefault(CHR1, []).append([a[0], a[1], 2])
                                else:
                                    sub_records[2].setdefault(CHR1, []).append([a[0], a[1], 2]) #INV
                                    #print("INV test:", a[0], a[1], 2, CHR1_Q_i1, CHR2_Q_i2, query_length)
                    elif len(set([CHR1_RL_i1, CHR2_RL_i2])) == 2:
                        idx = '|'.join([CHR1, CHR2])
                        if CHR1_RL_i1 == 1 and CHR1_G_i1 <= CHR2_G_i2: #R
                            records[4].setdefault(idx, []).append([CHR1_G_i1, CHR2_G_i2, 4, CHR1_RL_i1, CHR2_RL_i2])
                            outBnds.write('\t'.join(map(str, [idx,CHR1_G_i1, CHR2_G_i2, 4,CHR1_Q_i1,CHR2_Q_i2,CHR1_RL_i1,CHR2_RL_i2]))+'\n')
                        elif CHR1_G_i1 >= CHR2_G_i2:
                            records[4].setdefault(idx, []).append([CHR1_G_i1, CHR2_G_i2, 4, CHR1_RL_i1, CHR2_RL_i2])
                            outBnds.write('\t'.join(map(str, [idx,CHR1_G_i1, CHR2_G_i2, 4,CHR1_Q_i1,CHR2_Q_i2,CHR1_RL_i1,CHR2_RL_i2]))+'\n')
                    elif abs(CHR1_Q_i1 + CHR2_Q_i2 - query_length) <= 200:
                        idx = '|'.join([CHR1, CHR2])
                        records[4].setdefault(idx, []).append([CHR1_G_i1, CHR2_G_i2, 4, CHR1_RL_i1, CHR2_RL_i2])
                        outBnds.write('\t'.join(map(str, [idx,CHR1_G_i1, CHR2_G_i2, 4,CHR1_Q_i1,CHR2_Q_i2,CHR1_RL_i1,CHR2_RL_i2]))+'\n')


def splitSigs(chrom, svType, record, dbscan, sub=False):
    # Assuming clusterBND and cluster functions are defined elsewhere and optimized already.
    groups = []
    step = 300
    max_distance = 1500
    record_array = np.array(record, dtype=int)
    record_array = record_array[np.argsort(record_array[:, 0])]
    record_size = len(record_array)

    if record_size <= 1000:
        if svType == "BND":
            clusterBND(chrom, svType, record_array, dbscan, groups, sub)
        else:
            cluster(chrom, svType, record_array, dbscan, groups, sub)
        return groups

    start_index = 0
    while start_index < record_size:
        stop_index = min(start_index + step, record_size)

        if stop_index >= record_size or (np.linalg.norm(record_array[stop_index - 1, :2] - record_array[stop_index, :2]) > max_distance):
            segment = record_array[start_index:stop_index]
            start_index = stop_index
        else:
            end_index = np.argwhere(np.linalg.norm(record_array[stop_index - 1:record_size - 1, :2] - record_array[stop_index:record_size, :2], axis=1) > max_distance)
            if len(end_index) > 0:
                end_index = end_index[0, 0] + stop_index
                segment = record_array[start_index:end_index]
                start_index = end_index
            else:
                segment = record_array[start_index:]
                start_index = record_size

        if len(segment) > 2*step:
            segment = segment[:step]

        if svType == "BND":
            clusterBND(chrom, svType, segment, dbscan, groups, sub)
        else:
            cluster(chrom, svType, segment, dbscan, groups, sub)

    return groups

def cluster(chrom, svType, record, dbscan, groups, sub):
    clusters = dbscan.fit_predict(record[:,:2])
    labels = np.unique(clusters)
    for label in labels:
        cluster_points = record[clusters == label]
        if label != -1:
            pos1, pos2 = np.mean(cluster_points[:, :2], axis=0).astype(int)
            groups.append([chrom, pos1, pos2, svType, len(cluster_points)])
            #print("cluster test:", [chrom, pos1, pos2, svType, len(cluster_points)])
        elif not sub:
            # 优化点：避免for循环，用列表解析式代替
            groups.extend([[chrom, i[0], i[1], svType, 1] for i in cluster_points])

def clusterBND(chrom, svType, record, dbscan, groups, sub):
    clusters = dbscan.fit_predict(record[:,:2])
    labels = np.unique(clusters)
    for label in labels:
        cluster_points = record[clusters == label]
        if label != -1:
            pos1, pos2 = np.mean(cluster_points[:, :2], axis=0).astype(int)
            rl1 = pd.Series(cluster_points[:, 3]).value_counts().idxmax()
            rl2 = pd.Series(cluster_points[:, 4]).value_counts().idxmax()
            groups.append([chrom, pos1, pos2, svType, rl1, rl2])
        elif not sub:
            # 优化点：避免for循环，用列表解析式代替
            groups.extend([[chrom, i[0], i[1], svType, i[3], i[4]] for i in cluster_points])


def subMerge(groups):
    groups = iter(groups)
    head = next(groups)
    flag = True
    new_groups = []

    for i in groups:
        head_chrom = head[0]
        head_start = int(head[1])
        head_end = int(head[2])
        head_s = int(head[4])

        next_chrom = i[0]
        next_start = int(i[1])
        next_end = int(i[2])
        next_s = int(i[4])
        if head_chrom == next_chrom:
            if head_start <= next_start <= head_end or head_start <= next_end <= head_end:
                tmp = [head_start, head_end, next_start, next_end]
                new_groups.append([head_chrom, min(tmp), max(tmp), "INV", max([head_s, next_s])])
                head = i
                flag = False
            elif flag:
                new_groups.append(head)
                head = i
            else:
                flag = True
                head = i
        elif flag:
            new_groups.append(head)
            head = i
        else:
            flag = True
            head = i
    return  groups

def candidateSVLocation(bed_sorted_file, outPath, sample):
    startTime = time.time()
    min_quality =20
    outBnds = open(os.path.join(outPath, sample, sample+".outBnds.txt"),'w')

    with open(bed_sorted_file, 'rt') as f:
        while True:
            head = next(f).split('\t')
            if len(head) >5:
                break
        allDict = {}
        queryDict = {}
        reverseList = []
        records = [{},{},{},{},{}]
        sub_records = [{}, {}, {}, {}, {}]
        svTypeDict = {0:'DEL',1:'INS',2:'INV',3:'DUP',4:'BND'}
        dbscan = DBSCAN(eps=500, min_samples=2)
        for lx in f:
            line = lx.strip().split('\t')
            if head[9] != line[9]:
                if len(allDict) >0:
                    for chrom in allDict.keys():
                        sorted_index = sorted(range(len(allDict[chrom]['Q'])), key=lambda k: allDict[chrom]['Q'][k])
                        allDict[chrom]['G'] = [allDict[chrom]['G'][i] for i in sorted_index]
                        allDict[chrom]['Q'] = [allDict[chrom]['Q'][i] for i in sorted_index]
                        allDict[chrom]['RL'] = [allDict[chrom]['RL'][i] for i in sorted_index]
                    query_length = int(head[5])
                    #print("allDict:", allDict)
                    dealBND(allDict, queryDict, records, sub_records, query_length, Counter(reverseList), outBnds)
                allDict = {}
                queryDict = {}
                reverseList = []
                head = line
                continue
            elif int(line[6]) >= min_quality:
                if len(allDict) <1:
                    reverseList.append(int(head[7]))
                    reverseList.append(int(line[7]))
                    #def appendRecord(chrom, query_start, query_end, line9, recordHash, allDict, queryDict):
                    allDict, queryDict = appendRecord(head[0], int(head[3]), int(head[4]), int(head[7]), head[10:], hash(lx), allDict, queryDict)
                    allDict, queryDict = appendRecord(line[0], int(line[3]), int(line[4]), int(line[7]), line[10:], hash(lx), allDict, queryDict)
                else:
                    allDict, queryDict = appendRecord(line[0], int(line[3]), int(line[4]), int(line[7]), line[10:], hash(lx), allDict, queryDict)
                    reverseList.append(int(line[7]))

        #print("Time:", time.time()-startTime)
        startTime = time.time()
        #out_raw = open('out_raw.tab','w')
        with open(os.path.join(outPath, sample, sample+".out.tab"), 'w') as out:
            for i in range(len(records)):
                svType = svTypeDict[i]
                with open(os.path.join(outPath, sample, sample+"."+svType+".sigs"),'w') as outsvType:
                    for chrom in records[i].keys():
                        records[i][chrom] = splitSigs(chrom, svType, records[i][chrom], dbscan, sub=False)
                        for ii in  records[i][chrom]:
                            out.write('\t'.join(map(str, ii))+'\n')
                            outsvType.write('\t'.join(map(str, ii))+'\n')

                with open(os.path.join(outPath, sample, sample + ".sub." + svType + ".sigs"), 'w') as sub_outsvType:
                    for chrom in sub_records[i].keys():
                        sub_records[i][chrom] = splitSigs(chrom, svType, sub_records[i][chrom], dbscan, sub=True)
                        for ii in sub_records[i][chrom]:
                            out.write('\t'.join(map(str, ii)) + '\tsub\n')
                            sub_outsvType.write('\t'.join(map(str, ii)) + '\n')

    outBnds.close()
    print("Time:", time.time()-startTime)


if __name__ == "__main__":
    sample = "FY0_Father"
    out_path = "/home/lz/work_space/project/plotSV/plotSV5/callSV/tmp"
    tmp_path = os.path.join(out_path, "tmp")
    bedName = os.path.join(tmp_path, sample, sample + '.bed')
    bed_sorted_filie = bedName + ".readName.sorted"
    candidateSVLocation(bed_sorted_filie, tmp_path, sample)
