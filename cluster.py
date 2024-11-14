import os
import pickle
import time

import pysam
from sklearn.cluster import DBSCAN
from intervaltree import IntervalTree, Interval
import numpy as np
import multiprocessing
import traceback
import  pandas as pd
import hashlib
import warnings
from scipy.stats import nbinom, binom
from clusterSig import readSig, clusterSig
from sortedcontainers import SortedDict
from scipy.spatial import KDTree


def cluster1(sigs, dbscan, chrom, handle, sigT, fingerPrint, sigAllClusters, cand_SVs):
    sigs = np.array(sigs)
    _ = dbscan.fit_predict(sigs[:, :2].astype(int))

    labels = dbscan.labels_
    unique_labels = list(set(labels))
    unique_labels.sort()

    for k in unique_labels:
        if k == -1:
            continue
        tmp_index = (labels == k)
        cluster_xy = sigs[tmp_index][:, :2].astype(int)
        cluster = sigs[tmp_index][:, 2]

        if len(set([i.readName for i in cluster if i.mapQ >= 55])) <= 150:
            clusterSigObj = clusterSig(chrom, cluster_xy, cluster, handle, sigT, fingerPrint)
            fingerPrint += 1
            sigAllClusters[sigT].add(Interval(clusterSigObj.start_mode, clusterSigObj.stop_mode, clusterSigObj))
            cand_SVs[sigT] += 1#.append(clusterSigObj)

    return fingerPrint

def cluster_points1(points, min_distance=100, min_signals=2):
    clusters = []
    visited = set()
    # 创建 KD-Tree
    tree = KDTree(points)

    for i, point in enumerate(points):
        if i in visited:
            continue
        # 查找与当前点在 min_distance 范围内的点
        indices = tree.query_ball_point(point, r=min_distance)

        # 将这些点添加到当前聚类
        current_cluster = [i]
        visited.add(i)

        for j in indices:
            if j != i and j not in visited:
                current_cluster.append(j)
                visited.add(j)
                
        if len(current_cluster) >=2:
            clusters.append(current_cluster)

    return clusters

# def split_clusters(points, clusters, interval_start, interval_end):
#     for current_cluster in clusters:
#         if len(current_cluster[0]) >=2:
#             current_points = sorted(points[current_cluster], key=lambda x: x[2])#根据SV长度排序
#             head_index = 0
#             head_node = points[current_cluster[head_index]]
#             tmp_clusters = [current_cluster[head_index]]
#             for i in current_points:

#def merge_cluster(points, clusters, interval_start, interval_end)      
def merge_cluster_DI(clusters, points):
    
    pass        


def cluster_points_DI(points, chrom, interval_start, interval_end, sigT, min_distance, min_signals):
    clusters = []
    head_index = 0
    head_node = points[head_index]
    current_cluster = [head_index]
    for i, point in enumerate(points[1:]):
        i+=1
        dist =min([abs(head_node[0]-point[0]), abs(head_node[1]-point[1]), abs(head_node[0]-point[1]), abs(head_node[1]-point[0])])
        sv_cv = abs(head_node[2]-point[2])/min([head_node[2],point[2]])#长度波动水平
        overlap = min(head_node[1],point[1]) - max(head_node[0],point[0]) >0.5*max([head_node[2],point[2]])
        
        if (dist and overlap) or (dist < min_distance and sv_cv < 0.3 and abs(head_node[0]-point[0]) <min([head_node[2],point[2]])):
            current_cluster.append(i)
        else:  
            if len(current_cluster)>=2:
                std = np.mean([np.std(points[current_cluster][:,0]), np.std(points[current_cluster][:,1])])/np.mean(points[current_cluster][:,2])
                print("STD_RRR:",std/len(current_cluster))
                ratio = std/len(current_cluster)
                if  (ratio<0.05 and sigT=='I') or(ratio<0.15 and sigT=="D")  :
                    clusters.append(current_cluster)   
            elif len(current_cluster) == 1 and head_node[2]>500:                           
                clusters.append(current_cluster*2)
            current_cluster = [i]
            head_index = i  
            head_node = point
            
    if len(current_cluster)>=2:
        clusters.append(current_cluster)
                   
    for i in clusters:
        print('\n'.join(['\t'.join(map(str, ["QQQ_clusters", chrom, interval_start, interval_end, sigT]+j.tolist()+[abs(j[0]-j[1])])) for j in points[i]])+"\n###################################################")
        
    return clusters


def cluster_points_LR(points, chrom, interval_start, interval_end, sigT, min_distance, min_signals):
    clusters = []
    head_index = 0
    head_node = points[head_index]
    current_cluster = [head_index]
    for i, point in enumerate(points[1:]):
        i+=1
        dist =min([abs(head_node[0]-point[0]), abs(head_node[1]-point[1]), abs(head_node[0]-point[1]), abs(head_node[1]-point[0])])
        if dist < min_distance:
            current_cluster.append(i)
        else:  
            # if len(current_cluster) >=2:
            #     cv = np.std(points[current_cluster][:,-1])/np.mean(points[current_cluster][:,-1])
            #     if cv >= 0.1 : sorted(enumerate(a), key=lambda x: x[1])
            #         sort_current_cluster = sorted(enumerate(points[current_cluster]), key=lambda x: x[1][2])
            #         head_i = 0
            #         node_i = sort_current_cluster[0][1]
            #         tmp_cluster = [sort_current_cluster[0][0]]
            #         for j, node_j in sort_current_cluster[1:]:
            #             dist = node_i[0]+node_i[2]
            if len(current_cluster)>=2:
                clusters.append(current_cluster)    
            # if len(current_cluster) ==1:
            #                     clusters_index, pos_std, pos_mean, svlen_std, svlen_mean
            #     clusters.append([current_cluster, 0, points[current_cluster][:,-1][0])
            # else:
            #     clusters.append([current_cluster, np.std(points[current_cluster][:,-1]), np.mean(points[current_cluster][:,-1])])
            current_cluster = [i]
            head_index = i  
            head_node = point
            
    if len(current_cluster)>=2:
        clusters.append(current_cluster)
            
    # if len(current_cluster)>0:
    #     clusters.append(current_cluster)
    
    # if sigT in ['D', 'I']:
    #     n=0
    #     while len(clusters)>1 and n< len(clusters):
    #         head_node = clusters[n]
    #         next_node = clusters[n+1]
    #         if head_node[1]*len(head_node[0]) <
            
    #     for i in clusters[:-1]:
    #         j = clusters[i+1:]
    #         svlen_var_i = np.var(points[i][:,-1])
    #         svlen_var_j = np.var(points[j][:,-1])
    #         svlen_mean_i =  np.mean(points[i][:,-1]) 
    #         svlen_mean_j =  np.mean(points[j][:,-1]) 
    #         if abs(svlen_mean_i - svlen_mean_j)
                
    for i in clusters:
        print('\n'.join(['\t'.join(map(str, ["QQQ_clusters", chrom, interval_start, interval_end, sigT]+j.tolist()+[abs(j[0]-j[1])])) for j in points[i]])+"\n###################################################")
        
    return clusters


def cluster(interval_start, interval_end, sigs, min_distanc, min_signals, chrom, handle, sigT, fingerPrint, sigAllClusters, cand_SVs, bed):
    sorted_sigs = sorted(sigs, key=lambda x: x[0])
    sorted_sigs = np.array(sorted_sigs)
    start_time = time.time()
    try:
        if sigT in ['D', 'I']:
            cluster_index = cluster_points_DI(sorted_sigs[:, :3].astype(int), chrom, interval_start, interval_end, sigT, min_distanc, min_signals)
        else:
            cluster_index = cluster_points_LR(sorted_sigs[:, :2].astype(int), chrom, interval_start, interval_end, sigT, min_distanc, min_signals)
    except IndexError:
        return fingerPrint
    print("cluster_EEE",time.time() - start_time, bed)
    start_time = time.time()
    for tmp_index in cluster_index:
        cluster_xy = sorted_sigs[tmp_index][:, :2].astype(int)
        cluster = sorted_sigs[tmp_index][:, -1]
        start,stop = np.mean(cluster_xy[:,0]), np.mean(cluster_xy[:,1])
        depth = sum([1 for i in handle if i[0] <= start <= i[1] or i[0] <= stop <= i[1]])
        
        if len(set([i.readName for i in cluster if i.mapQ >= 45])) <= 200 and len(tmp_index)<=depth<=120:
            start_time = time.time()
            clusterSigObj = clusterSig(chrom, cluster_xy, cluster, handle, sigT, fingerPrint)
            print("cluster_WWW1",time.time() - start_time, len(list(handle)), len(cluster_xy), bed)
            fingerPrint += 1
            sigAllClusters[sigT].add(Interval(clusterSigObj.start_mode, clusterSigObj.stop_mode, clusterSigObj))
            cand_SVs[sigT] += 1#.append(clusterSigObj)
    print("cluster_WWW",time.time() - start_time, bed)
    return fingerPrint

# 插入区间的方法

def insert_interval(sorted_tree, interval):
    sorted_tree[interval.begin] = interval


def write_plk(file_plk, sigAllClusters):
    # 创建一个 SortedDict
    sorted_tree = SortedDict()
    # 将 IntervalTree 的区间插入到 SortedDict
    for sigT, interval_tree in sigAllClusters.items():
        for interval in interval_tree:
            insert_interval(sorted_tree, interval)
    
    if max([len(i) for i in sigAllClusters.values()]) >0:
        with open(file_plk, 'wb') as plkOut:
            pickle.dump(sigAllClusters, plkOut)
            pickle.dump(sorted_tree, plkOut)


def extraSigs(chrom, interval_start, interval_end, no_first, plk_name, bed, configer, bamfile, all_cand_svs, lock):
    # if configer.bam_cram:
    #     handle = pysam.AlignmentFile(bamfile, "rb", require_index=True)
    # else:
    #     handle = pysam.AlignmentFile(bamfile, "rc", require_index=True)

    sigAllClusters = {'D': IntervalTree(), 'R': IntervalTree(), 'L': IntervalTree(), 'I': IntervalTree()}
    sigTypes = {'1': 'L', '2': 'D', '3': 'I', '4': 'R', '5': 'A'}
    sigs = {'I':[],'D':[],'R':[],'L':[]}
    fingerPrint = 0
    num_nu = 0
    cand_SVs = {'D':0,'L':0,'I':0,'R':0}
    
    md5_hash = hashlib.md5()
    #handle = IntervalTree()
    handle = []
    start_time = time.time()
    for read in open(bed):
        md5_hash.update(read.strip().encode('utf-8'))
        read_hash = md5_hash.hexdigest()

        read = read.strip().split('\t')
        map_q = int(read[6])
        align_start = int(read[1])
        align_stop = int(read[2])
        handle.append([align_start, align_stop])
        #handle.addi(align_start, align_stop)
        #if no_first & (align_start <= interval_start <= align_stop):
        #   continue

        read_is_reverse = int(read[7])
        read_is_secondary = int(read[8])
        read_name = read[9]
        query_len = int(read[5])
        query_align_start = int(read[3])
        query_align_end = int(read[4])

        #print("read_name_read_hash:", read_name, read_hash)
        for sig in read[10:]:
            sig = sig.split(',')
            start = int(sig[0])
            stop = int(sig[1])
            svlen = abs(stop - start)
            sigType = sigTypes[sig[2]]

            if sigType in ['D', 'I']:
                if svlen < configer.min_sv_signal:
                    continue

                if sigType == "I":
                    read_sig_obj = readSig(start, stop, map_q, read_hash, read_is_reverse, read_is_secondary, read_name, 
                                         seq=sig[3], query_start=int(sig[4]), query_len=query_len )
                else:
                    read_sig_obj = readSig(start, stop, map_q, read_hash, read_is_reverse, read_is_secondary, read_name, 
                                           query_start=int(sig[3]), query_len=query_len)

                sigs[sigType].append([start, stop, svlen, read_sig_obj])
                print("extraSigs_WWW", chrom, interval_start, interval_end, sigType, start, stop, svlen)

            elif sigType == "L":
                read_sig_obj = readSig(stop, stop, map_q, read_hash, read_is_reverse, read_is_secondary, read_name, 
                                       query_start=int(sig[3]), query_len=query_len, sLR_other=align_stop, query_length=query_len)

                sigs[sigType].append([stop, stop, read_sig_obj])
                print("extraSigs_WWW", chrom, interval_start, interval_end, sigType, stop, stop, svlen)

            elif sigType == "R":
                read_sig_obj = readSig(start, start, map_q, read_hash, read_is_reverse, read_is_secondary, read_name,
                                       query_start=int(sig[3]), query_len=query_len, sLR_other=align_start, query_length=query_len)

                sigs[sigType].append([start, start, read_sig_obj])
                print("extraSigs_WWW", chrom, interval_start, interval_end, sigType, start, start, svlen)
    print("cluster_TTT",time.time() - start_time, bed)
    start_time = time.time()
    for sigT, sig in sigs.items():
        fingerPrint = cluster(interval_start, interval_end, sig, configer.cluster_param[sigT]['min_distanc'], configer.cluster_param[sigT]['min_signals'],
                                chrom, handle, sigT, fingerPrint, sigAllClusters, cand_SVs,bed)
    print("cluster_RRR",time.time() - start_time, bed)

#DBSCAN(eps=configer.cluster_param[sigT]['min_distanc'],
#min_samples=configer.cluster_param[sigT]['min_signals'])
    write_plk(plk_name, sigAllClusters)
    with lock:
        del cand_SVs['L']
        del cand_SVs['R']
        all_cand_svs[plk_name.split("/")[-1]] = cand_SVs
    #os.remove(bed)
    #os.remove(bed+".tbi")

#def run(chrom, interval_start, interval_end, no_first, plk_name, allSignals, bed, handle):
#    self.extraSigs(chrom, interval_start, interval_end, no_first, plk_name, bed, configer, handle)



if __name__ == "__main__":
    genomeFile = '/home/liuz/work_space/project/liuz/Genome/genome.fa'
    # bed = '/home/lz/work_space/Database/VISOR_test/2024_3_1/plotSV/callSV_new/SIM.bed.gz'
    bed = '/home/liuz/work_space/project/liuz/plotSV/plotSV/callSV_new/SIM.bed.gz'
    mutProces = 60
    genome = pysam.FastaFile(genomeFile)
    chromLengths = dict(zip(genome.references, genome.lengths))
    clu = cluster(bed, chromLengths)
    clu.run()
