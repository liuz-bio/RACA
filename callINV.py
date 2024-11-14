import os
import pickle
import time
import copy
import warnings
import pysam
import pandas as pd
import numpy as np
import random
import hashlib
import statistics
from dataclasses import dataclass
from intervaltree import IntervalTree, Interval
from clusterSig import readSig, clusterSig
from scipy.stats import ttest_ind
from sklearn.cluster import DBSCAN
from collections import Counter
from sortedcontainers import SortedDict


@dataclass(frozen=True)
class Configuration:
    genome: str
    win_size: int
    bam_files: list
    sample_ids: list
    samples_to_bams: dict
    min_map_quality: int
    min_sv_signal: int
    min_sv_size: int
    sv_types: list
    cluster_param: dict
    mut_proces: int
    chroms: list
    chrom_lengths: dict
    out_path: str
    tmp_path: str
    bam_cram: bool


class callINV:
    def __init__(self, out_file, configer, sample, mut_pos, shared_sig_trees, lock):
        self.INVs = []
        self.out_file = out_file
        self.configer = configer
        self.tree = IntervalTree()
        self.clusterObjBnd = []
        start_t1 = time.time()
        shared_sig_trees = list(shared_sig_trees)
        self.subSigTree = shared_sig_trees[1]["INV"]
        #self.subSigTrees = self.subSigRegion(sample)
        print("time_start1", time.time()-start_t1)
        self.sample = sample
        self.mut_pos = mut_pos
        self.lock = lock
        bedFileName = "/home/lz/work_space/project/plotSV/RACA/plotSV9/callSV/hg38.inv2.bed.gz"
        self.invbedfile = pysam.TabixFile(bedFileName) 
    
    def run(self):
        tSV = time.time()
        print('tSVX', 'INV', time.time()-tSV, len(self.INVs))
        #self.callINV(self.sample, self.mut_pos)
        print('tSV1', 'INV', time.time()-tSV, len(self.INVs))
        tSV = time.time()
        if len(self.INVs)>0:
            self.write_vcf()
        print('tSV2', 'INV', time.time()-tSV)
    
    def write_vcf(self):
        print("write_vcf INVs", len(self.INVs))
        with self.lock:
            with open(self.out_file, 'a') as out:
                out.write('\n'.join(['\t'.join(map(str, i)) for i in self.INVs])+'\n')
                out.flush()
    
    def mergeSigAllClusters(self, plk_files):
        sigAllClusters = {'L': None, 'R': None, 'D': None, 'I': None, }
        tmp = []
        for plk_file in plk_files:
            print("plk_file", plk_file)
            if os.path.exists(plk_file):
                with open(plk_file, 'rb') as plk:
                    tmp.append(pickle.load(plk))

        if len(tmp) >= 2:
            for sigT in ['L', 'R', 'D', 'I']:
                merged_tree = tmp[0][sigT]
                for trees in tmp[1:]:
                    merged_tree.update(trees[sigT])
                sigAllClusters[sigT] = merged_tree
            return sigAllClusters, True
        if len(tmp) == 1:
            return tmp[0], True
        else:
            return None, False
    
    def getSigAllClusters(self, chrom, n1, n2, sample):
        plk_files = [str(os.path.join(self.configer.tmp_path, sample,
                                          chrom + ".%d.%d.plk" % (
                                              n * self.configer.win_size + 1,
                                              (n + 1) * self.configer.win_size))) for n in [n1-1, n2]]
        tx = time.time()
        sigAllClusters, flag = self.mergeSigAllClusters(plk_files)
        print("getSigAllClusters1:", time.time() - tx)
        return sigAllClusters, flag
    
        # if n1 == n2:
        #     plk_files = [str(os.path.join(self.configer.tmp_path, sample,
        #                                   chrom + ".%d.%d.plk" % (
        #                                       n * self.configer.win_size + 1,
        #                                       (n + 1) * self.configer.win_size))) for n in [n1-1, n2]]
        #     tx = time.time()
        #     sigAllClusters, flag = self.mergeSigAllClusters(plk_files)
        #     print("getSigAllClusters1:", time.time() - tx)
        #     return sigAllClusters, flag
        #     #plk_file = str(os.path.join(self.configer.tmp_path, sample,
        #     #                            chrom + ".%d.%d.plk" % (
        #     #                                n1 * self.configer.win_size + 1, (n1 + 1) * self.configer.win_size)))
        #     #if os.path.exists(plk_file):
        #     #    with open(plk_file, 'rb') as plk:
        #     #        sigAllClusters = pickle.load(plk)
        #     #        return sigAllClusters, True
        #     #else:
        #     #    return None, False

        # else:
        #     plk_files = [str(os.path.join(self.configer.tmp_path, sample,
        #                                   chrom + ".%d.%d.plk" % (
        #                                       n * self.configer.win_size + 1,
        #                                       (n + 1) * self.configer.win_size))) for n in [n1, n2]]
        #     tx = time.time()
        #     sigAllClusters, flag = self.mergeSigAllClusters(plk_files)
        #     print("getSigAllClusters1:", time.time() - tx)
        #     return sigAllClusters, flag
    
    def subMerge(self, groups):
        # 预处理数据类型转换
        groups = sorted(groups, key=lambda x:  x[0])
        groups_iter = iter(groups)
        head = next(groups_iter)
        new_groups = []
        flag = True

        for i in groups_iter:
            head_start, head_end = head
            next_start, next_end = i

            if head_start <= next_start <= head_end or head_start <= next_end <= head_end:
                tmp = [head_start, head_end, next_start, next_end]
                new_groups.append([min(tmp), max(tmp)])
                flag = False
            elif flag:
                new_groups.append(head)
            else:
                flag = True

            head = i

        # 处理最后一个元素
        new_groups.append(head)

        return new_groups
    
    def clusterSub(self,chrom, svType, record, dbscan, subSigTree):
        clusters = dbscan.fit_predict(record[:, :2])
        labels = np.unique(clusters)
        for label in labels:
            cluster_points = record[clusters == label]
            if label != -1:
                pos1, pos2 = np.mean(cluster_points[:, :2], axis=0).astype(int)
                #subSigTree.addi(pos1, pos2)
                subSigTree.insert(pos1, pos2)
                #groups.append([chrom, pos1, pos2, svType, len(cluster_points)])
                print("cluster test:", [chrom, pos1, pos2, svType, len(cluster_points)])
            else:
                # 优化点：避免for循环，用列表解析式代替
                for i  in cluster_points:
                    #subSigTree.addi(i[0], i[1])
                    subSigTree.insert(i[0], i[1])
                #groups.extend([[chrom, i[0], i[1], svType, 1] for i in cluster_points])
    
    def splitChromSub(self,chrom, record, subSigTree, dbscan, svType):
        # Assuming clusterBND and cluster functions are defined elsewhere and optimized already.
        step = 300
        max_distance = 1500
        record_array = np.array(record, dtype=int)
        record_array = record_array[np.argsort(record_array[:, 0])]
        record_size = len(record_array)

        if record_size <= 1000:
            self.clusterSub(chrom, svType, record_array, dbscan, subSigTree)
            return

        start_index = 0
        while start_index < record_size:
            stop_index = min(start_index + step, record_size)

            if stop_index >= record_size or (
                    np.linalg.norm(record_array[stop_index - 1, :2] - record_array[stop_index, :2]) > max_distance):
                segment = record_array[start_index:stop_index]
                start_index = stop_index
            else:
                end_index = np.argwhere(np.linalg.norm(
                    record_array[stop_index - 1:record_size - 1, :2] - record_array[stop_index:record_size, :2],
                    axis=1) > max_distance)
                if len(end_index) > 0:
                    end_index = end_index[0, 0] + stop_index
                    segment = record_array[start_index:end_index]
                    start_index = end_index
                else:
                    segment = record_array[start_index:]
                    start_index = record_size

            if len(segment) > 2 * step:
                segment = segment[:step]

            self.clusterSub(chrom, svType, segment, dbscan, subSigTree)
    
    def subSigRegion(self, sample):
        #tmp/tmp/SIM/SIM.sub.INV.sigs
        sub_file = str(os.path.join(self.configer.tmp_path, sample, sample + ".sub.INV.sigs"))
        chrom_lines = {}
        for i in open(sub_file, 'r'):
            i = i.strip().split('\t')
            chrom_lines.setdefault(i[0], []).append([int(i[1]), int(i[2])]) #INV

        dbscan = DBSCAN(eps=500, min_samples=2)
        subSigTrees = {}
        out = open('ttttt.out','w')
        for chrom, lines in chrom_lines.items():
            subSigTree = IntervalTree()
            groups = self.subMerge(lines)
            self.splitChromSub(chrom, groups, subSigTree, dbscan, 'INV')
            subSigTrees[chrom] = subSigTree
            for i in groups:
                out.write('\t'.join([chrom]+list(map(str, i)))+"\n")
        out.close()

        return subSigTrees
    
    def checkSubSig_old(self, chrom, start, end):#从内往外
        try:
            overlap_nodes = self.subSigTree[chrom].overlap(start, end)
        except KeyError:
            return None, None, True
        overlap_nodes_size = len(overlap_nodes)
        if overlap_nodes_size < 1 or overlap_nodes_size > 10:
            return None, None, True

        # 提前计算 begins 和 ends 列表，并排序
        begins = sorted([i.begin - start for i in overlap_nodes], reverse=True)
        ends = sorted([i.end - end for i in overlap_nodes])

        # 计算最小距离
        min_dis = min([abs(i + ii) for i in begins for ii in ends])

        if min_dis < 100:  # 对称性
            for i in range(len(begins)):
                for ii in range(len(ends)):
                    if abs(begins[i] + ends[ii]) == min_dis:
                        print("checkSubSig:", chrom, start, end, begins[i] + start, ends[ii] + end, False)
                        return begins[i] + start, ends[ii] + end, False
        else:
            return None, None, True
    
    def process_min_dataccc(self, min_dat1, threshold, start, end, is_start): #当候选区域与补充候选区域为非中心对称时，同时read无法完全跨越INV, 则可能存在需要进行单边边界扩展的情况，但也需要注意误判
        k_i = [i for i in min_dat1.keys() if i <= threshold]
        k_l = len(k_i)
        if k_l > 1:
            
            if is_start:
                candidate_starts = [min_dat1[i][0] for i in k_i]
                print("candidate_starts:", k_i, start, end, candidate_starts)
                new_val = float('inf')
                for i in range(k_l):
                    new_val = min([min_dat1[k_i[i]][0], start])
                return new_val, end, False
            else:
                candidate_ends = [min_dat1[i][1] for i in k_i]
                print("candidate_ends:", k_i, start, end, candidate_ends)
                new_val = 0
                for i in range(k_l):
                    new_val = max([min_dat1[k_i[i]][1], end])

                return start, new_val, False
        elif k_l == 1:
            if is_start:
                print("candidate_starts1:", k_i, start, end, [min_dat1[i][0] for i in k_i])
                if min_dat1[k_i[0]][0] < start:
                    return min_dat1[k_i[0]][0], end, False
            else:
                print("candidate_ends1:", k_i, start, end, [min_dat1[i][1] for i in k_i])
                if min_dat1[k_i[0]][1] >end:
                    return start, min_dat1[k_i[0]][1], False
        return None, None, True
    
    def process_min_data(self, min_data_s, min_data_e, threshold, start, end): #当候选区域与补充候选区域为非中心对称时，同时read无法完全跨越INV, 则可能存在需要进行单边边界扩展的情况，但也需要注意误判
        k_starts = [i for i in min_data_s.keys() if i <= threshold and min_data_s[i][0] <start]
        k_ends = [i for i in min_data_e.keys() if i <= threshold and min_data_e[i][1] >end]

        new_start, new_end = start, end
        for i in k_starts:
            new_start = min([new_start, min_data_s[i][0]])

        for j in k_ends:
            new_end = max([new_end, min_data_e[j][1]])

        return new_start, new_end, False

    
    def checkSubSig(self, chrom, start, end):
        try:
            print("overlap_nodes1",  chrom, start, end)
            overlap_nodes = self.subSigTree[chrom].overlap(start, end)
        except KeyError:
            return None, None, True

        print("overlap_nodes", len(overlap_nodes))
        overlap_nodes_size = len(overlap_nodes)
        if overlap_nodes_size < 1 or overlap_nodes_size > 20:
            return None, None, True

        begins = [i.begin for i in overlap_nodes]
        ends = [i.end for i in overlap_nodes]
        
        min_data = {}
        min_data_s = {}
        min_data_e = {}
        for i in begins:
            for j in ends:
                j_e = j-end
                i_s  = i-start
                min_data[abs(i_s+j_e)+random.random()] = [i,j]
                min_data_s[abs(i_s)+random.random()] = [i,j]
                min_data_e[abs(j_e)+random.random()] = [i,j]

        print("min_data",chrom, start, end, min_data)  
        print("min_data_s",min_data_s)
        print("min_data_e",min_data_e)
        min_dis  = min(min_data.keys())
        if min_dis <100:# 扩展的start和end, 与扩展前的保持中心对称
            return min_data[min_dis][0], min_data[min_dis][1], False
        if min(min_data_s.keys()) < 100 and min(min_data_e.keys()) < 100:
            return self.process_min_data(min_data_s, min_data_e, 100, start, end)
        else:
            return None, None, True
        
        #if min_dis <100:# 扩展的start和end, 与扩展前的保持中心对称
        #    return min_data[min_dis][0], min_data[min_dis][1], False
        #else:
        #    if min(min_data_s.keys()) < 100 and min(min_data_e.keys()) < 100:
        #        return self.process_min_data(min_data_s, 100, start, end, is_start=False)
        #    elif min(min_data_e.keys()) < 100 and  min(min_data_s.keys()) < 100:
        #        return self.process_min_data(min_data_e, 100, start, end, is_start=True)
        #    else:
        #        return None, None, True
    
    def cluster(self, posX, dbscan, sigT, readNames, readHashs, mapQs, reverseReads, secondaryReads, obj_sv, nu):
        clusters = dbscan.fit_predict(posX.reshape(-1, 1))
        labels = np.unique(clusters)
        posRL = {}
        for label in labels:
            tmp = clusters == label
            cluster_points = posX[tmp]
            if label != -1:
                pos = np.mean(cluster_points).astype(int)
                posRL[pos] = [cluster_points, readNames[tmp], readHashs[tmp], mapQs[tmp], reverseReads[tmp], secondaryReads[tmp]]

        print("posRL:",posX, clusters, posRL)
        if len(posRL)>=2:
            vs = np.array([len(v[0]) for v in posRL.values()])
            ks = np.array([k for k in posRL.keys()])
            pos1, pos2 = ks[np.argsort(vs)][-2:].tolist()
            obj1 = copy.deepcopy(obj_sv)
            obj2 = copy.deepcopy(obj_sv)

            nu += 1
            obj1.start_mode = pos1
            obj1.stop_mode = pos1+1
            obj1.sigType = sigT
            obj1.readNames = set(posRL[pos1][1])
            obj1.readHashs = set(posRL[pos1][2])
            obj1.mapQs = posRL[pos1][3]
            obj1.support = len(posRL[pos1][0])
            obj1.fingerPrint = nu
            obj1.PosLocation = np.array(posRL[pos1][0])
            obj1.reverseReadsNu = np.sum([ i for i in posRL[pos1][4] if i])
            obj1.secondaryReadsNu = np.sum([ i for i in posRL[pos1][5] if i])

            nu+=1
            obj2.start_mode = pos2
            obj2.stop_mode = pos2 + 1
            obj2.sigType = sigT
            obj2.readNames = set(posRL[pos2][1])
            obj2.readHashs = set(posRL[pos2][2])
            obj2.mapQs = posRL[pos2][3]
            obj2.support = len(posRL[pos2][0])
            obj2.fingerPrint = nu
            obj2.PosLocation = np.array(posRL[pos2][0])
            obj1.reverseReadsNu = np.sum([ i for i in posRL[pos2][4] if i])
            obj1.secondaryReadsNu = np.sum([i for i in posRL[pos2][5] if i])
            return obj1, obj2, True
        else:
            return None, None, False
    
    def reBuildCluster(self, clusterObjL, clusterObjR, start, end):
        clusterTreeL, clusterTreeR = IntervalTree(), IntervalTree()
        dbscan = DBSCAN(eps=30, min_samples=2)
        L = np.array([i.start for i in clusterObjL.cluster])
        R = np.array([i.stop for i in  clusterObjR.cluster])
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings("error")
                t_statistic, p_value = ttest_ind(L, R)
        except RuntimeWarning:
            p_value = 1
            print("RuntimeWarning:")

        if p_value >0.05:
            print("p_value:", p_value, L, R)
            readNameLs = np.array([i.readName for i in clusterObjL.cluster])
            readHashLs = np.array([i.readHash for i in clusterObjL.cluster])
            mapQLs = np.array([i.mapQ for i in clusterObjL.cluster])
            reverseReadLs = np.array([i.read_is_reverse for i in clusterObjL.cluster])
            secondaryReadLs = np.array([i.read_is_secondary for i in clusterObjL.cluster])
            obj1, obj2, flag = self.cluster(L, dbscan, 'L', readNameLs, readHashLs, mapQLs, reverseReadLs, secondaryReadLs, clusterObjL, 0)
            if flag:
                #clusterTreeL.addi(obj1.start_mode, obj1.stop_mode, obj1)
                #clusterTreeL.addi(obj2.start_mode, obj2.stop_mode, obj2)
                clusterTreeL.insert(obj1.start_mode, obj1.stop_mode, obj1)
                clusterTreeL.insert(obj2.start_mode, obj2.stop_mode, obj2)
                readNameRs = np.array([i.readName for i in clusterObjR.cluster])
                readHashRs = np.array([i.readHash for i in clusterObjR.cluster])
                mapQRs = np.array([i.mapQ for i in clusterObjR.cluster])
                reverseReadRs = np.array([i.read_is_reverse for i in clusterObjR.cluster])
                secondaryReadRs = np.array([i.read_is_secondary for i in clusterObjR.cluster])
                obj1, obj2, flag = self.cluster(R, dbscan, 'R', readNameRs, readHashRs, mapQRs, reverseReadRs, secondaryReadRs, clusterObjR, 10)
                if flag:
                    #clusterTreeR.addi(obj1.start_mode, obj1.stop_mode, obj1)
                    #clusterTreeR.addi(obj2.start_mode, obj2.stop_mode, obj2)
                    clusterTreeR.insert(obj1.start_mode, obj1.stop_mode, obj1)
                    clusterTreeR.insert(obj2.start_mode, obj2.stop_mode, obj2)
                    posDeviation = 50
                    for clusterObjR_node in clusterTreeR.overlap(end - posDeviation, end + posDeviation):
                        for clusterObjL_node in clusterTreeL.overlap(start - posDeviation, start + posDeviation):
                            clusterObjL = clusterObjL_node.data
                            clusterObjR = clusterObjR_node.data
                            overlapReadNames = len(set(clusterObjL.readNames) & set(clusterObjR.readNames))
                            if overlapReadNames < 1:
                                return None, None, None, None, False
                            else:
                                return clusterObjL, clusterObjR, clusterTreeL, clusterTreeR, True
        return None, None, None, None, False
    
    def invRepeat(self, contig, start, end):
        center = (start+end)/2
        try:
            records = self.invbedfile.fetch(contig, start, end)
        except ValueError:
            records = []
        for record in records:
            line = record.strip().split("\t")
            print("invRepeatx:", contig, start, end, line, int(line[1]), int(line[2]))
            if abs(int(line[-1])-center) <int(line[5]): #对称中心偏差
                print("invRepeat:", contig, start, end, line, int(line[1]), int(line[2]))
                r_s = int(line[1])
                r_e = int(line[2])
                if (r_s <= start <= r_e) & (r_s <= end <= r_e):#断点位于重复片段内
                    return True
                elif (r_s-start+r_e-end)/2 <int(line[5]):
                    print("invRepeatZ:", contig, start, end, line, int(line[1]), int(line[2]))
                    return True
        return False

    
    def mutSigClusterDeal(self, contig, start, end, clusterTreeD, clusterTreeI, clusterTreeL, clusterTreeR):
        #elf.invRepeat(contig, start, end)
        #   return None, None, False

        posDeviation = 150
        svlen_inv = abs(end - start)
        #obj_Ds = [D_node.data for D_node in clusterTreeD.overlap(start - posDeviation, end + posDeviation)]
        #obj_Is = [I_node.data for I_node in clusterTreeI.overlap(start - posDeviation, end + posDeviation)]
        obj_Ds = [D_node.data for D_node in clusterTreeD.overlap(start - posDeviation, end + posDeviation) ] #inv附近存在小片段的del和ins信号需要被排除
        obj_Is = [I_node.data for I_node in clusterTreeI.overlap(start - posDeviation, end + posDeviation) ]



        print("mutSigClusterDeal Ds", start, end, [i.start_mode for i in obj_Ds])
        print("mutSigClusterDeal Is", start, end, [i.start_mode for i in obj_Is])
        D_l = len(obj_Ds)
        I_l = len(obj_Is)
        if D_l >10 or I_l  >10:###
            return None, None, True 
        print("obj_Ds", obj_Ds, D_l, [i.svlens for i in obj_Ds])
        print("obj_Is", obj_Is, I_l, [i.svlens for i in obj_Is])

        center = (start+end)/2
        if D_l > 0 and I_l > 0:
            D_lens_ratio = [0 for i in obj_Ds  if min([i.svlen_mode,svlen_inv])/max([i.svlen_mode,svlen_inv]) > 0.1]
            I_lens_ratio = [0 for i in obj_Is  if min([i.svlen_mode,svlen_inv])/max([i.svlen_mode,svlen_inv]) > 0.1]
            print("D_lens_ratio, I_lens_ratio:", D_lens_ratio, I_lens_ratio)
            if len(D_lens_ratio) >0 or len(I_lens_ratio) >0:
                return None, None, True

            tmp = {}
            for a in obj_Ds:
                for  b in obj_Is:
                    # 计算长度差的绝对值
                    diff = abs(a.svlen_mode - b.svlen_mode)
                    if diff < 50:
                        new_start = min([a.start_mode, b.start_mode])
                        new_end   = max([a.start_mode, b.start_mode])
                        tmp[center*(new_end-new_start)/abs(end-start)] = [new_start, new_end]
                        print("tmp2", center*(new_end-new_start)/abs(end-start), new_start, new_end)

            print("tmp:", tmp)
            if len(tmp) >20:
                return None, None, True
            elif len(tmp) > 0:
                if self.invRepeat(contig, start, end): #and min(D_l, I_l) / max(D_l, I_l) >= 0.8 
                    new_start, new_end =  tmp[min(tmp.keys())]
                    print("mutSigClusterDeal Pos", new_start, new_end, False)
                    return new_start, new_end, False
                else:
                    new_start, new_end =  tmp[min(tmp.keys())]
                    print("mutSigClusterDeal Pos1", new_start, new_end, False)
                    return new_start, new_end, False
            else:
                D_pos_bp =  [0 for i in obj_Ds  if min([abs(i.start_mode-start), abs(i.start_mode-end)])  >50]
                I_pos_bp =  [0 for i in obj_Is  if min([abs(i.start_mode-start), abs(i.start_mode-end)])  >50]
                if len(D_pos_bp) >0 or len(I_pos_bp) >0:#if len(D_lens_ratio) >0 or len(I_lens_ratio) >0:
                    return None, None, True
                else:   
                    return None, None, False

        elif D_l > 0:
            if D_l > 1:
                for obj_D in obj_Ds:
                    if obj_D.svlen_mode/svlen_inv >= 0.8: # INV 信号内存在与INV长度相似的DEL信号，该INV被过滤
                        return None, None, True

                svlen_Ds = [i.svlen_mode for i in obj_Ds]
                if max(svlen_Ds)-min(svlen_Ds) >30  or min(svlen_Ds) >50: #左右两侧DEL信号长度偏差较大则改INV信号可能来自单独的两个ins信号
                    return None, None, True
                else:
                    Ds1 = np.sort([start - i.start_mode for i in obj_Ds])[::-1]
                    Ds2 = np.sort([end - i.start_mode for i in obj_Ds])
                    if min(np.abs(Ds1 + Ds2)) < 100:
                        return None, None, False
                    else:
                        return None, None, True
            elif obj_Ds[0].svlen_mode >= 50:
                if min([obj_Ds[0].svlen_mode,svlen_inv])/max([obj_Ds[0].svlen_mode,svlen_inv]) >0.1:
                    return None, None, True
                else:
                    return None, None, False
            else:
                return None, None, False


        elif I_l > 0:

            if I_l > 1:
                for obj_I in obj_Is:
                    if obj_I.svlen_mode/svlen_inv >= 0.8: # INV 信号内存在与INV长度相似的INS信号，该INV被过滤
                        return None, None, True
                
                svlen_Is = [i.svlen_mode for i in obj_Is]
                if max(svlen_Is)-min(svlen_Is) >30  or min(svlen_Is) >50: #左右两侧INS信号长度偏差较大则改INV信号可能来自单独的两个ins信号
                    return None, None, True
                else:
                    Is1 = np.sort([start - i.start_mode for i in obj_Is])[::-1]
                    Is2 = np.sort([end - i.start_mode for i in obj_Is])
                    if min(np.abs(Is1 + Is2)) < 100:
                        return None, None, False
                    else:
                        return None, None, True

            elif obj_Is[0].svlen_mode >= 50:
                if min([obj_Is[0].svlen_mode,svlen_inv])/max([obj_Is[0].svlen_mode,svlen_inv]) >0.1:
                    return None, None, True
                else:
                    return None, None, False
            else:
                return None, None, False

        else:
            obj_LRs_center = [abs(node_L.data.start_mode-node_R.data.start_mode) for node_L in clusterTreeL.overlap(start - posDeviation, end + posDeviation)  for node_R in clusterTreeR.overlap(start - posDeviation, end + posDeviation)] #计算候选区域的所有左侧和右侧信号中心点， 然后计算其众数， 真实的INV的中心点和其众数接近

            obj_LRs_center_diff = [abs(obj_LRs_center[i]-obj_LRs_center[j]) for i in range(len(obj_LRs_center)-1) for j in range(1,len(obj_LRs_center)) if abs(obj_LRs_center[i]-obj_LRs_center[j]) <50]
            print("obj_LRs_center_diff:", obj_LRs_center_diff, obj_LRs_center)
            if len(obj_LRs_center_diff)>0:
                return None, None, False
            else:
                return None, None, True
    
    def call(self, contig, start, end, sigAllClusters):
        posDeviation = 150
        print("call:", contig, start, end, "INV")
        try:
            clusterTreeL = sigAllClusters['L']
        except:
            pass
        try:
            clusterTreeR = sigAllClusters['R']
        except:
            pass
        try:
            clusterTreeD = sigAllClusters['D']
        except:
            pass
        try:
            clusterTreeI = sigAllClusters['I']
        except:
            pass

        clusterR_nodes = clusterTreeR.overlap(end - posDeviation, end + posDeviation)
        clusterL_nodes = clusterTreeL.overlap(start - posDeviation, start + posDeviation)

        if len(clusterR_nodes)>0 and len(clusterL_nodes) >0:
            self.deal_R_L(contig, start, end, clusterTreeL, clusterTreeR, clusterTreeD, clusterTreeI, clusterL_nodes, clusterR_nodes, posDeviation)
        elif len(clusterR_nodes)>0: 
            clusterL_nodes = self.sing_sig_L_R(contig, start, clusterTreeL, "L")
            if clusterL_nodes != None:
                self.deal_R_L(contig, start, end, clusterTreeL, clusterTreeR, clusterTreeD, clusterTreeI, clusterL_nodes, clusterR_nodes, posDeviation)
            print("clusterR_nodes=0", clusterR_nodes, clusterL_nodes)
        elif len(clusterL_nodes) >0:
            clusterR_nodes = self.sing_sig_L_R(contig, end, clusterTreeR, "R")
            if clusterR_nodes != None:
                self.deal_R_L(contig, start, end, clusterTreeL, clusterTreeR, clusterTreeD, clusterTreeI, clusterL_nodes, clusterR_nodes, posDeviation)
            print("clusterL_nodes=0", clusterR_nodes, clusterL_nodes)
    
    def deal_R_L(self,contig, start, end, clusterTreeL, clusterTreeR, clusterTreeD, clusterTreeI, clusterL_nodes, clusterR_nodes, posDeviation):
        for clusterR_node in clusterR_nodes:
            for clusterL_node in clusterL_nodes:
                clusterObjL = clusterL_node.data
                clusterObjR = clusterR_node.data
                overlapReadNames = len(clusterObjL.readNames & clusterObjR.readNames)
                print("obj_sv", clusterObjL.stop_mode, clusterObjR.start_mode, [ i.start for i in clusterObjL.cluster],[ i.stop for i in clusterObjR.cluster])
                print("obj_sv1", clusterObjL.stop_mode, clusterObjR.start_mode, posDeviation, clusterObjL.support, clusterObjR.support, clusterObjL.reverseReadsNu, clusterObjR.reverseReadsNu, overlapReadNames, self.configer.min_sv_size)

                if overlapReadNames >= 1:
                    center_e1 = (clusterObjL.start_mode+clusterObjR.stop_mode)/2
                    center_e2 = (start+end)/2
                    print("INV center:", clusterObjR.start_mode, clusterObjL.stop_mode, center_e1, center_e2, abs(center_e1-center_e2))
                else:
                    continue

                if clusterObjL.stop_mode > clusterObjR.start_mode or np.abs(clusterObjL.stop_mode - clusterObjR.start_mode) <= self.configer.min_sv_size :
                    if np.abs(clusterObjL.stop_mode - clusterObjR.start_mode) <150:
                        print("obj_sv rebuild", clusterObjL.stop_mode, clusterObjR.start_mode)
                        clusterObjL, clusterObjR, TreeL, TreeR, fg = self.reBuildCluster(clusterObjL, clusterObjR, start, end)
                        print("obj_sv rebuild", fg)
                        if fg:
                            print("call elif:", contig, start, end)
                            self.mergeClusterSigObjGroupLS(posDeviation,
                                                            clusterObjL,
                                                            clusterObjR, TreeL,
                                                            TreeR,
                                                            overlapReadNames, clusterTreeD, clusterTreeI)
                else:
                    self.mergeClusterSigObjGroupLS(posDeviation, clusterObjL,
                                                clusterObjR, clusterTreeL, clusterTreeR,
                                                overlapReadNames, clusterTreeD, clusterTreeI)


    
    def sing_sig_L_R(self, contig, pos, clusterTreeLR, sigT):#补充获取无法聚类的R或L
        bedFileName = os.path.join(self.configer.tmp_path, self.sample, self.sample + ".bed.gz")
        with pysam.TabixFile(bedFileName) as bedfile:
            for record in bedfile.fetch(contig, pos-1, pos+1):
                read = record.strip().split("\t")
                if sigT == 'L':
                    new_pos = int(read[1])
                    query_start = int(read[3])
                else:
                    new_pos = int(read[2])
                    query_start = int(read[4])

                if abs(new_pos-pos) <100:
                    md5_hash = hashlib.md5()
                    md5_hash.update(record.strip().encode('utf-8'))
                    read_hash = md5_hash.hexdigest()
                    map_q = int(read[6])
                    read_is_reverse = int(read[7])
                    read_is_secondary = int(read[8])
                    read_name = read[9]
                    read_sig_obj = readSig(new_pos, new_pos, map_q, read_hash, read_is_reverse, read_is_secondary, read_name, query_start=query_start)
                    clusterSigObj = clusterSig(contig, np.array([[new_pos, new_pos]]), np.array([read_sig_obj]), None, sigT, 1234567)
                    node = Interval(clusterSigObj.start_mode, clusterSigObj.stop_mode, clusterSigObj)
                    #clusterTreeLR.add(Interval(clusterSigObj.start_mode, clusterSigObj.stop_mode, clusterSigObj))
                    clusterTreeLR.insert(clusterSigObj.start_mode, clusterSigObj.stop_mode, clusterSigObj)
                    return [node]
        return None
                    
            
    
    def callINV(self, sample, mut_pos):
        print("callINV123",mut_pos)
        mut_pos = iter(mut_pos)

        while True:
            chrom, start, end = next(mut_pos)
            head_n1 = int(start / self.configer.win_size)
            head_n2 = int(end / self.configer.win_size)
            sigAllClusters, flag = self.getSigAllClusters(chrom, head_n1, head_n2, sample)
            if flag:
                self.call(chrom, start, end, sigAllClusters)
                break


        for pos in mut_pos:
            print("pos:", pos)
            chrom, start, end = pos
            n1 = int(start / self.configer.win_size)
            n2 = int(end / self.configer.win_size)
            if head_n1==n1 and head_n2==n2:
                print("aaa:", head_n1, head_n2, flag, pos)
                self.call(chrom, start, end, sigAllClusters)
            else:
                sigAllClusters_tmp, flag = self.getSigAllClusters(chrom, n1, n2, sample)
                if flag:
                    head_n1, head_n2 = n1, n2
                    print("bbb:", head_n1, head_n2, flag, pos)
                    self.call(chrom, start, end, sigAllClusters_tmp)
                    sigAllClusters = sigAllClusters_tmp
                else:
                    print("ccc:", n1, n2, flag)
                    continue
    
    def min_abs_sum(self, no_target_Ls, no_target_Rs, clusterObjL_start, clusterObjR_start):
        if len(no_target_Ls) >5 or len(no_target_Rs) >5: ###
            return 1000
        A = np.array([i.start_mode - clusterObjL_start for i in no_target_Ls])
        B = np.array([i.start_mode - clusterObjR_start for i in no_target_Rs])
        print("min_abs_sum", A, B)
        print("min_abs_sum1", [i.start_mode for i in no_target_Ls], [i.start_mode for i in no_target_Rs])
    
        sum_A = A[:, np.newaxis] + B
        abs_sum = np.abs(sum_A)
    
        # 获取所有的绝对值和及其索引
        indices = np.dstack(np.unravel_index(np.argsort(abs_sum.ravel()), abs_sum.shape))[0]
    
        # 提取前3个最小值及其对应的索引
        if len(indices)>3:
            min_values = abs_sum[indices[:, 0], indices[:, 1]][:3]
            min_indices = indices[:3]
        else:
            min_values = abs_sum[indices[:, 0], indices[:, 1]]
            min_indices = indices

        for index, i in enumerate(min_indices):
            if no_target_Ls[i[0]].start_mode > no_target_Rs[i[1]].start_mode:
                return min_values[index]
        return 1000

       
    
    def mergeClusterSigObjGroupLS(self, posDeviation, clusterObjL, clusterObjR, clusterTreeL,
                                  clusterTreeR, overlapReadNames, clusterTreeD, clusterTreeI):
        pos_L = list(clusterObjL.PosLocation)
        pos_R = list(clusterObjR.PosLocation)
        pos_Ls = pd.Series(pos_L).value_counts()
        pos_Rs = pd.Series(pos_R).value_counts()

        obj_inv = copy.deepcopy(clusterObjR)
        overlapReadHashs = len(clusterObjL.readHashs & clusterObjR.readHashs)
        overlapReadHashName = [len(clusterObjL.readHashs | clusterObjR.readHashs), len(clusterObjL.readNames | clusterObjR.readNames)]
        obj_inv.readNamesRLNu = overlapReadNames  #
        obj_inv.readHashsRLNu = overlapReadHashs
        obj_inv.readNoHashRLNus = [clusterObjL.support, overlapReadHashs, clusterObjR.support,
                                            overlapReadNames, overlapReadHashName]

        start = time.time()

        svlen = pos_Rs.idxmax() - pos_Ls.idxmax()

        if svlen < 30:
            self.clusterObjBnd.append([clusterObjL, clusterObjR])
            # print("obj_inv.svlen_mode:", pos_Rs.idxmax(), pos_Ls.idxmax())
            return False, None

        svlens = [np.abs([pos_L[i] - pos_R[i]]) for i in range(np.min([len(pos_L), len(pos_R)]))]
        obj_inv.svlen_std = np.std(svlens)
        if obj_inv.svlen_std > 50:
            rangeSize = int(obj_inv.svlen_std)
        else:
            rangeSize = 50

        Rs = [i.data.support for i in clusterTreeR.overlap(Interval(np.min(pos_L) - rangeSize, np.min(pos_R) + rangeSize))
              if i.data.fingerPrint != clusterObjR.fingerPrint]
        Ls = [i.data.support for i in clusterTreeL.overlap(Interval(np.min(pos_L) - rangeSize, np.min(pos_R) + rangeSize))
              if i.data.fingerPrint != clusterObjL.fingerPrint]

        #除输入左侧和右侧断点外的其他右侧和左侧断点与其对应的断点应该来自不同的比对记录
        no_target_Rs = [i.data for i in clusterTreeR.overlap(Interval(pos_Ls.idxmax() - rangeSize, pos_Rs.idxmax() + rangeSize)) if (i.data.fingerPrint != clusterObjR.fingerPrint) and (len(i.data.readHashs & clusterObjL.readHashs)/len(i.data.readHashs) <0.7)]
        no_target_Ls = [i.data for i in clusterTreeL.overlap(Interval(pos_Ls.idxmax() - rangeSize, pos_Rs.idxmax() + rangeSize)) if (i.data.fingerPrint != clusterObjL.fingerPrint) and  (len(i.data.readHashs & clusterObjR.readHashs)/len(i.data.readHashs) <0.7)]
        #ax = [[i.data.start_mode, i.data.start_mode-clusterObjR.start_mode, len(i.data.readHashs & clusterObjL.readHashs)/len(i.data.readHashs), i.data.readHashs,  clusterObjL.readHashs] for i in clusterTreeR.overlap(Interval(pos_Ls.idxmax() - rangeSize, pos_Rs.idxmax() + rangeSize)) if (i.data.fingerPrint != clusterObjR.fingerPrint)]
        #bx = [[i.data.start_mode, i.data.start_mode-clusterObjL.start_mode, len(i.data.readHashs & clusterObjR.readHashs)/len(i.data.readHashs), i.data.readHashs,  clusterObjR.readHashs] for i in clusterTreeL.overlap(Interval(pos_Ls.idxmax() - rangeSize, pos_Rs.idxmax() + rangeSize)) if (i.data.fingerPrint != clusterObjL.fingerPrint)]
        #print("no_target_Rs_ax, no_target_Ls_bx:", ax, bx)
        Rs_Ls_start_mean = None
        Rs_start_mean = None
        Ls_Rs_start_mean = None
        Ls_start_mean = None

        if len(no_target_Rs) >0 and len(no_target_Ls) >0:
            if self.min_abs_sum(no_target_Ls, no_target_Rs, clusterObjL.start_mode, clusterObjR.start_mode) <= 300:
                start = pos_Ls.idxmax()
                end = pos_Rs.idxmax()
                start_t = time.time()
                ##############################################################################################
                sub_start, sub_end, flag1 = self.checkSubSig(obj_inv.chrom, start, end)
                print("checkSubSig1: ",sub_start, sub_end, flag1, obj_inv.chrom, start, end)
                if not flag1:
                    print("checkSubSig: ",sub_start, sub_end, flag1, obj_inv.chrom, start, end)
                    if sub_start <= start  and  end <= sub_end: #sub_start start end sub_end 扩展边界，需要进行对称性检查
                        #内部存在对称性，则从新定义边界
                        new_start, new_end, flag2 = self.mutSigClusterDeal(obj_inv.chrom, sub_start, sub_end, clusterTreeD, clusterTreeI, clusterTreeL, clusterTreeR)
                        if new_start == None and not flag2:
                            new_start, new_end = sub_start, sub_end
                    else:
                        new_start, new_end, flag2 = self.mutSigClusterDeal(obj_inv.chrom, start, end, clusterTreeD, clusterTreeI, clusterTreeL, clusterTreeR)
                else:
                    new_start, new_end, flag2 = self.mutSigClusterDeal(obj_inv.chrom, start, end, clusterTreeD, clusterTreeI, clusterTreeL, clusterTreeR)

                print("mutSigClusterDeal", new_start, new_end, flag2)
                if flag2:
                    return
                elif new_start != None:
                    obj_inv.start_mode = new_start
                    obj_inv.stop_mode = new_end
                else:
                    obj_inv.start_mode = pos_Ls.idxmax()
                    obj_inv.stop_mode = pos_Rs.idxmax()

                obj_inv.SVTYPE = 'INV'
                obj_inv.Rs = Rs
                obj_inv.Ls = Ls
                obj_inv.Rs_start_mean = Rs_start_mean
                obj_inv.Ls_start_mean = Ls_start_mean
                obj_inv.Rs_Ls_start_mean = Rs_Ls_start_mean
                obj_inv.Ls_Rs_start_mean = Ls_Rs_start_mean

                svlens = [np.abs([pos_L[i] - pos_R[i]]) for i in range(np.min([len(pos_L), len(pos_R)]))]
                #pos12 = [pos_Ls.idxmax(), pos_Rs.idxmax()]
                obj_inv.svlen_mode = svlen
                obj_inv.reverseReadsNus = [clusterObjL.reverseReadsNu, obj_inv.reverseReadsNu]
                obj_inv.secondaryReadsNus = [clusterObjL.secondaryReadsNu, obj_inv.secondaryReadsNu]
                mapQs = obj_inv.mapQs + clusterObjL.mapQs
                obj_inv.mapQ_mean = np.mean(mapQs)
                obj_inv.mapQ_std = np.std(mapQs)

                obj_inv.start_std = np.std(pos_L)
                obj_inv.stop_std = np.std(pos_R)
                obj_inv.svlen_std = np.std(svlens)
                #obj_inv.support = np.max([clusterObjL.support, clusterObjR.support])
                obj_inv.support = np.max([clusterObjL.support, clusterObjR.support,  overlapReadHashName[-1]])
                obj_inv.depth = np.max([clusterObjL.depth, clusterObjR.depth])
                #obj_inv.getDepth(np.min(pos_R + pos_L), np.max(pos_R + pos_L), bedGz)
                obj_inv.dr = obj_inv.depth - obj_inv.support
                # print("start2:", time.time() - start)
                start = time.time()
                # print("DV,DP:",obj_inv.support, obj_inv.depth)
                obj_inv.getGenotype_cuteSV()
                obj_inv.getGenotype_sniffles()
                # print("start3:", time.time() - start)


                try:
                    #self.tree.add(Interval(obj_inv.start_mode, obj_inv.stop_mode, obj_inv))
                    self.tree.insert(obj_inv.start_mode, obj_inv.stop_mode, obj_inv)
                except ValueError:
                    pass
                print("obj_inv3")
                flag, sv_record = self.svRecord(obj_inv)
                print("obj_inv4", flag, sv_record)
                if flag:
                    self.INVs.append(sv_record)
                    print("obj_inv5", flag, sv_record, len(self.INVs))
            else:
                print('None2', Rs_Ls_start_mean, Ls_Rs_start_mean, posDeviation)
        else:
            print('None1', Rs_Ls_start_mean, Ls_Rs_start_mean)

    
    def svRecord(self, obj_sv):
        flag = False
        SVid = "plotSV." + obj_sv.SVTYPE
        format = "GT1:GT2:DR:DV"
        af = obj_sv.support / obj_sv.depth
        reliability = "IMPRECISE"
        infos = ';'.join([reliability,
                          "SVTYPE=" + obj_sv.SVTYPE,
                          # "LNus="+str(obj_sv.LNus),
                          # "RNus="+str(obj_sv.RNus),
                          #"Ls_Rs_start_mean=" + str(obj_sv.Ls_Rs_start_mean),
                          #"Rs_Ls_start_mean=" + str(obj_sv.Rs_Ls_start_mean),
                          #"Ls_start_mean=" + str(obj_sv.Ls_start_mean),
                          #"Rs_start_mean=" + str(obj_sv.Rs_start_mean),
                          #"Ls=" + str(obj_sv.Ls),
                          #"Rs=" + str(obj_sv.Rs),
                          #"readNoHashRLNus=" + str(obj_sv.readNoHashRLNus),
                          #"readNamesRLNu=" + str(obj_sv.readNamesRLNu),
                          #"readHashsRLNu=" + str(obj_sv.readHashsRLNu),
                          "mapQ_mean=" + str(obj_sv.mapQ_mean),
                          "mapQ_std=" + str(obj_sv.mapQ_std),
                          "ReverseReadsNu=" + ','.join(map(str, obj_sv.reverseReadsNus)),
                          "SecondaryReadsNu=" + ','.join(map(str, obj_sv.secondaryReadsNus)),
                          "SVLEN=" + str(obj_sv.svlen_mode),
                          "END=" + str(obj_sv.stop_mode),
                          "SUPPORT=" + str(obj_sv.support),
                          "AF=" + str(round(af, 3)),
                          "STD_POS1=" + str(round(obj_sv.start_std, 3)),
                          "STD_POS2=" + str(round(obj_sv.stop_std, 3)),
                          "STD_SVLEN=" + str(round(obj_sv.svlen_std, 3)),
                          "RNAMES=NULL",
                          ])

        flag = self.filter(af, obj_sv)
        # if obj_sv.Ls_Rs_start_mean is None or obj_sv.Rs_Ls_start_mean is None:
        #    flag = False
        # elif abs(obj_sv.Ls_Rs_start_mean + obj_sv.Rs_Ls_start_mean) >100:
        #    flag = False
        return flag, [obj_sv.chrom, obj_sv.start_mode, SVid, "N", "<" + obj_sv.SVTYPE + ">", '60',
                      "PASS", infos, format, ':'.join(
                [obj_sv.gt1, obj_sv.gt2, str(obj_sv.dr), str(obj_sv.support)])]
    
    def filter(self, af, obj_sv):
        if af >= 0.001 and obj_sv.mapQ_mean>= 55:
            return True
        else:
            return False
        pass

if __name__ == '__main__':
    genome = "/home/lz/work_space/Database/VISOR_test/Genome_old/Genome.fa"
    #bamFiles = ["/public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/Pipelines/CCS/Sim.DEL.INS.DUP.INV.25x/winnowmap/bam/Sim.DEL.INS.DUP.INV.25x.sort.bam"]
    bamFiles = ["/home/lz/work_space/Database/VISOR_test/Chr1_two_new/T2/FY0/VISOR_LASeR_Father_DEL.DUP.INS.INV/sim.srt.bam"]
    #sampleIDs = ["SIM"]
    sampleIDs = ["FY0_Father"]
    sv_types = ["DEL", "INS", "INV", "DUP", "BND"]
    cluster_param = {'L': {'min_distanc': 150, 'min_signals': 2},
                     'R': {'min_distanc': 150, 'min_signals': 2},
                     'D': {'min_distanc': 300, 'min_signals': 2},
                     'I': {'min_distanc': 300, 'min_signals': 2},
                     'A': {'min_distanc': 50, 'min_signals': 5},
                     'E': {'min_distanc': 50, 'min_signals': 5}}
    chrom_lengths = dict([i.strip().split('\t')[:2] for i in open(genome + ".fai").read().strip().split('\n')])
    out_path = "/home/lz/work_space/project/plotSV/plotSV5/callSV/tmp"
    tmp_path = os.path.join(out_path, "tmp")
    os.makedirs(tmp_path, exist_ok=True)

    configer = Configuration(genome=genome,
                             win_size=1000000,
                             bam_files=bamFiles,
                             sample_ids=sampleIDs,
                             samples_to_bams=dict(zip(sampleIDs, bamFiles)),
                             min_map_quality=20,
                             min_sv_signal=25,
                             min_sv_size=30,
                             sv_types=sv_types,
                             cluster_param=cluster_param,
                             mut_proces=60,
                             chroms=list(chrom_lengths.keys()),
                             chrom_lengths=chrom_lengths,
                             out_path=out_path,
                             tmp_path=tmp_path,
                             bam_cram=True,
                             )

    sample = "SIM"
    mut_pos = 'aa'
    inv = callINV(configer, sample, mut_pos)
    pass
