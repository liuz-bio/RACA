import math
import os
import pickle
import time
import copy
import warnings

import pandas as pd
import numpy as np
from dataclasses import dataclass
from intervaltree import IntervalTree, Interval
from clusterSig import readSig, clusterSig
from scipy.stats import ttest_ind
from sklearn.cluster import DBSCAN
from collections import Counter
from Bio import Align
from itertools import combinations
from util import ref_sigs
import math
from sklearn.cluster import KMeans


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


class callINS:
    def __init__(self, out_file, configer, sample, part_plk_files, shared_sig_trees, data_inv_dup, lock):
        self.aligner = None
        self.INSs = []
        self.out_file = out_file
        self.configer = configer
        self.sample = sample
        self.part_plk_files = part_plk_files
        self.tree = IntervalTree()
        start_t1 = time.time()
        shared_sig_trees = list(shared_sig_trees)
        self.subSigTrees = shared_sig_trees[1]
        self.SigTrees = shared_sig_trees[0]
        self.lock = lock
        self.data_inv_dup = data_inv_dup
        #self.subSigTrees = self.subSigRegion(sample)
        print("time_start1", time.time()-start_t1)
    
    def run(self):
        tSV = time.time()
        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = 'global'  # 全局比对
        self.aligner.match_score = 2  # 匹配得分
        self.aligner.mismatch_score = -1  # 不匹配得分
        self.aligner.open_gap_score = -1  # 开放缺失得分
        self.aligner.extend_gap_score = -0.5  # 扩展缺失得分
        # self.aligner.mode = 'global'  # 全局比对
        # self.aligner.match_score = 1  # 匹配得分
        # self.aligner.mismatch_score = 0  # 不匹配得分
        # self.aligner.open_gap_score = 0  # 开放缺失得分
        # self.aligner.extend_gap_score = 0  # 扩展缺失得分
        self.callINS(self.sample, self.part_plk_files)
        self.pass_overlap()
        print('tSV1', 'INS', time.time()-tSV)
        tSV = time.time()
        if len(self.INSs)>0:
            self.write_vcf()
        print('tSV2', 'INS', time.time()-tSV)

    
    def write_vcf(self):
        print("write_vcf INSs")
        with self.lock:
            with open(self.out_file, 'a') as out:
                out.write('\n'.join(['\t'.join(map(str, i)) for i in self.INSs])+'\n')
                out.flush()
    
    def pass_overlap(self):
        num = len(self.INSs)
        #clusterSigObj.chrom, clusterSigObj.start_mode, clusterSigObj.stop_mode, clusterSigObj.svlen_mode, clusterSigObj.support
        INSs = [i[1] for i in self.INSs]
        print("pass_overlap3:", len(INSs))
        for i in range(num-1):
            a = self.INSs[i]
            for ii in range(i+1,num):
                b = self.INSs[ii]
                if a[0][0] == b[0][0]:
                    start_a, start_b = a[0][1], b[0][1]
                    stop_a, stop_b = a[0][1]+a[0][2], b[0][1]+b[0][2]
                    if (max([start_a, start_b]) <= min([stop_a,stop_b])) and (abs(start_a - stop_a) < 100 and  abs(a[0][2]-b[0][2]) <100):
                        try:
                            print("pass_overlap1:",a[1], b[1] )
                            if b[0][2] < a[0][2]:
                                INSs.remove(b[1])
                            else:
                                INSs.remove(a[1])
                        except ValueError:
                            continue
                    elif abs(start_b-start_a) <100 or abs(stop_a-stop_b) <100:
                        if b[0][2] > a[0][2]:
                            try:
                                print("pass_overlap22:",a[1], b[1] )
                                INSs.remove(a[1])
                            except ValueError:
                                continue
                        else:
                            try:
                                print("pass_overlap11:",a[1], b[1] )
                                INSs.remove(b[1])
                            except ValueError:
                                continue
        self.INSs = INSs
    
    def call(self, sample, sigAllClusters, task_start, task_end):
        try:
            clusterTreeI = sorted(list(sigAllClusters['I']))
            clusterTreeL = sigAllClusters['L']
            clusterTreeR = sigAllClusters['R']
        except KeyError:
            return

        clusterTreeI = clusterTreeI[task_start:task_end]
        print("clusterTreeI:", len(clusterTreeI))
        for obj_node in clusterTreeI:
            obj_ins = obj_node.data
            if obj_ins.support < self.configer.min_support:
                 continue
            flag1 = True  #
            
            #flag1 = self.filterIncludedINS(obj_ins, clusterTreeI, clusterTreeL, clusterTreeR)
            print("filterIncludedINS flag1", flag1, obj_ins.chrom, obj_ins.start_mode, obj_ins.stop_mode, obj_ins.svlens)
            obj_ins.median_svlen = np.median(obj_ins.svlens)

            flag2, SVrecord = self.svRecord(obj_ins, clusterTreeL, clusterTreeR)
            print("filterIncludedINS flag2", flag2, obj_ins.chrom, obj_ins.start_mode, obj_ins.stop_mode, obj_ins.svlens)
            if flag1 & flag2:
                self.INSs.append(SVrecord)
    
    def global_alignmentx(self, seqs, support, obj_ins):
        for i in range(support-1):
            for ii in range(i+1,support):
                alignment = self.aligner.align(seqs[i], seqs[ii])
                if alignment.score/(sum([len(seqs[i]), len(seqs[ii])])/2) < 0.6:
                    support = support -1
                    if support <2:
                        return False
                else:
                    support += 1

        if support <2:
            return  False

        sv_lens = [len(i) for i in seqs]
        if (support*obj_ins.svlen_mode <= 90) or (np.std(sv_lens) / np.mean(sv_lens) > 0.2):
            return False
        else:
            return True
    
    def calculate_kmer_frequency(self, sequence, k):
        kmer_count = Counter()
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            kmer_count[kmer] += 1
        total_kmers = sum(kmer_count.values())
        kmer_frequency = {kmer: count / total_kmers for kmer, count in kmer_count.items()}
        return kmer_frequency
    
    def compare_kmer_frequency(self, seq1, seq2, k):
        freq1 = self.calculate_kmer_frequency(seq1, k)
        freq2 = self.calculate_kmer_frequency(seq2, k)
        common_kmers = set(freq1.keys()) & set(freq2.keys())
        similarity = sum(min(freq1[kmer], freq2[kmer]) for kmer in common_kmers)
        return similarity
    
    def kmeans(self, sequences, svlen):
        svlens = np.array([len(i) for i in sequences])
        sequences = np.array(sequences)
        if np.std(svlens) >0:
            kms = KMeans(n_clusters=2)
            kms.fit(svlens.reshape(-1, 1))
            #mks_cl = Counter(kms.labels_)
            if abs(np.mean(svlens[kms.labels_==0]) - svlen) < abs(np.mean(svlens[kms.labels_==1]) - svlen):
                return sequences[kms.labels_==0].tolist()
            else:
                return sequences[kms.labels_==1].tolist()
        else:
            return sequences.tolist()
    
    def overall_similarity(self, sequences, k, ratio, svlen, chrom, start, end, flag):
        similarity_sum = 0
        svlens = [len(i) for i in sequences]
        cv = np.std(svlens)/np.mean(svlens)
        if cv >0.1:
            sequences = self.kmeans(sequences, svlen)
            if len(sequences) == 1:
                if not flag and cv < 0.3:
                    print("NO1 overall_similarity True:",  len(sequences), chrom, start, end, svlens)
                    return True
                else:
                    print("NO2 overall_similarity False:",  len(sequences), chrom, start, end, svlens)
                    return False
            
        pairs = combinations(sequences, 2)
        num = 0
        for seq1, seq2 in pairs:
            similarity = self.compare_kmer_frequency(seq1, seq2, k)
            similarity_sum += similarity
            num += 1
        if similarity_sum / len(sequences) >ratio: #math.log(svlen)*
            print("overall_similarity True:", similarity_sum / len(sequences), similarity_sum, len(sequences), num, chrom, start, end, svlens)
            return True
        else:
            print("overall_similarity False:", similarity_sum / len(sequences), similarity_sum, len(sequences), num, chrom, start, end, svlens)
            return False
    
    def filterIncludedINS(self, obj_ins, clusterTreeI, clusterTreeL, clusterTreeR):
        chrom, start, end = obj_ins.chrom, obj_ins.start_mode, obj_ins.stop_mode
        seqs = [i.seq for i in obj_ins.cluster]
        k = int(math.log(obj_ins.svlen_mode)*3)
        if obj_ins.support * obj_ins.svlen_mode < self.configer.min_sv_size*0.7*obj_ins.depth*0.5:#300: 25x测序深度
            ratio = 0.9
            flag = True
            return True
            #return self.overall_similarity(seqs, k, ratio, obj_ins.svlen_mode, chrom, start, end, flag) #21
        else:
            flag = False
            ratio = 0.005
            try:
                inv_tree = self.data_inv_dup["INV"][obj_ins.chrom]
            except KeyError:
                inv_tree = None

            try:
                dup_tree = self.data_inv_dup["DUP"][obj_ins.chrom]
            except KeyError:
                dup_tree = None

            #print("inv_tree_dup", inv_tree, dup_tree)
            if inv_tree is None and dup_tree is None:
                #return self.global_alignmentx([i.seq for i in obj_ins.cluster], obj_ins.support, obj_ins)     
                return True
                #return self.overall_similarity(seqs, k, ratio, obj_ins.svlen_mode, chrom, start, end, flag)#21
            elif inv_tree is None:
                if len(dup_tree.overlap(start, end)) > 0:
                    return False
                else:
                    return True
                    #return self.overall_similarity(seqs, k, ratio, obj_ins.svlen_mode, chrom, start, end, flag)#21
            elif dup_tree is None:
                if len(inv_tree.overlap(start, end)) > 0:
                    return False
                else:
                    return self.overall_similarity(seqs, k, ratio, obj_ins.svlen_mode, chrom, start, end, flag)#21
            elif len(inv_tree.overlap(start, end)) > 0 or len(dup_tree.overlap(start, end)) > 0:
                return False
            else:
                return True
                #return self.overall_similarity(seqs, k, ratio, obj_ins.svlen_mode, chrom, start, end, flag)#21
        # sv_lens = [abs(i.start - i.stop) for i in obj_ins.cluster]
        # if (obj_ins.support * obj_ins.svlen_mode <= 90) or (np.std(sv_lens) / np.mean(sv_lens) > 0.2):
        #     return False
        # else:
        #     return True
    
    def callINS(self, sample, part_plk_files):
        #for plk_file in part_plk_files:
        for plk_file, task_start, task_end in part_plk_files:
            if os.path.exists(plk_file):
                with open(plk_file, 'rb') as plk:
                    sigAllClusters = pickle.load(plk)
                    self.call(sample, sigAllClusters, task_start, task_end)
    
    def svRecord(self, clusterSigObj, clusterTreeL, clusterTreeR):
        Ls_readNames = [ii for i in clusterTreeL.overlap(
            Interval(clusterSigObj.start_mode - 100, clusterSigObj.start_mode + 100)) for ii in i.data.readNames]
        Rs_readNames = [ii for i in clusterTreeR.overlap(
            Interval(clusterSigObj.start_mode - 100, clusterSigObj.start_mode + 100)) for ii in i.data.readNames]
        clusterSigObj.support = clusterSigObj.support + len(set(Ls_readNames) | set(Rs_readNames))

        flag = False
        format = "GT1:GT2:DR:DV"
        af = clusterSigObj.support / clusterSigObj.depth
        reliability = "IMPRECISE"
        infos = ';'.join([reliability,
                          "SVTYPE=INS",
                          "MEDIAN_SVLEN=" + str(clusterSigObj.median_svlen),
                          "mapQ_mean=" + str(clusterSigObj.mapQ_mean),
                          "mapQ_std=" + str(clusterSigObj.mapQ_std),
                          "SVLEN=" + str(clusterSigObj.svlen_mode),
                          "END=" + str(clusterSigObj.start_mode + 1),
                          "SUPPORT=" + str(clusterSigObj.support),
                          "AF=" + str(round(af, 3)),
                          "STD_POS1=" + str(round(clusterSigObj.start_std, 3)),
                          "STD_POS2=" + str(round(clusterSigObj.stop_std, 3)),
                          "STD_SVLEN=" + str(round(clusterSigObj.svlen_std, 3)),
                          "READNAMES="+'|'.join(clusterSigObj.readNames),
                          "RNAMES=NULL",
                          ])

        flag = self.filter(af,clusterSigObj)
        #flag = True
        return flag, [[clusterSigObj.chrom, clusterSigObj.start_mode, clusterSigObj.svlen_mode, clusterSigObj.support],[clusterSigObj.chrom, clusterSigObj.start_mode, "plotSV.INS", "N", clusterSigObj.alt.upper(), '60',
                      "PASS", infos, format, ':'.join(
                [clusterSigObj.gt1, clusterSigObj.gt2, str(clusterSigObj.dr), str(clusterSigObj.support)])]]
    
    def filter(self, af, clusterSigObj):
        if af >= 0.001 and clusterSigObj.svlen_mode>=30:
            return True
        else:
            return False
        pass


class callINS_Signals:
    def __init__(self, out_file, configer, sample, part_plk_files, shared_sig_trees, data_inv_dup, lock):
        self.aligner = None
        self.INSs = []
        self.out_file = out_file
        self.configer = configer
        self.sample = sample
        self.part_plk_files = part_plk_files
        self.tree = IntervalTree()
        start_t1 = time.time()
        shared_sig_trees = list(shared_sig_trees)
        self.subSigTrees = shared_sig_trees[1]
        self.SigTrees = shared_sig_trees[0]
        self.lock = lock
        self.data_inv_dup = data_inv_dup
        #self.subSigTrees = self.subSigRegion(sample)
        print("time_start1", time.time()-start_t1)
    
    def getSigAllClusters(self, contig, n1, n2, sample):
        plk_files = [str(os.path.join(self.configer.tmp_path, sample,
                                          contig + ".%d.%d.plk" % (
                                              n * self.configer.win_size + 1,
                                              (n + 1) * self.configer.win_size))) for n in [n1-1, n2]]
        tx = time.time()
        sigAllClusters, flag = self.mergeSigAllClusters(plk_files)
        print("getSigAllClusters1:", time.time() - tx)
        return sigAllClusters, flag
    
    def call(self, contig, start, end, sigAllClusters):
        posDeviation = 150
        posX = 500
        print("call:", contig, start, end)
        clusterTreeL = sigAllClusters['L']
        clusterTreeR = sigAllClusters['R']
        clusterTreeD = sigAllClusters['D']
        clusterTreeI = sigAllClusters['I']

        clusterR_nodes = clusterTreeR.overlap(end - posX, end + posX)
        clusterL_nodes = clusterTreeL.overlap(start - posX, start + posX)
        objR_nodes_size = len(clusterR_nodes)
        objL_nodes_size = len(clusterL_nodes)
            
        if objR_nodes_size>0 and objL_nodes_size >0:
            self.deal_R_L(contig, start, end, clusterTreeL, clusterTreeR, clusterTreeD, clusterTreeI, clusterL_nodes, clusterR_nodes, posDeviation)
        elif objR_nodes_size>0: 
            clusterL_nodes = self.sing_sig_L_R(contig, start, clusterTreeL, "L")
            if clusterL_nodes != None:
                self.deal_R_L(contig, start, end, clusterTreeL, clusterTreeR, clusterTreeD, clusterTreeI, clusterL_nodes, clusterR_nodes, posDeviation)
            print("clusterR_nodes=0", clusterR_nodes, clusterL_nodes)
        elif objL_nodes_size >0:
            clusterR_nodes = self.sing_sig_L_R(contig, end, clusterTreeR, "R")
            if clusterR_nodes != None:
                self.deal_R_L(contig, start, end, clusterTreeL, clusterTreeR, clusterTreeD, clusterTreeI, clusterL_nodes, clusterR_nodes, posDeviation)
            print("clusterL_nodes=0", clusterR_nodes, clusterL_nodes)
    
    def deal_R_L(self,contig, start, end, clusterTreeL, clusterTreeR, clusterTreeD, clusterTreeI, clusterL_nodes, clusterR_nodes, posDeviation):
        idx = '|'.join(map(str, [contig, start, end]))
        a = [idx]
        for clusterSigObjR_node in clusterR_nodes:
            for clusterSigObjL_node in clusterL_nodes:
                #logger.info("BBB")
                clusterSigObjL = clusterSigObjL_node.data
                clusterSigObjR = clusterSigObjR_node.data

                if clusterSigObjL.recordNu > 2*clusterSigObjL.support or  clusterSigObjR.recordNu > 2*clusterSigObjR.support:
                    continue
                overlapReadNameLen = len(clusterSigObjL.readNames & clusterSigObjR.readNames)

                if abs(clusterSigObjL.stop_mode-clusterSigObjR.start_mode) >posDeviation and overlapReadNameLen>=2:
                    continue

                pos_L = list(clusterSigObjL.PosLocation)
                pos_R = list(clusterSigObjR.PosLocation)
                pos_Ls = pd.Series(pos_L).value_counts()
                pos_Rs = pd.Series(pos_R).value_counts()

                clusterSigObjIns = copy.deepcopy(clusterSigObjR)
                clusterSigObjIns.raw_signal = '|'.join(map(str,[contig, start_x, end_x]))

                seq_svlens = [abs(abs(clusterSigObjR.start_mode - clusterSigObjL.start_mode) - abs(clusterSigObjL.query_starts[i][0] - clusterSigObjR.query_starts[i][0])) if clusterSigObjR.start_mode <= clusterSigObjL.start_mode else abs(abs(clusterSigObjR.start_mode - clusterSigObjL.start_mode) + abs(clusterSigObjL.query_starts[i][0] - clusterSigObjR.query_starts[i][0])) for i in readNames ]
                sv_biase = [abs(clusterSigObjL.start_mode - clusterSigObjR.start_mode) for i in readNames]
                genome_svlen = abs(start-end)
                sv_svlen = int(np.mean(seq_svlens))
                # svlens = [sv_svlen, genome_svlen]
                # print("query_readName_LR2",sv_biase,min(svlens)/max(svlens), query_readName_LR)
                # #svlen_biase = np.mean(sv_biase)
                # svlen_biase = sv_biase
                # if min(svlens)/max(svlens) >0.5:
                #     if svlen_biase > 100 and len(query_readName_LR) >0:
                #         print("WWWW2:", contig, start, end, svlen_biase, query_readName_LR)
                #         new_start, new_end =self.dupRepeat(contig, start, end, svlen_biase, int(np.mean(seq_svlens)), readID_Pos)
                #     print("QQQQ2:", new_start, new_end, svlen_biase)
                # else:
                #     print("out5")
                #     return

                

                # print("clusterSigObj", clusterSigObjL.stop_mode,      , [ i.start for i in clusterSigObjL.cluster],[ i.stop for i in clusterSigObjR.cluster])
                # print("clusterSigObj1", clusterSigObjL.stop_mode, clusterSigObjR.start_mode, posDeviation, overlapReadNameLen, self.configer.min_sv_size)
            
                # size_L = sum([ 1 for i in clusterSigObjR.cluster if start - i.sLR_other >10])
                # size_R = sum([ 1 for i in clusterSigObjL.cluster if i.sLR_other - end >10])
                # b = sum([1 for i in clusterSigObjL.cluster if i.query_length > 1.5*abs(start-end)])+sum([1 for i in clusterSigObjR.cluster if i.query_length > 1.5*abs(start-end)])                
                # if b>=4:
                #     if size_L >0 or size_R>0:
                #         if idx in a:
                #             a.remove(idx)
                #             print("ERWERW:",contig, start, end, a, clusterSigObjL.start_mode, clusterSigObjR.start_mode, [i.sLR_other for i in clusterSigObjL.cluster ], [i.sLR_other for i in clusterSigObjR.cluster ])
                #         self.mergeClusterSigObjGroupLS(contig, clusterSigObjL,
                #                                         clusterSigObjR, clusterTreeL, clusterTreeR,
                #                                         overlapReadNameLen, clusterTreeD, clusterTreeI, start, end)
                # else:
                #     self.mergeClusterSigObjGroupLS(contig, clusterSigObjL,
                #                                         clusterSigObjR, clusterTreeL, clusterTreeR,
                #                                         overlapReadNameLen, clusterTreeD, clusterTreeI, start, end)
                #     print("ERWERW:",contig, start, end, a, clusterSigObjL.start_mode, clusterSigObjR.start_mode, [i.sLR_other for i in clusterSigObjL.cluster ], [i.sLR_other for i in clusterSigObjR.cluster ])
    
    def sing_sig_L_R(self, contig, pos, clusterTreeLR, sigT):#补充获取无法聚类的R或L
        #bedFileName = "/home/liuz/work_space/project/liuz/plotSV/plotSV8/callSV/tmp/tmp/SIM/SIM.bed.gz"
        bedFileName = os.path.join(self.configer.tmp_path, sample, sample + ".bed.gz")
        read_sig_objs = []
        sig_xy = []
        with pysam.TabixFile(bedFileName) as bedfile:
            for record in bedfile.fetch(contig, pos-1, pos+1):
                read = record.strip().split("\t")
                if sigT == 'L':
                    new_pos = int(read[1])
                    query_start = int(read[3])
                    sLR_other = int(read[2])
                else:
                    new_pos = int(read[2])
                    query_start = int(read[4])
                    sLR_other = int(read[1])

                if abs(new_pos-pos) <100:
                    md5_hash = hashlib.md5()
                    md5_hash.update(record.strip().encode('utf-8'))
                    read_hash = md5_hash.hexdigest()
                    map_q = int(read[6])
                    read_is_reverse = int(read[7])
                    read_is_secondary = int(read[8])
                    read_name = read[9]
                    query_len = int(read[5])
                    read_sig_obj = readSig(new_pos, new_pos, map_q, read_hash, 
                                           read_is_reverse, read_is_secondary, read_name, 
                                           query_start=query_start,  query_len=query_len, sLR_other=sLR_other, 
                                           query_length=query_len)
                    read_sig_objs.append(read_sig_obj)
                    sig_xy.append([new_pos, new_pos])
                    #clusterSigObj = clusterSig(contig, np.array([[new_pos, new_pos]]), np.array([read_sig_obj]), None, sigT, 1234567)
                    #node = Interval(clusterSigObj.start_mode, clusterSigObj.stop_mode, clusterSigObj)
                    #clusterTreeLR.add(Interval(clusterSigObj.start_mode, clusterSigObj.stop_mode, clusterSigObj))
                    #return [node]
        
        if len(read_sig_objs) >0:
            clusterSigObj = clusterSig(contig, np.array(sig_xy), np.array(read_sig_objs), None, sigT, 1234567)
            node = Interval(clusterSigObj.start_mode, clusterSigObj.stop_mode, clusterSigObj)
            #clusterTreeLR.add(Interval(clusterSigObj.start_mode, clusterSigObj.stop_mode, clusterSigObj))
            clusterTreeLR.insert(clusterSigObj.start_mode, clusterSigObj.stop_mode, clusterSigObj)
            return [node]
        else:
            return None

    
    def get_mut_pos(self,sample):
        with open(os.path.join(self.configer.tmp_path, sample, sample + '.' + "INS" + '.sigs'), 'r') as f:
            for line in f:
                contig, start, end = line[0], int(line[1]), int(line[2])
                if start == end:
                    end += 1
                elif start >end:
                    start, end = int(line[2]), int(line[1]) 
                yield contig, start, end
                
    
    def callINS(self, sample):
        mut_pos = self.get_mut_pos(sample)
        print(mut_pos)
        while True:
            contig, start, end = next(mut_pos)
            head_n1 = int(start / self.configer.win_size)
            head_n2 = int(end / self.configer.win_size)
            sigAllClusters, flag = self.getSigAllClusters(contig, head_n1, head_n2, sample)
            if flag:
                self.call(chrom, start, end, sigAllClusters)
                break

        for pos in mut_pos:
            print("pos:", pos)
            contig, start, end = pos
            n1 = int(start / self.configer.win_size)
            n2 = int(end / self.configer.win_size)
            if head_n1 == n1 and head_n2 == n2:
                print("aaa:", head_n1, head_n2, flag, pos)
                self.call(contig, start, end, sigAllClusters)
            else:
                sigAllClusters_tmp, flag = self.getSigAllClusters(contig, n1, n2, sample)
                if flag:
                    head_n1, head_n2 = n1, n2
                    print("bbb:", head_n1, head_n2, flag, pos)
                    self.call(contig, start, end, sigAllClusters_tmp)
                    sigAllClusters = sigAllClusters_tmp
                else:
                    print("ccc:", n1, n2, flag)
                    continue