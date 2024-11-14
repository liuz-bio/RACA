import os
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
from util import ref_sigs


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


class callDEL:
    def __init__(self, out_file, configer, sample, part_plk_files, shared_sig_trees, data_inv_dup, lock):
        self.DELs = []
        self.out_file = out_file
        self.sample = sample
        self.part_plk_files = part_plk_files
        self.configer = configer
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
        self.callDEL(self.sample, self.part_plk_files)
        print('tSV1', 'DEL', time.time()-tSV)
        tSV = time.time()
        self.pass_overlap()
        if len(self.DELs)>0:
            self.write_vcf()
        print('tSV2', 'DEL', time.time()-tSV)
    
    
    def pass_overlap(self):
        num = len(self.DELs)
        #clusterSigObj.chrom, clusterSigObj.start_mode, clusterSigObj.stop_mode, clusterSigObj.svlen_mode, clusterSigObj.support
        DELs = [i[1] for i in self.DELs]
        print("pass_overlap3:", len(DELs))
        for i in range(num-1):
            a = self.DELs[i]
            for ii in range(i+1,num):
                b = self.DELs[ii]
                if a[0][0] == b[0][0]:
                    start_a, start_b = a[0][1], b[0][1]
                    stop_a, stop_b = a[0][1]+a[0][2], b[0][1]+b[0][2]
                    if (max([start_a, start_b]) <= min([stop_a,stop_b])) and (abs(start_a - stop_a) < 100 and  abs(a[0][2]-b[0][2]) <100):
                        try:
                            print("pass_overlap1:",a[1], b[1] )
                            if b[0][2] < a[0][2]:
                                DELs.remove(b[1])
                            else:
                                DELs.remove(a[1])
                        except ValueError:
                            continue
                    elif abs(start_b-start_a) <100 or abs(stop_a-stop_b) <100:
                        if b[0][2] > a[0][2]:
                            try:
                                print("pass_overlap22:",a[1], b[1] )
                                DELs.remove(a[1])
                            except ValueError:
                                continue
                        else:
                            try:
                                print("pass_overlap11:",a[1], b[1] )
                                DELs.remove(b[1])
                            except ValueError:
                                continue
        self.DELs = DELs
       
    def write_vcf(self):
        print("write_vcf DELs")
        with self.lock:
            with open(self.out_file, 'a') as out:
                out.write('\n'.join(['\t'.join(map(str, i)) for i in self.DELs])+'\n')
                out.flush()
    
    def call_RL(self, clusterTreeD, clusterTreeI, clusterTreeL, clusterTreeR):
        clusterTreeL_nodes = sorted(list(clusterTreeL))
        clusterTreeR_nodes = sorted(list(clusterTreeR))
        for obj_node_R in clusterTreeR_nodes:
            for obj_node_L in clusterTreeL_nodes:
                obj_R = obj_node_R.data
                obj_L = obj_node_L.data
                start, stop = obj_R.start_mode, obj_L.start_mode
                
                if start >= stop:
                    continue
                overlap_reads = obj_L.readNames & obj_R.readNames
                if len(overlap_reads) >=2:
                    overlap_reads_size = len([0 for i in overlap_reads if abs(obj_L.query_starts[i][0]-obj_R.query_starts[i][0])< 0.1*abs(start-stop)])
                    if overlap_reads_size >=2:
                        size_D = len([0 for i in clusterTreeD.overlap(start+100, stop-100)])# if i.data.svlen_mode >30])
                        size_I = len([0 for i in clusterTreeI.overlap(start+100, stop-100)])# if i.data.svlen_mode >30])
                        size_L = len([0 for i in clusterTreeL.overlap(start+100, stop-100)])
                        size_R = len([0 for i in clusterTreeR.overlap(start+100, stop-100)])
                        size_all = size_D+size_I+size_L+size_R
                        if size_all <1:
                            obj_del = copy.deepcopy(obj_R)
                            obj_del.start_mode = start
                            obj_del.stop_mode = stop
                            obj_del.svlen_mode = abs(start-stop)
                            svlens = [ abs(i-j) for i in obj_L.PosLocation for j in obj_R.PosLocation]
                            obj_del.median_svlen = np.median(svlens)
                            
                            obj_del.svlen_std = np.std(svlens)
                            obj_del.SVTYPE = 'DEL'
                            obj_del.readNames = obj_L.readNames | obj_R.readNames

                            #clusterSigObjDel.reverseReadsNus = [obj_L.reverseReadsNu, obj_R.reverseReadsNu]
                            #clusterSigObjDel.secondaryReadsNus = [obj_L.secondaryReadsNu, obj_R.secondaryReadsNu]
                            mapQs = obj_R.mapQs + obj_L.mapQs
                            obj_del.mapQ_mean = np.mean(mapQs)
                            obj_del.mapQ_std = np.std(mapQs)
                            
                            obj_del.start_std = np.std(obj_L.PosLocation)
                            obj_del.stop_std = np.std(obj_R.PosLocation)
                        
                            obj_del.support = np.max([obj_L.support, obj_R.support,  len(obj_del.readNames)])
                            obj_del.depth = np.max([obj_L.depth, obj_R.depth])
                            obj_del.dr = obj_del.depth - obj_del.support

                            obj_del.getGenotype_cuteSV()
                            obj_del.getGenotype_sniffles()

                            #self.tree.add(Interval(clusterSigObjIns.start_mode, clusterSigObjIns.stop_mode, clusterSigObjIns))
                            print("clusterSigObjDel3")
                            flag, SVrecord = self.svRecord(obj_del, clusterTreeL, clusterTreeR)
                            print("clusterSigObjDel4", flag, SVrecord)
                            #logger.info("clusterSigObj6")
                            if flag:
                                print("clusterSigObjIns5", flag, SVrecord)
                                self.DELs.append(SVrecord)
                            #break
        
       
    def call(self, sample, sigAllClusters, task_start, task_end):
    #def call(self, sample, sigAllClusters, task_start, task_end):
        try:
            clusterTreeD = sorted(list(sigAllClusters['D']))
            
            clusterTreeL = sigAllClusters['L']
            clusterTreeR = sigAllClusters['R']
        except KeyError:
            return

        clusterTreeD = clusterTreeD[task_start:task_end]
        print("clusterTreeD:", len(clusterTreeD))
        for obj_node in clusterTreeD:
            obj_del = obj_node.data
            # if obj_del.support * obj_del.svlen_mode < 300:
            #     continue
            # flag1 = True  #
            flag1 = self.filterIncludedDEL(obj_del, clusterTreeD, clusterTreeL, clusterTreeR)
            print("filterIncludedDEL flag1", flag1, obj_del.chrom, obj_del.start_mode, obj_del.stop_mode)
            obj_del.median_svlen = np.median(obj_del.svlens)
            flag2, SVrecord = self.svRecord(obj_del, clusterTreeL, clusterTreeR)
            print("filterIncludedDEL flag2", flag2, obj_del.chrom, obj_del.start_mode, obj_del.stop_mode)
            if flag1 & flag2:
                self.DELs.append(SVrecord)
        
        try:
            clusterTreeI = sigAllClusters['I']
            clusterTreeD = sigAllClusters['D']
        except KeyError:
            return
        self.call_RL(clusterTreeD, clusterTreeI, clusterTreeL, clusterTreeR)
    
    def filterIncludedDEL(self, obj_del, clusterTreeD, clusterTreeL, clusterTreeR):
        chrom, start, end = obj_del.chrom, obj_del.start_mode, obj_del.stop_mode
        
        if obj_del.support * obj_del.svlen_mode < self.configer.min_sv_size*0.8*obj_del.depth*0.5:#300:
                return False
        else:
            try:
                inv_tree = self.data_inv_dup["INV"][obj_del.chrom]
            except KeyError:
                inv_tree = None

            try:
                dup_tree = self.data_inv_dup["DUP"][obj_del.chrom]
            except KeyError:
                dup_tree = None
            #print("inv_tree_dup", inv_tree, dup_tree)
            if inv_tree is None and dup_tree is None:
                return True
            elif inv_tree is None:
                if len(dup_tree.overlap(start, end)) > 0:
                    return False
                else:
                    return True
            elif dup_tree is None:
                if len(inv_tree.overlap(start, end)) > 0:
                    return False
                else:
                    return True
            elif len(inv_tree.overlap(start, end)) > 0 or len(dup_tree.overlap(start, end)) > 0:
                return False
            else:
                return True
    
    def callDEL(self, sample, part_plk_files):
        #for plk_file in part_plk_files:
        for plk_file, task_start, task_end in part_plk_files:
            if os.path.exists(plk_file):
                with open(plk_file, 'rb') as plk:
                    sigAllClusters = pickle.load(plk)
                    self.call(sample, sigAllClusters, task_start, task_end)
                #self.call(sample, sigAllClusters)
    
    def svRecord(self, clusterSigObj, clusterTreeL, clusterTreeR):
        Ls_readNames = [ii for i in clusterTreeL.overlap(
            Interval(clusterSigObj.start_mode - 25, clusterSigObj.stop_mode + 25)) for ii in i.data.readNames]
        Rs_readNames = [ii for i in clusterTreeR.overlap(
            Interval(clusterSigObj.start_mode - 25, clusterSigObj.stop_mode + 25)) for ii in i.data.readNames]
        clusterSigObj.support = len(set(Ls_readNames) | set(Rs_readNames) | set(clusterSigObj.readNames))

        flag = False
        format = "GT1:GT2:DR:DV"
        af = clusterSigObj.support / clusterSigObj.depth
        reliability = "IMPRECISE"
        infos = ';'.join([reliability,
                          "SVTYPE=DEL",
                          "MEDIAN_SVLEN=" + str(clusterSigObj.median_svlen),
                          "mapQ_mean=" + str(clusterSigObj.mapQ_mean),
                          "mapQ_std=" + str(clusterSigObj.mapQ_std),
                          "SVLEN=" + str(clusterSigObj.svlen_mode),
                          "END=" + str(clusterSigObj.stop_mode),
                          "SUPPORT=" + str(clusterSigObj.support),
                          "AF=" + str(round(af, 3)),
                          "STD_POS1=" + str(round(clusterSigObj.start_std, 3)),
                          "STD_POS2=" + str(round(clusterSigObj.stop_std, 3)),
                          "STD_SVLEN=" + str(round(clusterSigObj.svlen_std, 3)),
                          "READNAMES="+'|'.join(clusterSigObj.readNames),
                          "RNAMES=NULL",
                          ])

        flag = self.filter(af, clusterSigObj)
        #flag = True
        return flag, [[clusterSigObj.chrom, clusterSigObj.start_mode, clusterSigObj.svlen_mode, clusterSigObj.support], [clusterSigObj.chrom, clusterSigObj.start_mode, "plotSV.DEL", "N", "<DEL>", '60', "PASS", infos,
                      format, ':'.join(
                [clusterSigObj.gt1, clusterSigObj.gt2, str(clusterSigObj.dr), str(clusterSigObj.support)])]]

    
    def filter(self, af,clusterSigObj):
        if af >= 0.001 and clusterSigObj.svlen_mode>=30:
            return True
        else:
            return False
        pass