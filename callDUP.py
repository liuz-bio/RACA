import os
import pickle
import time
import copy
import pandas as pd
import numpy as np
from dataclasses import dataclass
from intervaltree import IntervalTree#, Interval
from clusterSig import readSig, clusterSig
from scipy.stats import ttest_ind
from sklearn.cluster import DBSCAN
import pysam
import pyalign
import hashlib
from Bio import Align
from loguru import logger
from collections import Counter
import gc
import tracemalloc
from scipy import stats

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

IntervalTreeadd
class callDUP:
    def __init__(self, out_file, configer, sample, mut_pos, shared_sig_trees, duplicate_records, lock):
        self.DUPs = []
        self.INSs = []
        self.configer = configer
        self.out_file = out_file
        #self.all_DUP = all_DUP
        #logger.info("all_DUP: "+str(id(self.all_DUP))+"|"+str(id(all_DUP)))
        self.lock = lock
        bed = str(os.path.join(configer.tmp_path, sample, sample+'.bed.gz'))
        self.genome_index = pysam.FastaFile(configer.genome)
        if self.configer.bam_cram:
            self.bam = pysam.AlignmentFile(self.configer.samples_to_bams[sample], "rb")
        else:
            self.bam = pysam.AlignmentFile(self.configer.samples_to_bams[sample], "rb")
        self.bed_index = pysam.TabixFile(bed)
        self.tree = IntervalTree()
        self.clusterObjBnd = []
        start_t1 = time.time()
        shared_sig_trees = list(shared_sig_trees)
        self.subSigTree = shared_sig_trees[1]["DUP"]
        self.duplicate_records = duplicate_records
        #self.subSigTrees = self.subSigRegion(sample)
        print("time_start1", time.time()-start_t1)
        tracemalloc.start()
        self.sample =sample
        self.mut_pos = mut_pos    
        bedFileName = "/home/lz/work_space/project/plotSV/RACA/plotSV9/callSV/hg38.dup2.bed.gz"
        self.dupbedfile = pysam.TabixFile(bedFileName) 

    
    def run(self):
        tSV = time.time()
        #self.callDUP(self.sample, self.mut_pos)
        print('tSV1', 'DUP', time.time()-tSV)
        tSV = time.time()
        self.pass_overlap_DUPs()
        self.pass_overlap_INSs()
        self.write_vcf()
        print('tSV2', 'DUP', time.time()-tSV)
 

    def write_vcf(self):
        print("write_vcf DUPs", self.out_file)
        with self.lock:
            with open(self.out_file, 'a') as out:
                out.write('\n'.join(['\t'.join(map(str, i)) for i in self.DUPs])+'\n'+'\n'.join(['\t'.join(map(str, i)) for i in self.INSs]))
                out.flush()
            
    
    def pass_overlap_DUPs(self):
        SVs = self.DUPs
        num = len(SVs)
        tmp_SVs = [i[1] for i in SVs]
        print("pass_overlap3:", len(SVs))
        for i in range(num-1):
            a = SVs[i]
            for ii in range(i+1,num):
                b = SVs[ii]
                if a[0][0] == b[0][0]:
                    if (a[0][1] <= b[0][1] and a[0][2]>=b[0][2]) and (abs(a[0][1] - b[0][1]) < 100 and  abs(a[0][2]-b[0][2]) <100):
                        try:
                            print("pass_overlap1:",a[1], b[1] )
                            tmp_SVs.remove(b[1])
                        except ValueError:
                            continue
                    elif b[0][1] <= a[0][1] and b[0][2]>=a[0][2] and (abs(a[0][1] - b[0][1]) < 100 and  abs(a[0][2]-b[0][2]) <100):
                        try:
                            print("pass_overlap2:",a[1], b[1] )
                            tmp_SVs.remove(a[1])
                        except ValueError:
                            continue
                    elif abs(b[0][1]-a[0][1]) <100 or abs(a[0][2]-b[0][2]) <100:
                        if b[0][3] > a[0][3]:
                            try:
                                print("pass_overlap22:",a[1], b[1] )
                                tmp_SVs.remove(a[1])
                            except ValueError:
                                continue
                        else:
                            try:
                                print("pass_overlap11:",a[1], b[1] )
                                tmp_SVs.remove(b[1])
                            except ValueError:
                                continue

        self.DUPs = tmp_SVs
        # with self.lock:
        #     for i in DUPs:
        #         self.all_DUP.append(i)
        #     logger.info("self.all_DUP: "+str(len(self.all_DUP)))
        print("pass_overlap4:", len(SVs))

    
    def pass_overlap_INSs(self):
        SVs = self.INSs
        num = len(SVs)
        tmp_SVs = [i[1] for i in SVs]
        print("pass_overlap3:", len(SVs))
        for i in range(num-1):
            a = SVs[i]
            for ii in range(i+1,num):
                b = SVs[ii]
                if a[0][0] == b[0][0]:
                    if abs(b[0][1]-a[0][1]) <100 or abs(a[0][2]-b[0][2]) <100:
                        if b[0][3] > a[0][3]:
                            try:
                                print("pass_overlap22:",a[1], b[1] )
                                tmp_SVs.remove(a[1])
                            except ValueError:
                                continue
                        else:
                            try:
                                print("pass_overlap11:",a[1], b[1] )
                                tmp_SVs.remove(b[1])
                            except ValueError:
                                continue

        self.INSs = tmp_SVs
        print("pass_overlap4:", len(SVs))
    
    
    def mergeSigAllClusters(self, plk_files):
        sigAllClusters = {'L': None, 'R': None, 'D': None, 'I': None, }
        tmp = []
        #logger.info("plk_files: "+str(plk_files))
        print(plk_files)
        for plk_file in plk_files:
            #logger.info("plk_file "+plk_file)
            if os.path.exists(plk_file):
                with open(plk_file, 'rb') as plk:
                    tmp.append(pickle.load(plk))
            #logger.info("CCC")

        if len(tmp) >= 2:
            for sigT in ['L', 'R', 'D', 'I']:
                merged_tree = tmp[0][sigT]
                for trees in tmp[1:]:
                    merged_tree.update(trees[sigT])
                    #merged_tree = merged_tree.union(trees[sigT])
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
    
    # def getSigAllClusters(self, chrom, n1, n2, sample):
    #     if n1 == n2:
    #         plk_files = [str(os.path.join(self.configer.tmp_path, sample,
    #                                       chrom + ".%d.%d.plk" % (
    #                                           n * self.configer.win_size + 1,
    #                                           (n + 1) * self.configer.win_size))) for n in [n1-1, n2]]
    #         tx = time.time()
    #         sigAllClusters, flag = self.mergeSigAllClusters(plk_files)
    #         print("getSigAllClusters2:", time.time() - tx)
    #         return sigAllClusters, flag
    #         # plk_file = str(os.path.join(self.configer.tmp_path, sample,
    #         #                             chrom + ".%d.%d.plk" % (
    #         #                             n1 * self.configer.win_size + 1, (n1 + 1) * self.configer.win_size)))
    #         # print(plk_file)
    #         # if os.path.exists(plk_file):
    #         #     with open(plk_file, 'rb') as plk:
    #         #         sigAllClusters = pickle.load(plk)
    #         #         return sigAllClusters, True
    #         # else:
    #         #     return None, False

    #     else:
    #         plk_files = [str(os.path.join(self.configer.tmp_path, sample,
    #                                               chrom + ".%d.%d.plk" % (
    #                                                   n * self.configer.win_size + 1,
    #                                                   (n + 1) * self.configer.win_size))) for n in [n1, n2]]
    #         tx = time.time()
    #         sigAllClusters, flag = self.mergeSigAllClusters(plk_files)
    #         print("getSigAllClusters2:", time.time()-tx)
    #         return sigAllClusters, flag

    
    def checkSubSigx(self, chrom, start, end):#从内往外
        overlap_nodes = self.subSigTree[chrom].overlap(start, end)
        overlap_nodes_size = len(overlap_nodes)
        if overlap_nodes_size<1 or overlap_nodes_size>10:
            return None, None, True

        begins = [i.begin-start for i in overlap_nodes]
        ends = [i.end-end for i in overlap_nodes]
        begins.sort(reverse=True)
        ends.sort()

        min_dis =  min([abs(i + ii) for i in begins for ii in ends])
        if min_dis < 100:  # 对称性
            #return start, end, False
            begins_size = len(begins)
            ends_size = len(ends)
            for i in range(begins_size):
                for ii in range(ends_size):
                    if abs(begins[i]+ends[ii]) == min_dis:
                        print("checkSubSig:", chrom, start, end, begins[i]+start, ends[ii]+end, False)
                        return begins[i]+start, ends[ii]+end, False
        else:
            return None, None, True

    
    def checkSubSig(self, chrom, start, end):#从内往外
        overlap_nodes = self.subSigTree[chrom].overlap(start, end)
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


    
    def global_alignment(self, chrom, start, end, readSeq):
        genomeSeq = self.genome_index.fetch(chrom, start, end)
        print("genomeSeq",genomeSeq)
        alignment = pyalign.global_alignment(genomeSeq.upper(), readSeq.upper(), gap_cost=0, eq=1, ne=-2)
        return alignment.score/len(readSeq)

    
    def global_alignmentx(self, chrom, start, end, readSeq):
        aligner = Align.PairwiseAligner()
        # 设置比对参数
        # aligner.mode = 'global'  # 全局比对
        # aligner.match_score = 1  # 匹配得分
        # aligner.mismatch_score = -1  # 不匹配得分
        # aligner.open_gap_score = -1  # 开放缺失得分
        # aligner.extend_gap_score = -0.5  # 扩展缺失得分

        aligner.mode = 'global'  # 全局比对
        aligner.match_score = 1  # 匹配得分
        aligner.mismatch_score = 0  # 不匹配得分
        aligner.open_gap_score = 0  # 开放缺失得分
        aligner.extend_gap_score = 0  # 扩展缺失得分

        genomeSeq = self.genome_index.fetch(chrom, start, end)
        alignment = aligner.align(genomeSeq.upper(), readSeq.upper())

        return alignment.score/len(readSeq)

    
    def lengthDI(self, contig, start, stop):
        data = [[],[],[],[]] #2 pos 3 pos 2 svlen 3 svlen
        readNames = []
        
        for i in self.bed_index.fetch(contig, start, stop): #重新获取候选DUP区域的D和I信号的长度
            i = i.strip().split('\t')
            readNames.append(i[9])
            for ii in i[10:]:

                ii = ii.strip().split(',')
                sig_start = int(ii[0])
                sig_end = int(ii[1])
                sig_size = abs(sig_end-sig_start)
                start_biase = start-sig_size
                stop_biase = stop+sig_size
                if ii[2] == '3'  and  start_biase <= sig_start <= stop_biase: #I
                    data[int(ii[2])].append(abs(sig_start- sig_end))
                    data[int(ii[2])-2].append([sig_start, sig_end])
                elif ii[2] == '2' and (start_biase <= sig_start <= stop_biase or start_biase <= sig_end <= stop_biase): #D
                    data[int(ii[2])].append(abs(sig_start- sig_end))
                    data[int(ii[2])-2].append([sig_start, sig_end])
                    
        return data, len(set(readNames))

    
    def mutSigClusterDeal_old(self, objL, objR, start, end, clusterTreeD, clusterTreeI):
        chrom = objL.chrom
        posDeviation = 150
        svlen = abs(end - start)
        obj_Ds = [D_node.data for D_node in clusterTreeD.overlap(start - posDeviation, end+posDeviation)]
        obj_Is = [I_node.data for I_node in clusterTreeI.overlap(start - posDeviation, end+posDeviation)]

        #ld = [len(objL.readNames & i.readNames) for i in obj_Ds]
        #rd = [len(objR.readNames & i.readNames) for i in obj_Ds]
        #li = [len(objL.readNames & i.readNames) for i in obj_Is]
        #ri = [len(objR.readNames & i.readNames) for i in obj_Is]

        print("mutSigClusterDeal Ds",start, end, [i.start_mode for i in obj_Ds], [i.svlen_mode for i in obj_Ds], abs(start-end))
        print("mutSigClusterDeal Is",start, end, [i.start_mode for i in obj_Is], [i.svlen_mode for i in obj_Is], abs(start-end))
        D_l = len(obj_Ds)
        I_l = len(obj_Is)
        print("obj_Ds", obj_Ds, D_l)
        print("obj_Is", obj_Is, I_l, )
        #logger.info("mutSigClusterDeal1 "+'\t'.join(map(str, [D_l, I_l])))
        data = self.lengthDI(chrom, start, end)
        
        
        xI = sorted(data['3'], reverse=True)
        xD = sorted(data['2'], reverse=True)
        xII = [i for i in data['3'] if abs(i - max(data['3'])) <= 100] #计算所有I长度与最大长度差大于100的长度
        xDD = [i for i in data['2'] if abs(i-max(data['2'])) <=100]
        print("xII:", start, end, xII)
        print("xDD:", start, end, xDD)

        if D_l > 0 and I_l > 0:
            Ds_svlen = np.array([i.svlen_mode for i in obj_Ds])
            Is_svlen = np.array([i.svlen_mode for i in obj_Is])
            D_svlen = np.array([min([i.svlen_mode, svlen]) / max([i.svlen_mode, svlen]) for i in obj_Ds])
            I_svlen = np.array([min([i.svlen_mode, svlen]) / max([i.svlen_mode, svlen]) for i in obj_Is])
            print("A_D_I:", Ds_svlen, Is_svlen, D_svlen, I_svlen)

            #logger.info("mutSigClusterDeal2 " + '\t'.join(map(str, D_svlen)) + '  |   ' + '\t'.join(map(str, I_svlen)))
  
            if max(D_svlen) < 0.1 and max(I_svlen) < 0.1:
                return None, None, False
            elif max(data['2']) > svlen*1.1 or max(data['3']) >svlen*1.1:
                 return None, None, True
            elif svlen*0.25<=np.max(xI[:len(xII)])<=svlen*0.75 and len(xII)>=1:
                return  None, None, True
            elif svlen*0.25<=np.max(xD[:len(xDD)])<=svlen*0.75 and len(xDD)>=1:
                return None, None, True
            elif  max(data['3']) >1.5*svlen  or max(data['2']) >1.5*svlen:
                return None, None, True
            elif max(D_svlen) >= 0.1 and max(Ds_svlen)>=100:
                max_D_obj = obj_Ds[np.argmax(D_svlen)]
                # mutSigClusterDeal None None True 16
                print("mutSigClusterDeal: D_svlen support", [i.support for i in obj_Ds], max_D_obj.support,
                      objL.support, objR.support)
                if len(max_D_obj.readNames & objL.readNames) >= 1 or len(max_D_obj.readNames & objR.readNames) >= 1:
                    return None, None, True
                else:
                    return None, None, False
            elif max(I_svlen) >= 0.1:
                max_I_obj = obj_Is[np.argmax(I_svlen)]

                score = self.global_alignmentx(chrom, start, end, max_I_obj.alt)
                #logger.info("mutSigClusterDeal score: " + str(score))
                print("mutSigClusterDeal score:", score, len(max_I_obj.alt))
                if score >= 0.5:
                    # if min([obj_I.svlen_mode, svlen_inv])/max([obj_I.svlen_mode, svlen_inv])>=0.85:
                    return None, None, False
                else:
                    return None, None, True
            else:
                return None, None, False

        elif D_l >0:
            D_svlen = np.array([i.svlen_mode for i in obj_Ds])
            print("A_D:", D_svlen)
            max_D_obj = obj_Ds[np.argmax(D_svlen)]
            # mutSigClusterDeal: D_svlen support [18] 18 16 16
            print("mutSigClusterDeal: D_svlen support", [i.support for i in obj_Ds], max_D_obj.support, objL.support,
                  objR.support)
            if ((svlen*0.25<=np.max(xD[:len(xDD)])) or (np.max(xD[:len(xDD)]) >= svlen )) and len(xDD)>=1:
                return None, None, False
            # elif max(data['2']) >1.5*svlen:
            #     return None, None, True
            elif max(D_svlen) >= 100:
                if len(max_D_obj.readNames & objL.readNames) >= 1 or len(max_D_obj.readNames & objR.readNames) >= 1:
                    return None, None, False
                else:
                    return None, None, True
            else:
                return None, None, False

        elif I_l>0:
            I_svlen = np.array([i.svlen_mode for i in obj_Is])
            print("A_I:", I_svlen)
            if (svlen*0.25<=np.max(xI[:len(xII)]) or np.max(xI[:len(xII)]) >= svlen) and len(xII)>=1:
                return None, None, False
            # elif  max(data['3']) >1.5*svlen:
            #     return None, None, True

            for obj_I in obj_Is:
                #if abs(obj_I.svlen_mode - svlen_inv) < 0.8 * max([obj_I.svlen_mode, svlen_inv]):
                if obj_I.svlen_mode >= 100:
                    score = self.global_alignmentx(chrom, start, end, obj_I.alt)
                    print("mutSigClusterDeal score:", score, len(obj_I.alt))
                    if score >=0.5:
                        return None, None, False
                    else:
                        return None, None, True
                else:
                    return None, None, False

        else:
            return None, None, False

    
    def mutSigClusterDeal(self, objL, objR, start, end, clusterTreeD, clusterTreeI, sv_biase, seq_svlens):
        chrom = objL.chrom
        posDeviation = 150
        #svlen = abs(end - start)
        seq_svlens1 = self.remove_outliers_iqr(pd.Series(seq_svlens))
        #svlen1 = int(np.mean(seq_svlens1))
        svlen1 = int(max(seq_svlens1))
        svlen2 = abs(objL.start_mode-objR.start_mode)
        obj_Ds = [D_node.data for D_node in clusterTreeD.overlap(start - posDeviation, end+posDeviation)]
        obj_Is = [I_node.data for I_node in clusterTreeI.overlap(start - posDeviation, end+posDeviation)]
        obj_Is_svlen = [i.svlen_mode for i in obj_Is]
        obj_Ds_svlen = [i.svlen_mode for i in obj_Ds]
        
        if objL.start_mode > objR.start_mode:
            return None, None, True
        
        if len(objL.readHashs & objR.readHashs)/max([len(objL.readHashs), len(objR.readHashs)]) >0.95: #避免INV类型信号的影响， 左侧断点和右侧断点来自同一条比对记录
            return None, None, True

        flag  = len([i for i in objL.readNames & objR.readNames if abs(objL.query_starts[i][0] - objR.query_starts[i][0]) <500 ]) >0
        

        data, map_recod_nu = self.lengthDI(chrom, start, end)
        d_start = [abs(start-i[0]) for i in data[0]]
        i_start = [abs(start-i[0]) for i in data[1]]
        if map_recod_nu>500: #比对记录数量远远超过测序深度
            return None, None, True


        xII_pos = [data[1][i] for i in range(len(data[3])) if abs(data[3][i] - max(data[3])) <= 100] #计算所有I长度与最大长度差大于100的长度
        xDD_pos = [data[0][i] for i in range(len(data[2]))  if abs(data[2][i]-max(data[2])) <=100]
        xII = [i for i in data[3] if abs(i-max(data[3]))<=100]
        xDD = [i for i in data[2] if abs(i-max(data[2]))<=100]
        
        xII_l = len(xII)
        xDD_l = len(xDD)
        print("xII:", start, end, xII)
        print("xDD:", start, end, xDD)
        
        
        if flag:
            if min([objL.support, objR.support]) < 2:
                print("HGHGHG:", objL.chrom, start, end, objL.support, objR.support, objL.recordNu, objR.recordNu, objL.readNames, objR.readNames)
                return None, None, True
            if xII_l>0:
                if max(xII) >1.5*svlen2:
                    return None, None, True
                if xDD_l>0:
                    if max(xDD) >svlen2*0.5:
                        return None, None, True
            if xDD_l>0:
                if max(xDD) >0.5*svlen2:
                    return None, None, True
            return None, None, False
        else: 
            if max([objL.support, objR.support])<=2*1.5: #
                print("HGHGHG0:", objL.chrom, start, end, objL.start_mode, objR.start_mode, objL.support, objR.support, objL.recordNu, objR.recordNu, objL.readNames, objR.readNames)
                return None, None, True
            elif xII_l>0:
                if len(obj_Is) >0 and len(obj_Ds) ==0:
                    if max(xII) >1.5*svlen2:
                        print("HGHGHG1:", objL.chrom, start, end, objL.start_mode, objR.start_mode, objL.support, objR.support, objL.recordNu, objR.recordNu, objL.readNames, objR.readNames)
                        return None, None, True
                    elif len(obj_Is) >=2:
                        try:#DUP内部由于重复导致错位匹配，产生多种INS信号
                            tmp_svlen, new_start, new_end = sorted([[abs(i.start_mode-j.start_mode), i.start_mode, j.start_mode]  for i in obj_Is[:-1] for j in obj_Is[1:] if (abs(i.start_mode-j.start_mode)  >abs(objL.start_mode-objR.start_mode)*0.2 and min([i.support,j.support])>3)], key=lambda x: x[0])[0]
                            return min([new_start, new_end]), max([new_start, new_end]), False
                        except IndexError:
                            if max(xII) > svlen1:
                                print("HGHGHG2:", objL.chrom, start, end, objL.start_mode, objR.start_mode, objL.support, objR.support, objL.recordNu, objR.recordNu, objL.readNames, objR.readNames)
                                return  None, None, True
                    elif abs(obj_Is[0].svlen_mode -svlen1) <100: #INS长度和DUP信号长度一致
                        return min([obj_Is[0].start_mode, start]), max([obj_Is[0].start_mode, end]), False
                    elif max(xII) >1.5*svlen2 :
                        print("HGHGHG3:", objL.chrom, start, end, objL.start_mode, objR.start_mode, objL.support, objR.support, objL.recordNu, objR.recordNu, objL.readNames, objR.readNames)
                        return None, None, True
                    elif 0.2*svlen2<=max(xII) <0.8*svlen2 :
                        return None, None, True
                elif len(obj_Is) >0 and xDD_l >0:#通过I D的长度差异区分INS和DUP, DUP的I,D信号通常较小
                    if max(xII) >svlen2 or max(xDD) >svlen2:
                        print("HGHGHG4:", objL.chrom, start, end, objL.start_mode, objR.start_mode, objL.support, objR.support, objL.recordNu, objR.recordNu, objL.readNames, objR.readNames)
                        return None, None, True
                    obj_I_Ds = obj_Is+obj_Ds
                    try:#DUP内部由于重复导致错位匹配，产生多种INS信号
                        tmp_svlen, new_start, new_end = sorted([[abs(i.start_mode-j.start_mode), i.start_mode, j.start_mode]  for i in obj_I_Ds[:-1] for j in obj_I_Ds[1:] if (abs(i.start_mode-j.start_mode)  >abs(objL.start_mode-objR.start_mode)*0.3 and min([i.support,j.support])>3)], key=lambda x: x[0])[0]
                        # if abs(new_start-new_end) < max([i.svlen_mode for i in obj_I_Ds]): #INS
                        #     return None, None, True
                        
                        return min([new_start, new_end]), max([new_start, new_end]), False
                    except IndexError:
                        if max(xII) > svlen1 or max(xDD) >svlen1:
                            print("HGHGHG5:", objL.chrom, start, end, objL.start_mode, objR.start_mode, objL.support, objR.support, objL.recordNu, objR.recordNu, objL.readNames, objR.readNames)
                            return  None, None, True
                        # elif 0.2*svlen2<=max(xII) <0.8*svlen2 :
                        #     print("HGHGHG6:", objL.chrom, start, end, objL.start_mode, objR.start_mode, objL.support, objR.support, objL.recordNu, objR.recordNu, objL.readNames, objR.readNames)
                        #     return None, None, True
                            
                    print("XQAZSW:"+str('|'.join(map(str, [i.start_mode for i in obj_Is])))+"MMM"+str('|'.join(map(str, [i.start_mode for i in obj_Ds]))), objL.chrom, objL.start_mode, objR.start_mode, sv_biase, seq_svlens, [i.svlen_mode for i in obj_Is], [i.svlen_mode for i in obj_Ds])
                    
                elif xDD_l>0:
                    if max(xDD) >svlen2 or max(xDD)>0.25*svlen2:
                        print("HGHGHG7:", objL.chrom, start, end, objL.start_mode, objR.start_mode, objL.support, objR.support, objL.recordNu, objR.recordNu, objL.readNames, objR.readNames)
                        return None, None, True
                elif max(xII) >svlen2 :
                    print("HGHGHG8:", objL.chrom, start, end, objL.start_mode, objR.start_mode, objL.support, objR.support, objL.recordNu, objR.recordNu, objL.readNames, objR.readNames)
                    return None, None, True
                elif max(xII) >svlen1:
                    print("HGHGHG9:", objL.chrom, start, end, objL.start_mode, objR.start_mode, objL.support, objR.support, objL.recordNu, objR.recordNu, objL.readNames, objR.readNames)
                    return None, None, True
             
            elif xDD_l >0 and len(obj_Ds)>0:
                if max(xDD) >svlen2 or max(xDD)>0.5*svlen2:
                    print("HGHGHG10:", objL.chrom, start, end, objL.start_mode, objR.start_mode, objL.support, objR.support, objL.recordNu, objR.recordNu, objL.readNames, objR.readNames)
                    return None, None, True
            elif min([objL.support, objR.support])<=2*1.5: #左右两侧信号虽然有共同的read 但read上的断点间距过大，同时如果信号较弱则误报的概率很大
                print("HGHGHG11:", objL.chrom, start, end, objL.start_mode, objR.start_mode, objL.support, objR.support, objL.recordNu, objR.recordNu, objL.readNames, objR.readNames)
                return None, None, True
                    
            # elif len(obj_Ds) >0:
                
            #     return None, None, True
            
            return None, None, False
    
    
    def unilateral(self, chrom, start, end, obj, clusterTreeD, clusterTreeI):
        k = 17
        seq = self.genome_index.fetch(chrom, start, end).upper()
        kmers = Counter([seq[i * k:(i + 1) * k] for i in range(len(seq) - k + 1)])

        ratio = (max(kmers.values())-1+k) / len(seq)
        print("kmers.values():", max(kmers.values()), ratio)
        if ratio > 0.80:
            clusterSigObjDup = copy.deepcopy(obj)
            print("unilateral: ratio, kmer-1+k, seqLen", ratio, max(kmers.values())-1+k, len(seq),"chrom|start|end", chrom, start, end, "recordNu|support|depth", obj.recordNu, obj.support, obj.depth)
            readNameRL = obj.readNames
            bx, flag = self.duplicate_records_deal(readNameRL)

            if flag:
                print("Done duplicate_records_deal")
                return


            clusterSigObjDup.bx = bx
            clusterSigObjDup.SVTYPE = "DUP"
            clusterSigObjDup.readNoHashRLNus = None
            clusterSigObjDup.readNamesRLNu = None  #
            clusterSigObjDup.readHashsRLNu = None
            clusterSigObjDup.Rs = None
            clusterSigObjDup.Ls = None
            clusterSigObjDup.Rs_start_mean = None
            clusterSigObjDup.Ls_start_mean = None
            clusterSigObjDup.Rs_Ls_start_mean = None
            clusterSigObjDup.Ls_Rs_start_mean = None
            pos_LR = (obj.PosLocation)

            clusterSigObjDup.svlen_mode = abs(end-start)
            clusterSigObjDup.reverseReadsNus = [obj.reverseReadsNu]
            clusterSigObjDup.secondaryReadsNus = [obj.secondaryReadsNu]
            mapQs = obj.mapQs
            clusterSigObjDup.mapQ_mean = np.mean(mapQs)
            clusterSigObjDup.mapQ_std = np.std(mapQs)
            clusterSigObjDup.start_mode = start
            clusterSigObjDup.stop_mode = end
            clusterSigObjDup.start_std = np.std(pos_LR)
            clusterSigObjDup.stop_std = np.std(pos_LR)
            clusterSigObjDup.svlen_std = 0
            clusterSigObjDup.dr = clusterSigObjDup.depth - clusterSigObjDup.support
            # print("start2:", time.time() - start)
            start_t = time.time()
            # print("DV,DP:",clusterSigObjDup.support, clusterSigObjDup.depth)
            clusterSigObjDup.getGenotype_cuteSV()
            clusterSigObjDup.getGenotype_sniffles()
            # print("start3:", time.time() - start)
            # logger.info("clusterSigObj5")
            print('start end',  start, end)
            self.tree.insert(start, end, clusterSigObjDup)
            #self.tree.addi(start, end, clusterSigObjDup)
            print("clusterSigObjInv3")
            flag, SVrecord = self.svRecord(clusterSigObjDup)
            print("clusterSigObjInv4", flag, SVrecord)
            # logger.info("clusterSigObj6")
            if flag:
                self.DUPs.append(SVrecord)
    
    def call(self, contig, start, end, sigAllClusters):
        posDeviation = 150
        posX = 200
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
            start_time = time.time()
            self.deal_R_L(contig, start, end, clusterTreeL, clusterTreeR, clusterTreeD, clusterTreeI, clusterL_nodes, clusterR_nodes, posDeviation)
            print("QQQ_deal_R_L1", time.time()-start_time)
        elif objR_nodes_size>0: 
            start_time = time.time()
            clusterL_nodes = self.sing_sig_L_R(contig, start, clusterTreeL, "L")
            print("QQQ1_deal_R_L", time.time()-start_time)
            if clusterL_nodes != None:
                start_time = time.time()
                self.deal_R_L(contig, start, end, clusterTreeL, clusterTreeR, clusterTreeD, clusterTreeI, clusterL_nodes, clusterR_nodes, posDeviation)
                print("QQQ_deal_R_L1", time.time()-start_time)
            print("clusterR_nodes=0", clusterR_nodes, clusterL_nodes)
        elif objL_nodes_size >0:
            start_time = time.time()
            clusterR_nodes = self.sing_sig_L_R(contig, end, clusterTreeR, "R")
            print("QQQ2_deal_R_L", time.time()-start_time)
            if clusterR_nodes != None:
                start_time = time.time()
                self.deal_R_L(contig, start, end, clusterTreeL, clusterTreeR, clusterTreeD, clusterTreeI, clusterL_nodes, clusterR_nodes, posDeviation)
                print("QQQ_deal_R_L2", time.time()-start_time)
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
                if overlapReadNameLen <1:
                    continue

                print("clusterSigObj", clusterSigObjL.stop_mode, clusterSigObjR.start_mode, [ i.start for i in clusterSigObjL.cluster],[ i.stop for i in clusterSigObjR.cluster])
                print("clusterSigObj1", clusterSigObjL.stop_mode, clusterSigObjR.start_mode, posDeviation, overlapReadNameLen, self.configer.min_sv_size)
            
                size_L = sum([ 1 for i in clusterSigObjR.cluster if start - i.sLR_other >10])
                size_R = sum([ 1 for i in clusterSigObjL.cluster if i.sLR_other - end >10])
                b = sum([1 for i in clusterSigObjL.cluster if i.query_length > 1.5*abs(start-end)])+sum([1 for i in clusterSigObjR.cluster if i.query_length > 1.5*abs(start-end)])                
                if b>=4:
                    if size_L >0 or size_R>0:
                        if idx in a:
                            a.remove(idx)
                            print("ERWERW:",contig, start, end, a, clusterSigObjL.start_mode, clusterSigObjR.start_mode, [i.sLR_other for i in clusterSigObjL.cluster ], [i.sLR_other for i in clusterSigObjR.cluster ])
                        self.mergeClusterSigObjGroupLS(contig, clusterSigObjL,
                                                        clusterSigObjR, clusterTreeL, clusterTreeR,
                                                        overlapReadNameLen, clusterTreeD, clusterTreeI, start, end)
                else:
                    self.mergeClusterSigObjGroupLS(contig, clusterSigObjL,
                                                        clusterSigObjR, clusterTreeL, clusterTreeR,
                                                        overlapReadNameLen, clusterTreeD, clusterTreeI, start, end)
                    print("ERWERW:",contig, start, end, a, clusterSigObjL.start_mode, clusterSigObjR.start_mode, [i.sLR_other for i in clusterSigObjL.cluster ], [i.sLR_other for i in clusterSigObjR.cluster ])
    
    def sing_sig_L_R(self, contig, pos, clusterTreeLR, sigT):#补充获取无法聚类的R或L
        #bedFileName = "/home/liuz/work_space/project/liuz/plotSV/plotSV8/callSV/tmp/tmp/SIM/SIM.bed.gz"
        bedFileName = os.path.join(self.configer.tmp_path, self.sample, self.sample + ".bed.gz")
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
    
    def callDUP(self, sample,  mut_pos):
        print(mut_pos)
        mut_pos = iter(mut_pos)
        print(mut_pos)
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
            if head_n1 == n1 and head_n2 == n2:
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
    
    def duplicate_records_deal(self, readNames):
        
        record_times = []
        for i in readNames:
            record_times.append(self.duplicate_records[i])

        a = [i+1 if i==0 else i for i in record_times]
        b = Counter(a)

        print("duplicate_records_deal:", readNames, a, record_times)
        #if (b[2]*2+b[3]*3)/all >=0.8:
        return b, False
    
    def dupRepeat_old(self, contig, start, end, seq_lens):
        seq_svlen = int(np.mean(seq_lens))
        dist = abs(start-end)
        max_dup_len = 0
        for record in self.dupbedfile.fetch(contig, start, end):
            line = record.strip().split("\t")
            dup_start = int(line[1])
            dup_end = int(line[2])
            dup_len = int(line[3])
            print("dupRepeatx:", contig, start, end, line, int(line[1]), int(line[2]))
            if dup_len >max_dup_len:
                max_dup_len = dup_len
            if (dup_end-dup_start) >= seq_svlen >= (dup_end-dup_start-2*dup_len):
                if (dup_start <= start <= dup_start+dup_len-5):
                    if  abs(dup_len-seq_svlen) <50:
                        new_start = start
                        new_end = end - dist
                        return new_start, new_end, "DUP"
                    else:
                        new_start = start
                        new_end = start+1
                        return new_start, new_end, "INS"
                elif (dup_end-dup_len+5 <= end <= dup_end):
                    if  abs(dup_len-seq_svlen) <50:
                        new_start = start+dist
                        new_end = end
                        return new_start, new_end, "DUP"
                    else:
                        new_start = end-1
                        new_end = end
                        return new_start, new_end, "INS"
            
        if abs(max_dup_len-seq_svlen) <50:
            new_start = start+int(seq_svlen/2)
            new_end = end-int(seq_svlen/2)
            return new_start, new_end, "DUP"
        else:
            new_start = int((start+end)/2)
            new_end = new_start+1
            return new_start, new_end, "INS"
    
    def get_genome_seq(self, contig, start, end):
        sequence = self.genome_index.fetch(reference=contig, start=start, end=end)
        return sequence.upper()
       
    def get_read_seq(self, read_seq, start, end):
        if start <0:
            start = 0 
        if end > len(read_seq):
            end = len(read_seq)
        return read_seq[start:end].upper()
    
    def calculate_seq_similarity(self, seq1, seq2, seq3):
        alig1 = pyalign.global_alignment(seq1, seq2, gap_cost=1.5, eq=1, ne=-1)
        alig2 = pyalign.global_alignment(seq1, seq3, gap_cost=1.5, eq=1, ne=-1)
        alig3 = pyalign.global_alignment(seq2, seq3, gap_cost=1.5, eq=1, ne=-1)
        return alig1.score*alig2.score*alig3.score/max((len(seq1),len(seq2),len(seq3)))**3
    
    def dupRepeat(self, contig, start, end, svlen_biase, svlen_raw, readID_Pos): 
        ratio_y_dict = []
        readIDs = readID_Pos.keys()
        for y_ratio in [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]:
            x_ratio = 1 - y_ratio
            L, L1 = start-int(y_ratio*svlen_biase), start
            R1, R = end, end+int(x_ratio*svlen_biase)
            genome_seq_R = self.get_genome_seq(contig, R1, R)
            genome_seq_L = self.get_genome_seq(contig, L, L1)
            if len(genome_seq_R) == 0 or len(genome_seq_L)==0:
                return start, end
            #RL_ratios = []
            reads = self.bam.fetch(contig, start, end)
            for read in reads:
                read_id = read.query_name.split("|")[0] 
                if read_id in readID_Pos:
                    l, r = readID_Pos[read_id]
                    r1,r2 = r-int(x_ratio*svlen_biase), r
                    r3, r4 = l+svlen_biase, l+int(x_ratio*svlen_biase)+svlen_raw
                    l1, l2 = l, l+int(y_ratio*svlen_biase)+1
                    l3, l4 = r-int(y_ratio*svlen_biase)-svlen_raw-1, r-svlen_biase
                    read_seq_R1 = self.get_read_seq(read.query_sequence, r1,r2)
                    read_seq_R2 = self.get_read_seq(read.query_sequence, r3, r4)
                    read_seq_L1 = self.get_read_seq(read.query_sequence, l1, l2)
                    read_seq_L2 = self.get_read_seq(read.query_sequence, l3, l4)
                    print("QSXtttB:", len(genome_seq_R), "|", len(read_seq_R1), "|", len(read_seq_R2), "||", len(genome_seq_L), "|",len(read_seq_L1), "|",len(read_seq_L2))
                    R_ratio = self.calculate_seq_similarity(genome_seq_R, read_seq_R1, read_seq_R2)
                    L_ratio = self.calculate_seq_similarity(genome_seq_L, read_seq_L1, read_seq_L2)
                    print("QSXtttA:",contig, start, end, y_ratio, L_ratio+ R_ratio, L_ratio*R_ratio, svlen_biase, L_ratio, R_ratio,len(genome_seq_L), len(read_seq_L1), len(read_seq_L2),len(genome_seq_R), len(read_seq_R1), len(read_seq_R2),"genome_L", L, L1, "ReadL1", l1, l2, "ReadL2", l3, l4, "ReadLen",len(read.query_sequence))
                    #RL_ratios.append(R_ratio+L_ratio)
                    ratio_y_dict.append([R_ratio+L_ratio,y_ratio])
            #ratio_y_dict[y_ratio] = RL_ratios
        sorted_ratio_y_dict = sorted(ratio_y_dict, key=lambda x: x[0])
        y_ratio = -1
        if len(sorted_ratio_y_dict) >=5:
            y_ratio = Counter([row[1] for row in sorted_ratio_y_dict[:5]]).most_common()[0][0]

        elif len(sorted_ratio_y_dict)>0:
            y_ratio = Counter([row[1] for row in sorted_ratio_y_dict]).most_common()[0][0]

        if y_ratio >=0:
            new_start, new_end = start-int(y_ratio*svlen_biase), end+int((1-y_ratio)*svlen_biase)
            return new_start, new_end
        else:
            return start, end
           
    
    def max_subarray(self, arr):
        max_sum = float('-inf')
        current_sum = 0
        start = 0
        max_start = 0
        max_length = 0

        for i in range(len(arr)):
            current_sum += arr[i]

            if current_sum > max_sum:
                max_sum = current_sum
                max_start = start
                max_length = i - start + 1

            if current_sum < 0:
                current_sum = 0
                start = i + 1

        return max_start, max_length
    
    def fetch_reads_from_region(self, contig, start, end, readNames, seq_lens, svlen_biase, query_LR):
        print("fetch_reads_from_region:", contig, start, end, readNames, seq_lens, svlen_biase, query_LR)
        # 打开BAM文件
        bam_file = "/home/liuz/work_space/project/liuz/plotSV/Sim.DEL.INS.DUP.INV.25x.sort.bam"
        if self.configer.bam_cram:
            handle = pysam.AlignmentFile(bam_file, "rb", require_index=True)
        else:
            handle = pysam.AlignmentFile(bam_file, "rc", require_index=True)
        # 提取比对到指定区域的reads
        reads = handle.fetch(contig, start, end)
        for read in reads:
                query_name = read.query_name.split('|')[0]
                if query_name in readNames:
                    query_seq = read.query_sequence
                    query_start = read.query_alignment_start
                    query_end = read.query_alignment_end
                    print("query_LR1", query_LR, len(query_seq))
                    read_seq1 = query_seq[query_LR[1]:].upper()
                    read_seq2 = query_seq[:query_LR[0]].upper()

                    ends = []
                    for seq_len in seq_lens:
                        genome_seq1 = self.genome_index.fetch(contig, read.reference_end, read.reference_end+seq_len).upper()
                        align = pyalign.semiglobal_alignment(read_seq1.upper(), genome_seq1.upper())
                        dat1 = list(self.max_subarray(align.s_to_t))

                        genome_seq2 = self.genome_index.fetch(contig, read.reference_start-seq_len, read.reference_start).upper()
                        align = pyalign.semiglobal_alignment(read_seq2.upper(), genome_seq2.upper())
                        dat2 = list(self.max_subarray(align.s_to_t))

                        ends.append(dat1[0]/(dat1[0]+dat2[0]))

                    #logger.info("ends: "+ str(ends))
                    
                    print("position: mean std cv", np.mean(ends), np.std(ends), np.std(ends)/np.mean(ends))
                    ends.remove(min(ends))
                    ends.remove(max(ends))
                    score_end = sum(ends)/len(ends)
                    print("position: score_end", score_end)
                    break

      
        # 关闭BAM文件
        handle.close()

        minSimScore = 0.6

        if score_end >minSimScore:
            return start, end-svlen_biase
        elif score_end <(1-minSimScore):
            return start+svlen_biase, end
        else:
            dist = int(svlen_biase/2)
            return start + dist, end - dist 

    
    def reverse_complement(self, dna_sequence):
        complement_map = str.maketrans('ATCGN', 'TAGCN')
        return dna_sequence.translate(complement_map)[::-1]
    
    def generate_kmers(self, sequence, k):
        # 直接使用生成器表达式生成kmer，减少内存占用
        return {sequence[i:i+k] for i in range(len(sequence) - k + 1)}

    def jaccard_similarity(self, set1, set2):
        intersection_size = len(set1 & set2)
        union_size = len(set1 | set2)
        return intersection_size / union_size if union_size != 0 else 0
    
    def calculate_similarity(self, seq1, seq2, k):
        kmer_set1 = self.generate_kmers(seq1, k)
        kmer_set2 = self.generate_kmers(seq2, k)
        return self.jaccard_similarity(kmer_set1, kmer_set2)
    
    def fetch_reads_from_region(self, contig, start, end, readNames, svlen_biase, query_LR, flag):
    # 打开BAM文件
        bam_file = "/home/liuz/work_space/project/liuz/plotSV/Sim.DEL.INS.DUP.INV.25x.sort.bam"
        if self.configer.bam_cram:
            handle = pysam.AlignmentFile(bam_file, "rb", require_index=True)
        else:
            handle = pysam.AlignmentFile(bam_file, "rc", require_index=True)
        reads = handle.fetch(contig, start, end)
        for read in reads:
            query_name = read.query_name.split('|')[0]
            if query_name in readNames:

                query_seq = read.query_sequence
                
                dat = []
                for end_ratio in [i/20 for i in range(1,20)]:
                    x = int(end_ratio*svlen_biase)
                    y = int((1-end_ratio)*svlen_biase)
                    dat1 = 0
                    dat2 = 0
                    genome_seq1 = self.genome_index.fetch(contig, read.reference_end-x, read.reference_end).upper()
                    genome_seq2 = self.genome_index.fetch(contig, read.reference_start, read.reference_start+y).upper()
                    svlen  = abs(end-start)-abs(query_LR[0]-query_LR[1])
                    if svlen <= 0:
                        print("svlen <0 :", contig, start, end, readNames, svlen_biase, query_LR, flag)
                        return start, end
                    e_f = query_LR[0]-query_LR[1]+2*svlen+y-x
                    z1 = query_LR[1]+e_f-x-svlen
                    z2 = query_LR[0]-e_f-y-svlen
                    if flag:
                        read_seq1 = self.reverse_complement(query_seq[z1:z1+x].upper())
                        read_seq2 = self.reverse_complement(query_seq[z2-y:z2].upper())
                        #read_seq1 = self.reverse_complement(query_seq[query_LR[1]-x:query_LR[1]].upper())
                        #read_seq2 = self.reverse_complement(query_seq[query_LR[0]:query_LR[0]+y].upper())
                    else:
                        #read_seq1 = query_seq[query_LR[1]-x:query_LR[1]].upper()
                        #read_seq2 = query_seq[query_LR[0]:query_LR[0]+y].upper()
                        read_seq1 = query_seq[z1:z1+x].upper()
                        read_seq2 = query_seq[z2-y:z2].upper()
                        
                    if len(read_seq1) >0 and len(genome_seq1) >0:
                        tmp = []
                        for k in range(5,37,2):
                            x = self.calculate_similarity(read_seq1, genome_seq1, k)
                            if x>0:
                                tmp.append(x)
                        if len(tmp) >0 or sum(tmp)>0:
                            dat1 =  np.mean(tmp)*np.mean(tmp)/(np.std(tmp)+1)
                    
                    if len(read_seq2) >0 and len(genome_seq2) >0:
                        tmp = []
                        for k in range(5,37,2):
                            x = self.calculate_similarity(read_seq2, genome_seq2, k)
                            if x >0:
                                tmp.append(x)
                        if len(tmp) >0 or sum(tmp)>0:
                            dat2 =   np.mean(tmp)*np.mean(tmp)/(np.std(tmp)+1)
                    
                    if dat1>0 and dat2 >0:
                        dat.append([abs(dat1-dat2),end_ratio])
                
                if len(dat) <= 0:
                    return start, end
                else:
                    print("out:", min(dat, key=lambda x: x[0]), flag)

                    endRatio = min(dat, key=lambda x: x[0])[1]
                    a = start+int(svlen_biase*(1-endRatio))
                    b = end-int(svlen_biase*endRatio)
                    return min([a,b]), max([a,b])
    
    def remove_outliers_iqr(self, data):
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR
        return data[(data >= lower_bound) & (data <= upper_bound)].to_list()
    
    def mergeClusterSigObjGroupLS(self, contig, clusterSigObjL, clusterSigObjR, sigClusterLs,
                                  sigClusterRs, overlapReadNameLen, clusterTreeD, clusterTreeI, start_x, end_x):

        pos_L = list(clusterSigObjL.PosLocation)
        pos_R = list(clusterSigObjR.PosLocation)
        pos_Ls = pd.Series(pos_L).value_counts()
        pos_Rs = pd.Series(pos_R).value_counts()

        clusterSigObjDup = copy.deepcopy(clusterSigObjR)
        clusterSigObjDup.raw_signal = '|'.join(map(str,[contig, start_x, end_x]))

        readNameRL = clusterSigObjL.readNames | clusterSigObjR.readNames
        bx, flag = self.duplicate_records_deal(readNameRL)

        if flag:
            print("Done duplicate_records_deal")
            return

        clusterSigObjDup.bx = bx

        overlapReadHashs = len(clusterSigObjL.readHashs & clusterSigObjR.readHashs)
        overlapReadHashName = [len(clusterSigObjL.readHashs | clusterSigObjR.readHashs), len(clusterSigObjL.readNames | clusterSigObjR.readNames)]
        print("overlapReadHashName:", start_x, end_x, clusterSigObjL.chrom, clusterSigObjL.start_mode, clusterSigObjR.chrom, clusterSigObjR.start_mode, clusterSigObjL.readHashs, clusterSigObjR.readHashs, clusterSigObjL.readHashs & clusterSigObjR.readHashs)
        print("verlapReadNameXX:",clusterSigObjL.readNames, clusterSigObjR.readNames,  clusterSigObjL.readNames & clusterSigObjR.readNames)
        clusterSigObjDup.readNamesRLNu = overlapReadNameLen  #
        clusterSigObjDup.readHashsRLNu = overlapReadHashs
        clusterSigObjDup.readNoHashRLNus = [clusterSigObjL.support, overlapReadHashs, clusterSigObjR.support,
                                            overlapReadNameLen, overlapReadHashName]


        svlen = pos_Rs.idxmax() - pos_Ls.idxmax()

        if svlen < 30 or svlen >300000:
            print("out1")
            self.clusterObjBnd.append([clusterSigObjL, clusterSigObjR])
            # print("clusterSigObjInv.svlen_mode:", pos_Rs.idxmax(), pos_Ls.idxmax())
            return

        svlens = [np.abs([pos_L[i] - pos_R[i]]) for i in range(np.min([len(pos_L), len(pos_R)]))]
        clusterSigObjDup.svlen_std = np.std(svlens)
        if clusterSigObjDup.svlen_std > 50:
            rangeSize = int(clusterSigObjDup.svlen_std)
        else:
            rangeSize = 50

        Rs = [i.data.support for i in sigClusterRs.overlap(Interval(np.min(pos_L) - rangeSize, np.min(pos_R) + rangeSize))
              if i.data.fingerPrint != clusterSigObjR.fingerPrint]
        Ls = [i.data.support for i in sigClusterLs.overlap(Interval(np.min(pos_L) - rangeSize, np.min(pos_R) + rangeSize))
              if i.data.fingerPrint != clusterSigObjL.fingerPrint]

        # if len(Rs) == 0 or len(Ls) == 0:
        #    return False, None

        

        ##############################################################################
        start = pos_Ls.idxmax()
        end = pos_Rs.idxmax()
        #logger.info("clusterSigObj2")
        # sub_start, sub_end, flag1 = self.checkSubSig(clusterSigObjDup.chrom, start, end)
        # #logger.info("clusterSigObj3")
        # if not flag1:
        #     if sub_start <= start and end <= sub_end:  # sub_start start end sub_end 扩展边界，需要进行对称性检查
        #         # 内部存在对称性，则从新定义边界
        #         new_start, new_end, flag2 = self.mutSigClusterDeal(clusterSigObjL, clusterSigObjR, sub_start, sub_end, clusterTreeD, clusterTreeI)
        #         if new_start == None and not flag2:
        #             new_start, new_end = sub_start, sub_end
        #     else:
        #         new_start, new_end, flag2 = self.mutSigClusterDeal(clusterSigObjL, clusterSigObjR, start, end, clusterTreeD, clusterTreeI)
        # else:
        readHashs_ratio = len(clusterSigObjL.readHashs & clusterSigObjR.readHashs)/max([len(clusterSigObjL.readHashs),  len(clusterSigObjR.readHashs)])
        if clusterSigObjR.start_mode >= clusterSigObjL.start_mode:
            #readNames = [i for i in clusterSigObjR.readNames & clusterSigObjL.readNames if abs(clusterSigObjL.query_starts[i][0] - clusterSigObjR.query_starts[i][0]) <500 ]
            readNames = clusterSigObjR.readNames & clusterSigObjL.readNames
        else:
            print('xxx1')
            return
        new_start, new_end = None, None
        if len(readNames) >0 :
            #seq_svlens = [abs(abs(clusterSigObjR.start_mode - clusterSigObjL.start_mode) - abs(clusterSigObjL.query_starts[i][0] - clusterSigObjR.query_starts[i][0]))  for i in readNames ]
            #seq_svlens = [abs(abs(clusterSigObjR.start_mode - clusterSigObjL.start_mode) - abs(clusterSigObjL.query_starts[i][0] - clusterSigObjR.query_starts[i][0])) if clusterSigObjL.query_starts[i][0] <= clusterSigObjR.query_starts[i][0] else abs(abs(clusterSigObjR.start_mode - clusterSigObjL.start_mode) + abs(clusterSigObjL.query_starts[i][0] - clusterSigObjR.query_starts[i][0])) for i in readNames ]
            seq_svlens = [abs(abs(clusterSigObjL.query_starts[i][0] - clusterSigObjR.query_starts[i][0])-abs(clusterSigObjR.start_mode - clusterSigObjL.start_mode)) if clusterSigObjL.query_starts[i][0] <= clusterSigObjR.query_starts[i][0] else abs(abs(clusterSigObjR.start_mode - clusterSigObjL.start_mode) + abs(clusterSigObjL.query_starts[i][0] - clusterSigObjR.query_starts[i][0])) for i in readNames ]
            sv_biase = [abs(clusterSigObjL.query_starts[i][0] - clusterSigObjR.query_starts[i][0]) for i in readNames]
            new_start, new_end, flag2 = self.mutSigClusterDeal(clusterSigObjL, clusterSigObjR, start, end, clusterTreeD, clusterTreeI, sv_biase, seq_svlens)
            print("new_start, new_end, flag2:",contig, start, end,  new_start, new_end, flag2)
        else:
            print('xxx2')
            return 
        
        if flag2:
            print("out2")
            return

        print("XSWQEDA1:", new_start, new_end, flag2)

        print("RRRR:", clusterSigObjR.chrom, start, end, [[clusterSigObjL.query_starts[i][0] , clusterSigObjR.query_starts[i][0] , clusterSigObjL.start_mode,  clusterSigObjR.start_mode] for i in clusterSigObjL.query_starts.keys() & clusterSigObjR.query_starts.keys()])
        readNames = [i for i in clusterSigObjL.query_starts.keys() & clusterSigObjR.query_starts.keys() if abs(clusterSigObjL.query_starts[i][0] - clusterSigObjR.query_starts[i][0]) <2000 and start < end]
        #readNames = [i for i in clusterSigObjL.query_starts.keys() & clusterSigObjR.query_starts.keys() if (clusterSigObjL.query_starts[i][0] - clusterSigObjR.query_starts[i][0]) > 0]
        print("XSWQEDA:",[ [clusterSigObjL.query_starts[i][0] , clusterSigObjR.query_starts[i][0] , clusterSigObjR.start_mode,  clusterSigObjL.start_mode] for i in readNames])
        #seq_svlens = [abs(abs(clusterSigObjR.start_mode - clusterSigObjL.start_mode) - abs(clusterSigObjL.query_starts[i][0] - clusterSigObjR.query_starts[i][0])) for i in readNames ]
        #seq_svlens = [abs(abs(clusterSigObjR.start_mode - clusterSigObjL.start_mode) - abs(clusterSigObjL.query_starts[i][0] - clusterSigObjR.query_starts[i][0])) if clusterSigObjL.query_starts[i][0] <= clusterSigObjR.query_starts[i][0] else abs(abs(clusterSigObjR.start_mode - clusterSigObjL.start_mode) + abs(clusterSigObjL.query_starts[i][0] - clusterSigObjR.query_starts[i][0])) for i in readNames ]
        seq_svlens = [abs(abs(clusterSigObjL.query_starts[i][0] - clusterSigObjR.query_starts[i][0])-abs(clusterSigObjR.start_mode - clusterSigObjL.start_mode)) if clusterSigObjL.query_starts[i][0] <= clusterSigObjR.query_starts[i][0] else abs(abs(clusterSigObjR.start_mode - clusterSigObjL.start_mode) + abs(clusterSigObjL.query_starts[i][0] - clusterSigObjR.query_starts[i][0])) for i in readNames ]
        print("RRRRLLLL:", contig, start, end, seq_svlens,readNames )
        if len(seq_svlens) >0:
            sv_biase = abs(clusterSigObjR.start_mode - clusterSigObjL.start_mode)-int(np.mean(seq_svlens))
            print("XSWQEDA seq_svlens:", contig, start, end, seq_svlens, sv_biase)
        else:
            print("XSWQEDA seq_svlens1:", contig, start, end, seq_svlens, sv_biase)
            return
        
        print("XSWQEDA query dist:",contig, start, end, [ [abs(clusterSigObjL.query_starts[i][0]-clusterSigObjR.query_starts[i][0])] for i in clusterSigObjL.query_starts.keys() & clusterSigObjR.query_starts.keys()] )
        
       

        #readNames = [i for i in clusterSigObjL.query_starts.keys() & clusterSigObjR.query_starts.keys() if (clusterSigObjL.query_starts[i][0] <  clusterSigObjR.query_starts[i][0] and abs(clusterSigObjL.query_starts[i][0] - clusterSigObjR.query_starts[i][0]) <100) ]
        seq_lens = [300, 400, 500, 600]
        
        print("XSWQEDA readName:", readNames, len(clusterSigObjR.readNames & clusterSigObjL.readNames))
        if len(seq_svlens) >0:
            if len(readNames)>3:
                readNames = readNames[:3]
            readID_Pos = {i:[clusterSigObjL.query_starts[i][0], clusterSigObjR.query_starts[i][0]] for i in readNames}

            query_readName_LR = [ [i, [clusterSigObjL.query_starts[i][0], clusterSigObjR.query_starts[i][0]], clusterSigObjR.query_starts[i][2], abs(clusterSigObjR.query_starts[i][0] - clusterSigObjL.query_starts[i][0]) ]for i in readNames if (clusterSigObjL.query_starts[i][0] > max(seq_lens)) and (clusterSigObjR.query_starts[i][0]+max(seq_lens)  < clusterSigObjR.query_starts[i][1]) and (clusterSigObjL.query_starts[i][0] <  clusterSigObjR.query_starts[i][0])]
            if max(seq_svlens) >= abs(start - end)*1.3:
                clusterSigObjDup.idx = '3'
                self.callINS(contig, start, end, seq_svlens, readID_Pos,  clusterSigObjDup, pos_Ls, pos_Rs, Ls, Rs, readHashs_ratio, clusterSigObjL, clusterSigObjR, overlapReadHashName, sv_biase, pos_L, pos_R, query_readName_LR)
                return 

            #[6475, 5928]
            clusterSigObjDup.SVTYPE = 'DUP'
            #query_readName_LR = [ [i, [clusterSigObjL.query_starts[i][0], clusterSigObjR.query_starts[i][0]], clusterSigObjR.query_starts[i][2], abs(clusterSigObjR.query_starts[i][0] - clusterSigObjL.query_starts[i][0]) ]for i in readNames if (clusterSigObjL.query_starts[i][0] > max(seq_lens)) and (clusterSigObjR.query_starts[i][0]+max(seq_lens)  < clusterSigObjR.query_starts[i][1]) and (clusterSigObjL.query_starts[i][0] <  clusterSigObjR.query_starts[i][0])]
            print("query_LR:", query_readName_LR)
            #seq_svlens = self.remove_outliers_iqr(pd.Series(seq_svlens))
            
            #clusterSigObjDup.svlen_mode = int(np.mean(seq_svlens))
            clusterSigObjDup.svlen_std = np.std(svlens)
            clusterSigObjDup.idx = '1'
        
            sv_svlen = int(np.mean(seq_svlens))
            if new_start != None:
                genome_svlen = abs(new_start-new_end)
                svlens = [sv_svlen, genome_svlen]
                print("query_readName_LR1",query_readName_LR)
                #svlen_biase = np.mean(sv_biase)
                svlen_biase = sv_biase
                if  min(svlens)/max(svlens) >0.5:
                    if svlen_biase > 1000 and len(query_readName_LR) >0:
                        #new_start, new_end = self.dupRepeat(contig, new_start, new_end, svlen_biase)
                        #new_start, new_end = self.dupRepeat(contig, start, end, svlen_biase, readID_Pos)
                        #new_start, new_end =  new_start, new_end
                        print("WWWW1:", contig, new_start, new_end, svlen_biase, int(np.mean(seq_svlens)), query_readName_LR)
                        #new_start, new_end = self.fetch_reads_from_region(contig, new_start, new_end, [query_readName_LR[0][0]], svlen_biase, query_readName_LR[0][1], query_readName_LR[0][2])
                    print("QQQQ1:", new_start, new_end, svlen_biase)
                else:
                    print("out4")
                    return
            elif len(query_readName_LR) >0:
                genome_svlen = abs(start-end)
                svlens = [sv_svlen, genome_svlen]
                print("query_readName_LR2",sv_biase,min(svlens)/max(svlens), query_readName_LR)
                #svlen_biase = np.mean(sv_biase)
                svlen_biase = sv_biase
                if min(svlens)/max(svlens) >0.5:
                    if svlen_biase > 100 and len(query_readName_LR) >0:
                        print("WWWW2:", contig, start, end, svlen_biase, query_readName_LR)
                        new_start, new_end =self.dupRepeat(contig, start, end, svlen_biase, int(np.mean(seq_svlens)), readID_Pos)
                    print("QQQQ2:", new_start, new_end, svlen_biase)
                else:
                    print("out5")
                    return
        elif len( clusterSigObjR.readNames & clusterSigObjL.readNames) >=2:
            clusterSigObjDup.svlen_std = np.std(svlens)
            clusterSigObjDup.SVTYPE = 'DUP'    
            clusterSigObjDup.idx = '2'
        else:
            return 
        #if readHashs_ratio >0.8:
        #    return 

        #logger.info("clusterSigObj4")
        print("mutSigClusterDeal", new_start, new_end, flag2, clusterSigObjDup.support)
        if new_start != None:
            clusterSigObjDup.start_mode = new_start
            clusterSigObjDup.stop_mode = new_end
            clusterSigObjDup.svlen_mode = abs(new_start-new_end)
        else:
            clusterSigObjDup.start_mode = pos_Ls.idxmax()
            clusterSigObjDup.stop_mode = pos_Rs.idxmax()
            clusterSigObjDup.svlen_mode = abs(pos_Ls.idxmax() - pos_Rs.idxmax())
        ##############################################################################
        clusterSigObjDup.Rs = Rs
        clusterSigObjDup.Ls = Ls
        clusterSigObjDup.readHashs_ratio = readHashs_ratio

        #svlens = [np.abs([pos_L[i] - pos_R[i]]) for i in range(np.min([len(pos_L), len(pos_R)]))]
        #pos12 = [pos_Ls.idxmax(), pos_Rs.idxmax()]
        
        clusterSigObjDup.reverseReadsNus = [clusterSigObjL.reverseReadsNu, clusterSigObjDup.reverseReadsNu]
        clusterSigObjDup.secondaryReadsNus = [clusterSigObjL.secondaryReadsNu, clusterSigObjDup.secondaryReadsNu]
        mapQs = clusterSigObjDup.mapQs + clusterSigObjL.mapQs
        clusterSigObjDup.mapQ_mean = np.mean(mapQs)
        clusterSigObjDup.mapQ_std = np.std(mapQs)
        
        # clusterSigObjDup.svlen_mode = svlen#
        # clusterSigObjDup.svlen_std = np.std(svlens)#
        # clusterSigObjDup.start_mode = pos_Ls.idxmax()#
        # clusterSigObjDup.stop_mode = pos_Rs.idxmax()#
        
        clusterSigObjDup.start_std = np.std(pos_L)
        clusterSigObjDup.stop_std = np.std(pos_R)
       
        #clusterSigObjDup.support = np.max([clusterSigObjL.support, clusterSigObjR.support])
        clusterSigObjDup.support = np.max([clusterSigObjL.support, clusterSigObjR.support,  overlapReadHashName[-1]])
        clusterSigObjDup.depth = np.max([clusterSigObjL.depth, clusterSigObjR.depth])
        clusterSigObjDup.dr = clusterSigObjDup.depth - clusterSigObjDup.support
        # print("start2:", time.time() - start)
        start = time.time()
        # print("DV,DP:",clusterSigObjDup.support, clusterSigObjDup.depth)
        clusterSigObjDup.getGenotype_cuteSV()
        clusterSigObjDup.getGenotype_sniffles()
        # print("start3:", time.time() - start)
        #logger.info("clusterSigObj5")
        self.tree.add(Interval(clusterSigObjDup.start_mode, clusterSigObjDup.stop_mode, clusterSigObjDup))
        print("clusterSigObjInv3")
        flag, SVrecord = self.svRecord(clusterSigObjDup)
        print("clusterSigObjInv4", flag, SVrecord)
        #logger.info("clusterSigObj6")
        if flag:
            print("clusterSigObjInv5", flag, SVrecord)
            self.DUPs.append(SVrecord)
    
    def callINS(self,contig, start, end, seq_svlens, readID_Pos, clusterSigObjDup,pos_Ls,pos_Rs, Ls, Rs, readHashs_ratio, clusterSigObjL, clusterSigObjR, overlapReadHashName, svlen_biase, pos_L, pos_R, query_readName_LR):
        clusterSigObjIns = clusterSigObjDup
        new_start, new_end = None, None
        if svlen_biase > 100 and len(query_readName_LR) >0:
            new_start, new_end =self.dupRepeat(contig, start, end, svlen_biase, int(np.mean(seq_svlens)), readID_Pos)
        
        if new_start != None:
            clusterSigObjIns.start_mode = new_start
            clusterSigObjIns.stop_mode = new_start+1
            clusterSigObjIns.svlen_mode = abs(new_start-new_end)
        else:
            clusterSigObjIns.start_mode = pos_Ls.idxmax()
            clusterSigObjIns.stop_mode = pos_Ls.idxmax()+1
            clusterSigObjIns.svlen_mode = abs(pos_Ls.idxmax() - pos_Rs.idxmax())
        ##############################################################################
        clusterSigObjIns.Rs = Rs
        clusterSigObjIns.Ls = Ls
        clusterSigObjIns.readHashs_ratio = readHashs_ratio
        clusterSigObjIns.SVTYPE = 'INS'
        
        clusterSigObjIns.reverseReadsNus = [clusterSigObjL.reverseReadsNu, clusterSigObjIns.reverseReadsNu]
        clusterSigObjIns.secondaryReadsNus = [clusterSigObjL.secondaryReadsNu, clusterSigObjIns.secondaryReadsNu]
        mapQs = clusterSigObjIns.mapQs + clusterSigObjL.mapQs
        clusterSigObjIns.mapQ_mean = np.mean(mapQs)
        clusterSigObjIns.mapQ_std = np.std(mapQs)
        
        clusterSigObjIns.start_std = np.std(pos_L)
        clusterSigObjIns.stop_std = np.std(pos_R)
       
        clusterSigObjIns.support = np.max([clusterSigObjL.support, clusterSigObjR.support,  overlapReadHashName[-1]])
        clusterSigObjIns.depth = np.max([clusterSigObjL.depth, clusterSigObjR.depth])
        clusterSigObjIns.dr = clusterSigObjIns.depth - clusterSigObjIns.support

        start = time.time()

        clusterSigObjIns.getGenotype_cuteSV()
        clusterSigObjIns.getGenotype_sniffles()

        #self.tree.add(Interval(clusterSigObjIns.start_mode, clusterSigObjIns.stop_mode, clusterSigObjIns))
        print("clusterSigObjIns3")
        flag, SVrecord = self.svRecord(clusterSigObjIns)
        print("clusterSigObjIns4", flag, SVrecord)
        #logger.info("clusterSigObj6")
        if flag:
            print("clusterSigObjIns5", flag, SVrecord)
            self.INSs.append(SVrecord)

    
    def svRecord(self, clusterSigObj):
        flag = False
        SVid = "plotSV." + clusterSigObj.SVTYPE
        format = "GT1:GT2:DR:DV"
        af = clusterSigObj.support / clusterSigObj.depth
        reliability = "IMPRECISE"
        infos = ';'.join([reliability,
                          "raw_signal="+clusterSigObj.raw_signal,
                          ";BX="+str(clusterSigObj.bx),
                          ";SVTYPE=" + clusterSigObj.SVTYPE,
                          ";IDX="+clusterSigObj.idx,
                          ";readHashs_ratio="+ str(clusterSigObj.readHashs_ratio),
                          # "LNus="+str(clusterSigObj.LNus),
                          # "RNus="+str(clusterSigObj.RNus),
                          ";readNoHashRLNus=" + str(clusterSigObj.readNoHashRLNus),
                          ";readNamesRLNu=" + str(clusterSigObj.readNamesRLNu),
                          ";readHashsRLNu=" + str(clusterSigObj.readHashsRLNu),
                          ";mapQ_mean=" + str(clusterSigObj.mapQ_mean),
                          ";mapQ_std=" + str(clusterSigObj.mapQ_std),
                          ";ReverseReadsNu=" + ','.join(map(str, clusterSigObj.reverseReadsNus)),
                          ";SecondaryReadsNu=" + ','.join(map(str, clusterSigObj.secondaryReadsNus)),
                          ";SVLEN=" + str(clusterSigObj.svlen_mode),
                          ";END=" + str(clusterSigObj.stop_mode),
                          ";SUPPORT=" + str(clusterSigObj.support),
                          ";AF=" + str(round(af, 3)),
                          ";STD_POS1=" + str(round(clusterSigObj.start_std, 3)),
                          ";STD_POS2=" + str(round(clusterSigObj.stop_std, 3)),
                          ";STD_SVLEN=" + str(round(clusterSigObj.svlen_std, 3)),
                          ";RNAMES=NULL",
                          ])

        flag = self.filter(af, clusterSigObj)
        # if clusterSigObj.Ls_Rs_start_mean is None or clusterSigObj.Rs_Ls_start_mean is None:
        #    flag = False
        # elif abs(clusterSigObj.Ls_Rs_start_mean + clusterSigObj.Rs_Ls_start_mean) >100:
        #    flag = False
        return flag, [[clusterSigObj.chrom, clusterSigObj.start_mode, clusterSigObj.stop_mode, clusterSigObj.support], [clusterSigObj.chrom, clusterSigObj.start_mode, SVid, "N", "<" + clusterSigObj.SVTYPE + ">", '60',
                      "PASS", infos, format, ':'.join(
                [clusterSigObj.gt1, clusterSigObj.gt2, str(clusterSigObj.dr), str(clusterSigObj.support)])]]
    
    def filter(self, af, clusterSigObj):
        if af >= 0.001:
            return True
        else:
            return False
        pass

if __name__ == '__main__':
    genome = "/public3/SVDataset_lz/backup/2023_11_14/Genome/hg38/genome.fa"
    bamFiles = ["/public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/Pipelines/CCS/Sim.DEL.INS.DUP.INV.25x/winnowmap/bam/Sim.DEL.INS.DUP.INV.25x.sort.bam"]
    sampleIDs = ["SIM"]
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
