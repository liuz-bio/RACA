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
import warnings
from scipy.stats import nbinom, binom

class readSig:
    def __init__(self, start, stop, mapQ, readHash, read_is_reverse, read_is_secondary, readName, seq=None):
        self.start = start
        self.stop  = stop
        self.read_is_reverse = read_is_reverse
        self.read_is_secondary = read_is_secondary
        self.readName = readName
        self.seq = seq
        self.mapQ = mapQ
        self.readHash = readHash


class clusterSig:
    def __init__(self, chrom, cluster_xy, cluster, bedGz, sigType, fingerPrint):
        self.chrom = chrom
        self.readNames = set([i.readName for i in cluster])
        self.readHashs = set([i.readHash for i in cluster])
        self.mapQs = [i.mapQ for i in cluster]
        self.support = len(self.readNames)
        self.cluster = cluster
        self.fingerPrint = fingerPrint
        self.sigType = sigType
        if sigType in ['D', 'I']:
            self.mapQ_mean = np.mean(self.mapQs)
            self.mapQ_std = np.std(self.mapQs)
            self.clusterStatistic(cluster_xy, cluster, bedGz, sigType)
            # print("cluster_xy_DIs:", list(cluster_xy))
        else:
            self.start_mode = int(np.mean(cluster_xy[:, 0]))
            self.stop_mode = self.start_mode + 1
            self.reverseReadsNu = np.sum([1 for i in cluster if i.read_is_reverse])
            self.secondaryReadsNu = np.sum([1 for i in cluster if i.read_is_secondary])
            self.PosLocation = cluster_xy[:, 0]

    def clusterStatistic(self, cluster_xy, cluster, bedGz, sigType):
        start = time.time()
        x = cluster_xy[:, 0]
        y = cluster_xy[:, 1]
        svlens = np.abs(x - y)
        pos_x = pd.Series(x).value_counts()
        pos_y = pd.Series(y).value_counts()

        if len(x) < 3:
            self.start_mode = np.max(x)
            x_index = np.where(x == self.start_mode)[0][0]
            self.stop_mode = y[x_index]
            self.svlen_mode = svlens[x_index]
            if sigType == 'I':
                self.alt = cluster[x_index].seq
            else:
                self.alt = "<DEL>"
        elif pos_x.max() > pos_y.max():
            self.start_mode = pos_x.idxmax()
            x_index = np.where(x == self.start_mode)[0][0]
            self.stop_mode = y[x_index]
            self.svlen_mode = svlens[x_index]
            if sigType == 'I':
                self.alt = cluster[x_index].seq
            else:
                self.alt = "<DEL>"
        else:
            self.stop_mode = pos_y.idxmax()
            y_index = np.where(y == self.stop_mode)[0][0]
            self.start_mode = x[y_index]
            self.svlen_mode = svlens[y_index]
            if sigType == 'I':
                self.alt = cluster[y_index].seq
            else:
                self.alt = "<DEL>"

        self.start_std = np.std(x)
        self.stop_std = np.std(y)
        self.svlen_std = np.std(svlens)

        # print("sigType:", sigType, self.start_mode, self.stop_mode, self.svlen_mode, self.readNames)
        # print("stat1:",time.time() - start)
        start = time.time()
        self.getDepth(min(x), max(y), bedGz)  # self.depth
        self.dr = self.depth - self.support  # self.dr
        # print("stat2:",time.time() - start)
        start = time.time()
        self.getGenotype_cuteSV()  # self.gt1
        # print("stat3:",time.time() - start)
        start = time.time()
        self.getGenotype_sniffles()  # self.gt2
        # print("stat4:",time.time() - start)

    def getDepth(self, start, stop, bedGz, chrom=None):
        if chrom is None:
            allSignals = bedGz.fetch(reference=self.chrom, start=start, end=stop)
        else:
            allSignals = bedGz.fetch(reference=chrom, start=start, end=stop)
        tmp_readNames = {}
        while True:
            try:
                tmp = next(allSignals).split('\t')
                if int(tmp[6])-500 >= 20:
                    tmp_readNames[tmp[9]] = 0
                # num += 1
            except StopIteration:
                break
        self.depth = len(tmp_readNames)

    def caculation_qual(self):
        pass

    def nomalized_count(self, depth, support):
        normalization_target = 100
        if depth > normalization_target:
            norm = normalization_target / float(depth)
            normalized_support = round(support * norm)
            normalized_depth = round(depth * norm)
            normalized_ref = normalized_depth - normalized_support
        else:
            normalized_support = support
            normalized_depth = depth
            normalized_ref = normalized_depth - normalized_support
        if normalized_ref < 0:
            print(depth, support)
        return normalized_support, normalized_depth, normalized_ref

    def getGenotype_sniffles(self):
        genotypes_p = {'0/0': 0.1, '0/1': 0.5, '1/1': 0.9}
        genotypes = ['0/0', '0/1', '1/1']
        # genotypes_p = {'0/0': 0.1, '0/1': 0.5, '1/1': 0.9}
        normalized_support, normalized_depth, normalized_ref = self.nomalized_count(self.depth, self.support)
        genotype_likelihoods = []
        for gt in genotypes:
            # q = nbinom.pmf(normalized_support, normalized_ref, genotypes_p[gt])
            # q = binomial_probability(normalized_support, normalized_depth, p)
            q = binom.pmf(normalized_support, normalized_depth, genotypes_p[gt])
            genotype_likelihoods.append(q)
        sum_likelihoods = sum(genotype_likelihoods)
        try:  # chr4 D [0.0, 0.0, 0.0] 0 2 2 0 -2
            with warnings.catch_warnings():
                warnings.filterwarnings("error")
                normalized_likelihoods = [q / sum_likelihoods for q in genotype_likelihoods]
                self.gt2 = genotypes[np.argmax(normalized_likelihoods)]
        except RuntimeWarning:
            #0 16 16 0 -16
            print("RuntimeWarning:", self.chrom, self.start_mode, self.stop_mode, self.readNames, self.sigType, genotype_likelihoods, self.depth, self.support,
                  normalized_support, normalized_depth, normalized_ref)
            self.gt2 = '0/0'

    def getGenotype_cuteSV(self):
        Prob = float(1 / 3)
        align_err = float(0.1)
        genotypes = ['0/0', '0/1', '1/1']

        normalized_support, normalized_depth, normalized_ref = self.nomalized_count(self.depth, self.support)

        ref_p = np.float64((1 - Prob) / 2 * pow((1 - align_err), normalized_ref) * pow(align_err, normalized_support))
        het_p = np.float64((1 - Prob) / 2 * pow(0.5, normalized_depth))
        # print(Prob,align_err,normalized_support,align_err, normalized_ref)
        var_p = np.float64(Prob * pow((1 - align_err), normalized_support) * pow(align_err, normalized_ref))

        genotypes_p = {'0/0': ref_p, '0/1': het_p, '1/1': var_p}

        # GL_P = [q / sum_likelihoods) for  q in genotype_likelihoods]
        # PL = [int(np.around(-10*log10(i))) for i in GL_P]
        # GQ = [int(-10*log10(GL_P[1] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[1]))]
        # QUAL = abs(np.around(-10*log10(GL_P[0]), 1))

        self.gt1 = genotypes[np.argmax([ref_p, het_p, var_p])]

class cluster:
    def __init__(self, bed, chromLengths):
        self.bed = bedle

        self.winSize = 10000000
        self.clusterParam = {'L': {'min_distanc': 150, 'min_signals': 2},
                             'R': {'min_distanc': 150, 'min_signals': 2},
                             'D': {'min_distanc': 300, 'min_signals': 2},
                             'I': {'min_distanc': 300, 'min_signals': 2},
                             'A': {'min_distanc': 50, 'min_signals': 5},
                             'E': {'min_distanc': 50, 'min_signals': 5}}
        self.chromLengths = chromLengths
        self.min_signal = 25
        self.sigTypes = {'1': 'L', '2': 'D', '3': 'I', '4': 'R', '5': 'A'}
        self.chromLengthGroups = self.splitChrom(chromLengths)

    def splitChrom(self, chromLengths):
        chromLengthGroups = {}
        for chrom, length in chromLengths.items():
            length = int(length)
            chromLengthGroups[chrom] = []
            if length < self.winSize:
                chromLengthGroups[chrom].append([1, length])
            else:
                multipleTimes = length // self.winSize
                for i in range(multipleTimes):
                    chromLengthGroups[chrom].append([i * self.winSize + 1, (i + 1) * self.winSize])
                if length % self.winSize > 0:
                    chromLengthGroups[chrom].append([multipleTimes * self.winSize + 1, length])
        return chromLengthGroups

    def cluster(self, sigs, dbscan, chrom, bedGz, sigT, fingerPrint, sigAllClusters):
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

            if len(set([i.readName for i in cluster])) <= 75:
                clusterSigObj = clusterSig(chrom, cluster_xy, cluster, bedGz, sigT, fingerPrint)
                if sigT in ['D', 'I'] and clusterSigObj.depth == 0:
                    continue
                fingerPrint += 1
                sigAllClusters[sigT].add(Interval(clusterSigObj.start_mode, clusterSigObj.stop_mode, clusterSigObj))
                #pickle.dump(clusterSigObj, plkOut)
                #clusterSigObj = None
                # sigClusters.add(Interval(clusterSigObj.start_mode, clusterSigObj.stop_mode, clusterSigObj))
        return fingerPrint

    def write_plk(self, file_plk, sigAllClusters):
        with open(file_plk, 'wb') as plkOut:
            pickle.dump(sigAllClusters, plkOut)
    def extraSigs(self, chrom, start, end, no_first, plk_name):
        sigAllClusters = {'D': IntervalTree(), 'R': IntervalTree(), 'L': IntervalTree(), 'I': IntervalTree()}
        bedGz = pysam.TabixFile(self.bed)
        try:
            allSignals = bedGz.fetch(reference=chrom, start=start, end=end)
        except ValueError:
            return sigAllClusters

        sigs = {}
        read_hash_nu = 0
        fingerPrint = 0
        num_nu = 0

        for read in allSignals:
            read_hash = read_hash_nu
            read_hash_nu += 1
            read = read.split('\t')
            map_q = int(read[6])
            align_start = int(read[1])
            align_stop = int(read[2])
            if no_first & (align_start<= end <= align_stop):
                continue

            read_is_reverse = int(read[7])
            read_is_secondary = int(read[8])
            read_name = read[9]

            for sig in read[10:]:
                sig = sig.split(',')
                start = int(sig[0])
                stop = int(sig[1])
                sigType = self.sigTypes[sig[2]]

                if sigType in ['D', 'I']:
                    if stop - start < self.min_signal:
                        continue

                    if sigType == "I":
                        read_sig_obj = readSig(start, stop, map_q, read_hash, read_is_reverse, read_is_secondary, read_name,
                                             seq=sig[3])
                    else:
                        read_sig_obj = readSig(start, stop, map_q, read_hash, read_is_reverse, read_is_secondary, read_name)
                    try:
                        sigs[sigType].append([start, stop, read_sig_obj])
                    except KeyError:
                        sigs[sigType] = [[start, stop, read_sig_obj]]
                elif sigType == "L":
                    read_sig_obj = readSig(stop, stop, map_q, read_hash, read_is_reverse, read_is_secondary, read_name)
                    try:
                        sigs[sigType].append([stop, stop, read_sig_obj])
                    except KeyError:
                        sigs[sigType] = [[stop, stop, read_sig_obj]]
                elif sigType == "R":
                    read_sig_obj = readSig(start, start, map_q, read_hash, read_is_reverse, read_is_secondary, read_name)
                    try:
                        sigs[sigType].append([start, start, read_sig_obj])
                    except KeyError:
                        sigs[sigType] = [[start, start, read_sig_obj]]

            if num_nu >= 200:
                for sigT, sig in sigs.items():
                    try:
                        dis = np.linalg.norm(np.array(sig[-1][:2]) - np.array(sig[-2][:2]))
                    except IndexError:
                        dis = 0

                    if dis >= 1000:
                        fingerPrint = self.cluster(sig[:-1], DBSCAN(eps=self.clusterParam[sigT]['min_distanc'],
                                                                    min_samples=self.clusterParam[sigT]['min_signals']),
                                                   chrom, bedGz, sigT, fingerPrint, sigAllClusters)
                        sigs[sigT] = [sig[-1]]
                    elif len(sig) >= 1500:
                        sigs[sigT] = []

                num_nu = max([len(i) for i in sigs.values()])
                #print("num_nu:", num_nu, chrom, [[i[-1][:2], i[-2][:2]] for i in sigs.values() if len(i) > 1],
                #      [np.linalg.norm(np.array(i[-1][:2]) - np.array(i[-2][:2])) for i in sigs.values() if len(i) > 1])
            else:
                num_nu += 1

        for sigT, sig in sigs.items():
            if len(sig) <= 1500:
                fingerPrint = self.cluster(sig, DBSCAN(eps=self.clusterParam[sigT]['min_distanc'],
                                                       min_samples=self.clusterParam[sigT]['min_signals']),
                                           chrom, bedGz, sigT, fingerPrint, sigAllClusters)

        self.write_plk(plk_name, sigAllClusters)
        #return sigAllClusters


    def handle_error(self, e):
        print(''.join(traceback.format_exception(type(e), e, e.__traceback__)))
        print("An error occurred in the child process:", e)

    def mergeChrom(self, chrom, tasks):
        sigAllClusters = {'D': IntervalTree(), 'R': IntervalTree(), 'L': IntervalTree(), 'I': IntervalTree()}
        for task in tasks:
            try:
                with open(task, 'rb') as in_put:
                    for sigT, sigTree in pickle.load(in_put).items():
                        for node in sigTree:
                            sigAllClusters[sigT].addi(node.begin, node.end, node.data)
                        #sigAllClusters[sigT] = IntervalTree(list(sigAllClusters[sigT]) + list(sigTree))
                os.remove(task)
            except FileNotFoundError:
                #print("FileNotFoundError:", task)
                continue

        for sigT, sigs in sigAllClusters.items():
            if len(sigs) > 0:
                with open(chrom + ".%s.plk" % sigT, 'wb') as plkOut:
                    #print("sigs:", sigs)
                    pickle.dump(sigs, plkOut)


    def run(self):
        self.mutProces = 20
        process_pool = multiprocessing.Pool(self.mutProces)
        tasks_dict = {}
        for chrom, datas in self.chromLengthGroups.items():
            no_first = False
            #print("datas:", datas)
            tasks_dict[chrom] = []
            for data in datas:
                #print("chrom:", chrom, data[0], data[1], no_first)
                plk_name = chrom+".%d.%d.plk"%(data[0], data[1])
                process_pool.apply_async(func=self.extraSigs,
                                                args=(chrom, data[0], data[1], no_first, plk_name, ),
                                                error_callback=self.handle_error)
                no_first = True
                tasks_dict[chrom].append(plk_name)

        process_pool.close()
        process_pool.join()

        process_pool = multiprocessing.Pool(self.mutProces)
        for chrom, tasks in tasks_dict.items():
            process_pool.apply_async(func=self.mergeChrom,
                                     args=(chrom, tasks,),
                                     error_callback=self.handle_error)
        process_pool.close()
        process_pool.join()



if __name__ == "__main__":
    genomeFile = '/home/liuz/work_space/project/liuz/Genome/genome.fa'
    # bed = '/home/lz/work_space/Database/VISOR_test/2024_3_1/plotSV/callSV_new/SIM.bed.gz'
    bed = '/home/liuz/work_space/project/liuz/plotSV/plotSV/callSV_new/SIM.bed.gz'
    mutProces = 20
    genome = pysam.FastaFile(genomeFile)
    chromLengths = dict(zip(genome.references, genome.lengths))
    clu = cluster(bed, chromLengths)
    clu.run()