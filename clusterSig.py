import numpy as np
import pandas as pd
import time
import warnings
from scipy.stats import nbinom, binom
from util import kmeans
# from shapely.geometry import LineString


class readSig:
    def __init__(self, start, stop, mapQ, readHash, read_is_reverse, read_is_secondary, readName,
                 seq=None,  query_start=None, query_len=None, sLR_other=None, query_length= None):
        self.start = start
        self.stop  = stop
        self.read_is_reverse = read_is_reverse
        self.read_is_secondary = read_is_secondary
        self.readName = readName
        self.seq = seq
        self.mapQ = mapQ
        self.readHash = readHash
        self.query_start = query_start
        self.query_len = query_len
        self.sLR_other = sLR_other #LS 或RS 另一侧的断点在参考基因组的位置
        self.query_length = query_length
        # self.read.query_alignment_start,
        # self.read.query_alignment_end,
        #read.query_length,
        #self.svlen = svlen


class clusterSig:
    def __init__(self, chrom, cluster_xy, cluster, handle, sigType, fingerPrint):
        self.chrom = chrom
        self.recordNu = len(cluster)
        self.readNames = set([i.readName for i in cluster])
        self.readHashs = set([i.readHash for i in cluster])
        
        self.query_starts = dict([[i.readName, [i.query_start, i.query_len, i.read_is_reverse]] for i in cluster])
        self.mapQs = [i.mapQ for i in cluster]
        self.support = len(self.readNames)
        self.cluster = cluster
        self.fingerPrint = fingerPrint
        self.sigType = sigType
        if sigType in ['D', 'I']:
            self.mapQ_mean = np.mean(self.mapQs)
            self.mapQ_std = np.std(self.mapQs)
            self.clusterStatistic(cluster_xy, cluster, handle, sigType)
            # print("cluster_xy_DIs:", list(cluster_xy))
        else:     
            self.start_mode = int(np.mean(cluster_xy[:, 0]))
            self.stop_mode = self.start_mode + 1
            self.reverseReadsNu = np.sum([1 for i in cluster if i.read_is_reverse])
            self.secondaryReadsNu = np.sum([1 for i in cluster if i.read_is_secondary])
            self.PosLocation = cluster_xy[:, 0]
            self.getDepth(self.start_mode, self.stop_mode, handle)
    
    def clusterStatistic(self, cluster_xy, cluster, handle, sigType):
        start = time.time()
        x = cluster_xy[:, 0]
        y = cluster_xy[:, 1]
        svlens = np.abs(x - y)
        # svlen_std = np.std(svlens)
        # svlen_mean = np.mean(svlens)
        # sub_svlens = []
        # if svlen_std/svlen_mean >=0.1:
        #     sub_svlens = kmeans(svlens)

        # if len(sub_svlens) >= 2:
        #     self.sub_support = True
        #     self.svlen_std = np.std(sub_svlens)
        # else:
        #     self.sub_support = False
        #     self.svlen_std = np.std(svlens)

        pos_x = pd.Series(x).value_counts()
        pos_y = pd.Series(y).value_counts()

        if len(x) < 3:
            self.start_mode = np.max(x)
            x_index = np.where(x == self.start_mode)[0][0]
            self.stop_mode = y[x_index]
            self.svlen_mode = max(svlens) #svlens[x_index]
            if sigType == 'I':
                self.alt = cluster[x_index].seq
            else:
                self.alt = "<DEL>"
        elif pos_x.max() > pos_y.max():
            self.start_mode = pos_x.idxmax()
            x_index = np.where(x == self.start_mode)[0][0]
            self.stop_mode = y[x_index]
            self.svlen_mode = int(np.mean(svlens)) #svlens[x_index]
            if sigType == 'I':
                self.alt = cluster[x_index].seq
            else:
                self.alt = "<DEL>"
        else:
            self.stop_mode = pos_y.idxmax()
            y_index = np.where(y == self.stop_mode)[0][0]
            self.start_mode = x[y_index]
            self.svlen_mode = int(np.mean(svlens)) #svlens[y_index]
            if sigType == 'I':
                self.alt = cluster[y_index].seq
            else:
                self.alt = "<DEL>"

        self.start_std = np.std(x)
        self.stop_std = np.std(y)
        self.svlen_std = np.std(svlens)
        self.svlens = svlens

        #print("sigType:", sigType, self.start_mode, self.stop_mode, self.svlen_mode, self.readNames)
        print("stat1:",time.time() - start)
        start = time.time()
        self.getDepth(min(x), max(y), handle)  # self.depth
        self.dr = self.depth - self.support  # self.dr
        # print("stat2:",time.time() - start)
        start = time.time()
        self.getGenotype_cuteSV()  # self.gt1
        # print("stat3:",time.time() - start)
        start = time.time()
        self.getGenotype_sniffles()  # self.gt2
        print("stat4:",time.time() - start)

    # def count_intersections(self, intervals, query):
    #     query_line = LineString([(query[0], 0), (query[1], 0)])
    #     count = 0
        
    #     for start, end in intervals:
    #         line = LineString([(start, 0), (end, 0)])
    #         if query_line.intersects(line):
    #             count += 1
        
    #     return count
    
    def getDepth(self, start, stop, handle, chrom=None):
        
        if handle is None:
            self.depth = 1
        else:
            start = start -100
            stop = stop + 100
            if start  <0:
                start = 1
                
            if chrom is None:
                start_time = time.time()
                #self.depth = len(handle.overlap(start, stop))
                self.depth = sum([1 for i in handle if i[0] <= start <= i[1] or i[0] <= stop <= i[1]])
                # self.depth = self.count_intersections(handle,(start, stop))
                if self.depth == 0:
                    self.depth = len(self.readNames)
                print("getDepthWWW:",self.depth, start, stop, time.time()-start_time)
                return
                #allSignals = handle.fetch(self.chrom, start, stop)
            else:
                #self.depth = len(handle.overlaps(start, stop))
                
                allSignals = handle.fetch(chrom, start, stop)
            tmp_readNames = {}

            for read in allSignals:
                #if read.mapping_quality < 20:
                #    continue
                tmp_readNames[read.query_name] = 0

            # while True:
            #     try:
            #         tmp = next(allSignals).split('\t')
            #         if int(tmp[6])>= 20:
            #             tmp_readNames[tmp[9]] = 0
            #         # num += 1
            #     except StopIteration:
            #         break
            self.depth = len(tmp_readNames)
            print("getDepthEEE:",self.depth)

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
            normalized_ref = 1
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
