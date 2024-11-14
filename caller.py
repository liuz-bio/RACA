import os
import subprocess
import sys
import pysam
import numpy as np
import pandas as pd
import pytz
import traceback
import pickle
import multiprocessing
from loguru import logger
import asyncio
from callSV import callSV
# from pathos.multiprocessing import ProcessingPool
#from candidateSVLocation import candidateSVLocation
from candidateSVLocation1 import candidateSVLocation
from cluster import extraSigs
# from pathos.multiprocessing import Manager
from scipy.stats import nbinom, binom
from math import log, log10
from datetime import datetime
from tzlocal import get_localzone
from sklearn.cluster import DBSCAN
from intervaltree import IntervalTree, Interval
from datetime import datetime

current_path = os.path.abspath(os.path.dirname(__file__))
print("current_path:", current_path)
current_path = current_path.replace("callSV_new", "")
sys.path.append(current_path)
print("current_path:", current_path)
# from callSV_new.callDels import callDels
# from callSV_new.callInss import callInss
# from callSV_new.callInvs import callInvs
# from callSV_new.callDups import callDups
# from callSV_new.callTras import callTras
from Bio import bgzf
import json
import time
import copy
import warnings
import random
import threading
from dataclasses import dataclass
from ExtractSig import extractSignals
from collections import Counter
from collections import defaultdict


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
    min_support: int
    sv_types: list
    cluster_param: dict
    mut_proces: int
    chroms: list
    chrom_lengths: dict
    out_path: str
    tmp_path: str
    bam_cram: bool
    overlap_size:int


class caller:
    def __init__(self, configer):
        self.configer = configer
        self.chromLengthGroups = self.splitChrom(configer.chrom_lengths)

    def handle_error(self, e):
        print(''.join(traceback.format_exception(type(e), e, e.__traceback__)))
        print("An error occurred in the child process:", e)

    async def write_file(self, out_f, input_file):
        try:
            with open(input_file, 'rb') as in_f:
                data = in_f.read()
                out_f.write(data)
        except FileNotFoundError:
            pass

    async def merge_files_async(self, input_files, output_file):
        with open(output_file, 'wb') as out_f:
            tasks = [self.write_file(out_f, input_file) for input_file in input_files]
            await asyncio.gather(*tasks)

    def extract_duplicate_records(self, duplicate_records, all_beds, readBedName):#提取多reads记录
        id_counter = Counter()
        for bedName in all_beds:
            for line in open(bedName, 'r'):
                line = line.strip().split()
                if len(line) ==10:
                    continue
                else:
                    id_counter[line[9]]+=1

        with open(readBedName, 'wb') as out:
            for bedName in all_beds:
                with open(bedName, 'rb') as file:
                    lines = []
                    for line in file:
                        record_id = line.decode().strip().split('\t')[9]
                        if id_counter[record_id] >= 2:
                            lines.append(line)
                        else:
                            del id_counter[record_id]
                    out.write(b''.join(lines))

        with open(duplicate_records, 'wb') as plkOut:
            pickle.dump(id_counter, plkOut)

    def mergeSignalBeds(self, all_beds, sample):#合并提取的所有信号，提取多reads记录
        time_start = time.time()
        time_start1 = time.time()
        bedName = os.path.join(self.configer.tmp_path, sample, sample + '.bed')

        all_beds = sorted(all_beds, key=lambda x: (x[0], x[1]))
        all_beds = [i[2] for i in all_beds if os.path.exists(i[2])]
        asyncio.run(self.merge_files_async(all_beds, bedName))

        readBedName = bedName + ".readName"
        duplicate_records = str(os.path.join(self.configer.tmp_path, sample, sample + ".duplicate_records.plk"))
        self.extract_duplicate_records(duplicate_records, all_beds, readBedName)
        print("duplicate_records time:", time.time() - time_start)
        time_start = time.time()

        task = "env LC_ALL=C sort -S 200M -k10V  --unique  --parallel=10 %s >%s.sorted" % (
            readBedName, readBedName)#--compress-program=pzstd 
        tast1 = "env LC_ALL=C sort -S 200M -k1,1 -k2,2n  --parallel=10 %s |bgzip  -@ 10 >%s.gz" % (
            bedName, bedName)#--compress-program=pzstd 
        p = subprocess.Popen(task, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        p1 = subprocess.Popen(tast1, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        exit_code = p.wait()
        exit_code = p1.wait()
        pysam.tabix_index("%s.gz"%bedName, force=True, seq_col=0, start_col=1, end_col=2, preset="bed")
        candidateSVLocation(bedName + ".readName.sorted", self.configer.tmp_path, sample) #计算候选SV
        print("time AAAALL", time.time()-time_start)
        print("time AAAALL ALL", time.time()-time_start1)

    def mergeSignalPlks(self, sample, chrom, plks):#将同一染色体不同区域的区间树进行合并
        sigAllClusters = {'D': IntervalTree(), 'R': IntervalTree(), 'L': IntervalTree(), 'I': IntervalTree()}

        tmp_plks = []
        for plk in plks:
            if os.path.exists(plk):
                tmp_plks.append(plk)
                with open(plk, 'rb') as in_put:
                    for sigT, sigTree in pickle.load(in_put).items():
                        for node in sigTree:
                            sigAllClusters[sigT].addi(node.begin, node.end, node.data)

        if max([len(i) for i in sigAllClusters.values()]) >0:
            out_plk = str(os.path.join(self.configer.tmp_path, sample, sample+'.'+chrom + ".plk"))
            with open(out_plk, 'wb') as plkOut:
                pickle.dump(sigAllClusters, plkOut)

        for plk in tmp_plks:
            if os.path.exists(plk):
                print("tmp_plks:", tmp_plks)
                os.remove(plk)


    def splitChrom(self, chrom_lengths):#染色体划分多个区域便于并行
        chromLengthGroups = {}
        for chrom, length in chrom_lengths.items():
            length = int(length)
            chromLengthGroups[chrom] = []
            if length < self.configer.win_size:
                chromLengthGroups[chrom].append([1, self.configer.win_size])
            else:
                multipleTimes = length // self.configer.win_size
                for i in range(multipleTimes):
                    chromLengthGroups[chrom].append([i * self.configer.win_size + 1, (i + 1) * self.configer.win_size])
                if length % self.configer.win_size > 0:
                    chromLengthGroups[chrom].append([multipleTimes * self.configer.win_size + 1, (multipleTimes+1) * self.configer.win_size])
        return chromLengthGroups

    def run(self):
        t1 = time.time()
        for sample, bam_file in self.configer.samples_to_bams.items():
            process_pool = multiprocessing.Pool(self.configer.mut_proces)
            all_beds = []
            self.chroms_plk = {}

            os.makedirs(str(os.path.join(self.configer.tmp_path, sample)), exist_ok=True)
            for chrom, lengths in self.chromLengthGroups.items():
                no_first = False
                self.chroms_plk[chrom] = []
                for length in lengths:
                    start, stop = length[0], length[1]
                    plk_name = str(os.path.join(self.configer.tmp_path, sample, chrom + ".%d.%d.plk" % (start, stop)))
                    self.chroms_plk[chrom].append(plk_name)
                    tmpBed = str(
                        os.path.join(self.configer.tmp_path, sample, '_'.join(map(str, [chrom, start, stop])) + ".bed"))
                    all_beds.append([chrom, start, tmpBed])
                    process_pool.apply_async(func=extractSignals,
                                             args=(configer, bam_file, chrom, start, stop, no_first, tmpBed,
                                                   plk_name,),
                                             error_callback=self.handle_error)
                    no_first = True
            process_pool.close()
            process_pool.join()

            logger.info(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()) + ' ' + "signals deal" + ': ' + ' finish...')
            print("all Time1:", time.time() - t1)
            t2 = time.time()
            process_pool = multiprocessing.Pool(self.configer.mut_proces)
            manager = multiprocessing.Manager()
            lock = manager.Lock()
            all_cand_svs = manager.dict()
            process_pool.apply_async(func=self.mergeSignalBeds, args=(all_beds, sample,), error_callback=self.handle_error)
            for chrom, lengths in self.chromLengthGroups.items():
                no_first = False

                # process_pool.apply_async(func=extraSigs, 
                #                         args=(chrom, lengths, sample, self.configer, bam_file, all_cand_svs, lock, ),
                #                         error_callback=self.handle_error)
                for length in lengths:
                    start, stop = length[0], length[1]
                    plk_name = str(os.path.join(self.configer.tmp_path, sample, chrom + ".%d.%d.plk" % (start, stop)))
                    tmpBed = str(os.path.join(self.configer.tmp_path, sample, '_'.join(map(str, [chrom, start, stop])) + ".bed"))                
                    if os.path.exists(tmpBed):
                        process_pool.apply_async(func=extraSigs,
                                                 args=(chrom, start, stop, no_first, plk_name, tmpBed, configer, bam_file, all_cand_svs, lock, ),
                                                 error_callback=self.handle_error)
                    no_first = True

            process_pool.close()
            process_pool.join()
                 

            # process_pool = multiprocessing.Pool(self.configer.mut_proces)
            # for chrom, plks in self.chroms_plk.items():
            #     process_pool.apply_async(func=self.mergeSignalPlks,
            #                              args=(sample, chrom, plks,),
            #                              error_callback=self.handle_error)
            #
            # process_pool.apply_async(func=self.mergeSignalBeds,
            #                          args=(chroms_beds, sample, ),
            #                          error_callback=self.handle_error)
            # process_pool.close()
            # process_pool.join()
            #self.mergeSignalBeds(all_beds, sample)
            logger.info(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()) + ' ' + "signals merged" + ': ' + ' finish...')
            print("all Time2:", time.time()- t2)
            with open('chroms_plk.all_cand_svs.plk', 'wb') as f:
                pickle.dump(self.chroms_plk, f)
                pickle.dump(dict(all_cand_svs), f)
            inv = callSV(configer, sample, dict(self.chroms_plk), dict(all_cand_svs))
            logger.info(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()) + ' ' + "call SV" + ': ' + ' finish...')
            #####################################################################################################################################################


# if __name__ == "__main__":
#     genome = "/home/lz/Data_sequence/2023_11_14/Genome/hg38/genome.fa"
#     bamFiles = ["/home/lz/Data_sequence/2023_11_14/working/Pipelines/CCS/HG002/minimap2/bam/HG002.sort.bam"]
#                #/home/lz/work_space/Database/VISOR_test/Chr1_two_new/T2/FY0/VISOR_LASeR_Father_DEL.DUP.INS.INV/sim.srt.bam
#     sampleIDs = ["HG002"]#, "FY0_Mother", "FY0_Son", "FY1_Father", "FY1_Mother", "FY1_Son", "FY2_Father", "FY2_Mother", "FY2_Son"]
#     sv_types = ["DEL", "INS", "INV", "DUP", "BND"]
#     cluster_param = {'L': {'min_distanc': 150, 'min_signals': 2},
#                      'R': {'min_distanc': 150, 'min_signals': 2},
#                      'D': {'min_distanc': 150, 'min_signals': 2},
#                      'I': {'min_distanc': 150, 'min_signals': 2},
#                      'A': {'min_distanc': 50, 'min_signals': 5},
#                      'E': {'min_distanc': 50, 'min_signals': 5}}
#     chrom_lengths = dict([i.strip().split('\t')[:2] for i in open(genome + ".fai").read().strip().split('\n')])
#     out_path = "/home/lz/work_space/project/plotSV/RACA/plotSV9/callSV/tmp"
#     tmp_path = os.path.join(out_path, "tmp")
#     os.makedirs(tmp_path, exist_ok=True)

#     configer = Configuration(genome=genome,
#                              win_size=1000000,
#                              bam_files=bamFiles,
#                              sample_ids=sampleIDs,
#                              samples_to_bams=dict(zip(sampleIDs, bamFiles)),
#                              min_map_quality=45,
#                              min_sv_signal=25,
#                              min_sv_size=30,
#                              min_support=2,
#                              sv_types=sv_types,
#                              cluster_param=cluster_param,
#                              mut_proces=20,
#                              chroms=list(chrom_lengths.keys()),
#                              chrom_lengths=chrom_lengths,
#                              out_path=out_path,
#                              tmp_path=tmp_path,
#                              bam_cram=True,
#                              overlap_size=2000,
#                              )

#     cal = caller(configer)
#     cal.run()
    
if __name__ == "__main__":
    #genome = "/home/lz/work_space/Database/VISOR_test/Genome_old/Genome.fa"
    genome = "/public3/SVDataset_lz/backup/2023_11_14/Genome/hg38/lra_CCS/genome.fa"
    bamFiles = ["/home/lz/work_space/Database/VISOR_test/Chr1_two_new/T2/FY2/VISOR_LASeR_Father_DEL.DUP.INS.INV/sim.srt.bam",
                "/home/lz/work_space/Database/VISOR_test/Chr1_two_new/T2/FY1/VISOR_LASeR_Mother_DEL.DUP.INS.INV/sim.srt.bam",
               "/public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/Pipelines/CCS/Sim.DEL.INS.DUP.INV.25x/minimap2/bam/Sim.DEL.INS.DUP.INV.25x.sort.bam"]
    bamFiles = ["/public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/Pipelines/CCS/Sim.DEL.INS.DUP.INV.25x/minimap2/bam/Sim.DEL.INS.DUP.INV.25x.sort.bam"]
    #sampleIDs = ["FY2_Father","FY1_Mother","SIM"]#, "FY0_Mother", "FY0_Son", "FY1_Father", "FY1_Mother", "FY1_Son", "FY2_Father", "FY2_Mother", "FY2_Son"]
    sampleIDs = ["SIM"]
    sv_types = ["DEL", "INS", "INV", "DUP", "BND"]
    cluster_param = {'L': {'min_distanc': 100, 'min_signals': 2},
                     'R': {'min_distanc': 100, 'min_signals': 2},
                     'D': {'min_distanc': 250, 'min_signals': 2},
                     'I': {'min_distanc': 250, 'min_signals': 2},
                     'A': {'min_distanc': 50, 'min_signals': 5},
                     'E': {'min_distanc': 50, 'min_signals': 5}}
    chrom_lengths = dict([i.strip().split('\t')[:2] for i in open(genome + ".fai").read().strip().split('\n')])
    out_path = "/home/lz/work_space/project/plotSV/RACA/plotSV9/callSV/tmp"
    tmp_path = os.path.join(out_path, "tmp")
    os.makedirs(tmp_path, exist_ok=True)

    configer = Configuration(genome=genome,
                             win_size=1000000,
                             bam_files=bamFiles,
                             sample_ids=sampleIDs,
                             samples_to_bams=dict(zip(sampleIDs, bamFiles)),
                             min_map_quality=45,
                             min_sv_signal=25,
                             min_sv_size=30,
                             min_support=2,
                             sv_types=sv_types,
                             cluster_param=cluster_param,
                             mut_proces=20,
                             chroms=list(chrom_lengths.keys()),
                             chrom_lengths=chrom_lengths,
                             out_path=out_path,
                             tmp_path=tmp_path,
                             bam_cram=True,
                             overlap_size=2000,
                             )

    cal = caller(configer)
    cal.run()
