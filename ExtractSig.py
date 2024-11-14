import gzip
import sys
import multiprocessing
from intervaltree import IntervalTree, Interval
import pysam
import pickle
import traceback
import os
#print("current_path")
#current_path = current_path.replace("callSV", "")
#current_path = current_path.replace("utils", "")
#print("当前路径为：%s ExtractSignal" % current_path)
parent_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_path)
from utils.candidateSVLocation import candidateSVLocation
from loguru import logger
import ctypes
import json
import time
import mmap
import subprocess
import numpy as np
from Bio import bgzf
from cluster import cluster


#提取SV信号（D,I,R,L）

def extractSignals_old(configer, bamfile, chrom, interval_start, interval_end, no_first, tmpBed, plk_name):
    min_map_quality = configer.min_map_quality
    lib = ctypes.CDLL('%s/golang/cigar.so'%parent_path)
    lib.DealWith.argtypes = [ctypes.c_char_p,
                             ctypes.POINTER(ctypes.c_int),
                             ctypes.c_char_p,
                             ctypes.POINTER(ctypes.c_int),
                             ]

    lib.DealWith.restype = ctypes.c_char_p

    if configer.bam_cram:
        handle = pysam.AlignmentFile(bamfile, "rb", require_index=True)
    else:
        handle = pysam.AlignmentFile(bamfile, "rc", require_index=True)


    allreads = handle.fetch(contig=chrom, start=interval_start, stop=interval_end)
    c_min_signal = ctypes.c_int(configer.min_sv_signal)

    allSignals = []
    start_time = time.time()
    for read in allreads:
        mapping_quality = read.mapping_quality

        if mapping_quality < min_map_quality:
            continue

        if no_first & (read.reference_start <= interval_start <= read.reference_end):
            continue

        is_reverse = '1' if read.is_reverse else '0'
        is_supplementary = '1' if read.is_supplementary else '0'
        query_name = read.query_name.split('|')[0]


        signals = lib.DealWith(read.cigarstring.encode(),
                               ctypes.c_int(read.reference_start),
                               read.query_sequence.encode(),
                               c_min_signal,
                              ).decode().strip("\t")

        allSignals.append([
            chrom,
            read.reference_start,
            read.reference_end,
            read.query_alignment_start,
            read.query_alignment_end,
            read.query_length,
            mapping_quality,
            is_reverse,
            is_supplementary,
            query_name,
            signals
        ])


    start_time2 = time.time()
    if len(allSignals) > 0:
        encoded_data = '\n'.join(['\t'.join(map(str, i)) for i in allSignals]) + "\n"
        encoded_data = encoded_data.encode('ascii')
        f = open(tmpBed, "wb")
        f.write(len(encoded_data)*b'\0')
        f.flush()
        f.close()
        with open(tmpBed, "r+b") as f:
            with mmap.mmap(f.fileno(), 0) as mm:
                mm.seek(0)
                mm.write(encoded_data)



def extractSignals(configer, bamfile, contig, interval_start, interval_end, no_first, tmpBed, plk_name):
    min_map_quality = configer.min_map_quality
    if configer.bam_cram:
        handle = pysam.AlignmentFile(bamfile, "rb", require_index=True)
    else:
        handle = pysam.AlignmentFile(bamfile, "rc", require_index=True)

    if interval_start >configer.overlap_size:
        interval_start -= configer.overlap_size
    interval_end   += configer.overlap_size
        
    allreads = handle.fetch(contig=contig, start=interval_start, stop=interval_end)
    allSignals = []
    for read in allreads:
        mapping_quality = read.mapping_quality

        if mapping_quality < 10 or read.query_length<1000:
             continue

        #if no_first & (read.reference_start <= interval_start <= read.reference_end):
        #   continue

        is_reverse = '1' if read.is_reverse else '0'
        is_supplementary = '1' if read.is_supplementary else '0'
        query_name = read.query_name.split('|')[0]
        
        cigar_tuples = read.cigartuples
        read_sequence = read.query_sequence

        ref_pos = read.reference_start
        read_pos = 0
        # 遍历 CIGAR 操作
        signals = []
        for cigar_op, length in cigar_tuples:
            if cigar_op == 0:  # M (Match)
                ref_pos += length
                read_pos += length
            elif cigar_op == 1:  # I (Insertion)
                if length >=configer.min_sv_signal:
                    signals.append(','.join(map(str, [ref_pos, ref_pos+length, 3, read_sequence[read_pos:read_pos+length], read_pos])))
                read_pos += length
            elif cigar_op == 2:  # D (Deletion)
                if length >=configer.min_sv_signal:
                    signals.append(','.join(map(str, [ref_pos, ref_pos+length, 2, read_pos])))
                ref_pos += length
            elif cigar_op == 4:  # S (Soft clip)
                if read_pos == 0:
                    signals.append(','.join(map(str, [ref_pos - length, ref_pos, 1, read_pos+length])))
                else:
                    signals.append(','.join(map(str, [ref_pos, ref_pos+ length, 4, read_pos])))
                read_pos += length
            elif cigar_op == 3:  # N (Skipped region from the reference)
                ref_pos += length
      
        allSignals.append([
            contig,
            read.reference_start,
            read.reference_end,
            read.query_alignment_start,
            read.query_alignment_end,
            read.query_length,
            mapping_quality,
            is_reverse,
            is_supplementary,
            query_name,
            '\t'.join(signals)
        ])


    start_time2 = time.time()
    if len(allSignals) > 0:
        encoded_data = '\n'.join(['\t'.join(map(str, i)) for i in allSignals]) + "\n"
        encoded_data = encoded_data.encode('ascii')
        f = open(tmpBed, "wb")
        f.write(len(encoded_data)*b'\0')
        f.flush()
        f.close()
        with open(tmpBed, "r+b") as f:
            with mmap.mmap(f.fileno(), 0) as mm:
                mm.seek(0)
                mm.write(encoded_data)

