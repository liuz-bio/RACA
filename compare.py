import logging
import argparse
import numpy as np
import time
import os
import sys
import multiprocessing
# 配置日志
logging.basicConfig(filename='app.log', level=logging.INFO)


class compareSV:
    def __init__(self, commVcf=None, baseVcf=None, base_h1=None, base_h2=None, rangeSize=1000, chroms=None, outPath="./out", svtypes=None):
        self.rangeSize = rangeSize
        self.chroms = chroms
        self.outPath = outPath
        self.svtypes = svtypes
        start = time.time()
        self.dealBaseComm(commVcf, baseVcf, base_h1, base_h2)
        print("start1:", time.time() - start)
        start = time.time()
        self.compare()
        print("start2:", time.time() - start)

    def process_chromosome(self, chrom, baseRecords):
        all_baseRecords = []
        all_commRecords = []
        tp_baseRecords = []
        tp_commRecords = []
        all_baseRecords.extend([ i[-1] for i in baseRecords])
        commRecords = self.comm.get(chrom)
        if commRecords != None:
            all_commRecords.extend([ i[-1] for i in commRecords])
            for commRecord in commRecords:
                for baseRecord in baseRecords:
                    basePos = sorted(baseRecord[1:3])
                    commPos = sorted(commRecord[1:3])
                    ratio = min([baseRecord[3], commRecord[3]])/max([baseRecord[3], commRecord[3]])
                    print("BBBB:", baseRecord, commRecord)
                    if  ((np.abs(basePos[0]-commPos[0]) <= self.rangeSize) or (np.abs(basePos[1]-commPos[1]) <= self.rangeSize)) and (baseRecord[4]==commRecord[4]) and (ratio>=0.7):
                        tp_commRecords.append(commRecord[-1])
                        print("BBBB:", baseRecord, commRecord)
                        break


            for baseRecord in baseRecords:
                for commRecord in commRecords:
                    basePos = sorted(baseRecord[1:3])
                    commPos = sorted(commRecord[1:3])
                    ratio = min([baseRecord[3], commRecord[3]])/max([baseRecord[3], commRecord[3]])
                    if  ((np.abs(basePos[0]-commPos[0]) <= self.rangeSize) or (np.abs(basePos[1]-commPos[1]) <= self.rangeSize)) and (baseRecord[4]==commRecord[4]) and (ratio>=0.7):
                        print("EEEE:", baseRecord, commRecord)
                        tp_baseRecords.append(baseRecord[-1])
                        break

        return tp_baseRecords, tp_commRecords, all_baseRecords, all_commRecords

    def compare(self):
        all_commRecords = []
        all_baseRecords = []
        tp_commRecords = []
        tp_baseRecords = []


        # Prepare arguments for multiprocessing
        tasks = [(chrom, baseRecords) for chrom, baseRecords in self.base.items() if len(self.chroms) == 0 or chrom in self.chroms]

        with multiprocessing.Pool() as pool:
            results = pool.starmap(self.process_chromosome, tasks)

        # Collect results from all processes
        for tp_base, tp_comm, all_base, all_comm in results:
            tp_commRecords.extend(tp_comm)
            tp_baseRecords.extend(tp_base)
            all_commRecords.extend(all_comm)
            all_baseRecords.extend(all_base)


        # Calculate false positives and metrics
        fp_commRecords = list(set(all_commRecords) - set(tp_commRecords))
        fp_baseRecords = list(set(all_baseRecords) - set(tp_baseRecords))
        
        tp_comm = len(tp_commRecords)
        tp_base = len(tp_baseRecords)
        fp_comm = len(fp_commRecords) 
        fp_base = len(fp_baseRecords)
        all_base = len(all_baseRecords)
        all_comm = len(all_commRecords)
        precision = tp_comm / all_comm if all_comm else 0
        recall = tp_base / all_base if all_comm else 0
        f1 = (2 * precision * recall / (recall + precision)) if (recall + precision) else 0

        infos = (f"fp_comm: {fp_comm}\n"
                f"tp_comm: {tp_comm}\n"
                f"all_comm: {len(all_commRecords)}\n"
                f"fp_base: {fp_base}\n"
                f"tp_base: {tp_base}\n"
                f"all_base: {all_base}\n"
                f"precision: {precision}\n"
                f"recall: {recall}\n"
                f"F1: {f1}\n")

        self.outSummary(tp_commRecords, fp_commRecords, tp_baseRecords, fp_baseRecords, infos)
        

    def outSummary(self, tp_commRecords, fp_commRecords, tp_baseRecords, fp_baseRecords, infos):
        os.makedirs(self.outPath, exist_ok=True)
        with open(os.path.join(self.outPath, "outSummary.txt"), "w") as out:
            out.write(infos)
        
        with open(os.path.join(self.outPath, "tp_commRecords.txt"), "w") as out:
            for lx in tp_commRecords:
                out.write(lx)
      
        with open(os.path.join(self.outPath, "fp_commRecords.txt"), "w") as out:
            for lx in fp_commRecords:
                out.write(lx)

        with open(os.path.join(self.outPath, "tp_baseRecords.txt"), "w") as out:
            for lx in tp_baseRecords:
                out.write(lx)

        with open(os.path.join(self.outPath, "fp_baseRecords.txt"), "w") as out:
            for lx in fp_baseRecords:
                out.write(lx)


    def dealBaseComm(self, commVcf, baseVcf, base_h1, base_h2):
        if commVcf == None:
            logging.error("Please enter a commVcf file!!!")
        else:
            self.comm = self.readVcf(commVcf)

        if baseVcf == None:
            if base_h1 == None:
                logging.error("Please enter a baseVcf or base_h1 file!!!")
            else:
                self.base = self.readSimH(base_h1, base_h2)
        else:
            self.base = self.readVcf(baseVcf)
       
    def readVcf(self, vcfFile):
        data = {}
        for line in open(vcfFile, 'r'):
            if "#" in line:
                continue
            
            lx = line.strip().split('\t')
            if len(lx) <5:
                continue
            chrom = lx[0]
            pos1  = int(lx[1])
            infos = dict(zip([i.split("=")[0] for i in  lx[7].split(";") if "=" in i], [i.split("=")[1] for i in lx[7].split(";") if "=" in i]))
            print(infos, self.svtypes)
            if infos['SVTYPE'].upper() not in self.svtypes:
                continue
            pos2  = int(infos['END'])
            svlen = abs(int(infos['SVLEN']))
            svtype = infos['SVTYPE']
            gt = lx[9].split(':')[0]
            if svtype == "INS":
                pos2= pos1+svlen
            print([chrom, pos1, pos2, svlen, svtype, gt, line])
            try:
                data[chrom].append([chrom, pos1, pos2, svlen, svtype, gt, line])
            except KeyError:
                data[chrom] = [[chrom, pos1, pos2, svlen, svtype, gt, line]]
        return data

    def readSimH(self, base_h1, base_h2):
        svtypeDict = {'deletion':'DEL', 'insertion':'INS', 'inversion':'INV', 'tandem duplication':'DUP'}
        gts = {}
        data = {}
        if base_h2 != None:
            h1 = set(open(base_h1, 'r').read().strip().split('\n'))
            h2 = set(open(base_h2, 'r').read().strip().split('\n'))
            h12 = (h1-h2) | (h2-h1)
            for lx in h12:
                gts[lx] = ['0/1']

        a = 0
        for line in open(base_h1, 'r'):
            a+=1
            lx = line.strip()
            try:
                gt = gts[lx]
            except KeyError:
                gt = '1/1'
            lx = lx.split('\t')
            chrom = lx[0]
            pos1 = int(lx[1])
            pos2 = int(lx[2])
            svtype = svtypeDict[lx[3]]
            svlen = abs(pos2 - pos1)
            if svtype == "INS":
                svlen = len(lx[4])
                pos2= pos1+svlen
            if svtype not in self.svtypes:
                continue
            try:
                data[chrom].append([chrom, pos1, pos2, svlen, svtype, gt, line])
            except KeyError:
                data[chrom] = [[chrom, pos1, pos2, svlen, svtype, gt, line]]
        print("AAA:", a, sum([len(data[chrom]) for chrom in data.keys()]))
        return data
            
         

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compare two SV datasets and calculate precision,recall, F1.')
    
    # 添加命令行参数
    parser.add_argument('-b', '--base', type=str, help='')
    parser.add_argument('-c', '--comm', type=str, required=True, help='')
    parser.add_argument('-s', '--simh', nargs='+', help='')
    parser.add_argument('-r', '--rangeSize', type=int, default=1000, help='')
    parser.add_argument('-t', '--svtypes', type=str, nargs='+', default=['DEL','INV','INS','DUP'], help='')
    parser.add_argument('--chroms', type=str, nargs='+', default=[],  help='')
    parser.add_argument('-o', '--outPath', type=str, default='./out', help='')
    
    args = parser.parse_args()
    
    if len(args.simh) == 0:
        compareSV(commVcf=args.comm, baseVcf=args.base, rangeSize=args.rangeSize, chroms=args.chroms, outPath=args.outPath, svtypes=args.svtypes)
    elif len(args.simh) == 1:
        compareSV(commVcf=args.comm, base_h1=args.simh[0], rangeSize=args.rangeSize, chroms=args.chroms, outPath=args.outPath, svtypes=args.svtypes)
    else:
        compareSV(commVcf=args.comm, base_h1=args.simh[0], base_h2=args.simh[1], rangeSize=args.rangeSize, chroms=args.chroms, outPath=args.outPath, svtypes=args.svtypes)
