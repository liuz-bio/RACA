import multiprocessing
import pickle
import os
import time
import traceback

import numpy as np
from dataclasses import dataclass
from sklearn.cluster import DBSCAN
from intervaltree import IntervalTree#, Interval
from clusterSig import readSig, clusterSig
from callINV import callINV
from callDUP import callDUP
from callINS import callINS
#from callINS import callINS
from callDEL import callDEL
import gc
from loguru import logger


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




class callSV:
    def __init__(self, configer, sample, chroms_plks, all_cand_svs):
        self.configer = configer
        print("callSV:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
        self.run(sample, chroms_plks, all_cand_svs)

    
    def write_vcf(self):
        print("write_vcf DUPs", self.out_file)
        with open(self.out_file, 'a') as out:
            out.write('\n'.join(['\t'.join(map(str, i)) for i in self.all_DUPs])+'\n')
            out.flush()
    
    def pass_overlap(self):
        num = len(self.all_DUPs)
        DUPs = [i[1] for i in self.all_DUPs]
        print("pass_overlap3:", len(DUPs))
        for i in range(num-1):
            a = self.all_DUPs[i]
            for ii in range(i+1,num):
                b = self.all_DUPs[ii]
                if a[0][0] == b[0][0]:
                    #print("pass_overlapdddd",a,b)
                    if abs(b[0][1]-a[0][1]) >100 or abs(a[0][2]-b[0][2]) >100 :
                        continue
                    elif (a[0][1] <= b[0][1] and a[0][2]>=b[0][2]) and (abs(a[0][1] - b[0][1]) < 100 and  abs(a[0][2]-b[0][2]) <100):
                        try:
                            print("pass_overlap1:",a[1], b[1] )
                            DUPs.remove(b[1])
                        except ValueError:
                            continue
                    elif b[0][1] <= a[0][1] and b[0][2]>=a[0][2] and (abs(a[0][1] - b[0][1]) < 100 and  abs(a[0][2]-b[0][2]) <100):
                        try:
                            print("pass_overlap2:",a[1], b[1] )
                            DUPs.remove(a[1])
                        except ValueError:
                            continue
                    elif abs(b[0][1]-a[0][1]) <100 or abs(a[0][2]-b[0][2]) <100:
                        if b[0][3] > a[0][3]:
                            try:
                                print("pass_overlap22:",a[1], b[1] )
                                DUPs.remove(a[1])
                            except ValueError:
                                continue
                        else:
                            try:
                                print("pass_overlap11:",a[1], b[1] )
                                DUPs.remove(b[1])
                            except ValueError:
                                continue
        self.all_DUPs = DUPs
    
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
                #print("cluster test:", [chrom, pos1, pos2, svType, len(cluster_points)])
            else:
                # 优化点：避免for循环，用列表解析式代替
                for i  in cluster_points:
                    subSigTree.insert(i[0], i[1])
                    #subSigTree.addi(i[0], i[1])
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

    #############################################################################################
    
    def candidateSV(self, sample):
        candSV = {}
        SigTrees = {}
        mun = 0
        for i in ["DEL", "DUP", "INS", "INV", "BND"]:
            candSV[i] = []
            SigTrees[i] = {}
            with open(os.path.join(self.configer.tmp_path, sample, sample + '.' + i + '.sigs'), 'r') as f:
                for line in f:
                    line = line .strip().split('\t')
                    start, end = int(line[1]), int(line[2])
                    if start == end:
                        end += 1
                    elif start >end:
                        start, end = int(line[2]), int(line[1]) #BND 临时 需要修改candaidateSVLocation.py
                        
                    if line[0] in SigTrees[i]:
                        SigTrees[i][line[0]] = IntervalTree()
                        SigTrees[i][line[0]].insert(start, end, int(line[4])
                    else:
                        SigTrees[i][line[0]].insert(start, end, int(line[4])
                    #SigTrees[i].setdefault(line[0], IntervalTree([start, end, int(line[4]))).addi(start, end, int(line[4]))
                    candSV[i].append([line[0], start, end])
                    if i in ["DUP", "INV"]:
                        mun += 1

        print("LLLLL", mun)
        return candSV, int(mun/self.configer.mut_proces)+1, SigTrees
    
    def subSigRegion(self, sample):
        subSigTrees = {}
        for sigT in ['DEL', 'DUP', 'INV', 'INS', 'BND']:
            subSigTrees[sigT] = {}
            sub_file = str(os.path.join(self.configer.tmp_path, sample, sample + ".sub.%s.sigs" % sigT))
            chrom_lines = {}
            for i in open(sub_file, 'r'):
                i = i.strip().split('\t')
                chrom_lines.setdefault(i[0], []).append([int(i[1]), int(i[2])])  # INV

            dbscan = DBSCAN(eps=500, min_samples=2)
            out = open('ttttt.%s.out'%sigT, 'w')
            for chrom, lines in chrom_lines.items():
                subSigTree = IntervalTree()
                groups = self.subMerge(lines)
                self.splitChromSub(chrom, groups, subSigTree, dbscan, sigT)
                subSigTrees[sigT][chrom] = subSigTree
                for i in groups:
                    out.write('\t'.join([chrom] + list(map(str, i))) + "\n")
            out.close()

        return subSigTrees
    
    def handle_error(self, e):
        print(''.join(traceback.format_exception(type(e), e, e.__traceback__)))
        print("An error occurred in the child process:", e)
    
    def call_del(self, *args ):
        try:
            callDEL(*args).run()
        except Exception as e:
            self.handle_error(e)
    
    def call_ins(self, *args ):
        try:
            callINS(*args).run()
        except Exception as e:
            self.handle_error(e)
    
    def call_inv(self, *args):
        print("call_inv")
        try:
            callINV(*args).run()
        except Exception as e:
            self.handle_error(e)
    
    def call_dup(self, *args):
        try:
            callDUP(*args).run()
        except Exception as e:
            self.handle_error(e)
    
    def run(self, sample, chroms_plks, all_cand_svs):
        #out = open('test.vcf', 'w')
        start_time = time.time()
        self.out_file = str(os.path.join(self.configer.out_path, sample+".vcf"))
        if os.path.exists(self.out_file):
            os.remove(self.out_file)
        candSV, p_size, SigTrees = self.candidateSV(sample)
        subSigTrees = self.subSigRegion(sample)
        duplicate_records = pickle.load(
            open(str(os.path.join(self.configer.tmp_path, sample, sample + ".duplicate_records.plk")), 'rb'))

        logger.info(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()) + ' ' + "duplicate_records deal" + ': ' + ' finish...')

        #process_pool = multiprocessing.Pool(self.configer.mut_proces)


        shared_sig_trees = [SigTrees, subSigTrees]
        process_pool = multiprocessing.Pool(self.configer.mut_proces)
        manager = multiprocessing.Manager()
        all_DUPs = manager.list()
        lock = manager.Lock()
        print("shared_sig_trees")
        print("Time1:", time.time()-start_time)
        start_time = time.time()
        tasks_to_submit = []
        p_size = 100
        for svT, candPos in candSV.items():
            candSV_size = len(candPos)

            if svT == "INV":
                if candSV_size < 2*p_size:
                    mut_pos = candPos
                    #print("INV", sample, candSV_size)
                    tasks_to_submit.append((self.call_inv, (self.out_file, self.configer, sample, mut_pos, shared_sig_trees, lock)))
                    # process_pool.apply_async(func=self.call_inv,
                    #                             args=(out_file, self.configer, sample, mut_pos, shared_sig_trees, lock,),
                    #                             error_callback=self.handle_error)
                elif candSV_size % p_size == 0:
                    gs = int(candSV_size / p_size)
                else:
                    gs = int(candSV_size / p_size) + 1

                for i in range(gs):
                    mut_pos = candPos[i*p_size:(i+1)*p_size]
                    #print("INV", sample, len(mut_pos))
                    tasks_to_submit.append((self.call_inv, (self.out_file, self.configer, sample, mut_pos, shared_sig_trees, lock)))
                    # process_pool.apply_async(func=self.call_inv,
                    #                         args=(out_file, self.configer, sample, mut_pos, shared_sig_trees, lock,),
                    #                         error_callback=self.handle_error)

            elif svT == "DUP":
                if candSV_size < 2*p_size:
                    mut_pos = candPos
                    #print("DUP", sample, candSV_size)
                    tasks_to_submit.append((self.call_dup, (self.out_file, self.configer, sample, mut_pos, shared_sig_trees, duplicate_records, lock)))
                    # process_pool.apply_async(func=self.call_dup,
                    #                         args=(out_file, self.configer, sample, mut_pos, shared_sig_trees, duplicate_records, lock, ),
                    #                         error_callback=self.handle_error)
                elif candSV_size % p_size == 0:
                    gs = int(candSV_size / p_size)
                else:
                    gs = int(candSV_size / p_size) + 1

                for i in range(gs):
                        mut_pos = candPos[i*p_size:(i+1)*p_size]
                        #print("DUP", sample, len(mut_pos))
                        tasks_to_submit.append((self.call_dup, (self.out_file, self.configer, sample, mut_pos, shared_sig_trees, duplicate_records, lock)))
                        # process_pool.apply_async(func=self.call_dup,
                        #                     args=(out_file, self.configer, sample, mut_pos, shared_sig_trees, duplicate_records, lock, ),
                        #                     error_callback=self.handle_error)

        logger.info(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()) + ' ' + "deal DUP INV group" + ': ' + ' finish...')
        for func, args in tasks_to_submit:
            process_pool.apply_async(func=func, args=args, error_callback=self.handle_error)
        process_pool.close()
        process_pool.join()

        #self.all_DUPs = all_DUPs
        #logger.info(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()) + ' ' + "deal DUP INV" + ': ' + ' finish...|'+str(len(all_DUPs)))

        #self.pass_overlap()
        #self.write_vcf()

        print("Time2:", time.time()-start_time)
        start_time = time.time()
        #logger.info(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()) + ' ' + "deal DUP INV" + ': ' + ' finish...|'+str(len(self.all_DUPs)))
        logger.info(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()) + ' ' + "deal DUP INV" + ': ' + ' finish...')
        #######################################################################################
        all_size_DEL_INS_svs = sum([ii for i in all_cand_svs.values() for ii in i.values()])
        print('all_size_DEL_INS_svs:', all_size_DEL_INS_svs)
        task_size = 100#int(all_size_DEL_INS_svs/self.configer.mut_proces)*3+1

        data_inv_dup = {'DUP':{},'INV':{}}
        for i in open(self.out_file, 'r'):
            i = i.strip().split('\t')
            if len(i)<8:
                continue
            lx = [ii for ii in i[7].split(';') if "=" in ii]
            infos = dict(zip([ii.split("=")[0] for ii in lx], [ii.split("=")[1] for ii in lx]))
            if infos["SVTYPE"] in ["DUP", "INV"]:
                try:
                    if i[0] in data_inv_dup[infos["SVTYPE"]]:
                        data_inv_dup[infos["SVTYPE"]][i[0]] = IntervalTree()
                        data_inv_dup[infos["SVTYPE"]][i[0]].insert(int(i[1]), int(infos["END"]))
                    else:
                        data_inv_dup[infos["SVTYPE"]][i[0]].insert(int(i[1]), int(infos["END"]))
                    #data_inv_dup[infos["SVTYPE"]].setdefault(i[0], IntervalTree([Interval(int(i[1]), int(infos["END"]))])).addi(int(i[1]), int(infos["END"]))
                except ValueError:
                    continue
        data_inv_dup = {'DUP':{},'INV':{}}
        print("Time3:", time.time()-start_time)
        start_time = time.time()


        process_pool = multiprocessing.Pool(self.configer.mut_proces)
        manager = multiprocessing.Manager()
        lock = manager.Lock()
        tmp_task_size = {'D': 0, 'I': 0}
        tmp_task = {'D': [], 'I': []}
        all_size = 0

        tasks_to_submit = []

        for plk_file, plk_data in all_cand_svs.items():
            plk_file = str(os.path.join(self.configer.tmp_path, sample, plk_file))

            for sigT, plk_size in plk_data.items():
                if tmp_task_size[sigT] + plk_size < task_size:
                    tmp_task[sigT].append([plk_file, 0, plk_size])
                    tmp_task_size[sigT] += plk_size
                elif tmp_task_size[sigT] >= task_size:
                    print('DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD')
                else:
                    nu = task_size - tmp_task_size[sigT]
                    tmp_task[sigT].append([plk_file, 0, nu])
                    all_size += sum([abs(i[1] - i[2]) for i in tmp_task[sigT]])
                    if sigT == 'D':
                        tasks_to_submit.append((self.call_del, (self.out_file, self.configer, sample, tmp_task[sigT], shared_sig_trees, data_inv_dup, lock)))
                    elif sigT == 'I':
                        tasks_to_submit.append((self.call_ins, (self.out_file, self.configer, sample, tmp_task[sigT], shared_sig_trees, data_inv_dup, lock)))
                    surplus = plk_size - nu
                    if surplus == 0:
                        tmp_task[sigT] = []
                        tmp_task_size[sigT] = 0
                    elif surplus < task_size:
                        tmp_task[sigT] = [[plk_file, nu, plk_size]]
                        tmp_task_size[sigT] = surplus
                    elif surplus >= task_size:
                        for i in range(int(surplus / task_size)):
                            all_size += task_size
                            task = [[plk_file, i * task_size + nu, (i + 1) * task_size + nu]]
                            if sigT == 'D':
                                tasks_to_submit.append((self.call_del, (self.out_file, self.configer, sample, task, shared_sig_trees, data_inv_dup, lock)))
                            else:
                                tasks_to_submit.append((self.call_ins, (self.out_file, self.configer, sample, task, shared_sig_trees, data_inv_dup, lock)))
                        tail = surplus % task_size
                        if tail > 0:
                            tmp_task[sigT] = [[plk_file, plk_size - tail, plk_size]]
                            tmp_task_size[sigT] = tail
                        else:
                            tmp_task[sigT] = []
                            tmp_task_size[sigT] = 0

        if len(tmp_task['D']) > 0:
            all_size += sum([abs(i[1] - i[2]) for i in tmp_task['D']])
            tasks_to_submit.append((self.call_del, (self.out_file, self.configer, sample, tmp_task['D'], shared_sig_trees, data_inv_dup, lock)))

        if len(tmp_task['I']) > 0:
            all_size += sum([abs(i[1] - i[2]) for i in tmp_task['I']])
            tasks_to_submit.append((self.call_ins, (self.out_file, self.configer, sample, tmp_task['I'], shared_sig_trees, data_inv_dup, lock)))
        logger.info(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()) + ' ' + "deal DEL INS group" + ': ' + ' finish...')
        # Batch submit tasks
        for func, args in tasks_to_submit:
            process_pool.apply_async(func=func, args=args, error_callback=self.handle_error)

        process_pool.close()
        process_pool.join()
        print("all_size:", all_size)
        print("Time4:", time.time()-start_time)
        start_time = time.time()
        logger.info(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()) + ' ' + "deal DEL INS" + ': ' + ' finish...')
        
        # for sigT, data1 in SigTrees.items():
        #     if sigT not in ['DUP','INV']:
        #         continue
        #     for chrom, tree in data1.items():
        #         for node in tree:
        #             data_inv_dup[sigT].setdefault(chrom, IntervalTree()).add(node)

        #print("data_inv_dup", data_inv_dup)
        ##########################################################
        # process_pool = multiprocessing.Pool(self.configer.mut_proces)
        # manager = multiprocessing.Manager()
        # lock = manager.Lock()

        # for _, plk_files in chroms_plks.items():
        #     step = 30
        #     nu = 0
        #     print("plk_files_process_pool_len", len(plk_files))
        #     while True:
        #         part_plk_files = plk_files[nu*step:(nu+1)*step]
        #         print("plk_files_process_pool", nu*step, (nu+1)*step, len(part_plk_files))
        #         if len(part_plk_files) >0:
        #             nu += 1
        #             # process_pool.apply_async(func=self.call_del,
        #             #                          args=(out_file, self.configer, sample, part_plk_files, shared_sig_trees, data_inv_dup, lock, ),
        #             #                          error_callback=self.handle_error)

        #             process_pool.apply_async(func=self.call_ins,
        #                                      args=(out_file, self.configer, sample, part_plk_files, shared_sig_trees, data_inv_dup, lock, ),
        #                                      error_callback=self.handle_error)
        #         else:
        #             break
        # process_pool.close()
        # process_pool.join()



        # process_pool = multiprocessing.Pool(self.configer.mut_proces)
        # manager = multiprocessing.Manager()
        # lock = manager.Lock()
        # tmp_task_size = {'D':0,'I':0}
        # tmp_task = {'D':[], 'I':[]}
        # all_size = 0
        # for plk_file, plk_data in all_cand_svs.items():
        #     plk_file = str(os.path.join(self.configer.tmp_path, sample, plk_file))

        #     for sigT, plk_size in plk_data.items():
        #         if tmp_task_size[sigT] + plk_size < task_size:
        #             tmp_task[sigT].append([plk_file, 0, plk_size])
        #             tmp_task_size[sigT] += plk_size
        #         elif tmp_task_size[sigT] >= task_size:
        #             print('DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD')
        #         else:
        #             nu = task_size - tmp_task_size[sigT]
        #             tmp_task[sigT].append([plk_file, 0, nu])
        #             if sigT == 'D':
        #                 all_size += sum([abs(i[1]-i[2]) for i in tmp_task[sigT]])
        #                 process_pool.apply_async(func=self.call_del,
        #                                            args=(out_file, self.configer, sample, tmp_task[sigT], shared_sig_trees, data_inv_dup, lock, ),
        #                                            error_callback=self.handle_error)
        #             elif sigT == 'I':
        #                 all_size += sum([abs(i[1]-i[2]) for i in tmp_task[sigT]])
        #                 process_pool.apply_async(func=self.call_ins,
        #                                             args=(out_file, self.configer, sample, tmp_task[sigT], shared_sig_trees, data_inv_dup, lock, ),
        #                                             error_callback=self.handle_error)
        #             surplus = plk_size - nu
        #             if surplus == 0:
        #                 tmp_task[sigT] = []
        #                 tmp_task_size[sigT] = 0
        #             elif surplus < task_size:
        #                 tmp_task[sigT] = [[plk_file, nu, plk_size]]
        #                 tmp_task_size[sigT] = surplus
        #             elif surplus >= task_size:
        #                 if sigT == 'D':
        #                     for i in range(int(surplus/task_size)):
        #                        all_size+=task_size
        #                        task = [[plk_file, i*task_size+nu, (i+1)*task_size+nu]]
        #                        process_pool.apply_async(func=self.call_del,
        #                                                 args=(out_file, self.configer, sample, task, shared_sig_trees, data_inv_dup, lock, ),
        #                                                 error_callback=self.handle_error)

        #                 else:
        #                     for i in range(int(surplus/task_size)):
        #                         all_size+=task_size
        #                         task = [[plk_file, i*task_size+nu, (i+1)*task_size+nu]]
        #                         process_pool.apply_async(func=self.call_ins,
        #                                                  args=(out_file, self.configer, sample, task, shared_sig_trees, data_inv_dup, lock, ),
        #                                                  error_callback=self.handle_error)
        #                 tail = surplus % task_size
        #                 if tail > 0:
        #                     tmp_task[sigT] = [[plk_file, plk_size-tail, plk_size]]
        #                     tmp_task_size[sigT] = tail
        #                 else:
        #                     tmp_task[sigT] = []
        #                     tmp_task_size[sigT] = 0
                    

        # if len(tmp_task['D']) >0:
        #     all_size += sum([abs(i[1]-i[2]) for i in tmp_task['D']])
        #     process_pool.apply_async(func=self.call_del,
        #                                                args=(out_file, self.configer, sample, tmp_task['D'], shared_sig_trees, data_inv_dup, lock, ),
        #                                                error_callback=self.handle_error)
            
        # if len(tmp_task['I']) >0:
        #     all_size += sum([abs(i[1]-i[2]) for i in tmp_task['I']])
        #     process_pool.apply_async(func=self.call_ins,
        #                                 args=(out_file, self.configer, sample, tmp_task['I'], shared_sig_trees, data_inv_dup, lock, ),
        #                                 error_callback=self.handle_error)
            
        # process_pool.close()
        # process_pool.join()
        # print("all_size:", all_size)




        


if __name__ == '__main__':
    #bamFiles = ["/public3/SVDataset_lz/backup/2023_11_14/Sim/bed/working/Pipelines/CCS/Sim.DEL.INS.DUP.INV.25x/winnowmap/bam/Sim.DEL.INS.DUP.INV.25x.sort.bam"]
    #sampleIDs = ["SIM"]
    genome = "/home/lz/work_space/Database/VISOR_test/Genome_old/Genome.fa"
    bamFiles = ["/home/lz/work_space/Database/VISOR_test/Chr1_two_new/T2/FY0/VISOR_LASeR_Father_DEL.DUP.INS.INV/sim.srt.bam"]
    sample = "FY0_Father"
    sampleIDs = ["FY0_Father"]
    sv_types = ["DEL", "INS", "INV", "DUP", "BND"]
    cluster_param = {'L': {'min_distanc': 150, 'min_signals': 2},
                     'R': {'min_distanc': 150, 'min_signals': 2},
                     'D': {'min_distanc': 300, 'min_signals': 2},
                     'I': {'min_distanc': 300, 'min_signals': 2},
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
                             min_map_quality=20,
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
                             )


    plks = os.popen("ls /home/lz/work_space/project/plotSV/RACA/plotSV9/callSV/tmp/tmp/FY0_Father/*.plk|grep -v duplicate_records").read().strip().split("\n")
    chroms_plks = {}
    for i in plks:
        ii = i.strip().split("/")
        chrom_tmp = ii[-1].split(".")[0]
        chroms_plks.setdefault(chrom_tmp,[]).append(i)

    with open('chroms_plk.all_cand_svs.plk', 'rb') as f:
        chroms_plk = pickle.load(f)
        all_cand_svs = pickle.load(f)

    inv = callSV(configer, sample, chroms_plk, all_cand_svs)
