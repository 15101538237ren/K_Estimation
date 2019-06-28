# -*- coding: utf-8 -*-

import os, re
import numpy as np
import scipy.io as sio
import multiprocessing as mp
import pandas as pd

NUMBER_OF_REPLICATES = [2, 3, 3, 3] # number of replicates at 0, 1, 4, 16
TIME_POINTS = ['0h', '1h', '4h', '16h']
CHROMOSOMES = [str(i) for i in range(1, 23)]
CHROMOSOMES.extend(['X', 'Y'])

def _workers_count():
    cpu_count = 1
    try:
        cpu_count = len(os.sched_getaffinity(0))
    except AttributeError:
        cpu_count = os.cpu_count()
    return cpu_count


def dividing_data_randomly(input_fp, out_fp1, out_fp2, sep="\s+",line_end = "\n", p =0.5):
    '''
    Dividing both #methylated reads (and #unmethy reads) by binomial distribution
    :param input_fp: input bed file path
    :param out_fp1: output file path for splitted 1
    :param out_fp2: output file path for splitted 2
    :return: None
    '''
    line_counter = 0
    ltws1 = [] # data list to write for split 1
    ltws2 = [] # data list to write for split 2

    with open(input_fp, "r") as input_file:
        print("dividing %s" % input_fp)
        line = input_file.readline()
        while line:
            line_counter += 1
            if line_counter % 500000 == 0:
                print("%d lines processed" % line_counter)
            line_contents = re.split(sep, line.strip(line_end))
            try:
                chr_i, start, end, methy_reads, unmethy_reads = line_contents

                methy_reads = int(methy_reads)

                methy_1 = int(np.random.binomial(methy_reads, p, 1)) #methylated reads for split 1
                methy_2 =  methy_reads - methy_1 # methylated reads for split 2

                unmethy_reads = int(unmethy_reads)
                unmethy_1 = int(np.random.binomial(unmethy_reads, p, 1))  # unmethylated reads for split 1
                unmethy_2 = unmethy_reads - unmethy_1  # unmethylated reads for split 2

                ltws1.append(
                    '\t'.join([chr_i, str(start), str(end), str(methy_1), str(unmethy_1)]))
                ltws2.append(
                    '\t'.join([chr_i, str(start), str(end), str(methy_2), str(unmethy_2)]))
            except ValueError as e:
                pass
            line = input_file.readline()

    with open(out_fp1, "w") as out_file:
        out_file.write((line_end).join(ltws1))
        out_file.write(line_end)
    with open(out_fp2, "w") as out_file:
        out_file.write((line_end).join(ltws2))
        out_file.write(line_end)
    print('Split %s successful' % input_fp)

def dividing_data_test():
    input_fp = '../DATA/Repli_BS/OTHER_DATA/0h_rep1_3_merged.bed'
    out_put_fp1 = '../DATA/Repli_BS/0h_rep1.bed'
    out_put_fp2 = '../DATA/Repli_BS/0h_rep2.bed'
    dividing_data_randomly(input_fp, out_put_fp1, out_put_fp2)

input_dir = '../DATA/Repli_BS/TMP'
output_dir = '../DATA/Repli_BS/MATLAB_DATA'
def prepare_matlab_format_data_for_mle(comb_index, combination_indexs, sep="\s+",line_end = "\n"):

    out_subdir = os.path.join(output_dir, str(comb_index))
    if not os.path.exists(out_subdir):
        os.makedirs(out_subdir)
    for chromosome in CHROMOSOMES:
        sites = []  # the list to store cpg sites
        mat_data = []
        for t_idx, time_point in enumerate(TIME_POINTS):
            repli_name = time_point + '_' + 'rep' + str(combination_indexs[t_idx])

            input_file_path = os.path.join(input_dir, repli_name, 'chr' + chromosome + '.bed')
            print("processing comb_idx: %d %s" % (comb_index, input_file_path))

            line_idx = 0
            with open(input_file_path, "r") as input_file:
                line = input_file.readline()
                while line:
                    line_contents = re.split(sep, line.strip(line_end))
                    chr_i, start, end, methy_reads, unmethy_reads = line_contents
                    methy_reads = int(methy_reads)
                    unmethy_reads = int(unmethy_reads)
                    if t_idx == 0:
                        sites.append(int(start))
                        data_arr = np.zeros((len(TIME_POINTS), 2))
                        mat_data.append(data_arr)

                    mat_data[line_idx][t_idx][0] = methy_reads
                    mat_data[line_idx][t_idx][1] = unmethy_reads
                    line = input_file.readline()
                    line_idx += 1
        mat_data = np.array(mat_data)
        MAT_DICT = {'AllDat': mat_data, 'sites': sites}
        sio.savemat(os.path.join(out_subdir, 'chr' + chromosome + '.mat'), MAT_DICT)
def merge_all_data_for_consistency_validation(chromosome, sep="\s+",line_end = "\n"):
    merged_data_dir = os.path.join(input_dir, 'MERGED_DATA')
    if not os.path.exists(merged_data_dir):
        os.makedirs(merged_data_dir)

    # methy_merged_data_fp = os.path.join(merged_data_dir, 'methy_reads_chr' + chromosome + '.bed')
    # unmethy_merged_data_fp = os.path.join(merged_data_dir,'unmethy_reads_chr' + chromosome + '.bed')

    merged_data = None
    NCOLS = 5 # site loc, 0h reads, 1h reads, 4h reads, 16h reads
    FILE_COUNTER = 0
    for t_idx, time_point in enumerate(TIME_POINTS):
        for repli_i in range(1, NUMBER_OF_REPLICATES[t_idx] + 1):
            repli_name = time_point + '_' + 'rep' + str(repli_i)

            input_file_path = os.path.join(input_dir, repli_name, 'chr' + chromosome + '.bed')
            print("processing %s" % (input_file_path))

            line_idx = 0
            with open(input_file_path, "r") as input_file:
                FILE_COUNTER += 1
                lines = [line for line in input_file]
                if FILE_COUNTER == 1:
                    NUMBER_OF_LINES = len(lines)
                    merged_data = np.zeros((NUMBER_OF_LINES, NCOLS, 2), dtype=np.int64)
                for line in lines:
                    if line_idx  % 500000 == 0:
                        print("%d lines processed" % line_idx)
                    line_contents = re.split(sep, line.strip(line_end))
                    _, start, _, methy_reads, unmethy_reads = line_contents
                    methy_reads = int(methy_reads)
                    unmethy_reads = int(unmethy_reads)

                    if FILE_COUNTER == 1:
                        merged_data[line_idx][0][0] = int(start)
                        merged_data[line_idx][0][1] = int(start)

                    merged_data[line_idx][t_idx + 1][0] += methy_reads
                    merged_data[line_idx][t_idx + 1][1] += unmethy_reads
                    line_idx += 1

    # np.savetxt(methy_merged_data_fp, merged_data[:, :, 0], fmt="%d\t%d\t%d\t%d\t%d", delimiter='\n')
    # np.savetxt(unmethy_merged_data_fp, merged_data[:, :, 1], fmt="%d\t%d\t%d\t%d\t%d", delimiter='\n')
    mat_data = merged_data[:, 1: NCOLS ,:]
    sites = list(merged_data[:, 0, 0])
    MAT_DICT = {'AllDat': mat_data, 'sites': sites}
    sio.savemat(os.path.join(merged_data_dir, 'chr' + chromosome + '.mat'), MAT_DICT)
def wrapper_for_merged_data_generation(row):
    chromosome= str(row[0])
    merge_all_data_for_consistency_validation(chromosome)
def call_for_parallel_merged_data_generation():
    num_workers = _workers_count()
    # Initiating pool
    print('Starting call_for_parallel_merged_data_generation with', num_workers, 'workers')
    parametersDF = pd.read_csv('CHROMOSOMES.txt', index_col=0, sep='\t',
                               lineterminator='\n')
    tups = parametersDF.itertuples(name=None)
    pool = mp.Pool(processes=num_workers)
    pool.map_async(wrapper_for_merged_data_generation, tups).get()

def wrapper(row):
    comb_index= row[0]
    combination_indexs = list(row[1:5])
    prepare_matlab_format_data_for_mle(comb_index, combination_indexs)
def combine_replicates():
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    idx = 1
    combination_index_fp = os.path.join(output_dir, 'INDEX_OF_COMBINATION.txt')
    header = '\t'.join(['idx', '0h', '1h', '4h', '16h'])
    ltws = [header]
    for i in range(1, NUMBER_OF_REPLICATES[0] + 1): # 0h
        for j in range(1, NUMBER_OF_REPLICATES[1] + 1): # 1h
            for k in range(1, NUMBER_OF_REPLICATES[2] + 1): # 4h
                for l in range(1, NUMBER_OF_REPLICATES[3] + 1): #16h
                    ltws.append('\t'.join([str(item) for item in [idx, i, j, k, l]]))
                    idx += 1

    with open(combination_index_fp, "w") as combination_index_file:
        content_to_write = '\n'.join(ltws)
        combination_index_file.write(content_to_write)
        combination_index_file.write('\n')

    num_workers = _workers_count()
    # Initiating pool
    print('Starting pool with', num_workers, 'workers')
    parametersDF = pd.read_csv(os.path.join(output_dir, 'INDEX_OF_COMBINATION.txt'), index_col=0, sep='\t', lineterminator='\n')
    tups = parametersDF.itertuples(name=None)
    pool = mp.Pool(processes=num_workers)
    pool.map_async(wrapper, tups).get()
if __name__ == "__main__":
    #combine_replicates()
    call_for_parallel_merged_data_generation()
