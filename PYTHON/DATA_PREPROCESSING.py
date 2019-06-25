# -*- coding: utf-8 -*-

import os, re
import numpy as np

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

                if unmethy_1 or methy_1:
                    # if not 0 unmethy and 0 methy: useless data, then append to output data list for split 1
                    ltws1.append(
                    '\t'.join([chr_i, str(start), str(end), str(methy_1), str(unmethy_1)]))
                if unmethy_2 or methy_2:
                    ltws2.append(
                    '\t'.join([chr_i, str(start), str(end), str(methy_2), str(unmethy_2)]))
            except ValueError as e:
                pass
            line = input_file.readline()

    with open(out_fp1, "w") as out_file:
        out_file.write((line_end).join(ltws1))
        out_file.write(line_end)
    with open(out_fp2, "w") as out_file:
        out_file.write((line_end).join(ltws1))
        out_file.write(line_end)
    print('Split %s successful' % input_fp)

if __name__ == "__main__":
    input_fp = '../DATA/Repli_BS/HUES64_WT_AllCells_1hBrdU_0Chase_nascent_rep1_3_merged.bed'
    out_put_fp1 = '../DATA/Repli_BS/OTHER_DATA/HUES64_WT_AllCells_1hBrdU_0Chase_nascent_random_split_rep1.bed'
    out_put_fp2 = '../DATA/Repli_BS/OTHER_DATA/HUES64_WT_AllCells_1hBrdU_0Chase_nascent_random_split_rep3.bed'
    dividing_data_randomly(input_fp, out_put_fp1, out_put_fp2)