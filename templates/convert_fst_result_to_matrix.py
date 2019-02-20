#!/usr/bin/env python2.7

import argparse,sys
import time as t

parser = argparse.ArgumentParser()
parser.add_argument("--fts_2by2_input", help="")
parser.add_argument("--fst_matrix_output", help="")
args = parser.parse_args()

def convert_fts_result_matrix(fts_2by2_input, fst_matrix_output):
    '''
    :param fts_2by2_input:
    :return fst_matrix_output:
    '''
    pops = []
    datas = {}
    fst_matrix = {}
    fst_matrix_out = open(fst_matrix_output, 'w')
    for line in open(fts_2by2_input):
        line = line.strip().split()
        pops.append(line[0])
        pops.append(line[1])
        datas[line[0]+'-'+line[1]] = str(float(line[2])*100)
    pops = sorted(list(set(pops)))
    for pop1 in pops:
        if pop1 not in fst_matrix:
            fst_matrix[pop1] = []
        for pop2 in pops:
            if pop1+'-'+pop2 in datas:
                fst_matrix[pop1].append(datas[pop1+'-'+pop2])
            elif pop2+'-'+pop1 in datas:
                fst_matrix[pop1].append(datas[pop2+'-'+pop1])
            else:
                fst_matrix[pop1].append('0')
    fst_matrix_out.writelines(','.join(['pop']+pops)+'\n')
    for pop in pops:
        fst_matrix_out.writelines(','.join([pop]+fst_matrix[pop])+'\n')
    fst_matrix_out.close()

convert_fts_result_matrix(args.fts_2by2_input, args.fst_matrix_output)