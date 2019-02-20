#!/usr/bin/env python2.7

import argparse,sys
import time

parser = argparse.ArgumentParser()
parser.add_argument("--fst_input", help="")
parser.add_argument("--fst_output", help="")
args = parser.parse_args()

def combine_weir_fst(fst_input, fst_output):
    '''
    :param fts_2by2_input:
    :return fst_matrix_output:
    '''
    fst_output_ = open(fst_output+".csv", 'w')
    fst_output_hD = open(fst_output+"_HighDiff.csv", 'w')
    fst_output_hD_1 = open(fst_output+"_HighDiff_pos.csv", 'w')
    datas = fst_input.split(' ')
    pops = []
    for data in datas:
        fst_data = data.split('__')
        pop1 = fst_data[0]
        pop2 = fst_data[1]
        fst = fst_data[2]
        pops.append(pop1+'.'+pop2+'.Fst')
        if datas.index(data) == 0:
            snps = [dat[0]+':'+dat[1] for dat in [ dat.split('\t') for dat in open(fst).readlines()[1:]]]
            fst_data_ = {key:[] for key in snps}
        pop_fst = [ dat.strip().split('\t') for dat in open(fst).readlines()[1:]]
        for dat in pop_fst:
            try:
                fst_data_[dat[0]+':'+dat[1]].append(str(abs(float(dat[2]))))
            except:
                print  dat
                fst_data_[dat[0]+':'+dat[1]].append('nan')
    fst_output_.writelines(','.join(['CHRM:POS']+pops)+'\n')
    fst_output_hD.writelines(','.join(['CHRM:POS']+pops)+'\n')
    i = 0
    for snp in snps:
        fst_output_.writelines(','.join([snp]+fst_data_[snp])+'\n')
        if any(i >= 0.5 for i in [float(t) for t in fst_data_[snp]]):
            fst_output_hD.writelines(','.join([snp]+fst_data_[snp])+'\n')
            fst_output_hD_1.writelines('\t'.join(snp.split(':'))+'\n')
    fst_output_.close()
    fst_output_hD.close()
    fst_output_hD_1.close()



combine_weir_fst(args.fst_input, args.fst_output)