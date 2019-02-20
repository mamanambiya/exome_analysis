#!/usr/bin/env python2.7
'''

'''
import argparse,sys,time

parser = argparse.ArgumentParser()
parser.add_argument("--frq_to_filter", default="${frq_to_filter}", help="")
parser.add_argument("--dataset_to_filter", default="${dataset_to_filter}", help="")
parser.add_argument("--dataset_frq_to_filter", default="${dataset_frq_to_filter}", help="")
parser.add_argument("--frq_file", default="${frq_file}", help="")
parser.add_argument("--frq_output", default="${frq_output}", help="")

args = parser.parse_args()

def filter_frq(dataset_to_filter, frq_to_filter, dataset_frq_to_filter, frq_file, frq_output):
    """
    :param bedFile:
    :param frq_file:
    :param frq_output:
    :return:
    """
    frq_to_filter_data = {}
    out = open(frq_output, 'w')
    for line in open(frq_to_filter):
        data = line.strip().split('\\t')
        chrm_pos = data[0]
        frq = data[1]
        if chrm_pos not in frq_to_filter_data:
            frq_to_filter_data[chrm_pos] = frq
    nline = 1
    frq_file_data = []
    frq_file_data_s = {}
    for line in open(frq_file):
        data = line.replace("#",'').strip().split('\\t')
        if nline == 1:
            out.writelines(data[0]+"\\t"+data[1]+"\\t"+dataset_to_filter+"_AF\\n")
        else:
            try:
                chrm_pos = data[0]
                frq = data[1]
                if chrm_pos in frq_to_filter_data:
                    if chrm_pos not in frq_file_data:
                        if len(frq.split(',')) == 1:
                            try:
                                if float(frq) <= 0.5:
                                    frq_file_data_s[chrm_pos] = frq
                                else:
                                    frq_file_data_s[chrm_pos] = str(float(frq)-0.5)
                                frq_file_data.append(chrm_pos)
                            except:
                                print line
            except:
                print line

        nline += 1
    for chrm_pos in frq_file_data:
        out.writelines(chrm_pos+"\\t"+frq_file_data_s[chrm_pos]+"\\t"+frq_to_filter_data[chrm_pos]+"\\n")
    out.close()


filter_frq(args.dataset_to_filter, args.frq_to_filter, args.dataset_frq_to_filter, args.frq_file, args.frq_output)
