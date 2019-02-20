#!/usr/bin/env python2.7
'''

'''
import argparse,sys,time

parser = argparse.ArgumentParser()
parser.add_argument("--mainTSV", default="${mainTSV}", help="Main TSV annotation file to check")
parser.add_argument("--otherTSVs", default="${otherTSVs}", help="One or many TSV annotation files, if many use ';' as separator")
parser.add_argument("--private_list_file", default="${private_list_file}", help="private sites annotation file")

args = parser.parse_args()

def intersection(main_tsv, tsv_list, private_list_file):
    """
    :param main_tsv:
    :param tsv_list:
    :param private_list:
    :return:
    """
    private_out = open(private_list_file, 'w')
    known_sites = set()
    for tsv in tsv_list.split(';'):
        for line in open(tsv):
            line = line.strip().split('\\t')
            if line[0][0] != '#':
                try:
                    if float(line[1]) > 0 :
                        known_sites.add(line[0])
                except:
                    print "Not valid MAF"

    ncount = 1
    for line in open(main_tsv):
        data = line.strip().split('\\t')
        ID = data[0].split(":")
        chrm = ID[0]
        pos = ID[1]
        if chrm != '#':
            if data[0] not in known_sites:
                private_out.writelines('\\t'.join([chrm, pos, pos])+"\\n")
        ncount += 1

    private_out.close()

intersection(args.mainTSV, args.otherTSVs, args.private_list_file)