#!/usr/bin/env python2.7

import argparse,sys
import time
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--reports_list", default="${reports_list}", help="")
parser.add_argument("--combine_report", default="${combine_report}", help="")
parser.add_argument("--combine_report_accum", default="${combine_report_accum}", help="")
parser.add_argument("--combine_report_summary", default="${combine_report_summary}", help="")
parser.add_argument("--title", default="${title}", help="")
args = parser.parse_args()

def combine_report_sample(report_list, combine_report, combine_report_accum, combine_report_summary, title="Type"):
    """
    :param report_list:
    :param combine_report:
    :return:
    """
    samples = ['Region', title]
    datas = {}
    datas_ = []
    for reports in report_list.split(','):
        reports = reports.split(':')
        sample = reports[0]
        report = reports[1]
        samples.append(sample)
        for line in open(report):
            if "Region" not in line:
                line = [it.strip() for it in line.strip().split(';')]
                label = line[0]
                data = line[1]
                if label not in datas_:
                    datas_.append(label)
                if label not in datas:
                    datas[label] = [data]
                else:
                    datas[label].append(data)

    output = open(combine_report, 'w')
    output.writelines(";".join(samples)+"\\n")
    output_accum = open(combine_report_accum, 'w')
    output_accum.writelines(";".join([samples[0]]+samples[2:])+"\\n")
    output_summary = open(combine_report_summary, 'w')
    output_summary.writelines(";".join(samples[0:2])+"\\n")
    for data in datas_:
        if data in datas:
            mean_val = round(float(np.mean([int(it) for it in datas[data]])),1)
            output.writelines(";".join([data, str(mean_val)]+datas[data])+"\\n")
            output_summary.writelines(";".join([data, str(mean_val)])+"\\n")
            region = '_'.join(data.strip().split(' '))
            if region in ["Total_Exome_SNV", "Exome_Singletons"]:
                output_accum.writelines(";".join([region]+datas[data])+"\\n")



combine_report_sample(args.reports_list, args.combine_report, args.combine_report_accum, args.combine_report_summary, args.title)