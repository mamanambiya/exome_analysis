#!/usr/bin/env python2.7

import argparse,sys
import time
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="AIBST", help="")
parser.add_argument("--reports", default="${reports}", help="")
parser.add_argument("--combine_report", default="AIBST_variants_report.csv", help="")
args = parser.parse_args()

def combine_report(dataset, reports, combine_report):
    """
    :param dataset:
    :param reports
    :param combine_report:
    :return:
    """
    reports = reports.split(',')
    datas = pd.read_csv(reports[0], sep=";")
    for report in reports[1:]:
        data = pd.read_csv(report, sep=";", usecols=[1])
        head = list(data.columns.values)[0]
        datas.insert(loc=datas.shape[1], column=head, value=data)
    datas.to_csv(combine_report)

combine_report(args.dataset, args.reports, args.combine_report)