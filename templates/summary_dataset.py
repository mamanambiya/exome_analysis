#!/usr/bin/env python2.7

import argparse,sys
import time

parser = argparse.ArgumentParser()
parser.add_argument("--vcf", help="")
parser.add_argument("--csv_report", help="")
args = parser.parse_args()

def summary(vcf, report):
    """
    :param vcf:
    :param report:
    :return:
    """
