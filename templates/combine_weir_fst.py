#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--weir_fst_files", default="${weir_fst_files}", help="")
parser.add_argument("--weir_fst_out", default="${weir_fst_out}", help="")
parser.add_argument("--fst_cutoff", default="${fst_cutoff}", help="")
args = parser.parse_args()


def combine_weir_fst(weir_fst_files, weir_fst_out, fst_cutoff=0.5):
    '''
    :param fts_2by2_input:
    :return fst_matrix_output:
    '''
    weirs = weir_fst_files.split(',')
    fst_cutoff = float(fst_cutoff)
    base = pd.read_csv(weirs[0], sep='\\s+', engine='python')
    base = base[base.iloc[:, 2] >= fst_cutoff].reset_index().drop(columns=[
        'index'])
    base = base[['CHROM', 'POS']]
    for weir in weirs[1:]:
        data = pd.read_csv(weir, sep='\\s+', engine='python')
        data = data[data.iloc[:, 2] >= fst_cutoff].reset_index().drop(columns=[
            'index'])
        data = data[['CHROM', 'POS']]
        base = pd.concat([base, data])

    base = base.groupby(by=['CHROM', 'POS']).size(
    ).reset_index().rename(columns={0: "Count"})
    base = base[base['Count'] >= 2].reset_index().drop(
        columns=['index']).sort_values(['CHROM', 'POS'])
    base.reset_index().drop(columns=['index']).drop_duplicates().to_csv(f'{weir_fst_out}.csv', index=False,
                                                                        header=True, sep='\\t')


combine_weir_fst(args.weir_fst_files, args.weir_fst_out, args.fst_cutoff)
