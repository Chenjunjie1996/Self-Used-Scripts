import pandas as pd
import sys
import os
import glob
from celescope.tools import utils

def parse_file(filter_contig, transc):
    filter_contig = glob.glob(f'{filter_contig}/06.summarize/filtered_contig_annotations.csv')[0]
    filter_contig = pd.read_csv(filter_contig)
    cell_nums = len(set(filter_contig['barcode'].tolist()))

    match_cell_barcodes, _match_cell_number = utils.read_barcode_file(transc)

    sgr_cbs = set(filter_contig['barcode'].tolist())
    df_sgr = pd.DataFrame(match_cell_barcodes, columns=['barcode'])
    df_match = pd.merge(df_sgr, filter_contig, on='barcode', how='inner')
    match_cells = len(set(df_match['barcode'].tolist()))

    return match_cells

if __name__ == '__main__':
    filter_contig = sys.argv[1]
    transc = sys.argv[2]
    match_cells = parse_file(filter_contig,transc)
    sys.stdout.write('Match Cell Number : ' + str(match_cells) +'\n')
