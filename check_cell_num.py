import pandas as pd
import sys
import os
import glob
from celescope.tools import utils

def parse_file(filter_contig, transc):
    """
    usage
    python /SGRNJ06/randd/USER/cjj/Script/Self-Used-Scripts/check_cell_num.py\
    /SGRNJ06/randd/PROJECT/IR_group/AccuraCode_TCR/2025sCircle_VDJ/20250717_LJ_Hum_GOTN/TCR2/GTBSA05038x2_TCR\
    /SGRNJ06/randd/PROJECT/IR_group/AccuraCode_TCR/2025sCircle_VDJ/20250717_LJ_Hum_GOTN/RNA/GTBSA05039x2_RNA
    """
    filter_contig = glob.glob(f'{filter_contig}/*/filtered_contig_annotations.csv')[0]
    filter_contig = pd.read_csv(filter_contig)
    cell_nums = len(set(filter_contig['barcode'].tolist()))

    match_cell_barcodes, _match_cell_number = utils.get_barcode_from_match_dir(transc)

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
