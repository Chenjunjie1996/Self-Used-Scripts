import os
import subprocess
import glob
import pandas as pd
import sys
from concurrent.futures import ProcessPoolExecutor
import xlwt
import os
import re
import celescope.tools.utils as utils


# HumanPrimaryCellAtlasData and MouseRNAseqData dataset
CELL_TYPE_DICT = {
    'TCR':['TCELLS', 'TCELL' 'NKTCELLS'],
    'BCR':['MATUREBCELL', 'PLASMACELLS', 'BCELLS', 'BCELL', 'PREBCELLCD34'],
    }
ROOT_DIR = os.path.dirname(__file__)
OUT_PATH = 'YingShe'


class Yinshe:
    def __init__(self, sample, rds, VDJ_file, outdir) -> None:
        self.sample = sample
        self.rds = rds
        self.VDJ_file = VDJ_file
        self.outdir = outdir
    
    @utils.add_log
    def run_mapping(self):
        cmd = (
            f'Rscript {ROOT_DIR}/yingshe.R '
            f'--sample {self.sample} '
            f'--rds {self.rds} '
            f'--VDJ {self.VDJ_file} '
            f'--outdir {self.outdir} '
        )
        subprocess.check_call(cmd, shell=True)
    
    def run(self):
        self.run_mapping()


def parse_mapfile(mapfile):

    df_mapfile = pd.read_csv(mapfile, sep='\t')
    path_list = df_mapfile['path']
    rds_list = df_mapfile['rds']
    sample_list = df_mapfile['sample']

    return path_list, rds_list, sample_list


def mapping(path, rds, sample, celltype):
    outdir = f'{OUT_PATH}/{sample}_{celltype}'
    utils.check_mkdir(outdir)
    
    if os.path.exists(f'{path}/05.match'): # cellranger
        try:
            VDJ_file = glob.glob(f'{path}/05.match/match_contigs.csv')[0] # 3.0, 4.0, 6.0
        except IndexError:
            VDJ_file = glob.glob(f'{path}/05.match/matched_contig_annotations.csv')[0] # 7.0
    
    elif os.path.exists(f'{path}/04.summarize'): # trust4
        VDJ_file = glob.glob(f'{path}/04.summarize/*_filtered_contig.csv')[0]
    
    elif os.path.exists(f'{path}/05.count_vdj'): # vdj-cdr3
        VDJ_file = glob.glob(f'{path}/05.count_vdj/*_cell_confident_count.tsv')[0]

    if os.path.exists(VDJ_file):
        runner = Yinshe(sample, rds, VDJ_file, outdir)
        runner.run()
    else:
        raise FileNotFoundError('please check vdj file')


def get_seqtype(path_list):
    run_shell = glob.glob(f'{path_list[0]}/../*.mapfile')[0]
    with open(run_shell) as f:
        for line in f.readlines():
            if 'BCR' in line:
                return 'BCR'
            return 'TCR'


def run_count(celltype):
    meta_list = glob.glob(f"./{OUT_PATH}/*_{celltype}/meta.csv")
    sample_list = [meta.split('/')[-2] for meta in meta_list]
    mapping_cell_type = CELL_TYPE_DICT[celltype]

    report = xlwt.Workbook()
    sheet = report.add_sheet('mapping count')

    row0 = ['number of mapping to T/B cells', 'number of T/B cells', 'percent']
    row1 = sample_list

    mapping_count, ZL_count = [], []
    for meta in meta_list:
        df = pd.read_csv(meta)

        if 'new_ident' in df.columns: # manual assign column name
            ident = 'new_ident'
        elif 'celltype' in df.columns: # singleR column name
            ident = 'celltype'

        df[ident] = df[ident].apply(lambda x: re.sub(r"[^a-zA-Z0-9]","", x).upper())
        df_count = df[df[ident].isin(mapping_cell_type)]

        mapping_count.append(int(df_count['Class'].value_counts()))
        ZL_count.append(int(df_count[ident].value_counts()))

    percent = [x / y for x, y in zip(mapping_count, ZL_count)]
    percent = [str(round(i,4)*100) + "%" for i in percent]

    for _i in range(len(row0)):
        sheet.write(0, _i+1, row0[_i])
    for _i in range(len(row1)):
        sheet.write(_i+1, 0, row1[_i])
        sheet.write(_i+1, 1, mapping_count[_i])
        sheet.write(_i+1, 2, ZL_count[_i])
        sheet.write(_i+1, 3, percent[_i])

    report.save(f'./{OUT_PATH}/mapping_count_{celltype}.xls')
    return report


def main():
    os.system(f'mkdir {OUT_PATH}')
    mapfile = sys.argv[1]
    path_list, rds_list, sample_list = parse_mapfile(mapfile)
    celltype = get_seqtype(path_list)

    with ProcessPoolExecutor(max_workers=10) as executor:
        for result in executor.map(mapping, path_list, rds_list, sample_list, [celltype]*len(path_list)):
            print(result, 'done')
    
    run_count(celltype)


if __name__ == '__main__':
    main()

