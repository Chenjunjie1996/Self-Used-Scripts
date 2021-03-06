import os
import subprocess
import glob
import pandas as pd
import sys
import argparse
from concurrent.futures import ProcessPoolExecutor
import xlwt

def run_mapping(rds, contig, sample, outdir):
    cmd = (
        f'Rscript /SGRNJ03/randd/cjj/Script/mappingManual/mappingManual.R '
        f'--rds {rds} '
        f'--VDJ {contig} '
        f'--sample {sample} '
        f'--outdir {outdir} '
    )
    subprocess.check_call(cmd, shell=True)


def parse_file(rds,contig_file):
    shell_file = '/'.join(os.path.abspath(rds).split('/')[:-2]) + "/" + "run.multSaps3.sh"
    rds_sample = []
    with open(shell_file, 'r') as shell:
        lines = (line.strip().strip('--ugname').strip().strip("\\").split(',') for line in shell if
                 line.strip().startswith('--ugname'))
        for line in lines:
            rds_sample = line

    contigs = contig_file.strip().split(",")
    contig_list = []
    for contig in contigs:
        if os.path.exists(f'{contig}/05.match/match_contigs.csv'):
            contig_file = glob.glob(f'{contig}/05.match/match_contigs.csv')[0]
            contig_list.append(contig_file)
        elif os.path.exists(f'{contig}/04.summarize'):
            contig_file = glob.glob(f'{contig}/04.summarize/*_filtered_contig.csv')[0]
            contig_list.append(contig_file)
        else:
            print("check contig file path")
            raise FileNotFoundError

    return rds_sample, contig_list

def parse_meta(mapping_cell_type):
    meta_list = glob.glob("./mappingManual/*_meta.csv")
    sample_list = ['_'.join(i.split('/')[-1].split('_')[:-1]) for i in meta_list]

    report = xlwt.Workbook()
    sheet = report.add_sheet('mapping count')

    row0 = ['number of mapping to T/B cells', 'number of T/B cells', 'percent']
    row1 = sample_list

    mapping_count, ZL_count = [], []
    for meta in meta_list:
        df = pd.read_csv(meta)
        try:
            df['new_ident'] = df['new_ident'].apply(lambda x: str(x).replace(" ", ""))
        except:
            pass
        df_count = df[df['new_ident'] == f'{mapping_cell_type}']
        _i,_j = df_count['Class'].value_counts(), df_count['new_ident'].value_counts()
        mapping_count.append(int(_i))
        ZL_count.append(int(_j))
    percent = [x / y for x, y in zip(mapping_count, ZL_count)]
    percent = [str(round(i,4)*100) + "%" for i in percent]

    for _i in range(len(row0)):
        sheet.write(0, _i+1, row0[_i])
    for _i in range(len(row1)):
        sheet.write(_i+1, 0, row1[_i])
        sheet.write(_i+1, 1, mapping_count[_i])
        sheet.write(_i+1, 2, ZL_count[_i])
        sheet.write(_i+1, 3, percent[_i])

    report.save('./mappingManual/mapping_count.xls')
    return report

def main():
    parser = argparse.ArgumentParser(description='contigfile')
    parser.add_argument('--contig_file', help='contig file', required=True)
    parser.add_argument('--rds', help='rds file', required=True)
    parser.add_argument('--mapping_cell_type', help='Tcells or Bcells', required=True)
    args = parser.parse_args()
    rds = args.rds
    outdir = "./mappingManual"
    os.system(f"mkdir {outdir}")
    contig_file = args.contig_file
    rds_sample, contig_list = parse_file(rds, contig_file)
    rds_list = [rds]* len(rds_sample)
    outdir_list = [outdir] * len(rds_sample)

    with ProcessPoolExecutor(max_workers=8) as executor:
        for result in executor.map(run_mapping, rds_list, contig_list, rds_sample, outdir_list):
            print(result, 'done')

    parse_meta(args.mapping_cell_type)


if __name__ == '__main__':
    main()

