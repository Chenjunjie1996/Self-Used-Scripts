import os
import subprocess
import glob
import pandas as pd
import sys
import argparse
from concurrent.futures import ProcessPoolExecutor
import xlwt
import os


def run_mapping(rds, contig, sample, outdir, contig_name):
    cmd = (
        f'Rscript /SGRNJ03/randd/cjj/Script/mappingManual_multi/mappingManual_multi.R '
        f'--rds {rds} '
        f'--VDJ {contig} '
        f'--sample {sample} '
        f'--outdir {outdir} '
        f'--contig_name {contig_name} '
    )
    subprocess.check_call(cmd, shell=True)


def parse_file(contig_file):
    contigs = contig_file.strip().split(",")
    contig_name_list = [os.path.abspath(i).split('/')[-1] for i in contigs]
    contig_list = []
    for contig in contigs:
        if os.path.exists(f'{contig}/05.match/match_contigs.csv'):
            contig_file = glob.glob(f'{contig}/05.match/match_contigs.csv')[0]
            contig_list.append(contig_file)
        elif os.path.exists(f'{contig}/04.summarize'):
            contig_file = glob.glob(f'{contig}/04.summarize/*_filtered_contig.csv')[0]
            contig_list.append(contig_file)
        elif os.path.exists(f'{contig}/03.summarize'):
            contig_file = glob.glob(f'{contig}/03.summarize/*_filtered_contig.csv')[0]
            contig_list.append(contig_file)
        elif os.path.exists(f'{contig}/05.count_vdj'):
            contig_file = glob.glob(f'{contig}/05.count_vdj/*_cell_confident_count.tsv')[0]
            contig_list.append(contig_file)
        else:
            print("check contig file path")
            raise FileNotFoundError

    return contig_list, contig_name_list


def parse_meta(mapping_cell_type, outdir):
    meta_list = glob.glob(f"./{outdir}/*_meta.csv")
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

    report.save(f'./{outdir}/mapping_count.xls')
    return report


def main():
    parser = argparse.ArgumentParser(description='contigfile')
    parser.add_argument('--contig_file', help='contig file', required=True)
    parser.add_argument('--rds', help='rds file', required=True)
    parser.add_argument('--mapping_cell_type', help='Tcells or Bcells', required=True)
    parser.add_argument('--sample_name', help='sample name in rds', required=True)
    args = parser.parse_args()
    contig_file = args.contig_file
    rds = args.rds
    sample_name = args.sample_name
    outdir = f"./mappingManual_{sample_name}"
    os.system(f"mkdir {outdir}")
    contig_list, contig_name_list = parse_file(contig_file)
    rds_list = [rds] * len(contig_list)
    outdir_list = [outdir] * len(contig_list)
    rds_sample = [sample_name] * len(contig_list)

    with ProcessPoolExecutor(max_workers=8) as executor:
        for result in executor.map(run_mapping, rds_list, contig_list, rds_sample, outdir_list, contig_name_list):
            print(result, 'done')

    parse_meta(args.mapping_cell_type, outdir)


if __name__ == '__main__':
    main()
