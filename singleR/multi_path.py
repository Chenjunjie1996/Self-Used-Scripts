"""
multi-path celescope directory
"""

import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import glob
import pandas as pd
import xlwt

import celescope.tools.utils as utils

ROOT_DIR = os.path.dirname(__file__)
OUT_ROOT = "singleR_annotation"


class NoMatrixError(Exception):
    pass


class singleR():
    def __init__(self, sample, outdir, matrix_file, species):
        self.sample = sample
        self.outdir = outdir
        self.matrix_file = matrix_file
        self.species = species
        self.meta = f'{self.outdir}/{self.sample}_metadata.tsv'

        # out
        self.rds = f'{self.outdir}/{self.sample}.rds'


    @utils.add_log
    def run_seurat(self):
        cmd = (
            f'Rscript {ROOT_DIR}/run_analysis.R '
            f'--sample {self.sample} '
            f'--outdir {self.outdir} '
            f'--matrix_file {self.matrix_file} '
            f'--mt_gene_list None '
            f'--save_rds True '
        )
        self.run_seurat.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)


    @utils.add_log
    def run_singleR(self):
        cmd = (
            f'Rscript {ROOT_DIR}/singleR.R '
            f'--species {self.species} '
            f'--rds {self.rds} '
            f'--sample {self.sample} '
            f'--outdir {self.outdir} '
        )
        subprocess.check_call(cmd, shell=True)
    
    def run_count(self):
        meta = pd.read_csv(self.meta,sep='\t')
        celltype_list = list(meta.celltype.unique())
        count_list,percent_list, cluster_list = [], [], []
        for _i in celltype_list:
            _count = meta[meta['celltype']==f'{_i}'].shape[0]
            count_list.append(_count)
            _percent = int(_count) / meta.shape[0]
            percent_list.append(_percent)
            _cluster = sorted(meta[meta['celltype']==f'{_i}'].seurat_clusters.unique().tolist())
            cluster_list.append(",".join([str(_j) for _j in _cluster]))
        percent_list = ['%.2f'%(i*100)  + "%" for i in percent_list]
        
        report = xlwt.Workbook()
        sheet = report.add_sheet('SingleR-Count')
        row0 = ["Cluster","Celltype","CellNumber","Percent"]
        for _i in range(len(row0)):
            sheet.write(0, _i, row0[_i])
        for _i in range(len(count_list)):
            sheet.write(_i+1, 0, cluster_list[_i])
            sheet.write(_i+1, 1, celltype_list[_i])
            sheet.write(_i+1, 2, count_list[_i])
            sheet.write(_i+1, 3, percent_list[_i])
        report.save(f'{self.outdir}/{self.sample}_count.xls')

    def run(self):
        self.run_seurat()
        self.run_singleR()
        self.run_count()


def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, sep='\t')
    path_list = df_mapfile['path']
    species_list = df_mapfile['species']
    return path_list, species_list


def run_single(sample, path, species):
    outdir = f'{OUT_ROOT}/{sample}'
    utils.check_mkdir(outdir)
    
    try:
        matrix_file = glob.glob(f'{path}/0*.count/{sample}_matrix_10X/')[0]
    except IndexError:
        matrix_file = glob.glob(f'{path}/0*.count/{sample}_filtered_feature_bc_matrix/')[0]
    else:
        raise NoMatrixError

    runner = singleR(sample, outdir, matrix_file, species)
    runner.run()


def main():
    os.system(f'mkdir {OUT_ROOT}')

    mapfile = sys.argv[1]
    path_list, species_list = parse_mapfile(mapfile)
    sample_list = [path.strip('/').split('/')[-1] for path in path_list]

    with ProcessPoolExecutor(max_workers=10) as executor:
        for result in executor.map(run_single, sample_list, path_list, species_list):
            print(result, 'done')


if __name__ == '__main__':
    main()
