import pandas as pd
import sys
import os
import xlwt
import subprocess

def run_plot(rds, UMI_file, sample, outdir):
    cmd = (
        f'Rscript /SGRNJ03/randd/cjj/Script/TCR_T/plot.R '
        f'--rds {rds} '
        f'--UMI_file {UMI_file} '
        f'--sample {sample} '
        f'--outdir {outdir}'
    )
    subprocess.check_call(cmd,shell=True)

def parse_file(outdir, sample):
    file = pd.read_csv(f'{outdir}/{sample}_meta.txt')

    report = xlwt.Workbook()
    sheet = report.add_sheet('TCRT-Count')

    row0 = [f"{sample}_cluster","CellNumber","Cell_With_TCR-T","Percent"]
    Clusters = sorted(list(file['seurat_clusters'].unique()))
    Clusters = [int(i) for i in Clusters]
    CellNumber_list,Cell_With_TCR_list = [],[]

    file_TCRT = file[file['UMI'] != 0]
    for cluster in Clusters:
        count = file[file['seurat_clusters'] == cluster].shape[0]
        CellNumber_list.append(count)
        countTCRT = file_TCRT[file_TCRT['seurat_clusters'] == cluster].shape[0]
        Cell_With_TCR_list.append(countTCRT)

    percent = [x / y for x, y in zip(Cell_With_TCR_list, CellNumber_list)]
    percent = ['%.2f'%(i*100)  + "%" for i in percent]

    for _i in range(len(row0)):
        sheet.write(0, _i, row0[_i])
    for _i in range(len(Clusters)):
        sheet.write(_i+1, 0, Clusters[_i])
        sheet.write(_i+1, 1, CellNumber_list[_i])
        sheet.write(_i+1, 2, Cell_With_TCR_list[_i])
        sheet.write(_i+1, 3, percent[_i])

    report.save(f'{outdir}/count.xls')
    return report


if __name__ == '__main__':
    sample, rds, UMI_file = sys.argv[1], sys.argv[2], sys.argv[3]
    outdir = os.getcwd() + "/" + sample
    os.system(f"mkdir {outdir}")
    run_plot(rds, UMI_file, sample, outdir)
    report = parse_file(outdir,sample)
