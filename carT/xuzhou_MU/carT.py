import os
import subprocess
import glob
import pandas as pd
import sys
import xlwt
import os

def run_plot(rds, UMI_tsne, sample, outdir):
    cmd = (
        f'Rscript /SGRNJ03/randd/cjj/Script/carT/xuzhou_MU/carT.R '
        f'--rds {rds} '
        f'--UMI_tsne {UMI_tsne} '
        f'--sample {sample} '
        f'--outdir {outdir} '
    )
    subprocess.check_call(cmd, shell=True)
    

def run_count(outdir, sample):
    df_meta = pd.read_csv(f'{outdir}/{sample}_metadata.tsv', sep='\t')

    report = xlwt.Workbook()
    sheet = report.add_sheet('CART-Count')
    row0 = [f"{sample}_celltype", "Cluster", "CellNumber", "Cell With CAR"]

    cell_type_list = list(df_meta['celltype'].unique()) #
    cluster_list, cell_number_list, car_cell_number_list = [], [], []

    for i in cell_type_list:
        cluster = sorted(df_meta[df_meta['celltype'] == f'{i}'].seurat_clusters.unique().tolist()) #
        cluster_list.append(",".join([str(int(x)) for x in cluster]))
        cell_number = df_meta[df_meta['celltype'] == f'{i}'].shape[0] # 
        cell_number_list.append(cell_number)
        car_cell_number = df_meta[ (df_meta['Class'] == 'positive') & (df_meta['celltype'] == f'{i}') ].shape[0] #
        car_cell_number_list.append(car_cell_number)

    for _i in range(len(row0)):
        sheet.write(0, _i, row0[_i])
    for _i in range(len(cell_type_list)):
        sheet.write(_i+1, 0, cell_type_list[_i])
        sheet.write(_i+1, 1, cluster_list[_i])
        sheet.write(_i+1, 2, cell_number_list[_i])
        sheet.write(_i+1, 3, car_cell_number_list[_i])

    report.save(f'{outdir}/{sample}_CART_count.xls')
    return report
    
    
def main():
    sample, rds, UMI_tsne = sys.argv[1], sys.argv[2], sys.argv[3]
    outdir = f"./CART_{sample}"
    os.system(f"mkdir {outdir}")
    run_plot(rds, UMI_tsne, sample, outdir)
    run_count(outdir, sample)


if __name__ == '__main__':
    main()






