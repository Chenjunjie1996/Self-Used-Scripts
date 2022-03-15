import pandas as pd
import sys
import os
import xlwt
import glob


def parse_file(singleR_file, tSNE_file):
    tSNE_file = pd.read_csv(tSNE_file, sep=',')
    tSNE_file.fillna(value = 0, inplace = True)
    tSNE_file = tSNE_file[tSNE_file['BZ_CART'] != 0]
    CAR_barcode = set(tSNE_file.barcode)

    df_meta = pd.read_csv(glob.glob(f'{singleR_file}/*_metadata.tsv')[0], sep='\t')
    sample = os.path.abspath(singleR_file).split('/')[-1]
    df_meta = df_meta.rename_axis('barcode').reset_index()
    df_meta['tag'] = df_meta['barcode'].apply(lambda x: 'positive' if x in CAR_barcode else 'negative')

    report = xlwt.Workbook()
    sheet = report.add_sheet('CART-Count')
    row0 = [f"{sample}_celltype", "Cluster", "CellNumber", "Cell With CAR", "Percent"]

    cell_type_list = list(df_meta['celltype'].unique())
    cluster_list, cell_number_list, car_cell_number_list, percent_list = [], [], [], []

    for i in cell_type_list:
        cluster = sorted(df_meta[df_meta['celltype'] == f'{i}'].seurat_clusters.unique().tolist())
        cluster_list.append(",".join([str(x) for x in cluster]))
        cell_number = df_meta[df_meta['celltype'] == f'{i}'].shape[0]
        cell_number_list.append(cell_number)
        car_cell_number = df_meta[ (df_meta['tag'] == 'positive') & (df_meta['celltype'] == f'{i}') ].shape[0]
        car_cell_number_list.append(car_cell_number)
    percent = [x / y for x, y in zip(car_cell_number_list, cell_number_list)]
    percent = ['%.2f'%(i*100)  + "%" for i in percent]

    for _i in range(len(row0)):
        sheet.write(0, _i, row0[_i])
    for _i in range(len(cell_type_list)):
        sheet.write(_i+1, 0, cell_type_list[_i])
        sheet.write(_i+1, 1, cluster_list[_i])
        sheet.write(_i+1, 2, cell_number_list[_i])
        sheet.write(_i+1, 3, car_cell_number_list[_i])
        sheet.write(_i+1, 4, percent[_i])

    report.save(f'./{sample}_CART_count.xls')
    return report


if __name__ == '__main__':
    singleR_file, tSNE_file = sys.argv[1], sys.argv[2]
    parse_file(singleR_file, tSNE_file)
