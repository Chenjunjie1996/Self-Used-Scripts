import pandas as pd
import sys
import glob
import xlwt
import re

def parse_meta():
    meta_list = glob.glob("./*_meta.csv")
    sample_list = ['_'.join(i.split('/')[-1].split('_')[:-1]) for i in meta_list]

    report = xlwt.Workbook()
    sheet = report.add_sheet('mapping count ')

    row0 = ['number of mapping to T/B cells', 'number of T/B cells', 'percent']
    row1 = sample_list
    mapping_count, ZL_count = [], []
    for meta in meta_list:
        df = pd.read_csv(meta)
        df['cluster'] = df['cluster'].apply(lambda x: re.sub(r"[^a-zA-Z0-9]","", str(x)).upper())
        df_count = df[df['cluster'].isin(Celltype)]
        _i,_j = df_count[df_count['Class']=='T/BCR'].shape[0], df_count.shape[0]
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

    report.save('./mapping_count.xls')
    return report

if __name__ == '__main__':
    mapping_type = sys.argv[1]
    if mapping_type == "T":
        Celltype = {'TCELLS', 'TCELL' 'NKTCELLS'}
    elif mapping_type == "B":
        Celltype = {'MATUREBCELL', 'PLASMACELLS', 'BCELLS', 'BCELL', 'PREBCELLCD34'}
    parse_meta()
