import pandas as pd
import sys
import glob
import xlwt

def parse_meta():
    meta_list = glob.glob("./mapping_res/*_meta.csv")
    sample_list = ['_'.join(i.split('/')[-1].split('_')[:-1]) for i in meta_list]

    report = xlwt.Workbook()
    sheet = report.add_sheet('mapping count ')

    row0 = ['number of mapping to T/B cells', 'number of T/B cells', 'percent']
    row1 = sample_list
    mapping_count, ZL_count = [], []
    for meta in meta_list:
        df = pd.read_csv(meta)
        
        df_count = df[df['CellTypes'].isin(Celltype)]
        _i,_j = df_count[df_count['Class']=='T/BCR'].shape[0], df_count['CellTypes'].shape[0]
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

    report.save('./mapping_res/mapping_count.xls')
    return report

if __name__ == '__main__':
    mapping_type = sys.argv[1]
    if mapping_type == "T":
        Celltype = {'T_cells','NKT_cells','T cells','NK T cells','Tcells'}
    elif mapping_type == "B":
        Celltype = {'Plasma_cells','B_cells','Mature_B_cell', 'Plasma cells', 'B cells','Bcells'}
    parse_meta()