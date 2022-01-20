import xlwt
import sys
import pandas as pd
from matplotlib_venn import venn2
from matplotlib import pyplot as plt
import os
import glob

# clonotypes result between cellranger and trust4


def parse_clontypes_file(Cellranger, Trust):
    cr_sample_name = os.path.basename(os.path.abspath(Cellranger))
    Cellranger = glob.glob(f'{Cellranger}/03.assemble/{cr_sample_name}/outs/clonotypes.csv')[0]
    df_cr = pd.read_csv(Cellranger,sep=',')
    Trust = glob.glob(f'{Trust}/04.summarize/clonotypes.csv')[0]
    df_tr = pd.read_csv(Trust,sep=',')
    # tr_sample_name = os.path.abspath(Trust)[-3]

    data_cr_list = df_cr['cdr3s_aa'].tolist()
    data_tr_list = df_tr['cdr3s_aa'].tolist()
    intersec_cr_tr = set(data_cr_list).intersection(set(data_tr_list))

    clonotypes_num = []
    intersec_clonotypes_num = []
    clonotypes_num.append(len(data_cr_list))
    clonotypes_num.append(len(data_tr_list))
    intersec_clonotypes_num.append((len(intersec_cr_tr)))

    df_cr_tr = df_cr[df_cr['cdr3s_aa'].isin(intersec_cr_tr)]
    df_tr_cr = df_tr[df_tr['cdr3s_aa'].isin(intersec_cr_tr)]
    df_cr_tr.to_csv('./df_cr_tr.txt', sep='\t', index=False)
    df_tr_cr.to_csv('./df_tr_cr.txt', sep='\t', index=False)

    data_cr_set = set(data_cr_list)
    data_tr_set = set(data_tr_list)
    subset_ = [data_cr_set, data_tr_set]
    venn2(subset_, set_labels=('Cellranger', 'Trust'), set_colors=('r', 'g'))
    plt.savefig('./Vnplot')

    return clonotypes_num, intersec_clonotypes_num


def write_to_sheet(clonotypes_num,intersec_clonotypes_num):
    report = xlwt.Workbook()
    sheet = report.add_sheet('clonotypes count')

    row0 = ['CR_data','TR_data']
    row1 = ['clonotypes_count', 'intersection_clonotypes_count']
    for data in range(1,3):
        sheet.write(data,0,row0[data-1])
        sheet.write(0,data,row1[data-1])
    for num in range(1,3):
        sheet.write(num,1,clonotypes_num[num-1])
    sheet.write(1,2,intersec_clonotypes_num[0])
    sheet.write(2,2,intersec_clonotypes_num[0])

    report.save('./clonotypes_count.xlsx')
    return report

if __name__ == '__main__':
    data_CR = sys.argv[1]
    data_TR = sys.argv[2]
    clonotypes_num, intersec_clonotypes_num = parse_clontypes_file(data_CR, data_TR)
    write_to_sheet(clonotypes_num, intersec_clonotypes_num)

