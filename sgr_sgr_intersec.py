import xlwt
import sys
import pandas as pd
from matplotlib_venn import venn2
from matplotlib import pyplot as plt
import os

def parse_clonotypes_file(data_SGR1, data_SGR2):
    SGR1 = data_SGR1.split('/')[5]
    SGR2 = data_SGR2.split('/')[5]
    data_SGR1 = pd.read_csv(data_SGR1,sep=',')
    data_SGR2 = pd.read_csv(data_SGR2,sep=',')
    data_SGR1_list = data_SGR1['cdr3s_aa'].tolist()
    data_SGR2_list = data_SGR2['cdr3s_aa'].tolist()
    intersec_SGR = set(data_SGR1_list).intersection(set(data_SGR2_list))

    clonotypes_num = []
    intersec_clonotypes_num = []
    clonotypes_num.append(len(data_SGR1_list))
    clonotypes_num.append(len(data_SGR2_list))
    intersec_clonotypes_num.append((len(intersec_SGR)))

    df_SGR1_SGR2 = data_SGR1[data_SGR1['cdr3s_aa'].isin(intersec_SGR)]
    df_SGR2_SGR1 = data_SGR2[data_SGR2['cdr3s_aa'].isin(intersec_SGR)]
    df_SGR1_SGR2.to_csv(f'./df_{SGR1}.txt', sep='\t', index=False)
    df_SGR2_SGR1.to_csv(f'./df_{SGR2}.txt', sep='\t', index=False)
    
    data_SGR1_set = set(data_SGR1_list)
    data_SGR2_set = set(data_SGR2_list)
    subset_ = [data_SGR1_set,data_SGR2_set]
    venn2(subset_, set_labels = (SGR1, SGR2), set_colors=('r', 'g'))
    plt.savefig('./Vnplot')
    
    return clonotypes_num, intersec_clonotypes_num, SGR1, SGR2


def write_to_sheet(clonotypes_num,intersec_clonotypes_num, SGR1, SGR2):
    report = xlwt.Workbook()
    sheet = report.add_sheet('clonotypes count')
    
    row0 = [SGR1,SGR2]
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

def main():
    data_SGR1 = sys.argv[1]
    data_SGR2 = sys.argv[2]
    clonotypes_num, intersec_clonotypes_num, SGR1, SGR2 = parse_clonotypes_file(data_SGR1, data_SGR2)
    write_to_sheet(clonotypes_num, intersec_clonotypes_num, SGR1, SGR2)

if __name__ == '__main__':
    main()