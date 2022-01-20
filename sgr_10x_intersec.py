import xlwt
import sys
import pandas as pd
from matplotlib_venn import venn2
from matplotlib import pyplot as plt

def parse_clonotypes_file(data_10X, data_SGR):
    data_10X = pd.read_csv(data_10X,sep=',')
    data_SGR = pd.read_csv(data_SGR,sep=',')
    data_10X_list = data_10X['cdr3s_nt'].tolist()
    data_SGR_list = data_SGR['cdr3s_nt'].tolist()
    intersec_SGR_10X = set(data_10X_list).intersection(set(data_SGR_list))

    clonotypes_num = []
    intersec_clonotypes_num = []
    clonotypes_num.append(len(data_10X_list))
    clonotypes_num.append(len(data_SGR_list))
    intersec_clonotypes_num.append((len(intersec_SGR_10X)))

    df_10X_sgr = data_10X[data_10X['cdr3s_nt'].isin(intersec_SGR_10X)]
    df_sgr_10X = data_SGR[data_SGR['cdr3s_nt'].isin(intersec_SGR_10X)]
    df_10X_sgr.to_csv('./df_10X_sgr.txt', sep='\t', index=False)
    df_sgr_10X.to_csv('./df_sgr_10X.txt', sep='\t', index=False)
    
    data_10X_set = set(data_10X_list)
    data_SGR_set = set(data_SGR_list)
    subset_ = [data_10X_set,data_SGR_set]
    venn2(subset_, set_labels = ('10X', 'SGR'), set_colors=('r', 'g'))
    plt.savefig('./Vnplot')
    
    return clonotypes_num, intersec_clonotypes_num


def write_to_sheet(clonotypes_num,intersec_clonotypes_num):
    report = xlwt.Workbook()
    sheet = report.add_sheet('clonotypes count')
    
    row0 = ['10X_data','SGR_data']
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
    data_10X = sys.argv[1]
    data_SGR = sys.argv[2]
    clonotypes_num, intersec_clonotypes_num = parse_clonotypes_file(data_10X, data_SGR)
    write_to_sheet(clonotypes_num, intersec_clonotypes_num)

if __name__ == '__main__':
    main()