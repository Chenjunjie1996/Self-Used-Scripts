import pandas as pd
import os
import glob
import xlwt
from collections import Counter

sample_name = []
l1,l2,l3,l4,l5 =[],[],[],[],[]

def parse_directory(directory):
    files = os.listdir(directory)
    for file in files:
        match_contig_file = f'{file}/03.assemble/match/match_contigs.csv'
        if os.path.exists(match_contig_file):
            sample_name.append(file)
    return sample_name


def count_chain(contig_file):
    for file in contig_file:
        data = pd.read_csv(file, sep=',')
        data = data[(data['productive']==True)&(data['full_length']==True)]
        cbs_nums = len(set(data.barcode.tolist()))
        three_len = 0
        three_len_except = 0
        kk_len = 0
        ll_len = 0
        chain_dict = data.groupby('barcode')['chain'].apply(lambda x:x.tolist()).to_dict()
        for i in chain_dict.values():
            if len(i) >= 3:
                three_len +=1
            if len(i) >= 3 and 'IGL' in i and 'IGK' in i:
                three_len_except += 1
        b = [Counter(i) for i in chain_dict.values()]
        for i in b:
            if i['IGL'] > 1:
                ll_len += 1
            elif i['IGK'] > 1:
                kk_len += 1
        l1.append(cbs_nums)
        l2.append(three_len)
        l3.append(three_len_except)
        l4.append(kk_len)
        l5.append(ll_len)
    return l1,l2,l3,l4,l5


def write_to_sheet():
    report = xlwt.Workbook()
    sheet = report.add_sheet('chains count')
    row0 = ['sample','All productive barcodes number','Barcodes with more than two chains','Barcodes with more than two chains except KK LL', 'Barcodes with KK chains', 'Barcodes with LL chains']
    for i in range(len(row0)):
        sheet.write(0,i,row0[i])
    for i in range(len(sample_name)):
        sheet.write(i+1,0,sample_name[i])
        sheet.write(i+1,1,l1[i])
        sheet.write(i+1,2,l2[i])
        sheet.write(i+1,3,l3[i])
        sheet.write(i+1,4,l4[i])
        sheet.write(i+1,5,l5[i])
    report.save('./chains_count.xlsx')

def main():
    directory = "./"
    contig_file = glob.glob(f'{directory}/*/03.assemble/match/match_contigs.csv')
    
    parse_directory(directory)
    count_chain(contig_file)
    write_to_sheet()
    
    
if __name__ == '__main__':

    main()