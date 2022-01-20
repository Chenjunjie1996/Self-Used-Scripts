import xlwt
import argparse
import glob
import sys
import os

report = xlwt.Workbook()

sample_name = []
region_log_list = []
annotation_stat_list = []

def parse_directory(directory):
    files = os.listdir(directory)
    for file in files:
        if os.path.isdir(file):
            annotation_stat = f'{file}/07.auto_annotation/stat.txt'
            region_log = f'{file}/03.STAR/{file}_region.log'
            sample_name.append(file)
            annotation_stat_list.append(annotation_stat)
            if os.path.exists(region_log):
                region_log_list.append(region_log)
            else:
                region_log = f'{file}/03.star/{file}_region.log'
                region_log_list.append(region_log)

    return sample_name, region_log_list, annotation_stat_list

def parse_directory_stat(directory):
    files = os.listdir(directory)
    for file in files:
        if os.path.isdir(file):
            annotation_stat = f'{file}/07.auto_annotation/stat.txt'
            region_log = f'{file}/03.STAR/stat.txt'
            sample_name.append(file)
            annotation_stat_list.append(annotation_stat)
            if os.path.exists(region_log):
                region_log_list.append(region_log)
            else:
                region_log = f'{file}/03.star/stat.txt'
                region_log_list.append(region_log)

    return sample_name, region_log_list, annotation_stat_list


def write_to_sheet1():
    sheet1 = report.add_sheet('Cell percentage count')
    row0 = ['sample', 'CellType', 'Percentage']
    for i in range(len(row0)):
        sheet1.write(0, i, row0[i])

    count_cell_list = []
    for i in range(len(sample_name)):
        sheet1.write(int(2 * i + 1), 0, sample_name[i])
        annotation_stat_file = annotation_stat_list[i]
        with open(annotation_stat_file, 'r') as f:
            for line in f.readlines():
                line = line.strip('\n').split(':')
                for chara in line:
                    count_cell_list.append(chara)

    for i in range(len(count_cell_list)):
        if i % 2 == 1:
            sheet1.write(int((i + 1) / 2), 2, count_cell_list[i])
        else:
            sheet1.write(int((i + 2) / 2), 1, count_cell_list[i])

    return sheet1


def write_to_sheet2():
    sheet2 = report.add_sheet('RIBOSOMAL COUNT')
    row0 = ['SAMPLE', 'RIBOSOMAL_BASES', 'PF_ALIGNED_BASES', 'Percentage']
    for i in range(len(row0)):
        sheet2.write(0, i, row0[i])  # create columns names

    RIBOSOMAL = []
    PF_ALIGNED_BASES = []

    for i in range(len(sample_name)):
        region_log_file = region_log_list[i]
        with open(region_log_file, 'r') as fh:
            for line in fh:
                if not line:
                    break
                if line.startswith('## METRICS CLASS'):
                    header = fh.readline().strip().split('\t')
                    data = fh.readline().strip().split('\t')
                    dict_log = dict(zip(header, data))
                    RIBOSOMAL.append(dict_log['RIBOSOMAL_BASES'])  # RIBOSOMAL_BASES, if not exist, test by CODING_BASES
                    PF_ALIGNED_BASES.append(dict_log['PF_ALIGNED_BASES'])  # PF_ALIGNED_BASES
                    break
    # count RIBOSOMAL BASES, write to sheet2
    row1 = [i for i in sample_name]
    row2 = [i for i in RIBOSOMAL]
    row3 = [i for i in PF_ALIGNED_BASES]
    row4 = []
    for i in range(len(row1)):  # count RIBOSOMAL_BASES / PF_ALIGNED_BASES
        row4.append(int(row2[i]) / int(row3[i]))

    for i in range(1, len(row1) + 1):
        sheet2.write(i, 0, row1[i - 1])
        sheet2.write(i, 1, row2[i - 1])
        sheet2.write(i, 2, row3[i - 1])
        sheet2.write(i, 3, row4[i - 1])

    return sheet2


def write_to_sheet2_stat():
    sheet2 = report.add_sheet('RIBOSOMAL COUNT')
    row0 = ['SAMPLE', 'Reads Mapped to rRNA']
    for i in range(len(row0)):
        sheet2.write(0, i, row0[i])  # create columns names

    count_ribosome_list = []
    Reads_Mapped_to_rRNA = []

    for i in range(len(sample_name)):
        region_log_file = region_log_list[i]
        with open(region_log_file, 'r') as fh:
            for line in fh:
                Reads_Mapped_to_rRNA.append(line.strip("\n").split(':'))
                dict_rRNA = dict(Reads_Mapped_to_rRNA)
            count_ribosome_list.append(dict_rRNA['Reads Mapped to rRNA'])

    # count RIBOSOMAL BASES, write to sheet2
    row1 = [i for i in sample_name]
    row2 = [i for i in count_ribosome_list]

    for i in range(1, len(row1) + 1):
        sheet2.write(i, 0, row1[i - 1])
        sheet2.write(i, 1, row2[i - 1])

    return sheet2


def main():
    directory = "./"
    method = sys.argv[1]

    if method == 'region':
        parse_directory(directory)
        write_to_sheet1()
        write_to_sheet2()
    elif method == 'stat':
        parse_directory_stat(directory)
        write_to_sheet1()
        write_to_sheet2_stat()
    
    report.save('./Report.xlsx')  # Save xlsx to current folder

if __name__ == '__main__':

    main()








