import pysam
import pandas as pd
import numpy as np
import argparse
import glob
from collections import defaultdict
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import os
from concurrent.futures import ProcessPoolExecutor


def parse_bam(out_path):
    try:
        contig_file = glob.glob(f'{out_path}//03.assemble/*/outs/filtered_contig_annotations.csv')[0]
        bam_file = glob.glob(f'{out_path}//03.assemble/*/outs/all_contig.bam')[0]
    except:
        contig_file = glob.glob(f'{out_path}//*/outs/filtered_contig_annotations.csv')[0]
        bam_file = glob.glob(f'{out_path}//*/outs/all_contig.bam')[0]

    contig = pd.read_csv(contig_file)
    filter_contig = contig[contig['productive'] == True]
    pro_dict = filter_contig.groupby('contig_id')['length'].apply(lambda x: x.tolist()).to_dict()
    pro_contig = list(pro_dict.keys())
    start_pos_list, end_pos_list = [], []
    bam = pysam.AlignmentFile(bam_file, "rb")
    read_limit = 0
    for read in bam:
        if str(read.reference_name) in pro_contig:
            read_limit += 1
            if read_limit <= 10000000:
                start_pos_list.append(read.reference_start)
                end_pos_list.append(read.reference_end)
            else:
                break
    bam.close()

    # start_index = max(start_pos_list) / 100
    end_index = max(end_pos_list) / 100
    start_list = [round(j/end_index) for j in start_pos_list]
    end_list = [round(j/end_index) for j in end_pos_list]
    count_dic = defaultdict(int)
    for i in range(len(start_list)):
        for j in range(start_list[i], end_list[i] + 1, 1):
            count_dic[j] += 1

    x_list = [i for i in range(101)]
    if len(count_dic.keys()) != len(x_list):
        difference_set = set(x_list) - set(count_dic.keys())
        for i in difference_set:
            count_dic[i] = 0
    sort_dic = sorted(count_dic.items(), key=lambda count_dic: count_dic[0])
    y_list = []
    for i in sort_dic:
        y_list.append(i[1])

    return y_list


def multi_run(out_SGR, out_10X):
    out_path_list = [out_SGR, out_10X]
    result = []
    with ProcessPoolExecutor(4) as pool:
        for res in pool.map(parse_bam, out_path_list):
            result.append(res)
    return result[0], result[1]


def make_plot(y_list_SGR, y_list_10X,  sample_name_SGR, sample_name_10X):

    x_list = [i for i in range(101)]
    plt.figure(figsize=(12, 9))
    plt.grid(linestyle="--")
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.plot(x_list, y_list_SGR, linewidth=1.5, color='red', label=f"SGR_{sample_name_SGR}")
    plt.plot(x_list, y_list_10X, linewidth=1.5, color='blue', label=f"10X_{sample_name_10X}")

    plt.title("Mapped Reads Distribution", fontsize=20)
    plt.xlabel("Length of Assembled chain (%)", fontsize=12)
    plt.ylabel("Read Count", fontsize=12)

    plt.tick_params(axis='both', labelsize=10)
    plt.legend(loc=1, numpoints=1)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=12, fontweight='bold')

    plt.savefig("./Mapped_reads_distribution/mapped_reads_distribution.png", bbox_inches='tight', dpi=300)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='reads distribution')
    parser.add_argument('--out_path1', help='SGR path', required=True)
    parser.add_argument('--out_path2', help='10X', required=True)
    args = parser.parse_args()
    out_SGR = args.out_path1
    out_10X = args.out_path2
    os.system(f"mkdir Mapped_reads_distribution")

    sample_name_SGR = os.path.abspath(f"{out_SGR}").split('/')[-1]
    sample_name_10X = os.path.abspath(f"{out_10X}").split('/')[-1]
    y_list_SGR, y_list_10X = multi_run(out_SGR, out_10X)
    make_plot(y_list_SGR, y_list_10X,  sample_name_SGR, sample_name_10X)
