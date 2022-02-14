import pysam
import pandas as pd
import argparse
import glob
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
import os


def parse_bam(out_path):
    contig_file = glob.glob(f'{out_path}//03.assemble/*/outs/filtered_contig_annotations.csv')[0]
    bam_file = glob.glob(f'{out_path}//03.assemble/*/outs/all_contig.bam')[0]

    contig = pd.read_csv(contig_file)
    filter_contig = contig[contig['productive'] == True]
    pro_dict = filter_contig.groupby('contig_id')['length'].apply(lambda x: x.tolist()).to_dict()
    pro_contig = list(pro_dict.keys())
    pos_list = []
    bam = pysam.AlignmentFile(bam_file, "rb")
    for read in bam:
        if str(read.reference_name) in pro_contig:
            pos_list.append(read.reference_start)
    bam.close()
    return pos_list


def make_plot(pos_list, sample_name):
    x_list = [i for i in range(101)]
    index = max(pos_list) / 100
    pos_count = [round(j/index) for j in pos_list]
    dic = Counter(pos_count)
    if len(dic.keys()) != len(x_list):
        difference_set = set(x_list) - set(dic.keys())
        for i in difference_set:
            dic[i] = 0
    sort_dic = sorted(dic.items(), key=lambda dic: dic[0])
    y_list = []
    for i in sort_dic:
        y_list.append(i[1])

    plt.figure(figsize=(12, 9))
    plt.grid(linestyle="--")
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.plot(x_list, y_list, linewidth=1.5, color='red', label=sample_name)

    plt.title("Mapped Reads Distribution", fontsize=20)
    plt.xlabel("Length of Assembled chain (%)", fontsize=12)
    plt.ylabel("Read Count", fontsize=12)

    plt.tick_params(axis='both', labelsize=10)
    plt.legend(loc=1, numpoints=1)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=12, fontweight='bold')

    plt.savefig(f"{sample_name}/mapped_reads_distribution.png", bbox_inches='tight', dpi=300)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='reads distribution')
    parser.add_argument('--out_path', help='file1', required=True)
    args = parser.parse_args()
    out_path = args.out_path
    sample_name = os.path.abspath(f"{out_path}").split('/')[-1]
    os.system(f"mkdir {sample_name}")
    pos_list = parse_bam(out_path)
    make_plot(pos_list, sample_name)
