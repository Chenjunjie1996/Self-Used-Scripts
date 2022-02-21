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
import math
from concurrent.futures import ProcessPoolExecutor


def parse_bam(out_path, sample_name, read_limit):
    try:
        contig_file = glob.glob(f'{out_path}//03.assemble/*/outs/filtered_contig_annotations.csv')[0]
        bam_file = glob.glob(f'{out_path}//03.assemble/*/outs/all_contig.bam')[0]
    except:
        contig_file = glob.glob(f'{out_path}//*/outs/filtered_contig_annotations.csv')[0]
        bam_file = glob.glob(f'{out_path}//*/outs/all_contig.bam')[0]

    contig = pd.read_csv(contig_file)
    productive_contig = contig[contig['productive'] == True]
    productive_contig = set(productive_contig['contig_id'])
    bam = pysam.AlignmentFile(bam_file, "rb")
    pos_dict = defaultdict(int)
    read_count = 0
    for read in bam:
        if read_count <= read_limit:
            if str(read.reference_name) in productive_contig:
                read_count += 1
                start, end = read.reference_start, read.reference_end
                ref_length = bam.get_reference_length(read.reference_name)
                start_percent, end_percent = math.floor(start / ref_length * 100), math.ceil(end / ref_length * 100)
                for index in range(start_percent, end_percent + 1):
                    pos_dict[index] += 1
        else:
            break
    bam.close()

    sort_dic = sorted(pos_dict.items(), key=lambda x: x[0])
    y_list = []
    for i in sort_dic:
        y_list.append(i[1])
    result_file = open(f"./Mapped_reads_distribution/{sample_name}.txt", "w")
    for pos, count in sort_dic:
        result_file.write(str(pos) + " : " + str(count))
        result_file.write("\n")
    result_file.close()
    return y_list


def multi_run(out_SGR, out_10X, sample_name_SGR, sample_name_10X, read_limit):
    out_path_list = [out_SGR, out_10X]
    sample_name = [sample_name_SGR, sample_name_10X]
    read_limit_list = [read_limit] * 2
    result = []
    with ProcessPoolExecutor(4) as pool:
        for res in pool.map(parse_bam, out_path_list, sample_name, read_limit_list):
            result.append(res)
    return result[0], result[1]


def make_plot(y_list_SGR, y_list_10X, sample_name_SGR, sample_name_10X):

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
    parser.add_argument('--read_limit', help= 'read_limit', required=True)
    args = parser.parse_args()
    out_SGR = args.out_path1
    out_10X = args.out_path2
    read_limit = int(args.read_limit)
    os.system(f"mkdir Mapped_reads_distribution")

    sample_name_SGR = os.path.abspath(f"{out_SGR}").split('/')[-1]
    sample_name_10X = os.path.abspath(f"{out_10X}").split('/')[-1]
    y_list_SGR, y_list_10X = multi_run(out_SGR, out_10X, sample_name_SGR, sample_name_10X, read_limit)
    make_plot(y_list_SGR, y_list_10X,  sample_name_SGR, sample_name_10X)
