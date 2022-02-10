import pandas as pd
import sys
import argparse
from concurrent.futures import ProcessPoolExecutor


def hamming_distance(string1,string2):
    distance = 0
    length, length2 = len(string1), len(string2)
    for i in range(length):
        if string1[i] != string2[i]:
            distance += 1
            if distance > 4:
                break
    return distance

def hamming_process(new_tagA1C, new_tagA1C_rev):
    res = set()
    for i in range(len(new_tagA1C)):
        for j in range(len(new_tagA1C_rev)):
            if hamming_distance(new_tagA1C[i], new_tagA1C_rev[j])<=4:
                break
            if j == len(new_tagA1C_rev) - 1:
                res.add(new_tagA1C[i])
    return res


def parse_file(file1, file2, file3, N):
    file1 = pd.read_csv(file1, sep='\t')
    file2 = pd.read_csv(file2, sep='\t')
    file3 = pd.read_csv(file3, sep='\t')
    file1.sort_values(by='umi_count',ascending=False,inplace=True)
    file2.sort_values(by='umi_count',ascending=False,inplace=True)
    file3.sort_values(by='umi_count',ascending=False,inplace=True)

    # top N tag
    tagA = set(file1.iloc[:N, ]['tag'].tolist())
    tagB = set(file2.iloc[:N, ]['tag'].tolist())
    tagAB = tagA.intersection(tagB)

    # tagA1
    tagA1 = tagA - tagAB

    # tagC
    tagC = set(file3.iloc[100:, ]['tag'].tolist())

    # tagA1C
    tagA1C = tagA1 | tagC

    # del consecutive same bases >= 3
    tagA1C_del_bp = set()
    for i in tagA1C:
        for j in range(1, 14):
            if i[j] == i[j - 1] and i[j] == i[j + 1]:
                tagA1C_del_bp.add(i)
                break
    new_tagA1C = tagA1C - tagA1C_del_bp

    # del hamming distance <= 4
    new_tagA1C = list(new_tagA1C)
    new_tagA1C = sorted(new_tagA1C)
    new_tagA1C_rev = new_tagA1C[::-1]
    index = len(new_tagA1C) // 8
    new_tagA1C_1, new_tagA1C_2, new_tagA1C_3, new_tagA1C_4 = new_tagA1C[:index], new_tagA1C[index:2 * index], new_tagA1C[2 * index:3 * index], new_tagA1C[3 * index:4 * index]
    new_tagA1C_5, new_tagA1C_6, new_tagA1C_7, new_tagA1C_8 = new_tagA1C[4 * index:5 * index], new_tagA1C[5 * index:6 * index], new_tagA1C[6 * index:7 * index], new_tagA1C[7 * index:]
    new_tagA1C_list = [new_tagA1C_1, new_tagA1C_2, new_tagA1C_3, new_tagA1C_4, new_tagA1C_5, new_tagA1C_6, new_tagA1C_7,new_tagA1C_8]
    result = []
    with ProcessPoolExecutor(8) as pool:
        for res in pool.map(hamming_process, new_tagA1C_list, [new_tagA1C_rev]*8):
            result.append(res)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='TAG')
    parser.add_argument('--file1', help='file1', required=True)
    parser.add_argument('--file2', help='file2', required=True)
    parser.add_argument('--file3', help='file3', required=True)
    parser.add_argument('--N', help='top N tag', required=True)
    args = parser.parse_args()
    file1 = args.file1
    file2 = args.file2
    file3 = args.file3
    N = int(args.N)
    result = parse_file(file1, file2, file3, N)
    outfile = open("./tag_seq.txt","w")
    for i in result:
        for j in i:
            outfile.write(j + "\n")
    outfile.close()
    # sys.stdout.write('REQUIRED TAG NUMBER : ' + len(res) + '\n')