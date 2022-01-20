from __future__ import print_function
import pandas as pd
from Bio.Seq import Seq
import pysam
import os
from Bio import pairwise2 as pw2
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
import scipy as sp
import sys
import glob

def reversed_compl(seq):
    return str(Seq(seq).reverse_complement())

def percent_match(seq1,seq2):
    global_align = pw2.align.globalxx(seq1,seq2)
    seq_length = min(len(seq1),len(seq2))
    matches = global_align[0][2]
    percent_match = (matches/seq_length) * 100
    return percent_match


def parse_reads(CellrangerPath, TrustPath):
    trust = glob.glob(f'{_TrustPath}/04.summarize/*_contig.csv')[0]
    cellranger = f'{_CellrangerPath}/03.assemble/all/filtered_contig_annotations.csv'

    df_10X = pd.read_csv(cellranger)
    cell_barcodes_10X = (set(df_10X['barcode'].tolist()))

    df_trust = pd.read_csv(trust, sep=',')
    df_trust = df_trust[(df_trust['chain'] == 'TRA') | (df_trust['chain'] == 'TRB')]
    df_trust = df_trust.sort_values(by='umis', ascending=False)
    df_trust.drop_duplicates(['barcode', 'chain'], inplace=True)

    Path = '/'.join(trust.split('/')[:-2])
    filterbc_rep = f'{Path}/03.assemble/assemble/barcoderepfl.tsv'
    trust_rep = f'{Path}/03.assemble/assemble/report.out'
    filterbc_rep = pd.read_csv(filterbc_rep, sep='\t')
    trust_rep = pd.read_csv(trust_rep, sep='\t')

    filterbc_rep = filterbc_rep.rename(columns={'#barcode': 'barcode'})
    filterbc_rep = filterbc_rep[(filterbc_rep['chain1'] != '*') | (filterbc_rep['chain2'] != '*')]
    filterbc = set(filterbc_rep['barcode'].tolist())

    trust_rep = trust_rep[trust_rep['cid_full_length'] >= 1]
    trust_rep = trust_rep[trust_rep['CDR3aa'] != 'out_of_frame']
    trust_rep = trust_rep.rename(columns={'cid': 'barcode', '#count': 'count'})
    trust_rep['barcode'] = trust_rep['barcode'].apply(lambda x: x.split('_')[0])
    trust_rep = trust_rep.sort_values(by='count', ascending=False)
    trust_rep = trust_rep[trust_rep['count'] >= 4]
    trust_bc = set(trust_rep['barcode'].tolist())

    df_trust = df_trust[df_trust['barcode'].isin(filterbc)]
    df_trust = df_trust[df_trust['barcode'].isin(trust_bc)]
    df_trust = df_trust[df_trust['productive'] == True]
    cell_barcodes = set(df_trust['barcode'].tolist())
    cell_barcodes = [reversed_compl(i) for i in cell_barcodes]

    intersec_cells = (set(cell_barcodes).intersection((cell_barcodes_10X)))

    return intersec_cells

def gene_fq(intersec_cells):
    os.system(f'mkdir -p ./Seq')
    full_len_10X = glob.glob(f'{_CellrangerPath}/03.assemble/*/outs/airr_rearrangement.tsv')[0]
    full_len_10X = pd.read_csv(full_len_10X, sep='\t')
    full_len_trust = glob.glob(f'{_TrustPath}/03.assemble/assemble/*_full_len.fa')[0]
    barcode_dict = pd.read_csv(f'{_CellrangerPath}/02.convert/barcode_correspond.txt', sep='\t',index_col=1)
    barcode_dict = barcode_dict.to_dict()['sgr']

    full_len_10X['cell_id'] = full_len_10X['cell_id'].apply(lambda x: reversed_compl(barcode_dict[x.split('-')[0]]))
    full_len_10X['sequence_id'] = full_len_10X['sequence_id'].apply(lambda x: reversed_compl(barcode_dict[x.split('-')[0]]) + '-' + x.split('-')[1])
    full_len_10X = full_len_10X[full_len_10X['cell_id'].isin(intersec_cells)]
    full_len_10X['v_call'] = (full_len_10X['v_call'].apply(lambda x: x[:3]))
    airr_A = full_len_10X[full_len_10X['v_call'] == 'TRA']
    airr_B = full_len_10X[full_len_10X['v_call'] == 'TRB']
    airr_A = airr_A.groupby('sequence_id')['sequence'].apply(lambda x: x.tolist()).to_dict()
    airr_B = airr_B.groupby('sequence_id')['sequence'].apply(lambda x: x.tolist()).to_dict()

    TRA_fq_10X = './Seq/TRA_10X.fq'
    TRB_fq_10X = './Seq/TRB_10X.fq'
    with open(TRA_fq_10X, 'w') as f:
        for key, values in airr_A.items():
            f.write('@' + key + '\n' + values[0] + '\n')
    with open(TRB_fq_10X, 'w') as f:
        for key, values in airr_B.items():
            f.write('@' + key + '\n' + values[0] + '\n')

    intersecRev_cells = [reversed_compl(i) for i in intersec_cells]
    TRA_fq_trust = './Seq/TRA_trust.fq'
    TRB_fq_trust = './Seq/TRB_trust.fq'
    trust_fq_TRA = open(TRA_fq_trust, 'w')
    trust_fq_TRB = open(TRB_fq_trust, 'w')
    with pysam.FastxFile(full_len_trust, 'r') as fa:
        for read in fa:
            name = read.name
            barcode = name.split('_')[0]
            seq = read.sequence
            comment = read.comment
            TR_type = comment.split(' ')[2][:3]
            if barcode in intersecRev_cells:
                if TR_type == 'TRA':
                    trust_fq_TRA.write('@' + reversed_compl(barcode) + '_' + name.split('_')[1] + '\n' + seq + '\n')
                elif TR_type == 'TRB':
                    trust_fq_TRB.write('@' + reversed_compl(barcode) + '_' + name.split('_')[1] + '\n' + seq + '\n')
    trust_fq_TRA.close()
    trust_fq_TRB.close()
    fq_list = [TRA_fq_10X, TRB_fq_10X, TRA_fq_trust, TRB_fq_trust]
    return fq_list

def parse_fq(fq_list):
    TRA_trust_dict = defaultdict(list)
    TRA_10X_dict = defaultdict(list)
    TRB_trust_dict = defaultdict(list)
    TRB_10X_dict = defaultdict(list)
    for fq in fq_list:
        fq_type = os.path.basename(fq)
        with pysam.FastxFile(fq, 'r') as f:
            for read in f:
                name = read.name
                if fq_type == 'TRA_10X.fq':
                    barcode = name.split('-')[0]
                    TRA_10X_dict[barcode].append(read.sequence)
                elif fq_type == 'TRB_10X.fq':
                    barcode = name.split('-')[0]
                    TRB_10X_dict[barcode].append(read.sequence)
                elif fq_type == 'TRA_trust.fq':
                    barcode = name.split('_')[0]
                    TRA_trust_dict[barcode].append(read.sequence)
                elif fq_type == 'TRB_trust.fq':
                    barcode = name.split('_')[0]
                    TRB_trust_dict[barcode].append(read.sequence)

    return TRA_10X_dict, TRB_10X_dict, TRA_trust_dict, TRB_trust_dict

def similarity_seq(TRA_10X_dict, TRB_10X_dict, TRA_trust_dict, TRB_trust_dict):
    similarity_TRA = defaultdict(list)
    similarity_TRB = defaultdict(list)
    intersec_TRA = set(TRA_10X_dict) & set(TRA_trust_dict)
    intersec_TRB = set(TRB_10X_dict) & set(TRB_trust_dict)

    for _TRA in TRA_trust_dict:
        if _TRA not in intersec_TRA:
            similarity_TRA[_TRA].append(0)
        else:
            if len(TRA_trust_dict[_TRA]) == 1 and len(TRA_10X_dict[_TRA]) == 1:
                correlation = percent_match(TRA_trust_dict[_TRA][0], TRA_10X_dict[_TRA][0])
                similarity_TRA[_TRA].append(correlation)
            elif len(TRA_trust_dict[_TRA]) > 1 or len(TRA_10X_dict[_TRA]) > 1:
                temp_list = []
                for _i in TRA_trust_dict[_TRA]:
                    for _j in TRA_10X_dict[_TRA]:
                        temp_list.append(percent_match(_i, _j))
                correlation = max(temp_list)
                similarity_TRA[_TRA].append(correlation)

    for _TRB in TRB_trust_dict:
        if _TRB not in intersec_TRB:
            similarity_TRB[_TRB].append(0)
        else:
            if len(TRB_trust_dict[_TRB]) == 1 and len(TRB_10X_dict[_TRB]) == 1:
                correlation = percent_match(TRB_trust_dict[_TRB][0], TRB_10X_dict[_TRB][0])
                similarity_TRB[_TRB].append(correlation)
            elif len(TRB_trust_dict[_TRB]) > 1 or len(TRB_10X_dict[_TRB]) > 1:
                temp_list = []
                for _i in TRB_trust_dict[_TRB]:
                    for _j in TRB_10X_dict[_TRB]:
                        temp_list.append(percent_match(_i, _j))
                correlation = max(temp_list)
                similarity_TRB[_TRB].append(correlation)

    return similarity_TRA, similarity_TRB

def Corplot(similarity_TRA, similarity_TRB):
    df_TRA = pd.DataFrame.from_dict(similarity_TRA,orient='index')
    df_TRA = df_TRA.rename_axis('barcode').reset_index()
    df_TRB = pd.DataFrame.from_dict(similarity_TRB,orient='index')
    df_TRB = df_TRB.rename_axis('barcode').reset_index()
    df_TRA = df_TRA.rename(columns={0:'similarity'})
    df_TRB = df_TRB.rename(columns={0:'similarity'})
    data = pd.merge(df_TRA, df_TRB, on='barcode',how='outer',suffixes=('_TRA','_TRB'))
    data.fillna(0, inplace=True)
    data.to_excel('./count.xlsx', index=False)
    data = data[['similarity_TRA','similarity_TRB']]
    sns.set(font_scale = 1.5, style = 'whitegrid')
    snsFig = sns.relplot(x="similarity_TRA", y="similarity_TRB", data=data)
    snsFig.set(xlim=[-0.1, 100.5])
    snsFig.set(ylim=[-0.1, 100.5])
    plt.xlabel("TRA Similarity")
    plt.ylabel("TRB Similarity")
    r = sp.stats.pearsonr(data['similarity_TRA'].tolist(), data['similarity_TRB'].tolist())
    plt.text(55, 30, s="cor=%.3f" % r[0])
    plt.plot([0, 100], [0, 100], ls="--", c="g")
    plt.title("Trust-4 Vs Cellranger")
    plt.savefig("./SimilarityCorrelation_zero.png",bbox_inches='tight',dpi=300)

'''
    snsFig.set(xlim=[60, 100.5])
    snsFig.set(ylim=[60, 100.5])
    plt.text(85, 65, s="cor=%.3f" % r[0])
    plt.plot([60, 100], [60, 100], ls="--", c="g")
    plt.savefig("./SimilarityCorrelation.png", bbox_inches='tight',dpi=300)
'''

if __name__ == '__main__':
    _CellrangerPath = sys.argv[1]
    _TrustPath = sys.argv[2]

    intersec_cells = parse_reads(_CellrangerPath, _TrustPath)
    fq_list = gene_fq(intersec_cells)
    TRA_10X_dict, TRB_10X_dict, TRA_trust_dict, TRB_trust_dict = parse_fq(fq_list)
    similarity_TRA, similarity_TRB = similarity_seq(TRA_10X_dict, TRB_10X_dict, TRA_trust_dict, TRB_trust_dict)
    Corplot(similarity_TRA, similarity_TRB)