import pandas as pd
import pysam
from Bio.Seq import  Seq
import glob
import os
import sys

def reversed_compl(seq):
    return str(Seq(seq).reverse_complement())

def parse_fa(sample_file):
    loc_file = open("./VDJ_seq_loc/loc.csv","w")
    full_len_fa = glob.glob(f'{sample_file}/03.assemble/assemble/*_full_len.fa')[0]
    with pysam.FastxFile(full_len_fa) as fa:
        for read in fa:
            name = read.name
            barcode = name.split('_')[0]
            attrs = read.comment.split(' ')
            v_gene = attrs[2]
            d_gene = attrs[3]
            j_gene = attrs[4]
            cdr3 = attrs[8]
            chain = attrs[2][:3]
            sequence = read.sequence
            string = '\t'.join([barcode, name, chain, v_gene, d_gene, j_gene, cdr3, sequence])
            loc_file.write(f'{string}\n')
    loc_file.close()

    df = pd.read_csv("./VDJ_seq_loc/loc.csv",sep='\t')
    sample_name = os.path.basename(os.path.abspath(sample_file))
    filtered_contig = glob.glob(f'{sample_file}/04.summarize/{sample_name}_filtered_contig.csv')[0]
    filtered_contig = pd.read_csv(filtered_contig)
    filtered_contig_id = set(filtered_contig['contig_id'].tolist())
    df.columns = ['barcode', 'contig_id', 'chain', 'V-region', 'D-region', 'J-region', 'CDR3', 'Seq']
    df['barcode'] = df['barcode'].apply(lambda x: reversed_compl(x))
    df['contig_id'] = df['contig_id'].apply(lambda x: reversed_compl(x.split('_')[0]) + '_' + x.split('_')[1])
    df = df[df['contig_id'].isin(filtered_contig_id)]

    df['V-region'] = df['V-region'].apply(
        lambda x: x.split(':')[1].replace('(', '').replace(')', '') if not x == '*' else '*')
    df['D-region'] = df['D-region'].apply(
        lambda x: x.split(':')[1].replace('(', '').replace(')', '') if not x == '*' else '*')
    df['J-region'] = df['J-region'].apply(
        lambda x: x.split(':')[1].replace('(', '').replace(')', '') if not x == '*' else '*')
    df['CDR3_loc'] = df['CDR3'].apply(lambda x: x.split(':')[0].replace('CDR3', '').replace('(', '').replace(')', ''))

    df['V-start'] = df['V-region'].apply(lambda x: int(x.split('-')[0]) if not x == '*' else '*')
    df['V-end'] = df['V-region'].apply(lambda x: int(x.split('-')[1]) if not x == '*' else '*')
    df['D-start'] = df['D-region'].apply(lambda x: int(x.split('-')[0]) if not x == '*' else '*')
    df['D-end'] = df['D-region'].apply(lambda x: int(x.split('-')[1]) if not x == '*' else '*')
    df['J-start'] = df['J-region'].apply(lambda x: int(x.split('-')[0]) if not x == '*' else '*')
    df['J-end'] = df['J-region'].apply(lambda x: int(x.split('-')[1]) if not x == '*' else '*')
    df['CDR3-start'] = df['CDR3_loc'].apply(lambda x: int(x.split('-')[0]) if not x == '*' else '*')
    df['CDR3-end'] = df['CDR3_loc'].apply(lambda x: int(x.split('-')[1]) if not x == '*' else '*')

    df['V-start'] = df['CDR3-start']
    df['J-end'] = df['CDR3-end']

    V_start_index_list = df['V-start'].tolist()
    V_end_index_list = [i + 1 for i in df['V-end'].tolist()]
    D_start_index_list = df['D-start'].tolist()
    D_end_index_list = df['D-end'].tolist()
    D_end_index_list = [i + 1 if i != '*' else '*' for i in D_end_index_list]
    J_start_index_list = df['J-start'].tolist()
    J_end_index_list = [i + 1 for i in df['J-end'].tolist()]

    Seq_list = df['Seq'].tolist()
    CDR3_V_Seq, CDR3_D_Seq, CDR3_J_Seq = [], [], []

    for i in range(len(Seq_list)):
        V_start_index, V_end_index = V_start_index_list[i], V_end_index_list[i]
        J_start_index, J_end_index = J_start_index_list[i], J_end_index_list[i]
        CDR3_V_Seq.append(Seq_list[i][V_start_index:V_end_index])
        CDR3_J_Seq.append(Seq_list[i][J_start_index:J_end_index])

    for i in range(len(Seq_list)):
        D_start_index, D_end_index = D_start_index_list[i], D_end_index_list[i]
        if D_start_index != '*':
            CDR3_D_Seq.append(Seq_list[i][D_start_index:D_end_index])
        else:
            CDR3_D_Seq.append('*')

    df['CDR3_V_seq'] = CDR3_V_Seq
    df['CDR3_D_seq'] = CDR3_D_Seq
    df['CDR3_J_seq'] = CDR3_J_Seq

    df_final = df[['barcode', 'contig_id', 'chain', 'CDR3_V_seq', 'CDR3_D_seq', 'CDR3_J_seq', 'CDR3']]
    df_final['CDR3'] = df_final['CDR3'].apply(lambda x: x.split('=')[1])

    df_final.to_csv('./VDJ_seq_loc/CDR3_VDJ_seq.csv',
                    sep='\t', index=False)


if __name__ == '__main__':
    sample_file = sys.argv[1]
    os.system(f"mkdir VDJ_seq_loc")
    parse_fa(sample_file)
