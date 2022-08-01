import pandas as pd
import pysam
import sys
import copy
import glob
import matplotlib.pyplot as plt
from Bio.Seq import Seq

"""
columns
Clonotype(clonotype_id|V-gene|D-gene|J-gene|CDR3|C-gene|)       Cell_Number     Cgene   Sequence
"""
def get_anno_seq(cr_out):
    airr = pd.read_csv(f'{cr_out}/airr_rearrangement.tsv',sep='\t')
    contig = pd.read_csv(f'{cr_out}/filtered_contig_annotations.csv')
    clonotype = pd.read_csv(f'{cr_out}/clonotypes.csv')
    
    contig = contig[contig['productive']==True]
    contig_for_plot = copy.deepcopy(contig)
    contig = contig[['contig_id','raw_clonotype_id','cdr3','reads']]
    contig.rename(columns={'contig_id':'sequence_id'},inplace=True)
    df_merge = pd.merge(airr,contig,on='sequence_id',how='inner')
    del df_merge['clone_id']
    df_merge.rename(columns={'raw_clonotype_id':'clonotype_id'},inplace=True)
    
    df_merge_f = pd.merge(df_merge,clonotype,on='clonotype_id')
    df_merge_f['chain'] = df_merge_f['v_call'].apply(lambda x:x[:3])
    df_merge_f['id_num'] = df_merge_f['clonotype_id'].apply(lambda x: int(x[9:]))
    df_merge_f.sort_values(['id_num','chain','reads'],ascending=[True,True,False],inplace=True)
    
    df_merge_f.drop_duplicates(['clonotype_id','chain'],inplace=True)
    df_merge_f = df_merge_f.where(df_merge_f.notnull(),'None')
    clonotype_id_list = list(df_merge_f.clonotype_id.unique())
    
    outfile = open('./anno_seq_aa.fasta','w')
    outfile.write("Clonotype(clonotype_id|V-gene|D-gene|J-gene|CDR3|C-gene|)\tCell_Number\tCgene\tSequence\n")
    for i in clonotype_id_list:
        df_tmp = df_merge_f[df_merge_f['clonotype_id']==i]
        outfile.write(f'{list(df_tmp.clonotype_id)[0]}|')
        for i in range(df_tmp.shape[0]):
            outfile.write(f'{list(df_tmp.v_call)[i]}|{list(df_tmp.d_call)[i]}|{list(df_tmp.j_call)[i]}|{list(df_tmp.cdr3)[i]}|{list(df_tmp.c_call)[i]}|')
        outfile.write(f'\t{list(df_tmp.frequency)[0]}\t')
        outfile.write('|'.join(list(df_tmp.c_call)))
        outfile.write('\n')
        for i in range(df_tmp.shape[0]):
            outfile.write(f'ChianType:{list(df_tmp.chain)[i]}\n{list(df_tmp.sequence_aa)[i]}\n')
    outfile.close()


    outfile = open('./anno_seq_nt.fasta','w')
    outfile.write("Clonotype(clonotype_id|V-gene|D-gene|J-gene|CDR3|C-gene|)\tCell_Number\tCgene\tSequence\n")
    for i in clonotype_id_list:
        df_tmp = df_merge_f[df_merge_f['clonotype_id']==i]
        outfile.write(f'{list(df_tmp.clonotype_id)[0]}|')
        for i in range(df_tmp.shape[0]):
            outfile.write(f'{list(df_tmp.v_call)[i]}|{list(df_tmp.d_call)[i]}|{list(df_tmp.j_call)[i]}|{list(df_tmp.cdr3)[i]}|{list(df_tmp.c_call)[i]}|')
        outfile.write(f'\t{list(df_tmp.frequency)[0]}\t')
        outfile.write('|'.join(list(df_tmp.c_call)))
        outfile.write('\n')
        for i in range(df_tmp.shape[0]):
            seq_nt = list(df_tmp.sequence)[i]
            seq_aa = list(df_tmp.sequence_aa)[i]
            for j in range(len(seq_nt)):
                if str(Seq(seq_nt[j:]).translate()) == seq_aa:
                    seq_nt_write = seq_nt[j:]
                    if len(seq_nt_write) % 3 == 1:
                        seq_nt_write = seq_nt[j:-1]
                    elif len(seq_nt_write) % 3 == 2:
                        seq_nt_write = seq_nt[j:-2]
                    outfile.write(f'ChianType:{list(df_tmp.chain)[i]}\n{seq_nt_write}\n')
                    break


    df_h = contig_for_plot[contig_for_plot['chain']=='IGH']
    tmp_dict = df_h.groupby('c_gene')['raw_clonotype_id'].count().to_dict()
    try:
        del(tmp_dict['None'])
    except KeyError:
        pass

    tmp_dict,reverse_tmp_dict = sorted(tmp_dict.items(),key=lambda x:x[1]), sorted(tmp_dict.items(),key=lambda x:x[1], reverse=True)
    
    label, value = [], []
    if len(tmp_dict) % 2 == 0:
        for i in range(len(tmp_dict)//2):
            label.append(tmp_dict[i][0])
            label.append(reverse_tmp_dict[i][0])
            value.append(tmp_dict[i][1])
            value.append(reverse_tmp_dict[i][1])
    else:
        for i in range(len(tmp_dict)//2):
            label.append(tmp_dict[i][0])
            label.append(reverse_tmp_dict[i][0])
            value.append(tmp_dict[i][1])
            value.append(reverse_tmp_dict[i][1])
        label.append(tmp_dict[len(tmp_dict)//2][0])
        value.append(tmp_dict[len(tmp_dict)//2][1])
    
    plt.pie(value, labels=label,autopct='%.2f%%',pctdistance=0.8,explode=[0.1]*len(label))
    plt.title('IG types')
    plt.savefig('./IGH.png', dpi=280)


if __name__ == '__main__':
    cr_out = sys.argv[1]
    get_anno_seq(cr_out)
    
    
