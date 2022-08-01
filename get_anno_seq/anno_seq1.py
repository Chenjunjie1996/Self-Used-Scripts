# -*- coding:utf-8 -*-

import pandas as pd
import pysam
import pandas as pd
import sys
from collections import defaultdict
from Bio.Seq import Seq

"""
columns
Clonotype	Frequency	IGH_cdr3	IGH_cdr3_nt	contig_id_IGH	IGH_aa	IGH_nt	IGLorK_cdr3	IGLorK_cdr3_nt	contig_id_IGLorKs	IGLorK_aa	IGLorK_nt
"""
def get_anno_seq(cr_out):    
    airr = pd.read_csv(f'{cr_out}/airr_rearrangement.tsv',sep='\t')
    contig = pd.read_csv(f'{cr_out}/filtered_contig_annotations.csv')
    clonotype = pd.read_csv(f'{cr_out}/clonotypes.csv')
    
    # cr 4.0
    if str(airr.clone_id[0]) =='nan':
        clonotype['idnum'] = clonotype['clonotype_id'].apply(lambda x: int(x[9:]))
        clonotype.sort_values(by='idnum', inplace=True)
        del clonotype['idnum']
    del clonotype['proportion']
    clonotype['cdr3s_aa'] = clonotype['cdr3s_aa'].apply(lambda x:x.split(';'))
    clonotype['cdr3s_nt'] = clonotype['cdr3s_nt'].apply(lambda x:x.split(';'))
    
    aa_list = clonotype['cdr3s_aa'].tolist()
    aa_dict = {}
    for i,j in enumerate(aa_list,1):
        aa_dict[i] = j
    
    IGH_cdr3 = dict()
    IGL_cdr3 = dict()
    for key,value in aa_dict.items():
        if len(value)==2: 
            for i in value:
                if i.startswith('IGH'):
                    IGH_cdr3['clonotype'+str(key)]= (i.split(':')[-1])
                else:
                    IGL_cdr3['clonotype'+str(key)]= (i.split(':')[-1])
    df_IGH_cdr3 = pd.DataFrame(list(IGH_cdr3.items()),columns=['clonotype_id', 'IGH_cdr3'])
    df_IGL_cdr3 = pd.DataFrame(list(IGL_cdr3.items()),columns=['clonotype_id', 'IGLorK_cdr3'])

    df_merge = pd.merge(clonotype,df_IGH_cdr3,on='clonotype_id',how='inner')
    df_merge = pd.merge(df_merge, df_IGL_cdr3,on='clonotype_id',how='inner')
    del df_merge['cdr3s_aa']
    
    nt_list = clonotype['cdr3s_nt'].tolist()
    nt_dict = {}
    for i,j in enumerate(nt_list,1):
        nt_dict[i] = j
    IGH_cdr3_nt = dict()
    IGL_cdr3_nt = dict()
    for key,value in nt_dict.items():
        for i in value:
            if i.startswith('IGH'):
                IGH_cdr3_nt['clonotype'+str(key)]=(i.split(':')[-1])
            else:
                IGL_cdr3_nt['clonotype'+str(key)]=(i.split(':')[-1])

    df_IGH_cdr3_nt = pd.DataFrame(list(IGH_cdr3_nt.items()),columns=['clonotype_id', 'IGH_cdr3_nt'])
    df_IGL_cdr3_nt = pd.DataFrame(list(IGL_cdr3_nt.items()),columns=['clonotype_id', 'IGLorK_cdr3_nt'])
    df_merge = pd.merge(df_merge,df_IGH_cdr3_nt,on='clonotype_id',how='inner')
    df_merge = pd.merge(df_merge,df_IGL_cdr3_nt,on='clonotype_id',how='inner')
    del df_merge['cdr3s_nt']
    
    df_merge['Clonotype'] = df_merge['IGH_cdr3'] + '_' + df_merge['IGLorK_cdr3']
    df_merge =df_merge[['clonotype_id','Clonotype','frequency','IGH_cdr3','IGH_cdr3_nt','IGLorK_cdr3','IGLorK_cdr3_nt']]

    airr = airr[['clone_id','sequence_id','sequence','sequence_aa']]
    airr.rename(columns={'sequence_id':'contig_id'},inplace=True)
    
    # cr4.0
    contig = contig[contig['productive']==True]
    if str(airr.clone_id[0]) =='nan':
        airr = pd.merge(airr,contig[['contig_id','raw_clonotype_id']],on='contig_id',how='inner')
        airr.clone_id = airr.raw_clonotype_id
        del airr['raw_clonotype_id']
        
    contig.sort_values('reads',ascending=False,inplace=True)
    contig= contig[['contig_id','chain']]
    df_merge1 = pd.merge(airr,contig,on='contig_id')
    df_merge1.rename(columns={'clone_id':'clonotype_id'},inplace=True)
    df_merge1['idnum'] = df_merge1['clonotype_id'].apply(lambda x: int(x[9:]))
    df_merge1.sort_values(['idnum','chain'],inplace=True)
    del df_merge1['idnum']
    
    target_contigid = df_merge['clonotype_id'].tolist()
    df_merge1 = df_merge1[df_merge1['clonotype_id'].isin(target_contigid)]
    df_merge1.drop_duplicates(['clonotype_id','chain'],inplace=True)
    df_merge1_IGH = df_merge1[df_merge1['chain']=='IGH']
    df_merge1_IGL = df_merge1[df_merge1['chain']!='IGH']
    del df_merge1_IGH['chain']
    del df_merge1_IGL['chain']
    df_merge1_IGH.rename(columns={'contig_id':'contig_id_IGH'},inplace=True)
    df_merge1_IGL.rename(columns={'contig_id':'contig_id_IGLorKs'},inplace=True)
    df_merge1_IGH.rename(columns={'sequence':'IGH_nt'},inplace=True)
    df_merge1_IGL.rename(columns={'sequence':'IGLorK_nt'},inplace=True)
    df_merge1_IGH.rename(columns={'sequence_aa':'IGH_aa'},inplace=True)
    df_merge1_IGL.rename(columns={'sequence_aa':'IGLorK_aa'},inplace=True)
    
    df_merge_f = pd.merge(df_merge,df_merge1_IGH,on='clonotype_id')
    df_merge_f = pd.merge(df_merge_f,df_merge1_IGL,on='clonotype_id')
    df_merge_f.rename(columns={'frequency':'Frequency'},inplace=True)
    df_merge_f = df_merge_f[['clonotype_id','Clonotype','Frequency','IGH_cdr3','IGH_cdr3_nt','contig_id_IGH','IGH_aa','IGH_nt',
                            'IGLorK_cdr3','IGLorK_cdr3_nt','contig_id_IGLorKs','IGLorK_aa','IGLorK_nt']]
    df_merge_f.to_csv('./filtered_contig_annotations_pair_productive_with_sequence.tsv',sep='\t')


if __name__ == '__main__':
    cr_out = sys.argv[1]
    get_anno_seq(cr_out)


