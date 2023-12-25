import pysam
import pandas as pd
import glob


def main():
    clonotypes = glob.glob("./*/04.summarize/clonotypes.csv")
    annos = glob.glob("./*/04.summarize/filtered_contig_annotations.csv")
    fastas = glob.glob("./*/04.summarize/filtered_contig.fasta")
    
    for i in range(len(clonotypes)):
        sample = clonotypes[i].split('/')[-3]
    
        cl = pd.read_csv(clonotypes[i])
        cl.cdr3s_aa = cl.cdr3s_aa.apply(lambda x: x.split(';'))
        cl = cl[cl['cdr3s_aa'].apply(lambda x: len(x)) == 2]
        cl = cl[cl['cdr3s_aa'].apply(lambda x: x[0]).str.contains('IGH')]
        cl = cl[~cl['cdr3s_aa'].apply(lambda x: x[1]).str.contains('IGH')]
        cl['cdr3s_aa'] = cl['cdr3s_aa'].apply(lambda x: ';'.join(x))
        cl = cl[['clonotype_id', 'cdr3s_aa', 'cdr3s_nt']]
        
        anno = pd.read_csv(annos[i])
        anno = anno[['barcode' ,'contig_id', "raw_clonotype_id", "chain", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3"]]
    
        s1, s2 = [], []
        with pysam.FastxFile(fastas[i]) as f:
            for read in f:
                s1.append(read.name)
                s2.append(read.sequence)
        fl_seq_dict = dict(zip(s1, s2))
    
        anno['fl_seq'] = anno['contig_id'].apply(lambda x: fl_seq_dict[x])
        anno['t'] = anno['raw_clonotype_id'].apply(lambda x: x[9:])
        anno['t'] = anno['t'].astype(int)
        anno = anno.sort_values(['t','barcode'])
        del anno['t']
        anno = anno.rename(columns = {'raw_clonotype_id': 'clonotype_id'})
    
        df_merge = pd.merge(cl, anno, on='clonotype_id', how="inner")
        df_merge.to_csv(f"./{sample}.tsv", sep='\t', index=False)


if __name__ == '__main__':
    main()