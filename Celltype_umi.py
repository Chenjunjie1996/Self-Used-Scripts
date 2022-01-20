import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys
import subprocess

def run_rds(rds):
    cmd = (
        f'Rscript /SGRNJ03/randd/cjj/Script/extratmeta.R '
        f'--rds {rds} '
        f'--outdir {outdir}'
    )
    subprocess.check_call(cmd, shell=True)

def parse_file(meta,contig_file):
    sample_name = os.path.basename(os.path.abspath(contig_file)).split("_")[0]
    contig = pd.read_csv(contig_file)

    meta.rename(columns={'Unnamed: 0': 'barcode'}, inplace=True)
    meta = meta[meta['barcode'].apply(lambda x: x.split('_')[0] == f'{sample_name}')]
    meta['barcode'] = meta['barcode'].apply(lambda x: x.split('_')[1])
    meta = meta[['barcode', 'new_ident']]
    pdmerge = pd.merge(contig, meta, on='barcode', how='left')
    
    plt.xticks(rotation=60)
    snsFig = sns.violinplot(x="umis", y="new_ident", data=pdmerge,
                            linewidth=1, 
                            width=0.8, 
                            palette='muted',        
                            inner='box'
                            )
    snsFig = snsFig.get_figure()
    snsFig.savefig("./Vlnplot.png", bbox_inches='tight', dpi=300)

if __name__ == '__main__':
    outdir = "./"
    rds = sys.argv[1]
    contig_file = sys.argv[2]

    run_rds(rds)
    meta = pd.read_csv("./meta.csv", sep=',')
    parse_file(meta, contig_file)

