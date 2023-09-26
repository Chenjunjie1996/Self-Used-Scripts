from celescope.tools.matrix import CountMatrix
import pandas as pd
import numpy as np
import argparse
import glob
import os
import functools


def get_metrics_dict(mtx_path, doublet_threshold):
    mtx = CountMatrix.from_matrix_dir(mtx_path)
    m = mtx.get_matrix()
    f = mtx.get_features()

    human = []
    mouse = []
    for i,gn in enumerate(f.gene_name):
        if gn == gn.upper():
            human.append(i)
        else:
            mouse.append(i)
    umi_h = m.tocsr()[human,:].sum(axis=0)
    umi_m = m.tocsr()[mouse,:].sum(axis=0)
    df_m = pd.DataFrame(umi_m.T,columns=['mouse'])
    df_h = pd.DataFrame(umi_h.T,columns=['human'])
    df = pd.merge(df_m,df_h,left_index=True,right_index=True)
    df['umi_sum'] = df.sum(axis=1)
    df['human_fraction'] = df['human'] / df['umi_sum'] * 100
    df['mouse_fraction'] = df['mouse'] / df['umi_sum'] * 100
    df['identity'] = 'doublet'
    df.loc[df['human_fraction'] < doublet_threshold,'identity'] = 'mouse'
    df.loc[df['mouse_fraction'] < doublet_threshold,'identity'] = 'human'
    df.loc[df['identity'] == 'human','ambient_fraction'] = df[df['identity'] == 'human'].mouse_fraction
    df.loc[df['identity'] == 'mouse','ambient_fraction'] = df[df['identity'] == 'mouse'].human_fraction

    n_cell = df.shape[0]
    n_doublet = df[df['identity'] == 'doublet'].shape[0]
    doublet_fraction = round(n_doublet / n_cell * 100,2)
    dict = {
        'n_cell': n_cell,
        'n_human_cell': df[df['identity'] == 'human'].shape[0],
        'n_mouse_cell': df[df['identity'] == 'mouse'].shape[0],
        'doublet_umi_fraction_threshold': doublet_threshold,
        'n_doublet': n_doublet,
        'doublet_cell_fraction': doublet_fraction,
        'median_ambient_umi_fraction': round(df['ambient_fraction'].median(),2),
        'mean_ambient_umi_fraction': round(df['ambient_fraction'].mean(),2),
        'median_umi_human': np.median(df[df['identity'] == 'human'].human) ,
        'median_umi_mouse': np.median(df[df['identity'] == 'mouse'].mouse)
    }

    return dict

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-c','--celescope_dir', type = str, required = True)
    parser.add_argument('--doublet_threshold', type=float, default=25.0)
    parser.add_argument('--out_prefix', type=str, default='human_mouse')
    args = parser.parse_args()

    cur_dir = os.getcwd()
    dfs = []
    for d in args.celescope_dir.split(','):
        os.chdir(d)
        paths = glob.glob(f'*/*count/*filtered_feature_bc_matrix')
        for path in paths:
            sample = path.split('/')[0]
            dict = get_metrics_dict(path, args.doublet_threshold)
            dfs.append(pd.DataFrame.from_dict(dict,orient='index',columns=[sample]))

    df = functools.reduce(lambda left,right: pd.merge(left,right,left_index=True,right_index=True), dfs)
    os.chdir(cur_dir)
    df.to_csv(f'{args.out_prefix}.tsv',sep='\t')


if __name__ == '__main__':
    main()