import sys
import pandas as pd


def get_pairs_file(annotation_file, type):
    df = pd.read_csv(annotation_file)
    df['productive'].replace({'True': True, 'None': False}, inplace=True)
    df_productive = df[df['productive'] == True]

    if type == "BCR":
        df_heavy = df_productive[(df_productive['chain'] == "IGH")] 
        df_light = df_productive[(df_productive['chain'] == "IGL") | (df_productive['chain'] == "IGK")]
    else:
        df_heavy = df_productive[(df_productive['chain'] == "TRA")] 
        df_light = df_productive[(df_productive['chain'] == "TRB")]

    df_pair = pd.merge(df_heavy, df_light, on='barcode', how='inner')

    return df_pair


def get_pairs_clonotypes(df_match, out_file):

    df_match['chain_cdr3aa'] = df_match[['chain', 'cdr3']].apply(':'.join, axis=1)

    match_clonotypes = open(out_file, 'w')
    match_clonotypes.write('barcode\tcdr3s_aa\n')
    for cb in set(df_match.barcode):
        temp = df_match[df_match['barcode']==cb].sort_values(by='chain', ascending=True)
        chain_pair = ';'.join(temp['chain_cdr3aa'].tolist())
        match_clonotypes.write(f'{cb}\t{chain_pair}\n')
    match_clonotypes.close()

    df_match_clonetypes = pd.read_csv(out_file, sep='\t', index_col=None)
    df_match_clonetypes = df_match_clonetypes.groupby('cdr3s_aa', as_index=False).agg({'barcode': 'count'})
    df_match_clonetypes.rename(columns={'barcode': 'frequency'}, inplace=True)
    sum_f = df_match_clonetypes['frequency'].sum()
    df_match_clonetypes['proportion'] = df_match_clonetypes['frequency'].apply(lambda x: x/sum_f)
    df_match_clonetypes['clonotype_id'] = [f'clonotype{i}' for i in range(1, df_match_clonetypes.shape[0]+1)]
    df_match_clonetypes = df_match_clonetypes.reindex(columns=['clonotype_id', 'cdr3s_aa', 'frequency', 'proportion'])
    df_match_clonetypes.sort_values(by='frequency', ascending=False, inplace=True)
    df_match_clonetypes.to_csv(out_file, sep=',', index=False)


def main():
    annotation_file, type = sys.argv[1], sys.argv[2]
    out_annotation, out_pair_clonotypes = "productive_contig_annotations.csv", "out_pair_clonotypes"

    df_pair = get_pairs_file(annotation_file, type)
    df_pair.to_csv(out_annotation, sep=',', index=False)

    get_pairs_clonotypes(df_pair, out_pair_clonotypes)


if __name__ == '__main__':
    main()
