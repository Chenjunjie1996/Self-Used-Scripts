# 计算双链的细胞数
import sys
import glob
import os
import pandas as pd


def main(input_dir, seqtype):
    annotation_file = glob.glob(f"{input_dir}/05.match/matched_productive_contig_annotations.csv")
    
    if not annotation_file:
        raise FileNotFoundError("File matched_productive_contig_annotations.csv not found")
    
    df = pd.read_csv(annotation_file[0])
    
    if seqtype == "BCR":
        df_chain_heavy = df[(df['chain'] == 'IGH')] 
        df_chain_light = df[(df['chain'] == 'IGL') | (df['chain'] =='IGK')]
    else:
        df_chain_heavy = df[df['chain'] == 'TRA']
        df_chain_light = df[df['chain'] == 'TRB']
    
    df_chain_heavy = df_chain_heavy.groupby("barcode").filter(lambda x: (len(x) == 1))
    df_chain_light = df_chain_light.groupby("barcode").filter(lambda x: (len(x) == 1))
    heavy_light_pair = len( set(df_chain_heavy.barcode) & set(df_chain_light.barcode) )
    
    match_cells = len(set(df.barcode))
    one_chain_cells = df.groupby("barcode").filter(lambda x: (len(x) == 1))
    one_chain_cells = len(set(one_chain_cells.barcode))
    three_chain_cells = df.groupby("barcode").filter(lambda x: (len(x) >= 3))
    three_chain_cells = len(set(three_chain_cells.barcode))
    sample_name = os.path.abspath(input_dir).split('/')[-1]
    
    print(
        f"Sample : {sample_name}\n"
        f"Total match cells : {match_cells}\n"
        f"Cells with single chain : {one_chain_cells}\n"
        f"Cells with more than 3 chain : {three_chain_cells}\n"
        f"Cells with single heavy-light pair chain : {heavy_light_pair}"
    )


if __name__ == '__main__':
    input_dir, seqtype = sys.argv[1], sys.argv[2]
    main(input_dir, seqtype)