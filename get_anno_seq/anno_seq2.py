import pandas as pd
import json
import sys
import glob

"""
cell_id	clone_id	sequence_id	chain	sequence	sequence_aa
"""
def get_anno_seq(cr_out, seqtype):
    with open(f"{cr_out}/02.convert/barcode_convert.json", 'r') as f:
        tenX_sgr = json.load(f)
    
    airr = pd.read_csv(glob.glob(f"{cr_out}/03.assemble/*/outs/airr_rearrangement.tsv", sep='\t')[0])
    df = pd.read_csv(glob.glob(f"{cr_out}/03.assemble/*/outs/filtered_contig_annotations.csv")[0])
    
    bc_chain_dict = df.groupby("barcode")["chain"].apply(set).to_dict()
    if seqtype == "TCR":
        pair_chain_dict = {key: value for key, value in bc_chain_dict.items() if len(value)==2}
    else:
        pair_chain_dict = {key: value for key, value in bc_chain_dict.items() if len(value)>=2 and "IGH" in value}
    df = df[df["barcode"].isin(pair_chain_dict)]
    df_pair_chain = df.groupby("barcode").filter(lambda x: (len(x) == 2))
    
    airr = airr[airr["cell_id"].isin(set(df_pair_chain.barcode))]
    airr["cell_id"] = airr["cell_id"].apply(lambda x: tenX_sgr[x.split('-')[0]])
    airr["sequence_id"] = airr["sequence_id"].apply(lambda x: tenX_sgr[x.split('-')[0]]+'_'+x.split('_')[1]+'_'+x.split('_')[2])
    airr["chain"] = airr["j_call"].apply(lambda x: x[0:3])
    airr = airr[["cell_id","clone_id","sequence_id","chain", "sequence","sequence_aa"]]
    
    airr["temp"] = airr["clone_id"].apply(lambda x: int(x[9:]))
    airr = airr.sort_values(["temp","cell_id", "chain"])
    del airr["temp"]
    
    airr.to_csv(f"./full_len_seq.csv", sep=',', index=False)

if __name__ == '__main__':
    cr_out = sys.argv[1]
    seqtype = sys.argv[2]
    get_anno_seq(cr_out, seqtype)