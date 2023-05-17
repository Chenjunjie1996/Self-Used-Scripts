import glob
import os
import sys
import pandas as pd
import scanpy as sc


def split_barcode(rna_path):
    h5ad = glob.glob(f"{rna_path}/06.analysis/*.h5ad")[0]
    marker = glob.glob(f"{rna_path}/06.analysis/*_markers.tsv")[0]
    
    if not h5ad:
        raise FileNotFoundError("h5ad file not found")

    if not marker:
        raise FileNotFoundError("marker file not found")
    
    sample_name = h5ad.split('/')[-1].split('.')[0]
    
    h5ad = sc.read_h5ad(h5ad)
    marker = pd.read_csv(marker, sep='\t')
    
    marker = marker.sort_values(['cluster', 'avg_log2FC'], ascending=[True, False])
    marker = marker.groupby("cluster", as_index=False).head(20)
    marker = marker.groupby('cluster')['gene'].apply(lambda x: x.tolist()).to_dict()
    
    hs_cluster, mu_cluster = [], []
    
    for k, v in marker.items():
        count = 0
        for gene in v:
            if gene[1:].isupper():
                count += 1
            elif gene[1:].islower():
                count -= 1
    
        if count > 0:
            hs_cluster.append(str(k - 1))
        elif count < 0:
            mu_cluster.append(str(k - 1))
    
    hs_barcode = h5ad.obs[h5ad.obs["cluster"].isin(hs_cluster)]
    mu_barcode = h5ad.obs[h5ad.obs["cluster"].isin(mu_cluster)]
    hs_barcode = list(hs_barcode.index)
    mu_barcode = list(mu_barcode.index)
    
    print(
        f"sample name: {sample_name}\n"
        f"hs_barcode: {hs_cluster} ; {len(hs_barcode)}\n"
        f"mu_barcode: {mu_cluster} ; {len(mu_barcode)}\n"
    )

    for species in ["hs", "mu"]:
        out_dir = f"./{sample_name}_{species}/05.count/{sample_name}_matrix_10X"
        if not os.path.exists(out_dir):
            os.system(f"mkdir -p {out_dir}")
    
        if species == "hs":
            barcodes = hs_barcode
        else:
            barcodes = mu_barcode

        with open(f"{out_dir}/barcodes.tsv", 'w') as f:
            for i in barcodes:
                f.write(f"{i}\n")


if __name__ == '__main__':
    split_barcode(rna_path = sys.argv[1])
    
    