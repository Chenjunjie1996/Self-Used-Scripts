import pandas as pd
import os
import glob


def main():
    
    dir_name = "split_annotations"
    if not os.path.exists(dir_name):
        os.system(f"mkdir -p {dir_name}")

    barcode_sample = pd.read_csv('barcode_sample.xls',sep='\t')
    barcode_sample_dict = barcode_sample.set_index('barcode')['sample'].to_dict()
    
    file = glob.glob("*/05.count_vdj/*_filtered_annotations.csv")[0]
    df = pd.read_csv(file)
    index_set = set(df["barcode"])
    
    for index in index_set:
        df1 = df[df["barcode"]==index]
        df1 = df1.sort_values("umis", ascending=False)
        df1 = df1.drop_duplicates(["chain", "cdr3"])
        df1["raw_clonotype_id"] =  [f"clonotype{i}" for i in range(1, len(df1)+1)]
        df1["barcode"] = df1["sequence_id"].apply(lambda x: x.split(":")[0])
        sample = barcode_sample_dict[df1.iloc[1]["barcode"]]
        df1["barcode"] = ["{}_{}".format(barcode, i + 1) for barcode, i in zip(df1["barcode"], range(len(df1)))]
        
        df1.to_csv(f"{dir_name}/{sample}_filtered_annotations.csv", sep=",", index=False)

if __name__ == '__main__':
    main()