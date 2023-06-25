import pandas as pd 
import glob


def count_multi_chain():
    df_list = glob.glob("./*/03.assemble/*/outs/all_contig_annotations.csv")
    sample_list = [i.split('/')[-3] for i in df_list]
    
    out_file = open("./multi_chains_cells.txt", 'w')
    for i in range(len(df_list)):
        df_tmp = pd.read_csv(df_list[i])
        df_productive = df_tmp[df_tmp["productive"]==True]
        df_cell = df_productive[df_productive["is_cell"]==True]
        total_cells = len(set(df_cell.barcode))
        print(sample_list[i])
        print(f"assembled cells : {total_cells}")

    
        df_not_cell = df_productive[df_productive["is_cell"]==False]
        df_not_cell = df_not_cell.groupby("barcode").filter(lambda x: (len(x) > 4))
        multi_chains = len(set(df_not_cell.barcode))
        print(f"multi productive chains: {multi_chains}\n")
    
        out_file.write(f"{sample_list[i]}\nassembled cells : {total_cells}\nmulti productive chains: {multi_chains}\n\n")
    out_file.close()
    

if __name__ == '__main__':
    count_multi_chain()