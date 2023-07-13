import glob
import pandas as pd


class Count:
    """
    Count types, frequencies, diversity of different index.
    """
    def __init__(self):
        self.df_clonotypes = glob.glob("*/05.count_vdj/*_clonetypes.csv")
    
    def __call__(self):
        
        for df in self.df_clonotypes:
            out_file = f"{df.split('/')[0]}_info.txt"
            df = pd.read_csv(df)
            
            df = df.groupby("Index").agg({'Frequency':'sum','ClonotypeID':'count','Diversity':'min'})\
                .reset_index().rename(columns={'Frequency':'Count','ClonotypeID':'Types'})\
                    .sort_values("Diversity", ascending=False)\
                        .to_csv(out_file, sep='\t', index=False)


if __name__ == '__main__':
    Count()()