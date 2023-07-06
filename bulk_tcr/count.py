import glob
import pandas as pd


class Count:
    """
    Count types, frequencies, diversity of different index.
    """
    def __init__(self):
        self.df_clonotypes = glob.glob("*/05.count_vdj/*_clonetypes.csv")[0]
        self.out_file = "clonotypes_info.txt"
    
    def __call__(self):
        df = pd.read_csv(self.df_clonotypes)
        
        df.groupby("Index").agg({'Frequency':'sum','ClonotypeID':'count','Diversity':'min'})\
            .reset_index().rename(columns={'Frequency':'Count','ClonotypeID':'Types'})\
                .sort_values("Diversity", ascending=False)\
                    .to_csv(self.out_file, sep='\t', index=False)


if __name__ == '__main__':
    Count()()