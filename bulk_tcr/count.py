import glob
import pandas as pd


class Count:
    """
    Count types, frequencies, diversity of different index.
    """
    def __init__(self):
        self.clonotypes = glob.glob("*/05.count_vdj/*_clonetypes.csv")
        self.out_file = "clonotypes_info.txt"
        
    def __call__(self):
        
        df_final = pd.DataFrame(columns=["Sample", "Index", "Count", "Types", "Diversity"])
        
        for file in self.clonotypes:
            df = pd.read_csv(file)
            df = df.groupby("Index").agg({'Frequency':'sum','ClonotypeID':'count','Diversity':'min'})\
                .reset_index().rename(columns={'Frequency':'Count','ClonotypeID':'Types'})\
                    .sort_values("Diversity", ascending=False)
            df["Sample"] = file.split('/')[0]
            df_final = pd.concat([df_final, df])

        df_final.to_csv(self.out_file, sep='\t', index=False)


if __name__ == '__main__':
    Count()()