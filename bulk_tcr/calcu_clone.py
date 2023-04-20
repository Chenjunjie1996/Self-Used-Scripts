import glob
import pandas as pd


class Calcu_clone:
    """
    Calculate counts and kinds of clonotypes
    """
    def __init__(self):
        self.clonotypes_list = glob.glob("*/05.count_vdj/*_clonetypes.csv")
        self.clonotypes_list = sorted(
            self.clonotypes_list, key=lambda x: int(x.split('/')[-1].split('_')[0].split("index")[-1].rstrip('K'))
        )
        self.sample_list = [i.split('/')[-1].split('_')[0] for i in self.clonotypes_list]

        self.out_file = "clonotypes_info.txt"
    
    def __call__(self):
        type_list, count_list = [], []
        
        for i in self.clonotypes_list:
            df_tmp = pd.read_csv(i)
            type_list.append(df_tmp.shape[0])
            count_list.append(sum(df_tmp.Frequency))
        
        df = pd.DataFrame({"sample":self.sample_list , "types":type_list, "count":count_list})
        df.to_csv(self.out_file, sep='\t', index=False)


if __name__ == '__main__':
    Calcu_clone()()