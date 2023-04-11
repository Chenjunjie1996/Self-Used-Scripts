import pandas as pd
import argparse
import plotly.express as px
import glob
import os
import seaborn as sns
import matplotlib.pyplot as plt


def check_mkdir(dir_name):
    """if dir_name is not exist, make one"""
    if not os.path.exists(dir_name):
        os.system(f"mkdir -p {dir_name}")
        

class Plot:
    
    def __init__(self, args):
        self.path = args.path
        self.df_list = glob.glob(f"{self.path}/*/05.count*/*_clonetypes.csv")
        self.sample_list = [i.split('/')[-1].split('_')[0] for i in self.df_list]

        self.df_list = sorted(self.df_list, key=lambda x: int(x.split('/')[-1].split('_')[0].split("index")[-1].rstrip('K')))
        self.sample_list = sorted(self.sample_list, key=lambda x: int(x.split('index')[-1].rstrip('K')))
        self.sample_list = [f"index_{i.split('index')[-1].rstrip('K')}" for i in self.sample_list]
        

    def __call__(self):
        check_mkdir(f"{self.path}/plot_index")
        self.imshow()
        self.barplot()

    def imshow(self):
        number_list = [pd.read_csv(df).shape[0] for df in self.df_list]
        df = pd.DataFrame({"sample":self.sample_list, "clonotypes": number_list})
        data = [df["clonotypes"].tolist()[i:i+8] for i in range(0, len(df), 8)]
        
        fig = px.imshow(data, text_auto=True, aspect="auto", title="96 Index", labels=dict(color="Clonotypes number"))
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)
        fig.write_image(f"{self.path}/plot_index/heatmap.pdf")

    def barplot(self):
        # 读入96个clonotype.csv文件并合并成一个数据框
        df_all = []
        for i in range(len(self.df_list)):
            df = pd.read_csv(self.df_list[i])
            df["Sample"] = self.sample_list[i]
            df_all.append(df)
        
        df_merge = pd.concat(df_all)
        pivot_table = df_merge.pivot_table(index="ClonotypeID", columns="Sample", values="Frequency")
        pivot_table = pivot_table.reindex(columns=self.sample_list)

        # 绘制热力图
        f, ax = plt.subplots(figsize=(30, 10))
        ax = sns.heatmap(pivot_table, cmap="RdBu_r", center=100, vmax=1000)
        plt.savefig(f"{self.path}/plot_index/barplot.pdf", dpi=300, bbox_inches="tight")
        #plt.show()
        
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="split index for bulk vdj")
    parser.add_argument("--path", help="analysis dir", required=True)
    args = parser.parse_args()
    Plot(args)()