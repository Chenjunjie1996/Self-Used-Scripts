import pandas as pd
from pandas import Series
import glob
import plotly
import plotly.express as px  # import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns
import matplotlib.pyplot as plt

import plotly.io as pio
#pio.templates.default = 'plotly_white'


def calc_corr(a,b):
    s1 = Series(a)
    s2 = Series(b)
    return round(s1.corr(s2), 3)


# human tcr
df = pd.read_csv("/SGRNJ06/randd/USER/cjj/celedev/vdj_bulk/20230727plot/9/9.csv")
df = df.dropna(axis=0)
df


df["Raw_Reads"] = df["Raw_Reads"].str.replace(',','').astype(int)
df["UMI_Counts"] = df["UMI_Counts"].str.replace(',','').astype(int)
df["Clonotype_Diversity"] = df["Clonotype_Diversity"].astype(str).str.replace(',','').astype(float)
df["UMIs_Mapped_Confidently_To_VJ_Gene"] = df["UMIs_Mapped_Confidently_To_VJ_Gene"].apply(lambda x: x.split('(')[-1]).str.replace(')','').str.replace("%",'').astype(float)
df["UMIs_Mapped_To_TRA"] = df["UMIs_Mapped_To_TRA"].apply(lambda x: x.split('(')[-1]).str.replace(')','').str.replace("%",'').astype(float)
df["UMIs_Mapped_To_TRB"] = df["UMIs_Mapped_To_TRB"].apply(lambda x: x.split('(')[-1]).str.replace(')','').str.replace("%",'').astype(float)
df["species"] = "TCR"
df["index"] = df.index
sample = "Mus_0620Spleen"
OUT_PATH = "/SGRNJ06/randd/USER/cjj/celedev/vdj_bulk/20230727plot/9"


# 1 reads-umi 相关性散点图
fig = px.scatter(
    df, # scatter绘制散点图
    x='Raw_Reads',    #  x轴
    y='UMI_Counts',  # y轴
    #color = 'UMI_Counts',
    labels = calc_corr(df["Raw_Reads"], df["UMI_Counts"]),
    width=1680, height=1024
)
fig.update_traces(marker_size=10)
fig.update_layout(title={"text":f'{sample} (Correlation={calc_corr(df["Raw_Reads"], df["UMI_Counts"])})',
                         'y':0.98, 'x':0.5, 'xanchor': 'center', 'yanchor': 'top'},
                  xaxis_title = "Raw_Reads", yaxis_title = "UMI_Counts", font=dict(size=18,color="Black"))

fig.write_image(f"{OUT_PATH}/scatter.png", scale=4)
#fig.show()


# 2 UMI mapping to AnyVdj/TRA/TRB 箱线图
fig = go.Figure()
fig.add_trace(go.Box(y=list(df.UMIs_Mapped_To_TRA),name="UMIs_Mapped_To_TRA"))
fig.add_trace(go.Box(y=list(df.UMIs_Mapped_To_TRB),name="UMIs_Mapped_To_TRB"))
fig.add_trace(go.Box(y=list(df.UMIs_Mapped_Confidently_To_VJ_Gene),name="UMIs_Mapped_Confidently_To_VJ_Gene"))
fig.update_layout(title={"text":f'{sample}',
                         'y':0.98, 'x':0.40, 'xanchor': 'center', 'yanchor': 'top'},
                          yaxis_title = "Mapping Percent(%)", font=dict(size=18,color="Black"))

fig.write_image(f"{OUT_PATH}/boxplot.png", scale=4, width=1680, height=1024)
#fig.show()


# 3 reads-diversity 相关性散点图
fig = px.scatter(
    df, # scatter绘制散点图
    x='Raw_Reads',    #  x轴
    y='Clonotype_Diversity',  # y轴
    #color = 'Clonotype_Diversity',
    labels = calc_corr(df["Raw_Reads"], df["Clonotype_Diversity"]),
    width=1680, height=1024
)
fig.update_traces(marker_size=10)
fig.update_layout(title={"text":f'{sample} Diversity (Correlation={calc_corr(df["Raw_Reads"], df["Clonotype_Diversity"])})',
                         'y':0.98, 'x':0.5, 'xanchor': 'center', 'yanchor': 'top'},
                  xaxis_title = "Raw_Reads", yaxis_title = "Clonotype_Diversity", font=dict(size=18,color="Black"))

fig.write_image(f"{OUT_PATH}/diversity_read_scatter.png", scale=4)
#fig.show()


# 4 umi-diversity 相关性散点图
fig = px.scatter(
    df, # scatter绘制散点图
    x='UMI_Counts',    #  x轴
    y='Clonotype_Diversity',  # y轴
    #color = 'Clonotype_Diversity',
    labels = calc_corr(df["UMI_Counts"], df["Clonotype_Diversity"]),
    width=1680, height=1024
)
fig.update_traces(marker_size=10)
fig.update_layout(title={"text":f'{sample} Diversity (Correlation={calc_corr(df["UMI_Counts"], df["Clonotype_Diversity"])})',
                         'y':0.98, 'x':0.5, 'xanchor': 'center', 'yanchor': 'top'},
                  xaxis_title = "UMI_Counts", yaxis_title = "Clonotype_Diversity", font=dict(size=18,color="Black"))

fig.write_image(f"{OUT_PATH}/diversity_umi_scatter.png", scale=4)
#fig.show()


# 5. sample diversity柱状图
fig = px.bar(df, x="Sample", y="Clonotype_Diversity", color="Sample")
fig.update_layout(title={"text":f'{sample}',
                         'y':0.98, 'x':0.45, 'xanchor': 'center', 'yanchor': 'top'},
                          yaxis_title = "Clonotype_Diversity", font=dict(size=18,color="Black"))
fig.write_image(f"{OUT_PATH}/Diversity_barplot.png", scale=4, width=1680, height=1024)
#fig.show()