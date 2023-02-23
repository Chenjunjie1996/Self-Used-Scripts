import pandas as pd
import plotly.express as px

file = "path/sample-matched_clonotypes.csv"
df = pd.read_csv(file)
sample = '-'.join(file.split('/')[-1].split('-')[:2])

fig = px.pie(df.head(10), values='frequency', names='cdr3s_aa', title=f'{sample}-T Top10 Clonotypes')

fig.update_traces(hoverinfo='label+percent', textinfo='percent', textfont_size=20,
                  marker=dict(line=dict(color='#000000', width=2)))

fig.show()