import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.gridspec as gridspec


from matplotlib import cm
colormap = cm.bwr
colormap.set_bad(color='lightgray')


plt.rc('xtick', labelsize=8) 
plt.rc('ytick', labelsize=8) 


df = pd.read_csv("results/heatmap_input.csv")
df.dropna(subset=[c for c in df.columns if("GO:" in c)], inplace=True)


gocols = [c for c in df.columns if ("GO" in c)]

cntr=1
for col in gocols:
    df[col].replace(1,cntr,inplace=True)
    df[col].replace(0,np.nan,inplace=True)
    cntr+=1

df.sort_values(by=[c for c in df.columns if ("GO:" in c)],inplace=True)
df.sort_values(by="cholesterol biosynthetic process GO:0006695",inplace=True)

df.index = df['gene'].values
df.drop("gene",inplace=True,axis=1)

go_df = df[[c for c in df.columns if ("GO:" in c)]]
cyt = df[[c for c in df.columns if ("_soluble" in c)]]
nuc = df[[c for c in df.columns if ("_insoluble1" in c)]]
pm = df[[c for c in df.columns if ("_insoluble2" in c)]]
rna = df[[c for c in df.columns if ("_rna" in c)]]


fig = plt.figure(figsize=(8,14))

ax = fig.add_axes([0.15,0.1,.2,.8])
ax2 = fig.add_axes([0.36,0.1,.15,.8])
ax3 = fig.add_axes([0.52,0.1,.15,.8])
ax4 = fig.add_axes([0.68,0.1,.15,.8])
ax5 = fig.add_axes([0.84,0.1,.067,.8])



sns.heatmap(go_df, cmap="tab20", ax=ax, square=False, cbar=False, xticklabels=True, yticklabels=True, linewidths=0.1, linecolor='black')

for _, spine in ax.spines.items():
    spine.set_visible(True)


cyt.replace(0, np.nan, inplace=True)
nuc.replace(0, np.nan, inplace=True)
pm.replace(0, np.nan, inplace=True)
rna.replace(0, np.nan, inplace=True)

sns.heatmap(cyt, vmin=-1.01, vmax=1.01, cmap=colormap, ax=ax2, square=False, cbar=False, xticklabels=True, yticklabels=False, linewidths=0.05, linecolor='black')
sns.heatmap(nuc, vmin=-1.01, vmax=1.01, cmap=colormap, ax=ax3, square=False, cbar=False, xticklabels=True, yticklabels=False, linewidths=0.05, linecolor='black')
sns.heatmap(pm, vmin=-1.01, vmax=1.01, cmap=colormap, ax=ax4, square=False, cbar=False, xticklabels=True, yticklabels=False, linewidths=0.05, linecolor='black')
sns.heatmap(rna, vmin=-1.01, vmax=1.01, cmap=colormap, ax=ax5, square=False, cbar=True, cbar_kws={"ticks":[-1,-0.5,0,0.5,1]}, xticklabels=True, yticklabels=False, linewidths=0.05, linecolor='black')
#fig.colorbar(ax5.collections[0], ax=ax2,location="right", use_gridspec=False, pad=0.2)

ax.tick_params(axis='x',rotation=90)
axes = [ax2,ax3,ax4,ax5]
for a in axes:
    a.yaxis.set_ticklabels([])
    a.set_yticks([])
    a.tick_params(axis='x',rotation=90)
    for _, spine in a.spines.items():
        spine.set_visible(True)

plt.savefig("results/enrichment_figure_20190414_gobio2018_lfc1.0.svg", bbox_inches="tight")

