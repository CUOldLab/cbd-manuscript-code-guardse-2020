import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(font_scale=.6)

def prune(df, frac):
    df = df[["gene"]+[c for c in df.columns if (c.endswith("_Corr"))]+["q_values"]]
    #rename columns
    df.columns = ["gene"]+[c+"_"+frac for c in df.columns[1:]]
    return df



df_cyt = pd.read_csv("/home/jjacob/R/edge/proteome_diff_analysis/20180827/results/for_kerri/CYT_13tp.txt")
df_nuc = pd.read_csv("/home/jjacob/R/edge/proteome_diff_analysis/20180827/results/for_kerri/NUC_13tp.txt")
df_pm = pd.read_csv("/home/jjacob/R/edge/proteome_diff_analysis/20180827/results/for_kerri/PM_13tp.txt")

df_cyt = prune(df_cyt, "soluble")
df_nuc = prune(df_nuc, "insoluble1")
df_pm = prune(df_pm, "insoluble2")

merged = df_cyt.merge(df_nuc, how='outer', left_on="gene", right_on="gene")
merged = merged.merge(df_pm, how='outer', left_on="gene", right_on="gene")

merged.to_csv("merged_proteins.csv", index=False)
merged = merged[(merged["q_values_soluble"]<0.01) | (merged["q_values_insoluble1"]<0.01) | (merged["q_values_insoluble2"]<0.01)]
merged.replace(-0.0001, np.nan, inplace=True)
merged.index = merged["gene"]
merged.to_csv("merged_significant_proteins_1perc_fdr.csv", index=False)
merged.drop(["gene", "q_values_soluble", "q_values_insoluble1", "q_values_insoluble2"], axis=1, inplace=True)
merged.replace([-np.inf,np.inf,np.nan], 0.0, inplace=True)

ax = sns.clustermap(merged, col_cluster=False, cmap='seismic', square=True, figsize=(15,60), vmin=-1, vmax=1)
ax.ax_row_dendrogram.set_visible(False)

plt.tight_layout()
plt.savefig("CBD_proteome_clustergram.svg")
plt.close()








