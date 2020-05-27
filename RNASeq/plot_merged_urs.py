import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns





times = ["3h","6h","12h","24h"]
df = pd.read_csv("significant_URs_counts.csv",sep="\t")

new_ids = []
new_colors = []
for idx,row in df.iterrows():
	print row
	split=row["ID"].split("_")
	new_ids.append(split[0]+"_"+split[2][0])
	if row["geneset_dir"]=="up":
		new_colors.append("b")
	else:
		new_colors.append("r")
fig, ax = plt.subplots(figsize=(5,12))
df.index = new_ids
df.drop(["ID","geneset_dir"],axis=1,inplace=True)
df.columns = times
df.replace(np.nan,0,inplace=True)
ax = sns.heatmap(df, linewidths=1,annot=True,vmin=-8, vmax=8, cbar=False)
plt.yticks(fontsize=16,rotation=0) 
plt.xticks(fontsize=16,rotation=45) 
plt.xlabel("Gene Set Count",fontsize=18)
plt.ylabel("Upstream Regulator",fontsize=18)
plt.tight_layout()


[t.set_color(i) for (i,t) in
 zip(new_colors,ax.yaxis.get_ticklabels())]
 
plt.savefig("figures/IPA_reduced_counts.svg",dpi=1000)
