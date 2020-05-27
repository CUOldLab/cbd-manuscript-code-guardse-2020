import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns




def combine(paths,times):
	dfs = []; cntr = 0
	for path in paths:
		df = pd.read_csv(path,sep="\t")
		df = df[(df["ID"].str.contains("Inhibited")) | (df["ID"].str.contains("Activated"))]
		df.index=df["ID"].values
		df = df[["Gene_Count"]]
		df.columns = ["Gene_Set_Count"+times[cntr]]
		dfs.append(df)
		cntr+=1
	merged = pd.concat(dfs,axis=1)
	zcols = [c for c in merged.columns if ("Count") in c]
	merged = merged[merged[zcols].apply(np.abs).apply(np.max, axis=1)>=2]
	merged = merged[merged.isnull().sum(axis=1)<2]
	merged.insert(0,"ID",merged.index.values)
	return(merged)
	




directory = "merged/"
times = ["3h","6h","12h","24h"]
combined_up = combine(["merged/3h_up.csv","merged/6h_up.csv","merged/12h_up.csv","merged/24h_up.csv"],times)
combined_down = combine(["merged/3h_down.csv","merged/6h_down.csv","merged/12h_down.csv","merged/24h_down.csv"],times)
combined_up.to_csv("results/up_counts.txt",sep="\t",index=False)
combined_down.to_csv("results/down_counts.txt",sep="\t",index=False)

