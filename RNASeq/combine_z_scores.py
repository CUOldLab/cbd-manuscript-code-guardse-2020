import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns



def return_overlapping_geneset(paths,times):
	dfs = []; cntr = 0
	gene2set = {}
	for path in paths:
		df = pd.read_csv(path,sep="\t")
		for dx,row in df.iterrows():
			UR = row["ID"]
			gene_set = row["Overlapping DS genes"].split(",")

			if UR in gene2set:
				gset = gene2set[UR]
				gset += gene_set
				gset = list(set(gene_set))
				gene2set[UR] = gset

			else:
				gene2set[UR] = gene_set
	gset_strings=[]
	for val in gene2set.values():
		gstr=""
		for g in val:
			gstr+=g.strip()+","
		gstr=gstr[0:len(gstr)-1]
		gset_strings.append(gstr)

	dfout = pd.DataFrame(gene2set.keys(),columns=["ID"])
	dfout["gene_set"] = gset_strings
	dfout=dfout[(dfout["ID"].str.contains("Inhibited"))|(dfout["ID"].str.contains("Activated"))]
	return(dfout)		



def combine(paths,times):
	dfs = []; cntr = 0
	overlapping_gene_set = []
	for path in paths:
		df = pd.read_csv(path,sep="\t")
		df.index=df["ID"].values
		df = df[["Mean Z-score"]]
		df.columns = ["Mean Z-score_"+times[cntr]]
		dfs.append(df)
		cntr+=1
	merged = pd.concat(dfs,axis=1)
	zcols = [c for c in merged.columns if ("Z-score") in c]
	merged = merged[merged[zcols].apply(np.abs).apply(np.max, axis=1)>=2]
	merged = merged[merged.isnull().sum(axis=1)<2]
	merged.insert(0,"ID",merged.index.values)
	return(merged)
	




directory = "merged/"
times = ["3h","6h","12h","24h"]
gene_set_up = return_overlapping_geneset(["merged/3h_up.csv","merged/6h_up.csv","merged/12h_up.csv","merged/24h_up.csv"],times)
gene_set_dn = return_overlapping_geneset(["merged/3h_down.csv","merged/6h_down.csv","merged/12h_down.csv","merged/24h_down.csv"],times)
gene_set_up.to_csv("results/gene_set_up.txt",sep="\t",index=False)
gene_set_dn.to_csv("results/gene_set_down.txt",sep="\t",index=False)
#combined_up = combine(["merged/3h_up.csv","merged/6h_up.csv","merged/12h_up.csv","merged/24h_up.csv"],times)
#combined_down = combine(["merged/3h_down.csv","merged/6h_down.csv","merged/12h_down.csv","merged/24h_down.csv"],times)


#combined_up.to_csv("results/up_zscores.txt",sep="\t",index=False)
#combined_down.to_csv("results/down_zscores.txt",sep="\t",index=False)

