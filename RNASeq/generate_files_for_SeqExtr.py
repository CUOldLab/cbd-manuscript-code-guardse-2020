import pandas as pd


combined = pd.read_csv("results/down_zscores.txt",sep="\t")
df = pd.read_csv("results/gene_set_down.txt",sep="\t")
df_ids = df["ID"].values
short_ids=[]
for i in df_ids:
	short_ids.append(i.split("_")[0])
df["short_name"] = short_ids


ids = combined["ID"].values
short_merged_ids=[]
for i in ids:
	short_merged_ids.append(i.split("_")[0])




for idx,row in df.iterrows():
	if row["short_name"] in short_merged_ids:
		UR = row["short_name"]
		gene_set = row["gene_set"].split(",")
		genes=[]
		for g in gene_set:
			genes.append(g)
		dfout = pd.DataFrame(genes,columns=["genes in set"])
		dfout.to_csv("motif_data/down/"+UR+".txt",sep="\t",index=False,header=False)


