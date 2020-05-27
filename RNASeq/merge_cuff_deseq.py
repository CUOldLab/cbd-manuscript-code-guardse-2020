import pandas as pd
import numpy as np
import os



def merge(df1,df2):
	if "Flags" in df1.columns:
		df1.drop(["Flags"],axis=1,inplace=True)
	if "Flags" in df2.columns:
		df2.drop(["Flags"],axis=1,inplace=True)

	uni_index = df1["Upstream Regulator"]+"_"+df1["Molecule Type"]+"_"+df1["Predicted Activation State"]
	df1.index = uni_index.values
	df1.drop(["Upstream Regulator","Expr Log Ratio","Molecule Type","Predicted Activation State","p-value of overlap","Mechanistic Network"],axis=1,inplace=True)
	df1.columns = ["Activation Z-score cuff","Target molecules in dataset cuff"]

	uni_index = df2["Upstream Regulator"]+"_"+df2["Molecule Type"]+"_"+df2["Predicted Activation State"]
	df2.index = uni_index.values
	df2.drop(["Upstream Regulator","Expr Log Ratio","Molecule Type","Predicted Activation State","p-value of overlap","Mechanistic Network"],axis=1,inplace=True)
	df2.columns = ["Activation Z-score deseq","Target molecules in dataset deseq"]

	merged = pd.concat([df1,df2],axis=1)
	merged.insert(0,"ID",merged.index.values)
	merged['Activation Z-score cuff']=pd.to_numeric(merged['Activation Z-score cuff'],errors='coerce')
	merged['Activation Z-score deseq']=pd.to_numeric(merged['Activation Z-score deseq'],errors='coerce')
	return merged
	
	

def overlap_merge(df):
	

	df = df[~df["Activation Z-score cuff"].isnull()]
	df = df[~df["Activation Z-score deseq"].isnull()]
	
	if df.shape[0] == 0:
		return df
	inter = []
	gene_count = []
	for idx,row in df.iterrows():
		tlist1 = row["Target molecules in dataset cuff"].split(",")
		tlist2 = row["Target molecules in dataset deseq"].split(",")
		intersection = np.intersect1d(tlist1,tlist2)
		gene_count.append(len(intersection))
		
		
		ident_str=""
		for ident in intersection:
			ident_str = ident_str + ident + ","
		ident_str = ident_str[0:len(ident_str)-1]
		inter.append(ident_str)

	
	df.drop(["Target molecules in dataset cuff","Target molecules in dataset deseq"],axis=1,inplace=True)
	df["Overlapping DS genes"] = inter
	df["Gene_Count"] = gene_count
	df["Mean Z-score"] = df[["Activation Z-score deseq","Activation Z-score cuff"]].apply(np.mean,axis=1)
	df.drop(["Activation Z-score cuff","Activation Z-score deseq"],axis=1,inplace=True)
	return df


directory = "data/"
times = ["3h","6h","12h","24h"]


for t in times:
	cuff_down = pd.read_csv(directory+t+"_CuffDiff_down_UpstreamRegs.txt",sep="\t")
	cuff_up = pd.read_csv(directory+t+"_CuffDiff_up_UpstreamRegs.txt",sep="\t")
	deseq_down = pd.read_csv(directory+t+"_DEseq_down_UpstreamRegs.txt",sep="\t")
	deseq_up = pd.read_csv(directory+t+"_DEseq_up_UpstreamRegs.txt",sep="\t")

	down = merge(cuff_down,deseq_down)
	down = overlap_merge(down)
	down.to_csv("merged/"+t+"_down.csv",sep="\t",index=False)
	

	up = merge(cuff_up,deseq_up)
	up = overlap_merge(up)
	up.to_csv("merged/"+t+"_up.csv",sep="\t",index=False)
	
	

