import pandas as pd
import numpy as np


sig_rna = pd.read_csv("/home/jjacob/python/SPARTA/testperiod_2017/rna_heatmap/merged_significant_transcripts_1perc_fdr.csv")
sig_prot = pd.read_csv("/home/jjacob/python/SPARTA/testperiod_2017/proteome_heatmap/merged_significant_proteins_1perc_fdr.csv")

#cut gene ids on lfc
sig_rna.replace([np.inf,-np.inf], np.nan, inplace=True)

sig_rna["max_exp"] = np.nanmax(np.abs(sig_rna[[c for c in sig_rna.columns if ("log2(fold_change)_" in c)]].values), axis=1)
sig_rna = sig_rna[sig_rna["max_exp"] > 1.0]

#start by putting all ids into a set
all_ids = list(set(list(sig_prot["gene"].values) + [v.split("_")[0] for v in sig_rna["unique_id"].values]))




#merge all_data
all_rna = pd.read_csv("/home/jjacob/python/SPARTA/testperiod_2017/rna_heatmap/merged_cuffdiff.csv", sep="\t")
all_rna.insert(0,"gene",[v.split("_")[0] for v in all_rna["unique_id"].values])
all_rna["min_qval"] = np.nanmin(all_rna[[c for c in all_rna.columns if ("q_value" in c)]].values, axis=1)
all_rna.sort_values(by="min_qval")
all_rna.replace([-np.inf, np.inf], np.nan, inplace=True)
all_rna = all_rna[["gene"] + [c for c in all_rna.columns if ("log2(fold_change)" in c)]]
all_rna.drop_duplicates(subset=["gene"], inplace=True)
all_rna.index = all_rna["gene"].values
print(["gene"]+ [c+"_rna" for c in all_rna.columns[1:]])
all_rna.columns = ["gene"]+ [c+"_rna" for c in all_rna.columns[1:]]
all_rna.drop("gene", axis=1, inplace=True)

all_prot = pd.read_csv("/home/jjacob/python/SPARTA/testperiod_2017/proteome_heatmap/merged_proteins.csv")
all_prot.replace(-0.0001, np.nan, inplace=True)
all_prot["min_qval"] = np.nanmin(all_prot[[c for c in all_prot.columns if ("q_values" in c)]].values, axis=1)
all_prot.drop_duplicates(subset=["gene"], inplace=True)
all_prot.replace([-np.inf, np.inf], np.nan, inplace=True)
all_prot = all_prot[["gene"] + [c for c in all_prot.columns if ("TRT" in c)]]
all_prot.index = all_prot["gene"].values
all_prot.drop("gene", axis=1, inplace=True)

#mergify
merged = all_prot.merge(all_rna, how="outer", left_index=True, right_index=True)
merged.insert(0,"gene",merged.index.values)
merged=merged[merged["gene"].isin(all_ids)]


#and merge with go terms
go_terms = pd.read_csv("results/revigo_graph.csv")
go_terms.index = go_terms["go annotation"].values #transpose the df
go_terms = go_terms.transpose()
go_terms.insert(0,"gene",go_terms.index.values)
go_terms.drop("go annotation",axis=0,inplace=True)

merged = go_terms.merge(merged, how="outer", left_on="gene", right_on="gene")
merged.dropna(subset = [c for c in merged.columns if (("TRT" in c) | ("log2" in c))], how="all", inplace=True)


merged.to_csv("results/heatmap_input.csv", index=False)









