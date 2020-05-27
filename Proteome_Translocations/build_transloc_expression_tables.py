import pandas as pd
import numpy as np

df_trans = pd.read_csv("output/transloc_table.csv")
expressions = pd.read_csv("output/reduced_expressions.csv")
cut = -0.8
exp_cols = [c for c in expressions.columns if(("Corr" in c)&("Err" not in c))]
times=["0","10","20","40","80","120","180","360","540","720","900","1080","1440","2880","5760"]

#pm_cyt
pm_cyt = df_trans[["gene","pm_cyt"]]
pm_cyt = pm_cyt[pm_cyt["pm_cyt"]<cut]
exp_series = []
for gene in pm_cyt["gene"].values:
   pm_exp = list(expressions[(expressions["fraction"]=="pm")&(expressions["gene"]==gene)][exp_cols].values[0])
   cyt_exp = list(expressions[(expressions["fraction"]=="cyt")&(expressions["gene"]==gene)][exp_cols].values[0])
   exp_series.append(pm_exp+cyt_exp)

dfout = pd.DataFrame(np.array(exp_series),columns = [t+"_pm" for t in times]+[t+"_cyt" for t in times])
dfout.insert(0,"gene",pm_cyt["gene"].values)
dfout.to_csv("output/pm_cyt.csv",sep="\t",index=False)


#pm_nuc
pm_nuc = df_trans[["gene","pm_nuc"]]
pm_nuc = pm_nuc[pm_nuc["pm_nuc"]<cut]
exp_series = []
for gene in pm_nuc["gene"].values:
   pm_exp = list(expressions[(expressions["fraction"]=="pm")&(expressions["gene"]==gene)][exp_cols].values[0])
   nuc_exp = list(expressions[(expressions["fraction"]=="nuc")&(expressions["gene"]==gene)][exp_cols].values[0])
   exp_series.append(pm_exp+nuc_exp)

dfout = pd.DataFrame(np.array(exp_series),columns = [t+"_pm" for t in times]+[t+"_nuc" for t in times])
dfout.insert(0,"gene",pm_nuc["gene"].values)
dfout.to_csv("output/pm_nuc.csv",sep="\t",index=False)


#cyt_nuc
cyt_nuc = df_trans[["gene","cyt_nuc"]]
cyt_nuc = cyt_nuc[cyt_nuc["cyt_nuc"]<cut]
exp_series = []
for gene in cyt_nuc["gene"].values:
   cyt_exp = list(expressions[(expressions["fraction"]=="cyt")&(expressions["gene"]==gene)][exp_cols].values[0])
   nuc_exp = list(expressions[(expressions["fraction"]=="nuc")&(expressions["gene"]==gene)][exp_cols].values[0])
   exp_series.append(cyt_exp+nuc_exp)

dfout = pd.DataFrame(np.array(exp_series),columns = [t+"_cyt" for t in times]+[t+"_nuc" for t in times])
dfout.insert(0,"gene",cyt_nuc["gene"].values)
dfout.to_csv("output/cyt_nuc.csv",sep="\t",index=False)
