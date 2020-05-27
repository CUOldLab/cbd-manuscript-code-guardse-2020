import pandas as pd
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from sklearn import preprocessing
from sklearn.decomposition import PCA


def return_exp_var(values,thresh):
   #reduce to N PCA components until 95% explained variance
   for i in range(2,15):
      model = PCA(n_components=i)
      transformed = model.fit_transform(values)
      exp_var = np.sum(model.explained_variance_ratio_)
      if exp_var>0.9:
         return(i)
   return(values.shape[1])


def transform(values_scaled,n_comp):
   model = PCA(n_components=n_comp)
   transformed = model.fit_transform(values)
   return(transformed)

def insert_fraction(df,name):
   df.insert(0,"fraction",[name for i in range(df.shape[0])])
   return(df)

def drop_missing(df):
   df.replace(-0.0001,np.nan,inplace=True)
   df.dropna(axis=0,inplace=True)
   return(df)


def find_translocs(df,pc_cols):
   pm_sub = df[df["fraction"]=="pm"]
   cyt_sub = df[df["fraction"]=="cyt"]
   nuc_sub = df[df["fraction"]=="nuc"]
   correls=[]
   genes=[]

   for idx,row in df.iterrows():
      gene = row["gene"]
      if gene in genes:
         continue
      genes.append(gene)
      pm_short = pm_sub[pm_sub["gene"]==gene][pc_cols].values
      cyt_short = cyt_sub[cyt_sub["gene"]==gene][pc_cols].values
      nuc_short = nuc_sub[nuc_sub["gene"]==gene][pc_cols].values
      pm_cyt = 0;pm_nuc = 0;cyt_nuc = 0;

      if ((pm_short.shape[0]==1)&(cyt_short.shape[0]==1)):
         pm_cyt = pearsonr(pm_short[0,:],cyt_short[0,:])[0]
      if ((pm_short.shape[0]==1)&(nuc_short.shape[0]==1)):
         pm_nuc = pearsonr(pm_short[0,:],nuc_short[0,:])[0]
      if ((cyt_short.shape[0]==1)&(nuc_short.shape[0]==1)):
         cyt_nuc = pearsonr(cyt_short[0,:],nuc_short[0,:])[0]
      correls.append([pm_cyt,pm_nuc,cyt_nuc])

   dfout = pd.DataFrame(np.array(correls),columns=["pm_cyt","pm_nuc","cyt_nuc"])
   dfout.insert(0,"gene",genes)
   return(dfout)
      



#try it first on the CYT fraction
cyt = pd.read_csv('/home/jjacob/R/edge/proteome_diff_analysis/20180223/results/CYT/raw/merged.txt')
pm = pd.read_csv('/home/jjacob/R/edge/proteome_diff_analysis/20180223/results/PM/raw/merged.txt')
nuc = pd.read_csv('/home/jjacob/R/edge/proteome_diff_analysis/20180223/results/NUC/raw/merged.txt')

cyt = insert_fraction(cyt,"cyt")
nuc = insert_fraction(nuc,"nuc")
pm = insert_fraction(pm,"pm")

cyt = drop_missing(cyt)
nuc = drop_missing(nuc)
pm = drop_missing(pm)


sigprots_cyt = list(cyt[cyt["q_values"]<0.01]["gene"].values)
sigprots_pm = list(pm[pm["q_values"]<0.01]["gene"].values)
sigprots_nuc = list(nuc[nuc["q_values"]<0.01]["gene"].values)
siggenes = list(set(sigprots_cyt+sigprots_pm+sigprots_nuc))

merged = pd.concat([cyt,nuc,pm],axis=0)


##reduce noise contribution using PCA.  then look for correl/anti-correl
exp_cols = [c for c in merged.columns if(("Corr" in c)&("Err" not in c))]
values = merged[exp_cols].values
values_scaled = preprocessing.scale(values)
n_comp = return_exp_var(values,0.95)
transformed = transform(values_scaled,n_comp)
pc_cols = []
for i in range(transformed.shape[1]):
   merged["pc"+str(i)] = transformed[:,i]
   pc_cols.append("pc"+str(i))

merged = merged[merged["gene"].isin(siggenes)]
merged.to_csv("output/reduced_expressions.csv",index=False)
transloc_table = find_translocs(merged,pc_cols)
transloc_table.to_csv("output/transloc_table.csv",index=False)










