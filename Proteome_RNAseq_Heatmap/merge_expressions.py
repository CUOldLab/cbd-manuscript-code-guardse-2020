import pandas as pd
import numpy as np

exp_df = pd.read_csv('data/cufflinks_diffexp.txt',sep="\t")
exp_df.index = exp_df['gene_id']+"&&"+exp_df['chromosome']



edge_df = pd.read_csv('timeseries_significance_qvalues.csv')
edge_df.index = edge_df['proteins.gene_id']+"&&"+edge_df['proteins.chromosome']
edge_df = edge_df[edge_df.columns[2:]]


merged = pd.concat([exp_df,edge_df],axis=1)
merged = merged[merged.isnull().sum(axis=1)==0]

expcols = ['3H_LFC','6H_LFC','12H_LFC','24H_LFC']
merged['max_expression'] = merged[expcols].apply(np.abs).apply(np.max,axis=1).values
merged = merged[merged['max_expression']>1.0]
merged = merged[merged['q_vals']<0.05]
merged.to_csv('final_results.txt',sep="\t",index=False)
