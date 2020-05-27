import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

edge_res = pd.read_csv('corrected.txt')
edge_res = edge_res[edge_res['bh_corrected']<0.05]
cuff_exp = pd.read_csv('data/diffexp.txt',sep="\t")
expcols = [c for c in cuff_exp.columns if("_diffexp" in c)]
stdcols = [c for c in cuff_exp.columns if("_stdev" in c)]
max_expressions = []
expression_vals = []
stdev_vals = []

for idx,row in edge_res.iterrows():
	short = cuff_exp[(cuff_exp['tracking_id'] == row['proteins.tracking_id']) & (cuff_exp['locus'] == row['proteins.locus'])]
	expression_max = short[expcols].apply(np.abs).apply(np.max,axis=1).values[0]
	max_expressions.append(expression_max)
	expression_vals.append(short[expcols].values[0])
	stdev_vals.append(short[stdcols].values[0])

	'''
	if expression_max>0.5:
		plt.errorbar(range(4),short[expcols].values[0],yerr=short[stdcols].values[0])
		plt.title(row['proteins.tracking_id'],fontsize=16)
		plt.xlabel("TIME",fontsize=14)
		plt.ylabel("Log2(Ratio)",fontsize=14)
		plt.xticks(range(4), ['3H','6H','12H','24H'],fontsize=14)
		plt.xticks(fontsize=14)
		plt.yticks(fontsize=14)
		plt.tight_layout()
		plt.savefig('figures/'+row['proteins.tracking_id']+".png")
		plt.clf()
	'''

expression_vals = np.array(expression_vals)
stdev_vals = np.array(stdev_vals)

for i in range(len(expcols)):
	edge_res[expcols[i]] = expression_vals[:,i]

for i in range(len(stdcols)):
	edge_res[stdcols[i]] = stdev_vals[:,i]

edge_res['diffexp_max'] = max_expressions
edge_res = edge_res[edge_res['diffexp_max']>0.5]
edge_res.to_csv('corrected_reduced.txt',sep="\t",index=False)
