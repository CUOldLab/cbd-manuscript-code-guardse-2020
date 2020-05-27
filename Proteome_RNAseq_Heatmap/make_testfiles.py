import pandas as pd
import numpy as np
import re
import sys



def run():
	root_out = 'edge_inp/'
	df = pd.read_csv('data/cufflinks_merged.txt',sep="\t")
	df.replace(0,np.nan,inplace=True)
	df = df[df.isnull().sum(axis=1)==0]
	names = df[['gene_id','chromosome']]
	names.to_csv('data/names.txt',sep="\t")
	df = df[df.columns[2:]].apply(np.log2)
	
	df.to_csv(root_out+'Expr.csv',index=False,header=False)


	tgroup=[]
	for c in df.columns:
		if "VEH" in c:
			tgroup.append("VEH")
		else:
			tgroup.append("CBD")

	pd.DataFrame(tgroup,columns=["class"]).to_csv(root_out+'class.csv',index=False)



	times=[]
	for c in df.columns:
		time = c.split("_")[1]
		r = re.compile("([0-9]+)([a-zA-Z]+)")
		m = r.match(time)
		times.append(float(m.group(1)))
	pd.DataFrame(times,columns=["time"]).to_csv(root_out+'time.csv',index=False)


	indices=[]
	cnt=0;cindex=1
	for c in df.columns:
		indices.append(c.split("_")[2])

	pd.DataFrame(indices,columns=["individual"]).to_csv(root_out+'indices.csv',index=False)


run()


	



