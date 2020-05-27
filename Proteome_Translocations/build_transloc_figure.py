import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import preprocessing
from sklearn.decomposition import PCA



def append_annotations(annot,all_expressions):
	colors=[]
	for idx,row in all_expressions.iterrows():
		_human = row["Protein"]
		annotation = annot[annot["Entry name"]==_human]
		comp_desc = str(annotation["Subcellular location [CC]"])	
		if("Cell membrane" in comp_desc):
			colors.append("green")
		#elif("Endoplasmic reticulum" in comp_desc):
		#	colors.append("red")
		else:
			colors.append("black")
	return(colors)



short_list = pd.read_csv("output/translocators.csv")
all_expressions = pd.read_csv("output/all_expressions.csv")
annot = pd.read_csv("annotations/human.tab",sep="\t")
colors = append_annotations(annot,all_expressions)
all_expressions["colors"] = colors
all_expressions = all_expressions[all_expressions['colors']!='black']


exp_cols = [c for c in all_expressions.columns if(("Corr" in c)&("Err" not in c))]
data = all_expressions[exp_cols]
data = preprocessing.scale(data)

model = PCA(n_components=2)
transformed = model.fit_transform(data)
plt.scatter(transformed[:,0],transformed[:,1],alpha=1,s=10,c=all_expressions["colors"])
plt.xlim([-10,10])
plt.ylim([-10,10])
plt.show()
