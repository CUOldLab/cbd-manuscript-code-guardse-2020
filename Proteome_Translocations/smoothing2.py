import pandas as pd
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr




def scale(nuc,pm):
	scaled_nuc = (nuc-np.min(nuc))/(np.max(nuc)-np.min(nuc))
	scaled_pm = (pm-np.min(pm))/(np.max(pm)-np.min(pm))
	return(scaled_nuc,scaled_pm)

def funfit(nuc):
	for i in range(len(nuc)-1): #startpoint
		for j in range(1,len(nuc)): #endpoint
			plt.plot([i,j],[nuc[i],nuc[j]],color='gray',alpha=0.1)
	


df = pd.read_csv("hk1_exp.csv")
nuc = df[[c for c in df.columns if ("nuc" in c)]].values[0]
pm = df[[c for c in df.columns if ("pm" in c)]].values[0]

sample_data_x = np.linspace(0,2*3.14,50)
vals = np.sin(sample_data_x)
print(vals)

funfit(vals)
plt.show()
nuc,pm = scale(nuc,pm)

'''
smoothed_nuc = gaussian_filter1d(nuc,2)
smoothed_pm = gaussian_filter1d(pm,2)


basal_correlation = spearmanr(nuc,pm)[0]
smoothed_correlation = spearmanr(smoothed_nuc,smoothed_pm)[0]
print(basal_correlation,smoothed_correlation)

plt.plot(nuc)
plt.plot(smoothed_nuc)
plt.plot(smoothed_pm)
plt.plot(pm)
plt.show()
'''
