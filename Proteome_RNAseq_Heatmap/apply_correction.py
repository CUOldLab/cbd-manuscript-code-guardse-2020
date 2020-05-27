import pandas as pd
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests

df = pd.read_csv('timeseries_significance_qvalues.csv',sep=",")
corrected = multipletests(df['p_values'].values,alpha=0.05,method='bonferroni')[1]
df['bh_corrected'] = corrected
df.to_csv('corrected.txt',index=False)
