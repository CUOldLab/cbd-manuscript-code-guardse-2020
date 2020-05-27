import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pylab
from scipy.optimize import curve_fit
import seaborn as sns
from pylab import rcParams
rcParams['figure.figsize'] = 8, 7


df = pd.read_csv('SKNBE2_sensors_20171002_corrected_timepoint_study.csv')
full_df = pd.read_csv('SKNBE2_sensors_20171002_corrected.csv')

dose = 13
sensors = ['glucose','Cyto_Ca2','pm_charge','AMP_kinase','Lactate','Glutamine','Adam17_Protease','ER_Ca2','PKD_Kinase','ATP','Membrane_Potential','Erk_Kinase','Pyruvate','mTOR']



for s in sensors:
	col = s+"_"+str(dose)
	plt.plot(full_df['Time'].values,full_df[col].values,lw=2,label="Normalized FRET ratio")
	plt.scatter(df['Time'].values,df[col].values,s=60,color='red',alpha=0.5,label="Selected Time Points for Sampling")
	plt.legend(loc=2,prop={'size':16})
	plt.title(s+" ; Dose = 13uM",fontsize=18)
	plt.xlabel("Time (hours)",fontsize=18)
	plt.ylabel("Normalized FRET Ratio",fontsize=18)
	plt.savefig('./selected_timepoint_figures_20171002/'+s+".png")
	plt.clf()
		
