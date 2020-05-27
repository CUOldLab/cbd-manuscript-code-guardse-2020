import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pylab
from scipy.optimize import curve_fit
import seaborn as sns


df = pd.read_csv('SKNBE2_sensors_20171002_corrected.csv')
doses = np.array([100,66.7,44.4,29.6,19.8,13,8.8,5.9,3.9,2.6,1.73,1.16,0.771,0.514,0.343,0])
doses = list(doses[range(len(doses)-1,-1,-1)])
sensors = ['glucose','Cyto_Ca2','pm_charge','AMP_kinase','Lactate','Glutamine','Adam17_Protease','ER_Ca2','PKD_Kinase','ATP','Membrane_Potential','Erk_Kinase','Pyruvate','mTOR']



def plotTimeByDose():
	for s in sensors:
		for d in doses:
			plt.plot(df['Time'].values,df[s+"_"+str(d)].values,label=str(d)+"uM",lw=3,alpha=0.5)
		plt.title("Sensor = "+str(s))
		plt.legend()
		plt.show()


def plotTimeBySensor(sensors,doses,df):
	for d in doses:
		for s in sensors:
			plt.plot(df['Time'].values,df[s+"_"+str(d)].values,label=s,lw=3,alpha=0.5)
		plt.title("Dose = "+str(d)+"uM")
		plt.legend()
		plt.show()


def reject_outliers(data, m=1):
	data = np.array(data)
   	return data[abs(data - np.mean(data)) < m * np.std(data)]


def plotPercResp(sensor,doses,df):
	def sigmoid(x, x0, k):
		y = 1 / (1 + np.exp(-k*(x-x0)))
		return y

	times = []
	EC50s=[]
	for time in df["Time"].values:
		activity_perc = []
		short = df[df['Time']==time]
		dat = np.flip(short[[c for c in short.columns if(sensor in c)]].values[0], axis=0)
		#plt.plot(dat,lw=3,alpha=0.5,label="raw ratio",color='blue')
	
		for val in dat:
			val = val-np.min(dat)
			activity_perc.append(val/(np.max(dat)-np.min(dat)))
		
		try:
			popt, pcov = curve_fit(sigmoid, doses, activity_perc)
		except:
			continue

		if(popt[0]<0) | (popt[0]>np.max(doses)):
			continue
		
		times.append(time)
		EC50s.append(popt[0])

		x = np.linspace(np.min(doses), np.max(doses), 50)
		y = sigmoid(x, *popt)
		pylab.axvline(popt[0],color='red')
		pylab.plot(doses, activity_perc, 'o', label='data')
		pylab.plot(x,y, label='fit')
		pylab.ylim(0, 1.05)
		pylab.legend(loc='best')
		pylab.title(sensor+", EC50 = "+str(popt[0])+", Time = "+str(time), fontsize=18)
		pylab.xlabel("Dose, uM", fontsize=18)
		pylab.ylabel("Percentage Response", fontsize=18)
		#pylab.savefig('')
		pylab.savefig("sigmoid_fits/"+sensor+"_time_"+str(round(time,2))+".png")	
		pylab.clf()

	plt.plot(times,EC50s,lw=3,alpha=0.6)
	plt.xlabel("Time",fontsize=18)
	plt.ylabel("EC50",fontsize=18)
	plt.title("EC50 @ Time for: "+sensor,fontsize=18)
	plt.savefig("sigmoid_summary_figures/"+sensor+".png")
	plt.clf()
	return EC50s


all_EC50s = []
for s in sensors:
	EC50s = plotPercResp(s,doses,df)
	all_EC50s = all_EC50s+EC50s

#all_EC50s = reject_outliers(all_EC50s)
plt.hist(all_EC50s,bins=20,histtype='step',color='blue',lw=2,label="Distribution of all EC50s")
#plt.hist(reject_outliers(all_EC50s),bins=20,histtype='step',color='red',lw=2,label="Distribution of all EC50s, outliers removed, cutoff= 1 std")
plt.legend(prop={'size':16})
plt.title("Median EC50: "+str(np.median(all_EC50s)),fontsize=18)
plt.show()







