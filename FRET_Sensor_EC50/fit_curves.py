import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv('Fret_data_with_replicates.csv')
df = df[df['Time']>=2]
sensors = ['glucose','Cytosolic_Ca','PM_charge','AMP_kinase','Lactate','Glutamine','Adam17_protease','ER_Ca','PDK_Kinase','ATP','Membrane_potential','Erk_kinase','Pyruvate','mTOR']
times = df['Time'].values
dose = '19.8'



for sens in sensors:

	pref = sens+"_"
	plt.scatter(times,df[[pref+'19.8_R1']].values,color='blue',s=60)
	plt.scatter(times,df[[pref+'19.8_R2']].values,color='blue',s=60)
	plt.scatter(times,df[[pref+'0_R1']].values,color='orange',s=60)
	plt.scatter(times,df[[pref+'0_R2']].values,color='orange',s=60)


	x = times
	y1 = np.mean(df[[pref+'19.8_R1',pref+'19.8_R2']].values,axis=1)
	y2 = np.mean(df[[pref+'0_R1',pref+'0_R2']].values,axis=1)


	treatment_fit = np.polynomial.Polynomial.fit(x, y1, 3)
	control_fit = np.polynomial.Polynomial.fit(x, y2, 3)


	fx = np.linspace(0,np.max(times),100)
	fy_treat = treatment_fit(fx)
	fy_ctrl = control_fit(fx)

	roots = (treatment_fit - control_fit).roots()
	putative_root=0
	for r in roots:
		if (r>0)&(r<15):
			putative_root=r
			break


	plt.plot(fx,fy_treat,color='black')
	plt.plot(fx,fy_ctrl,color='black')

	plt.axvline(putative_root,label="Min_Effect_Time: "+str(round(putative_root,2))+" hours",color='red',lw=3)
	plt.title(sens+" ; Dose : 19.8uM",fontsize=16)
	plt.legend(prop={'size':16})
	plt.xlabel("Time (Hours)",fontsize=16)
	plt.ylabel("Response",fontsize=16)
	plt.savefig('curve_fits_intersection/'+sens+".png")
	plt.clf()
	#plt.show()

