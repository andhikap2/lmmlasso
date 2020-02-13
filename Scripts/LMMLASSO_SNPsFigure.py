'''
Making figures to show the effects of including SNPs
	Keeping the number of alphas tested low to keep results consistent with what you have in the tables

'''

import sys
#import limix.vardec.vardec as va 
import scipy.linalg as LA
from sklearn.linear_model import Lasso
import numpy as np 
from sklearn.metrics import mean_squared_error
import math
from sklearn.model_selection import KFold 
import gc
import scipy as sp
import os
import csv
import time
import pandas as pd

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
plantings=np.load("Plantings_NEW.npy")

import matplotlib
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
K=np.load('K_NEWER.npy')


os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")
snps=np.load('SNPs_0.1.npy') #If x were given as just snps then they would be fixed effects right
snp05=np.load('SNPs_0.5.npy')
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/New_Microclim")

env0=pd.read_csv('Microclimate_Daily_threshold_0_0.csv',sep=',')





environment0=np.array(env0)




X0=np.concatenate((snps,environment0),axis=1)
X250=np.concatenate((snp05,environment0),axis=1)


assert K.shape[0] == K.shape[1], 'K MATRIX IS NOT SYMMETRICAL'
assert K.shape[0] == env0.shape[0], 'NO. OF INDIVIDUALS DOES NOT MATCH K MATRIX DIMENSIONS'
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
y=np.load("DTB_NEW.npy")
assert y.shape[0] == X.shape[0], 'NO. OF OBSERVATIONS DOES NOT MATCH NO. OF INDIVIDUALS'

'''Setting LASSO parameters'''

alphas = 2.**(sp.linspace(-10,10,10)) #list of alphas to test
n_splits=10
N = X.shape[0]
kf = KFold(n_splits,shuffle=True,random_state=None)
n_alphas = len(alphas)
MSE_train = sp.zeros((n_splits,n_alphas))
MSE_test  = sp.zeros((n_splits,n_alphas))
W_nonzero = sp.zeros((n_splits,n_alphas))
rsquared  = sp.zeros((n_splits,n_alphas))
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")
import lmmlasso


'''Using SNPs'''

lasso = lmmlasso.LmmLasso(warm_start=True,fit_intercept=False,tol=0.5) #note the tolerance value
MSE_train,MSE_test,W_nonzero, rsquared = lmmlasso.runCrossValidation(lasso,X0,y,alphas,n_splits=10,K=K,verbose=True)
MSE_train_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_train.mean(axis=0))).derivative(n=2)#Interpolating the values for alphas within -2 to 10...
MSE_test_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_test.mean(axis=0))).derivative(n=2)
alphas_inter = 2.**(sp.linspace(-10,10,100))
idx_train = sp.argmin(MSE_train_inter(alphas_inter)) 
idx_test = sp.argmin(MSE_test_inter(alphas_inter))
alpha_cv = (float(alphas_inter[idx_train])+float(alphas_inter[idx_test]))/2
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
kf12 = KFold(n_splits,shuffle=True,random_state=12) #12 because the random state is 12
#Model fitting step
N = X0.shape[0]
#kf = KFold(n_splits,shuffle=True,random_state=12)
MSE_train_final = sp.zeros((n_splits,))
MSE_test_final  = sp.zeros((n_splits,))
W_nonzero_final = sp.zeros((n_splits,))
rsquared_final  = sp.zeros((n_splits,)) 
lasso.set_params(alpha=alpha_cv)


for train_index,test_index in kf12.split(X0):
	X_train, X_test, y_train, y_test, plantings_train, plantings_test= X0[train_index], X0[test_index], y[train_index], y[test_index], plantings[train_index], plantings[test_index] #Split into training and testing data sets
	K_train = K[train_index][:,train_index]
	K_test  = K[test_index][:,train_index]

lasso = lasso.fit(X_train,y_train,K=K_train)
weights = lasso.coef_
Y_hat = lasso.predict(X_test, K_test)
pd_Yhat=pd.DataFrame(Y_hat)
pd_Yhat.reset_index(drop=True,inplace=True)
pd_ytest=pd.DataFrame(y_test)
pd_plantings=pd.DataFrame(plantings_test)
for_plotting0=pd.concat([pd_plantings,pd_ytest,pd_Yhat],axis=1)
for_plotting0.columns=("Planting","Y_test","Y_hat")


'''Using SNPs at threshold of 0.5'''


N = X250.shape[0]
#kf = KFold(n_splits,shuffle=True,random_state=12)
MSE_train_final = sp.zeros((n_splits,))
MSE_test_final  = sp.zeros((n_splits,))
W_nonzero_final = sp.zeros((n_splits,))
rsquared_final  = sp.zeros((n_splits,)) 
lasso.set_params(alpha=alpha_cv)


for train_index,test_index in kf12.split(X250):
	X_train, X_test, y_train, y_test, plantings_train, plantings_test= X250[train_index], X250[test_index], y[train_index], y[test_index], plantings[train_index], plantings[test_index] #Split into training and testing data sets
	K_train = K[train_index][:,train_index]
	K_test  = K[test_index][:,train_index]

lasso = lasso.fit(X_train,y_train,K=K_train)
weights = lasso.coef_
Y_hat = lasso.predict(X_test, K_test)
pd_Yhat=pd.DataFrame(Y_hat)
pd_Yhat.reset_index(drop=True,inplace=True)
pd_ytest=pd.DataFrame(y_test)
for_plotting250=pd.concat([pd_ytest,pd_Yhat],axis=1)
for_plotting250.columns=("Y_test250","Y_hat250")





for_plottingfinal=pd.concat([for_plotting0,for_plotting250],axis=1)



HalleFall2006=for_plottingfinal.loc[for_plottingfinal['Planting']=="HalleFall2006"]
NorwichSummer2006=for_plottingfinal.loc[for_plottingfinal['Planting']=="NorwichSummer2006"]
NorwichSummer2007=for_plottingfinal.loc[for_plottingfinal['Planting']=="NorwichSummer2007"]
NorwichSpring2007=for_plottingfinal.loc[for_plottingfinal['Planting']=="NorwichSpring2007"]
NorwichFall2006=for_plottingfinal.loc[for_plottingfinal['Planting']=="NorwichFall2006"]
OuluFall2007=for_plottingfinal.loc[for_plottingfinal['Planting']=="OuluFall2007"]
ValenciaFall2006=for_plottingfinal.loc[for_plottingfinal['Planting']=="ValenciaFall2006"]




import matplotlib.pyplot as plt
matplotlib.use('agg')
plt.figure()


ranges=range(int(np.amax(y_test)))
#ranges=range(250)
plt.plot(ranges,ranges)
v=str(alpha_cv)


plt.xlabel('Observed DTB')
plt.ylabel('Predicted DTB')

#plotting with SNPs

plt.plot(HalleFall2006.iloc[:,1],HalleFall2006.iloc[:,2],color='#e881ab',marker='o',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichSummer2006.iloc[:,1],NorwichSummer2006.iloc[:,2],color='#e881ab',marker='X',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichSummer2007.iloc[:,1],NorwichSummer2007.iloc[:,2],color='#e881ab',marker='p',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichSpring2007.iloc[:,1],NorwichSpring2007.iloc[:,2],color='#e881ab',marker='v',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichFall2006.iloc[:,1],NorwichFall2006.iloc[:,2],color='#e881ab',marker='d',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(OuluFall2007.iloc[:,1],OuluFall2007.iloc[:,2],color='#e881ab',marker='P',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(ValenciaFall2006.iloc[:,1],ValenciaFall2006.iloc[:,2],color='#e881ab',marker='*',linestyle='None',label='',mew=1,mec='black',ms=7.0)



#plotting at SNP = 0.5



plt.plot(HalleFall2006.iloc[:,1],HalleFall2006.iloc[:,4],color='#aae881',marker='o',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichSummer2006.iloc[:,1],NorwichSummer2006.iloc[:,4],color='#aae881',marker='X',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichSummer2007.iloc[:,1],NorwichSummer2007.iloc[:,4],color='#aae881',marker='p',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichSpring2007.iloc[:,1],NorwichSpring2007.iloc[:,4],color='#aae881',marker='v',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichFall2006.iloc[:,1],NorwichFall2006.iloc[:,4],color='#aae881',marker='d',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(OuluFall2007.iloc[:,1],OuluFall2007.iloc[:,4],color='#aae881',marker='P',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(ValenciaFall2006.iloc[:,1],ValenciaFall2006.iloc[:,4],color='#aae881',marker='*',linestyle='None',label='',mew=1,mec='black',ms=7.0)


from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='#e881ab', lw=2),
                Line2D([0], [0], color='#aae881', lw=2),
              ]

plt.legend((custom_lines),('r2 = 0.1','r2 = 0.5'),loc='upper left',frameon=True,fontsize='small')

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/LMMLASSO/Figures")


plt.savefig('PredvsObs_DifferentSNPThresholds.png')



alpha_cv_np=np.array(alpha_cv)
c=np.array(sum(weights!=0))
W_nonzero=np.array(c)
rmse_test=math.sqrt(np.mean(MSE_test_final))
r2_model=np.mean(rsquared_final)
results=np.array((alpha_cv_np,W_nonzero,rmse_test,r2_model))














