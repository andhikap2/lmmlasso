#Scipy Kendall rank calculation
#daily minmax model only (No par no ptt no tt no daylength
import sys
import limix.vardec.vardec as va 
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
import limix.vardec.vardec as va 

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
plantings=np.load("Plantings_NEW.npy")

import matplotlib
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
K=np.load('K_NEWER.npy') #Load the preferred K-matrix
#K=np.load("K_MATRIX_PLINKIBS.npy")
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")
snps=np.load('SNPs_0.1.npy') #If x were given as just snps then they would be fixed effects right
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/New_Microclim")
env=pd.read_csv('Microclimate_minmaxDaily_threshold_0_0.csv',sep=',')
logical=env.columns.str.startswith(('PAR','TT','PTT','daylength') )
env=env.iloc[:,~logical] #get columns that don't start with PAR
environment=np.array(env)
X=np.concatenate((snps,environment),axis=1)
assert K.shape[0] == K.shape[1], 'K MATRIX IS NOT SYMMETRICAL'
assert K.shape[0] == env.shape[0], 'NO. OF INDIVIDUALS DOES NOT MATCH K MATRIX DIMENSIONS'
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
y=np.load("DTB_NEW.npy")
assert y.shape[0] == X.shape[0], 'NO. OF OBSERVATIONS DOES NOT MATCH NO. OF INDIVIDUALS'
#Running LMMLASSO
alphas = 2.**(sp.linspace(-2,10,10)) #list of alphas to test
n_splits=10
N = X.shape[0]
kf = KFold(n_splits,shuffle=True,random_state=None)
n_alphas = len(alphas)
MSE_train = sp.zeros((n_splits,n_alphas))
MSE_test  = sp.zeros((n_splits,n_alphas))
W_nonzero = sp.zeros((n_splits,n_alphas))
rsquared  = sp.zeros((n_splits,n_alphas))
kf.get_n_splits(X)
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
import lmmlasso
lasso = lmmlasso.LmmLasso(warm_start=True,fit_intercept=False,tol=0.5) #note the tolerance value
MSE_train,MSE_test,W_nonzero, rsquared = lmmlasso.runCrossValidation(lasso,X,y,alphas,n_splits=10,K=K,verbose=True)
MSE_train_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_train.mean(axis=0))).derivative(n=2)#Interpolating the values for alphas within -2 to 10...
MSE_test_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_test.mean(axis=0))).derivative(n=2)
alphas_inter = 2.**(sp.linspace(2,10,100))
idx_train = sp.argmin(MSE_train_inter(alphas_inter)) 
idx_test = sp.argmin(MSE_test_inter(alphas_inter))
alpha_cv = (float(alphas_inter[idx_train])+float(alphas_inter[idx_test]))/2
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
kf12 = KFold(n_splits,shuffle=True,random_state=12) #12 because the random state is 12
kf12.get_n_splits(X)
#Model fitting step
N = X.shape[0]
kf = KFold(n_splits,shuffle=True,random_state=12)
MSE_train_final = sp.zeros((n_splits,))
MSE_test_final  = sp.zeros((n_splits,))
W_nonzero_final = sp.zeros((n_splits,))
rsquared_final  = sp.zeros((n_splits,)) 
kendall_final = sp.zeros((n_splits,))
kf.get_n_splits(X)
lasso.set_params(alpha=alpha_cv)

ifold = 0 #what is ifold
for train_index,test_index in kf.split(X):
    print(('running fold %d'%ifold)) 
    X_train, X_test, y_train, y_test= X[train_index], X[test_index], y[train_index], y[test_index]
    K_train = K[train_index][:,train_index]
    K_test  = K[test_index][:,train_index]
    model_fit=lasso.fit(X_train,y_train,K=K_train) #Fit the model to the training set
    ytrain_star=model_fit.predict(X_train,K_train)
    ytest_star=model_fit.predict(X_test,K_test)
    MSE_train_final[ifold]=mean_squared_error(ytrain_star,y_train)
    MSE_test_final[ifold]=mean_squared_error(ytest_star,y_test)
    W_nonzero_final[ifold]=sp.sum(model_fit.coef_!=0)
    rsquared_final[ifold]=(sp.stats.pearsonr(ytest_star,y_test)[0])**2 
    kendall[ifold]=sp.stats.kendalltau(ytest_star,y_test)[0]       
    ifold +=1

for train_index,test_index in kf12.split(X):
	X_train, X_test, y_train, y_test, plantings_train, plantings_test= X[train_index], X[test_index], y[train_index], y[test_index], plantings[train_index], plantings[test_index] #Split into training and testing data sets
	K_train = K[train_index][:,train_index]
	K_test  = K[test_index][:,train_index]

lasso.set_params(alpha=alpha_cv)
lasso = lasso.fit(X_train,y_train,K=K_train)
weights = lasso.coef_
Y_hat = lasso.predict(X_test, K_test)
pd_Yhat=pd.DataFrame(Y_hat)
pd_Yhat.reset_index(drop=True,inplace=True)
pd_ytest=pd.DataFrame(y_test)
pd_plantings=pd.DataFrame(plantings_test)
for_plotting=pd.concat([pd_plantings,pd_ytest,pd_Yhat],axis=1)
for_plotting.columns=("Planting","Y_test","Y_hat")

HalleFall2006=for_plotting.loc[for_plotting['Planting']=="HalleFall2006"]
NorwichSummer2006=for_plotting.loc[for_plotting['Planting']=="NorwichSummer2006"]
NorwichSummer2007=for_plotting.loc[for_plotting['Planting']=="NorwichSummer2007"]
NorwichSpring2007=for_plotting.loc[for_plotting['Planting']=="NorwichSpring2007"]
NorwichFall2006=for_plotting.loc[for_plotting['Planting']=="NorwichFall2006"]
OuluFall2007=for_plotting.loc[for_plotting['Planting']=="OuluFall2007"]
ValenciaFall2006=for_plotting.loc[for_plotting['Planting']=="ValenciaFall2006"]

import matplotlib.pyplot as plt
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/LMMLASSO/Figures")

plt.figure()


ranges=range(int(np.amax(y_test)))
plt.plot(ranges,ranges)
v=str(alpha_cv)

plt.xlabel('Observed DTB')
plt.ylabel('Predicted DTB')
plt.plot(HalleFall2006.iloc[:,1],HalleFall2006.iloc[:,2],color='#8f0000',marker='o',linestyle='None',label='HalleFall2006')
plt.plot(NorwichSummer2006.iloc[:,1],NorwichSummer2006.iloc[:,2],color='#c54259',marker='X',linestyle='None',label='NorwichSummer2006')
plt.plot(NorwichSummer2007.iloc[:,1],NorwichSummer2007.iloc[:,2],color='#e881ab',marker='p',linestyle='None',label='NorwichSummer2007')
plt.plot(NorwichSpring2007.iloc[:,1],NorwichSpring2007.iloc[:,2],color='#ffc2f3',marker='v',linestyle='None',label='NorwichSpring2007')
plt.plot(NorwichFall2006.iloc[:,1],NorwichFall2006.iloc[:,2],color='#d890e7',marker='d',linestyle='None',label='NorwichFall2006')
plt.plot(OuluFall2007.iloc[:,1],OuluFall2007.iloc[:,2],color='#9f66e1',marker='P',linestyle='None',label='OuluFall2007')
plt.plot(ValenciaFall2006.iloc[:,1],ValenciaFall2006.iloc[:,2],color='#3a47de',marker='*',linestyle='None',label='ValenciaFall2006')


plt.legend(loc='upper left',frameon=True,fontsize='small')
plt.savefig('PredvsObs_SNPs_0.1_minmaxDaily.png')


alpha_cv_np=np.array(alpha_cv)
c=np.array(sum(weights!=0))
W_nonzero=np.array(c)
rmse_test=math.sqrt(np.mean(MSE_test_final))
r2_model=np.mean(rsquared_final)
results=np.array((alpha_cv_np,W_nonzero,rmse_test,r2_model))