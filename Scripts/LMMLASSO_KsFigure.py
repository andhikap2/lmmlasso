#Making figures to show the effects of using different K matrices

'''
test the following k matrices

emmaxibs
ldak+0.5
ldak-0.5
emmaxbn
limix

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
K=np.load('K_NEWER.npy') #limix
Kldakpos=np.load('K_MATRIX_LDAKPOS0.50.npy')
Kldakneg=np.load('K_MATRIX_LDAKNEG0.50.npy')
Kemmaxbn=np.load('K_MATRIX_EMMAXBN.npy')
Kemmaxibs=np.load('K_MATRIX_EMMAXIBS.npy')


os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")
snps=np.load('SNPs_0.1.npy') #If x were given as just snps then they would be fixed effects right
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/New_Microclim")

env0=pd.read_csv('Microclimate_Daily_threshold_0_0.csv',sep=',')
environment0=np.array(env0)




X0=np.concatenate((snps,environment0),axis=1)

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


#K limix

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
N = X.shape[0]
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

#Kldakpos


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
N = X.shape[0]
#kf = KFold(n_splits,shuffle=True,random_state=12)
MSE_train_final = sp.zeros((n_splits,))
MSE_test_final  = sp.zeros((n_splits,))
W_nonzero_final = sp.zeros((n_splits,))
rsquared_final  = sp.zeros((n_splits,)) 
lasso.set_params(alpha=alpha_cv)


for train_index,test_index in kf12.split(X0):
	X_train, X_test, y_train, y_test, plantings_train, plantings_test= X0[train_index], X0[test_index], y[train_index], y[test_index], plantings[train_index], plantings[test_index] #Split into training and testing data sets
	K_train = Kldakpos[train_index][:,train_index]
	K_test  = Kldakpos[test_index][:,train_index]

lasso = lasso.fit(X_train,y_train,K=K_train)
weights = lasso.coef_
Y_hat = lasso.predict(X_test, K_test)
pd_Yhat=pd.DataFrame(Y_hat)
pd_Yhat.reset_index(drop=True,inplace=True)
pd_ytest=pd.DataFrame(y_test)
pd_plantings=pd.DataFrame(plantings_test)
for_plottingldakpos=pd.concat([pd_plantings,pd_ytest,pd_Yhat],axis=1)
for_plottingldakpos.columns=("Planting","Y_test","Y_hat")


#Kldakneg


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
N = X.shape[0]
#kf = KFold(n_splits,shuffle=True,random_state=12)
MSE_train_final = sp.zeros((n_splits,))
MSE_test_final  = sp.zeros((n_splits,))
W_nonzero_final = sp.zeros((n_splits,))
rsquared_final  = sp.zeros((n_splits,)) 
lasso.set_params(alpha=alpha_cv)


for train_index,test_index in kf12.split(X0):
	X_train, X_test, y_train, y_test, plantings_train, plantings_test= X0[train_index], X0[test_index], y[train_index], y[test_index], plantings[train_index], plantings[test_index] #Split into training and testing data sets
	K_train = Kldakneg[train_index][:,train_index]
	K_test  = Kldakneg[test_index][:,train_index]

lasso = lasso.fit(X_train,y_train,K=K_train)
weights = lasso.coef_
Y_hat = lasso.predict(X_test, K_test)
pd_Yhat=pd.DataFrame(Y_hat)
pd_Yhat.reset_index(drop=True,inplace=True)
pd_ytest=pd.DataFrame(y_test)
pd_plantings=pd.DataFrame(plantings_test)
for_plottingldakneg=pd.concat([pd_plantings,pd_ytest,pd_Yhat],axis=1)
for_plottingldakneg.columns=("Planting","Y_test","Y_hat")



#Kemmaxbn

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
N = X.shape[0]
#kf = KFold(n_splits,shuffle=True,random_state=12)
MSE_train_final = sp.zeros((n_splits,))
MSE_test_final  = sp.zeros((n_splits,))
W_nonzero_final = sp.zeros((n_splits,))
rsquared_final  = sp.zeros((n_splits,)) 
lasso.set_params(alpha=alpha_cv)


for train_index,test_index in kf12.split(X0):
	X_train, X_test, y_train, y_test, plantings_train, plantings_test= X0[train_index], X0[test_index], y[train_index], y[test_index], plantings[train_index], plantings[test_index] #Split into training and testing data sets
	K_train = Kemmaxbn[train_index][:,train_index]
	K_test  = Kemmaxbn[test_index][:,train_index]

lasso = lasso.fit(X_train,y_train,K=K_train)
weights = lasso.coef_
Y_hat = lasso.predict(X_test, K_test)
pd_Yhat=pd.DataFrame(Y_hat)
pd_Yhat.reset_index(drop=True,inplace=True)
pd_ytest=pd.DataFrame(y_test)
pd_plantings=pd.DataFrame(plantings_test)
for_plottingemmaxbn=pd.concat([pd_plantings,pd_ytest,pd_Yhat],axis=1)
for_plottingemmaxbn.columns=("Planting","Y_test","Y_hat")

#Kemmaxibs


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
N = X.shape[0]
#kf = KFold(n_splits,shuffle=True,random_state=12)
MSE_train_final = sp.zeros((n_splits,))
MSE_test_final  = sp.zeros((n_splits,))
W_nonzero_final = sp.zeros((n_splits,))
rsquared_final  = sp.zeros((n_splits,)) 
lasso.set_params(alpha=alpha_cv)


for train_index,test_index in kf12.split(X0):
	X_train, X_test, y_train, y_test, plantings_train, plantings_test= X0[train_index], X0[test_index], y[train_index], y[test_index], plantings[train_index], plantings[test_index] #Split into training and testing data sets
	K_train = Kemmaxibs[train_index][:,train_index]
	K_test  = Kemmaxibs[test_index][:,train_index]

lasso = lasso.fit(X_train,y_train,K=K_train)
weights = lasso.coef_
Y_hat = lasso.predict(X_test, K_test)
pd_Yhat=pd.DataFrame(Y_hat)
pd_Yhat.reset_index(drop=True,inplace=True)
pd_ytest=pd.DataFrame(y_test)
pd_plantings=pd.DataFrame(plantings_test)
for_plottingemmaxibs=pd.concat([pd_plantings,pd_ytest,pd_Yhat],axis=1)
for_plottingemmaxibs.columns=("Planting","Y_test","Y_hat")

for_plottingfinal=pd.concat([for_plotting0,for_plottingemmaxbn,for_plottingemmaxibs,for_plottingldakneg,for_plottingldakpos],axis=1)
for_plottingfinal=for_plottingfinal.iloc[:,[0,1,2,4,5,7,8,10,11,13,14]]


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

#plotting klimix

plt.plot(HalleFall2006.iloc[:,1],HalleFall2006.iloc[:,2],color='#e881ab',marker='o',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichSummer2006.iloc[:,1],NorwichSummer2006.iloc[:,2],color='#e881ab',marker='X',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichSummer2007.iloc[:,1],NorwichSummer2007.iloc[:,2],color='#e881ab',marker='p',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichSpring2007.iloc[:,1],NorwichSpring2007.iloc[:,2],color='#e881ab',marker='v',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichFall2006.iloc[:,1],NorwichFall2006.iloc[:,2],color='#e881ab',marker='d',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(OuluFall2007.iloc[:,1],OuluFall2007.iloc[:,2],color='#e881ab',marker='P',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(ValenciaFall2006.iloc[:,1],ValenciaFall2006.iloc[:,2],color='#e881ab',marker='*',linestyle='None',label='',mew=1,mec='black',ms=7.0)




#plotting kemmaxbn


plt.plot(HalleFall2006.iloc[:,3],HalleFall2006.iloc[:,4],color='#aae881',marker='o',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichSummer2006.iloc[:,3],NorwichSummer2006.iloc[:,4],color='#aae881',marker='X',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichSummer2007.iloc[:,3],NorwichSummer2007.iloc[:,4],color='#aae881',marker='p',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichSpring2007.iloc[:,3],NorwichSpring2007.iloc[:,4],color='#aae881',marker='v',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichFall2006.iloc[:,3],NorwichFall2006.iloc[:,4],color='#aae881',marker='d',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(OuluFall2007.iloc[:,3],OuluFall2007.iloc[:,4],color='#aae881',marker='P',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(ValenciaFall2006.iloc[:,3],ValenciaFall2006.iloc[:,4],color='#aae881',marker='*',linestyle='None',label='',mew=1,mec='black',ms=7.0)




#plotting kemmaxibs


plt.plot(HalleFall2006.iloc[:,5],HalleFall2006.iloc[:,6],color='#002a32',marker='o',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichSummer2006.iloc[:,5],NorwichSummer2006.iloc[:,6],color='#002a32',marker='X',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichSummer2007.iloc[:,5],NorwichSummer2007.iloc[:,6],color='#002a32',marker='p',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichSpring2007.iloc[:,5],NorwichSpring2007.iloc[:,6],color='#002a32',marker='v',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichFall2006.iloc[:,5],NorwichFall2006.iloc[:,6],color='#002a32',marker='d',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(OuluFall2007.iloc[:,5],OuluFall2007.iloc[:,6],color='#002a32',marker='P',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(ValenciaFall2006.iloc[:,5],ValenciaFall2006.iloc[:,6],color='#002a32',marker='*',linestyle='None',label='',mew=1,mec='black',ms=7.0)


#plotting kldakneg


plt.plot(HalleFall2006.iloc[:,7],HalleFall2006.iloc[:,8],color='#ffc6ac',marker='o',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichSummer2006.iloc[:,7],NorwichSummer2006.iloc[:,8],color='#ffc6ac',marker='X',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichSummer2007.iloc[:,7],NorwichSummer2007.iloc[:,8],color='#ffc6ac',marker='p',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichSpring2007.iloc[:,7],NorwichSpring2007.iloc[:,8],color='#ffc6ac',marker='v',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichFall2006.iloc[:,7],NorwichFall2006.iloc[:,8],color='#ffc6ac',marker='d',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(OuluFall2007.iloc[:,7],OuluFall2007.iloc[:,8],color='#ffc6ac',marker='P',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(ValenciaFall2006.iloc[:,7],ValenciaFall2006.iloc[:,8],color='#ffc6ac',marker='*',linestyle='None',label='',mew=1,mec='black',ms=7.0)


#plotting kldakpos


plt.plot(HalleFall2006.iloc[:,9],HalleFall2006.iloc[:,10],color='#f40076',marker='o',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichSummer2006.iloc[:,9],NorwichSummer2006.iloc[:,10],color='#f40076',marker='X',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichSummer2007.iloc[:,9],NorwichSummer2007.iloc[:,10],color='#f40076',marker='p',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichSpring2007.iloc[:,9],NorwichSpring2007.iloc[:,10],color='#f40076',marker='v',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(NorwichFall2006.iloc[:,9],NorwichFall2006.iloc[:,10],color='#f40076',marker='d',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(OuluFall2007.iloc[:,9],OuluFall2007.iloc[:,10],color='#f40076',marker='P',linestyle='None',label='',mew=1,mec='black',ms=7.0)
plt.plot(ValenciaFall2006.iloc[:,9],ValenciaFall2006.iloc[:,10],color='#f40076',marker='*',linestyle='None',label='',mew=1,mec='black',ms=7.0)





from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='#e881ab', lw=2),
                Line2D([0], [0], color='#aae881', lw=2),
                Line2D([0], [0], color='#002a32', lw=2),
                Line2D([0], [0], color='#ffc6ac', lw=2),
                Line2D([0], [0], color='#f40076', lw=2)
              ]

plt.legend((custom_lines),('Limix','EMMAX-BN','EMMAX-IBS','LDAK -0.5','LDAK +0.5'),loc='upper left',frameon=True,fontsize='small')
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/LMMLASSO/Figures")

plt.savefig('PredvsObs_Ks.png')
