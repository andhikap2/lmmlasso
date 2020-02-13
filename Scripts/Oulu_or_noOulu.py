#!/usr/bin/env python3
#No Oulu or yes Oulu


'''APAP'''

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

r2=str(0.1)
Th=str(0)
N=str(24)

if(N==720):
	b="Monthly"
elif (N==168):
	b="Weekly"
elif (N==str(24)):
	b="Daily"
elif (N==1):
	b="Hourly"
else:
	print("RECONSIDER PERIOR CHOICE")

import matplotlib
matplotlib.use('Agg')

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
K=np.load('K_NEWER.npy')
#K=np.load('K_NEW.npy')

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")

snpfilename=["SNPs_",r2,".npy"]
snpfilename="".join(snpfilename)
snps=np.load(snpfilename) #If x were given as just snps then they would be fixed effects right

envfilename=["Microclimate_",b,"_threshold_",Th,".csv"]
envfilename="".join(envfilename)
env=pd.read_csv(envfilename,sep=',')
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

##idx_train=sp.argmin
alpha_cv = (float(alphas_inter[idx_train])+float(alphas_inter[idx_test]))/2

import pylab as pl
pl.figure(figsize=[20,4])
pls = pl.subplot(1,3,1)
pls.plot(sp.log2(alphas),MSE_train.mean(axis=0),linewidth=2)
pl.axvline(sp.log2(alpha_cv),color='r')
pl.xlabel('log alpha')
pl.ylabel('training error')
pl.grid(True)

pls = pl.subplot(1,3,2)
pls.plot(sp.log2(alphas),MSE_test.mean(axis=0),linewidth=2)
pl.axvline(sp.log2(alpha_cv),color='r')
pl.xlabel('log alpha')
pl.ylabel('test error')
pl.grid(True)

pls = pl.subplot(1,3,3)
pls.plot(sp.log2(alphas),W_nonzero.mean(axis=0),linewidth=2)
pl.axvline(sp.log2(alpha_cv),color='r')
pl.xlabel('log alpha')
pl.ylabel('number of nonzero coefficients')
pl.grid(True)

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Figures")
#pl.savefig("Alphas_SNPs_0.1_Microclimate_Daily_Temp_0_APAP.png") #save the figure. Example file name: Alphas_SNPs_0.95_Microclimate_Monthly.png
pl.savefig("Alphas_SNPs_0.1_Microclimate_Daily_Temp_0_APAP_KNEWER.png")


#Now we make the predictions
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
    rsquared_final[ifold]=sp.stats.pearsonr(ytest_star,y_test)[0]        
    ifold +=1




#Generate a figure 

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

plt.figure()


ranges=range(int(np.amax(y_test)))
plt.plot(ranges,ranges)
v=str(alpha_cv)

plot_title=['LMMLASSO_SNPs_',r2,'_Microclimate_',b,'APAP']
plot_title="".join(plot_title)
plt.title(plot_title)
plt.xlabel('DTB_observed')
plt.ylabel('DTB_predicted')
plt.plot(HalleFall2006.iloc[:,1],HalleFall2006.iloc[:,2],color='#8f0000',marker='o',linestyle='None',label='HalleFall2006')
plt.plot(NorwichSummer2006.iloc[:,1],NorwichSummer2006.iloc[:,2],color='#c54259',marker='X',linestyle='None',label='NorwichSummer2006')
plt.plot(NorwichSummer2007.iloc[:,1],NorwichSummer2007.iloc[:,2],color='#e881ab',marker='p',linestyle='None',label='NorwichSummer2007')
plt.plot(NorwichSpring2007.iloc[:,1],NorwichSpring2007.iloc[:,2],color='#ffc2f3',marker='v',linestyle='None',label='NorwichSpring2007')
plt.plot(NorwichFall2006.iloc[:,1],NorwichFall2006.iloc[:,2],color='#d890e7',marker='d',linestyle='None',label='NorwichFall2006')
plt.plot(OuluFall2007.iloc[:,1],OuluFall2007.iloc[:,2],color='#9f66e1',marker='P',linestyle='None',label='OuluFall2007')
plt.plot(ValenciaFall2006.iloc[:,1],ValenciaFall2006.iloc[:,2],color='#3a47de',marker='*',linestyle='None',label='ValenciaFall2006')
plt.figtext(0,2,rmse_string)
plt.figtext(0,1.5,r2_string)
plt.legend(loc='lower right')

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Figures")
#plt.savefig("PredvsObs_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+"_APAP.png")
plt.savefig("PredvsObs_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+"_APAP_KNEWER.png")
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/LMMLASSO_RESULTS")
#np.save("W_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+"_APAP.npy",weights) #save the coefficients
np.save("W_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+"_APAP_KNEWER.npy",weights) #save the coefficients

alpha_cv_np=np.array(alpha_cv)
c=np.array(sum(weights!=0))
W_nonzero=np.array(c)


rmse_test=math.sqrt(np.mean(MSE_test_final))
r2_model=np.mean(rsquared_final)




results=np.array((alpha_cv_np,W_nonzero,rmse_test,r2_model))
#np.savetxt("Parameters_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+"_APAP.csv",results)
np.savetxt("Parameters_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+"_APAP_KNEWER.csv",results)

  

#APNO
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")





#Making a noOulu dataset
idx=pd.read_csv('Multitrait_noOulu_indices_offset.csv',header=None) #indices for non-Oulu plantings
idx=idx.iloc[1:,1] #second column contains the indices for subsetting
idx=np.array(idx)
idx=idx.astype(int)


plantings_NO=plantings[idx]
K_NO=K[idx,:]
K_NO=K_NO[:,idx]
snps_NO=snps[idx,:]
environment_NO=environment[idx,:]
X_NO=np.concatenate((snps_NO,environment_NO),axis=1)

y_NO=y[idx]

#CV

N = X_NO.shape[0]
kf = KFold(n_splits,shuffle=True,random_state=12)
MSE_train_final = sp.zeros((n_splits,))
MSE_test_final  = sp.zeros((n_splits,))
W_nonzero_final = sp.zeros((n_splits,))
rsquared_final  = sp.zeros((n_splits,)) 
kf.get_n_splits(X)
lasso.set_params(alpha=alpha_cv)

#################################################### The model should be fitted on ALL PLANTINGS THEN FITED ON NO...

ifold = 0 
for train_index,test_index in kf.split(X):
	print(('running fold %d'%ifold)) 
	X_train, X_test, y_train, y_test, plantings_train, plantings_test= X[train_index], X[test_index], y[train_index], y[test_index], plantings[train_index], plantings[test_index]
	K_train = K[train_index][:,train_index]
	K_test  = K[test_index][:,train_index]
	model_fit=lasso.fit(X_train,y_train,K=K_train) #Fit the model to the training set which contains all plantings
	NO_logical=plantings_test!='OuluFall2007' #get idx for non-Oulu plantings in the test test
	NO_idx=np.nonzero(NO_logical)
	NO_idx=np.array(NO_idx)
	NO_idx=NO_idx.flat #flatten the tuple
	X_test_NO=X_test[NO_idx] #subset the testing set to exclude Oulu
	y_test_NO=y_test[NO_idx]
	K_test_NO=K[NO_idx][:,train_index]
	ytrain_star=model_fit.predict(X_train,K_train)
	ytest_star=model_fit.predict(X_test_NO,K_test_NO) #predict on the noOulu testing set
	MSE_train_final[ifold]=mean_squared_error(ytrain_star,y_train)
	MSE_test_final[ifold]=mean_squared_error(ytest_star,y_test_NO) #get the MSE 
	W_nonzero_final[ifold]=sp.sum(model_fit.coef_!=0)
	rsquared_final[ifold]=sp.stats.pearsonr(ytest_star,y_test_NO)[0] #get r2 value
	ifold +=1

	#Plotting...
rmse=math.sqrt(np.mean(MSE_test_final))
mean_rsquared=np.mean(rsquared_final)



for train_index,test_index in kf12.split(X):
	X_train, X_test, y_train, y_test, plantings_train, plantings_test= X[train_index], X[test_index], y[train_index], y[test_index], plantings[train_index], plantings[test_index]
	NO_logical=plantings_test!='OuluFall2007' #get idx for non-Oulu plantings in the test test
	NO_idx=np.nonzero(NO_logical)
	NO_idx=np.array(NO_idx)
	NO_idx=NO_idx.flat #flatten the tuple
	plantings_test_NO=plantings_test[NO_idx]
	X_test_NO=X_test[NO_idx]
	y_test_NO=y_test[NO_idx]
	K_train = K[train_index][:,train_index]
	K_test_NO=K[NO_idx][:,train_index]

lasso.set_params(alpha=alpha_cv)
lasso = lasso.fit(X_train,y_train,K=K_train)
weights = lasso.coef_
Y_hat = lasso.predict(X_test_NO, K_test_NO)
pd_Yhat=pd.DataFrame(Y_hat)
pd_Yhat.reset_index(drop=True,inplace=True)
pd_ytest=pd.DataFrame(y_test_NO)
pd_plantings=pd.DataFrame(plantings_test_NO) #Plantings test no
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

plt.figure()


ranges=range(int(np.amax(y_test)))
plt.plot(ranges,ranges)
v=str(alpha_cv)

plot_title=['LMMLASSO_SNPs_',r2,'_Microclimate_',b,'APNO']
plot_title="".join(plot_title)
plt.title(plot_title)
plt.xlabel('DTB_observed')
plt.ylabel('DTB_predicted')
plt.plot(HalleFall2006.iloc[:,1],HalleFall2006.iloc[:,2],color='#8f0000',marker='o',linestyle='None',label='HalleFall2006')
plt.plot(NorwichSummer2006.iloc[:,1],NorwichSummer2006.iloc[:,2],color='#c54259',marker='X',linestyle='None',label='NorwichSummer2006')
plt.plot(NorwichSummer2007.iloc[:,1],NorwichSummer2007.iloc[:,2],color='#e881ab',marker='p',linestyle='None',label='NorwichSummer2007')
plt.plot(NorwichSpring2007.iloc[:,1],NorwichSpring2007.iloc[:,2],color='#ffc2f3',marker='v',linestyle='None',label='NorwichSpring2007')
plt.plot(NorwichFall2006.iloc[:,1],NorwichFall2006.iloc[:,2],color='#d890e7',marker='d',linestyle='None',label='NorwichFall2006')
plt.plot(OuluFall2007.iloc[:,1],OuluFall2007.iloc[:,2],color='#9f66e1',marker='P',linestyle='None',label='OuluFall2007')
plt.plot(ValenciaFall2006.iloc[:,1],ValenciaFall2006.iloc[:,2],color='#3a47de',marker='*',linestyle='None',label='ValenciaFall2006')
plt.figtext(0,2,rmse_string)
plt.figtext(0,1.5,r2_string)
plt.legend(loc='lower right')

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Figures")
plt.savefig("PredvsObs_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+"_APNO.png")

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/LMMLASSO_RESULTS")
np.save("W_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+"_APNO.npy",weights) #save the coefficients

alpha_cv_np=np.array(alpha_cv)
c=np.array(sum(weights!=0))
W_nonzero=np.array(c)


rmse_test=math.sqrt(np.mean(MSE_test_final))
r2_model=np.mean(rsquared_final)



results=np.array((alpha_cv_np,W_nonzero,rmse_test,r2_model))
np.savetxt("Parameters_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+"_APNO.csv",results)



#NONO



alphas = 2.**(sp.linspace(-2,10,10)) #list of alphas to test
n_splits=10
N = X.shape[0]
kf = KFold(n_splits,shuffle=True,random_state=None)
n_alphas = len(alphas)
MSE_train = sp.zeros((n_splits,n_alphas))
MSE_test  = sp.zeros((n_splits,n_alphas))
W_nonzero = sp.zeros((n_splits,n_alphas))
rsquared  = sp.zeros((n_splits,n_alphas))
kf.get_n_splits(X_NO)

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")

import lmmlasso

lasso = lmmlasso.LmmLasso(warm_start=True,fit_intercept=False,tol=0.5) #note the tolerance value

MSE_train,MSE_test,W_nonzero, rsquared = lmmlasso.runCrossValidation(lasso,X_NO,y_NO,alphas,n_splits=10,K=K_NO,verbose=True)
MSE_train_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_train.mean(axis=0))).derivative(n=2)#Interpolating the values for alphas within -2 to 10...
MSE_test_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_test.mean(axis=0))).derivative(n=2)

alphas_inter = 2.**(sp.linspace(2,10,100))
idx_train = sp.argmin(MSE_train_inter(alphas_inter)) 
idx_test = sp.argmin(MSE_test_inter(alphas_inter))

##idx_train=sp.argmin
alpha_cv = (float(alphas_inter[idx_train])+float(alphas_inter[idx_test]))/2

import pylab as pl
pl.figure(figsize=[20,4])
pls = pl.subplot(1,3,1)
pls.plot(sp.log2(alphas),MSE_train.mean(axis=0),linewidth=2)
pl.axvline(sp.log2(alpha_cv),color='r')
pl.xlabel('log alpha')
pl.ylabel('training error')
pl.grid(True)

pls = pl.subplot(1,3,2)
pls.plot(sp.log2(alphas),MSE_test.mean(axis=0),linewidth=2)
pl.axvline(sp.log2(alpha_cv),color='r')
pl.xlabel('log alpha')
pl.ylabel('test error')
pl.grid(True)

pls = pl.subplot(1,3,3)
pls.plot(sp.log2(alphas),W_nonzero.mean(axis=0),linewidth=2)
pl.axvline(sp.log2(alpha_cv),color='r')
pl.xlabel('log alpha')
pl.ylabel('number of nonzero coefficients')
pl.grid(True)

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Figures")
pl.savefig("Alphas_SNPs_0.1_Microclimate_Daily_Temp_0_NONO.png") #save the figure. Example file name: Alphas_SNPs_0.95_Microclimate_Monthly.png


#Now we make the predictions
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
kf.get_n_splits(X)
lasso.set_params(alpha=alpha_cv)


ifold = 0 #what is ifold
for train_index,test_index in kf.split(X_NO):
    print(('running fold %d'%ifold)) 
    X_train, X_test, y_train, y_test= X_NO[train_index], X_NO[test_index], y_NO[train_index], y_NO[test_index]
    K_train = K_NO[train_index][:,train_index]
    K_test  = K_NO[test_index][:,train_index]
    model_fit=lasso.fit(X_train,y_train,K=K_train) #Fit the model to the training set
    ytrain_star=model_fit.predict(X_train,K_train)
    ytest_star=model_fit.predict(X_test,K_test)
    MSE_train_final[ifold]=mean_squared_error(ytrain_star,y_train)
    MSE_test_final[ifold]=mean_squared_error(ytest_star,y_test)
    W_nonzero_final[ifold]=sp.sum(model_fit.coef_!=0)
    rsquared_final[ifold]=sp.stats.pearsonr(ytest_star,y_test)[0]        
    ifold +=1

rmse_test_final=math.sqrt(np.mean(MSE_test_final))
mrsquared_final=np.mean(rsquared_final)

#Generate a figure 

for train_index,test_index in kf12.split(X):
	X_train, X_test, y_train, y_test, plantings_train, plantings_test= X_NO[train_index], X_NO[test_index], y_NO[train_index], y_NO[test_index], plantings_NO[train_index], plantings_NO[test_index] #Split into training and testing data sets
	K_train = K_NO[train_index][:,train_index]
	K_test  = K_NO[test_index][:,train_index]

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

plt.figure()


ranges=range(int(np.amax(y_test)))
plt.plot(ranges,ranges)
v=str(alpha_cv)

plot_title=['LMMLASSO_SNPs_',r2,'_Microclimate_',b,'APAP']
plot_title="".join(plot_title)
plt.title(plot_title)
plt.xlabel('DTB_observed')
plt.ylabel('DTB_predicted')
plt.plot(HalleFall2006.iloc[:,1],HalleFall2006.iloc[:,2],color='#8f0000',marker='o',linestyle='None',label='HalleFall2006')
plt.plot(NorwichSummer2006.iloc[:,1],NorwichSummer2006.iloc[:,2],color='#c54259',marker='X',linestyle='None',label='NorwichSummer2006')
plt.plot(NorwichSummer2007.iloc[:,1],NorwichSummer2007.iloc[:,2],color='#e881ab',marker='p',linestyle='None',label='NorwichSummer2007')
plt.plot(NorwichSpring2007.iloc[:,1],NorwichSpring2007.iloc[:,2],color='#ffc2f3',marker='v',linestyle='None',label='NorwichSpring2007')
plt.plot(NorwichFall2006.iloc[:,1],NorwichFall2006.iloc[:,2],color='#d890e7',marker='d',linestyle='None',label='NorwichFall2006')
##plt.plot(OuluFall2007.iloc[:,1],OuluFall2007.iloc[:,2],color='#9f66e1',marker='P',linestyle='None',label='OuluFall2007')##
plt.plot(ValenciaFall2006.iloc[:,1],ValenciaFall2006.iloc[:,2],color='#3a47de',marker='*',linestyle='None',label='ValenciaFall2006')
plt.figtext(0,2,rmse_string)
plt.figtext(0,1.5,r2_string)
plt.legend(loc='lower right')

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Figures")
plt.savefig("PredvsObs_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+"_NONO.png")

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/LMMLASSO_RESULTS")
np.save("W_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+"_NONO.npy",weights) #save the coefficients

alpha_cv_np=np.array(alpha_cv)
c=np.array(sum(weights!=0))
W_nonzero=np.array(c)



rmse_test=math.sqrt(np.mean(MSE_test_final))
r2_model=np.mean(rsquared_final)



results=np.array((alpha_cv_np,W_nonzero,rmse_test,r2_model))
np.savetxt("Parameters_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+"_NONO.csv",results)



#NOAP##


os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
K=np.load('K_NEW.npy')
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")
snpfilename=["SNPs_",r2,".npy"]
snpfilename="".join(snpfilename)
snps=np.load(snpfilename) #If x were given as just snps then they would be fixed effects right
envfilename=["Microclimate_",b,"_threshold_",Th,".csv"]
envfilename="".join(envfilename)
env=pd.read_csv(envfilename,sep=',')
environment=np.array(env)
X=np.concatenate((snps,environment),axis=1)
assert K.shape[0] == K.shape[1], 'K MATRIX IS NOT SYMMETRICAL'
assert K.shape[0] == env.shape[0], 'NO. OF INDIVIDUALS DOES NOT MATCH K MATRIX DIMENSIONS'
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
y=np.load("DTB_NEW.npy")


#CV

N = X.shape[0]
kf = KFold(n_splits,shuffle=True,random_state=12)
MSE_train_final = sp.zeros((n_splits,))
MSE_test_final  = sp.zeros((n_splits,))
W_nonzero_final = sp.zeros((n_splits,))
rsquared_final  = sp.zeros((n_splits,)) 
kf.get_n_splits(X)
lasso.set_params(alpha=alpha_cv) #Same alpha_cv as NONO




ifold = 0 
for train_index,test_index in kf.split(X):
	print(('running fold %d'%ifold)) 
	X_train, X_test, y_train, y_test, plantings_train, plantings_test= X[train_index], X[test_index], y[train_index], y[test_index], plantings[train_index], plantings[test_index]
	NO_logical=plantings_train!='OuluFall2007' #get idx for non-Oulu plantings in the training set
	NO_idx=np.nonzero(NO_logical)
	NO_idx=np.array(NO_idx)
	NO_idx=NO_idx.flat #flatten the tuple
	X_train_NO=X_train[NO_idx]
	y_train_NO=y_train[NO_idx]
	K_train_NO=K[NO_idx][:,NO_idx]
	K_test  = K[test_index][:,NO_idx]
	model_fit=lasso.fit(X_train_NO,y_train_NO,K=K_train_NO) #Fit the model to the training set
	ytrain_star=model_fit.predict(X_train_NO,K_train_NO)
	ytest_star=model_fit.predict(X_test,K_test)
	MSE_train_final[ifold]=mean_squared_error(ytrain_star,y_train_NO)
	MSE_test_final[ifold]=mean_squared_error(ytest_star,y_test)
	W_nonzero_final[ifold]=sp.sum(model_fit.coef_!=0)
	rsquared_final[ifold]=sp.stats.pearsonr(ytest_star,y_test)[0]        
	ifold +=1


rmse=math.sqrt(np.mean(MSE_test_final))
mean_rsquared=np.mean(rsquared_final)



#Plotting...


for train_index,test_index in kf12.split(X):
	X_train, X_test, y_train, y_test, plantings_train, plantings_test= X[train_index], X[test_index], y[train_index], y[test_index], plantings[train_index], plantings[test_index]
	NO_logical=plantings_train!='OuluFall2007' #get idx for non-Oulu plantings in the train set
	NO_idx=np.nonzero(NO_logical)
	NO_idx=np.array(NO_idx)
	NO_idx=NO_idx.flat #flatten the tuple
	plantings_train_NO=plantings_train[NO_idx]
	X_train_NO=X_train[NO_idx]
	y_train_NO=y_train[NO_idx]
	K_train_NO = K[NO_idx][:,NO_idx]
	K_test=K[test_index][:,NO_idx]


lasso.set_params(alpha=alpha_cv)
lasso = lasso.fit(X_train_NO,y_train_NO,K=K_train_NO)
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

plt.figure()


ranges=range(int(np.amax(y_test)))
plt.plot(ranges,ranges)
v=str(alpha_cv)

plot_title=['LMMLASSO_SNPs_',r2,'_Microclimate_',b,'NOAP']
plot_title="".join(plot_title)
plt.title(plot_title)
plt.xlabel('DTB_observed')
plt.ylabel('DTB_predicted')
plt.plot(HalleFall2006.iloc[:,1],HalleFall2006.iloc[:,2],color='#8f0000',marker='o',linestyle='None',label='HalleFall2006')
plt.plot(NorwichSummer2006.iloc[:,1],NorwichSummer2006.iloc[:,2],color='#c54259',marker='X',linestyle='None',label='NorwichSummer2006')
plt.plot(NorwichSummer2007.iloc[:,1],NorwichSummer2007.iloc[:,2],color='#e881ab',marker='p',linestyle='None',label='NorwichSummer2007')
plt.plot(NorwichSpring2007.iloc[:,1],NorwichSpring2007.iloc[:,2],color='#ffc2f3',marker='v',linestyle='None',label='NorwichSpring2007')
plt.plot(NorwichFall2006.iloc[:,1],NorwichFall2006.iloc[:,2],color='#d890e7',marker='d',linestyle='None',label='NorwichFall2006')
plt.plot(OuluFall2007.iloc[:,1],OuluFall2007.iloc[:,2],color='#9f66e1',marker='P',linestyle='None',label='OuluFall2007')
plt.plot(ValenciaFall2006.iloc[:,1],ValenciaFall2006.iloc[:,2],color='#3a47de',marker='*',linestyle='None',label='ValenciaFall2006')
plt.figtext(0,2,rmse_string)
plt.figtext(0,1.5,r2_string)
plt.legend(loc='lower right')

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Figures")
plt.savefig("PredvsObs_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+"_NOAP.png")

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/LMMLASSO_RESULTS")
np.save("W_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+"_NOAP.npy",weights) #save the coefficients

alpha_cv_np=np.array(alpha_cv)
c=np.array(sum(weights!=0))
W_nonzero=np.array(c)


rmse_test=math.sqrt(np.mean(MSE_test_final))
r2_model=np.mean(rsquared_final)



results=np.array((alpha_cv_np,W_nonzero,rmse_test,r2_model))
np.savetxt("Parameters_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+"_NOAP.csv",results)


