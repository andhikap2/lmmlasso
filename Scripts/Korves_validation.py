#Korves2006 external validation
#training the minmax model to get the model object [training on whole dataset]

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

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
plantings=np.load("Plantings_NEW.npy")

import matplotlib
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
K=np.load('K_MATRIX_EMMAXIBS.npy') #Load the preferred K-matrix
#K=np.load("K_MATRIX_PLINKIBS.npy")
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/New_Microclim")
#snps=np.load('SNPs_0.1.npy') #If x were given as just snps then they would be fixed effects right
env=pd.read_csv('Microclimate_minmaxDaily_threshold_0_0.csv',sep=',')
logical=env.columns.str.startswith(('PAR','TT','PTT','daylength') )
env=env.iloc[:,~logical] #get columns that don't start with PAR
environment=np.array(env)
#X=np.concatenate((snps,environment),axis=1)
X=environment
assert K.shape[0] == K.shape[1], 'K MATRIX IS NOT SYMMETRICAL'
assert K.shape[0] == env.shape[0], 'NO. OF INDIVIDUALS DOES NOT MATCH K MATRIX DIMENSIONS'
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
y=np.load("DTB_NEW.npy")
assert y.shape[0] == X.shape[0], 'NO. OF OBSERVATIONS DOES NOT MATCH NO. OF INDIVIDUALS'
#Running LMMLASSO
alphas = 2.**(sp.linspace(-10,10,10)) #list of alphas to test
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
alphas_inter = 2.**(sp.linspace(-10,10,100))
idx_train = sp.argmin(MSE_train_inter(alphas_inter)) 
idx_test = sp.argmin(MSE_test_inter(alphas_inter))
alpha_cv = (float(alphas_inter[idx_train])+float(alphas_inter[idx_test]))/2
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
kf12 = KFold(n_splits,shuffle=True,random_state=12) #12 because the random state is 12
kf12.get_n_splits(X)

#Model fitting 

lasso.set_params(alpha=alpha_cv) #alpha_cv=4.0
lasso=lasso.fit(X,y,K=K)

#Fast model fitting version
# alpha_cv=4.0
# lasso=lmmlasso.LmmLasso(warm_start=True,fit_intercept=False,tol=0.05)
# lasso=lasso.fit(X,y,K=K)



import os
os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Daymet")

korves_env=pd.read_csv("korves_env_2002-2003.csv",sep=',')

falldays=korves_env.iloc[300:(300+203),:]
#springdays=korves_env.iloc[(365+90):(365+90+203),:]

tmin=falldays.iloc[:,3]
tmin=np.array(tmin)

# tmin=springdays.iloc[:,3]
# tmin=np.array(tmin)

# tmax=springdays.iloc[:,2]
# tmax=np.array(tmax)


tmax=falldays.iloc[:,2]
tmax=np.array(tmax)

predictors=np.concatenate((tmin,tmax))
predictors=np.reshape(predictors,(-1,406))

os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Daymet")

korves_data=pd.read_csv("korves_accessiondata.csv",sep='\\')
korves_ecotypes=korves_data.iloc[:,1]
korves_ecotypes=korves_ecotypes.dropna() 


os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
ids=np.load('Ecotype_ids_NEW_ORDERED.npy')

logical2=korves_data.loc[:,'Ecotype_id'].isin(ids)
korves_data_filtered=korves_data.loc[logical2,:] #data only for those ecotypes that are in the k2029 dataset
filtered_ecotypes=korves_data_filtered.loc[:,'Ecotype_id']

korves_dtb=korves_data_filtered.iloc[:,12]


k_matrix=[]






for i in filtered_ecotypes:
	logical=ids==i
	if sum(logical)==0:
		break
	else:
		index=np.where(logical)
		eco_k=K[index[0][0],]
		k_matrix.append(eco_k)

k_matrix=np.array(k_matrix)

predictors2=np.repeat(predictors,k_matrix.shape[0],axis=0)
predictions=lasso.predict(predictors2,k_matrix)


import matplotlib.pyplot as plt 
plt.xlabel('Median DTB')
plt.ylabel('Predicted DTB')
ranges=range(0,250)
plt.plot(ranges,ranges)
plt.plot(korves_dtb,predictions,marker='X',color='#e881ab',linestyle='None')
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/LMMLASSO/Figures_Corrected")
plt.savefig('PredvObs_Korves2006.png')



#calculating rmse and r2
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from scipy.stats import kendalltau
from math import sqrt

import scipy


rmse = sqrt(mean_squared_error(korves_dtb, predictions))
kendall=kendalltau(korves_dtb,predictions)

def rsquared(x, y):
    """ Return R^2 where x and y are array-like."""
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    return r_value**2


rsquared(korves_dtb,predictions)
