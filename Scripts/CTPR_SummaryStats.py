#Making CTPRL summary statistics file
import numpy as np 
import os 
import pandas as pd

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")
snps=np.load("SNPs_0.1.npy")
snp_ids=pd.read_csv("SNPs_0.1.prune.in",header=None)

#Calculating minor allele frequencies

from collections import Counter

frequencies=np.zeros((86760,1),dtype=float)
for i in range(snps.shape[1]):
    count = Counter(snps[:, i])
    count_0=count[0]
    freq0=count[0]/(count[0]+count[1]) #frequency of the zero allele
    frequencies[i,0]=freq0


logical=frequencies[:,0]>=0.5
logical_index=np.nonzero(logical)[0]

for i in logical_index:
	count=Counter(snps[:,i])
	count_0=count[0]
	freq0=1-(count[0]/(count[0]+count[1]))
	frequencies[i,0]=freq0

sum(frequencies[:,0]>0.5)==0 #check...



#calculating marker effects for each phenotype (?
'''how would you do this... in the lmmlasso framework you used the kinships to account for those effects... but in this case if you did that all the effects would be 0
can CTPRL account for random effects and non-genetic covariates?

> just try getting the betas using regular lasso first I guess, and get standard error estimates by getting betas across the 10 folds :)
> the order of phenotypes is DTB - LfL - SeedNum
> regress the phenotype on a single SNP + K-matrix, get the beta and standard error estimates of that
'''

#getting dtb betas including standard errors 

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
K=np.load("K_MATRIX_EMMAXIBS.npy")
y=np.load("DTB_NEW.npy")
# os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/New_Microclim")
# env=pd.read_csv('Microclimate_Hourly_threshold_7_0.csv',sep=',') #the best model overall 
# environment=np.array(env)

import lmmlasso
import scipy as sp 
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
import scipy.stats

lasso = lmmlasso.LmmLasso(warm_start=True,fit_intercept=False,tol=0.5) #note the tolerance value

#Recalculate alpha for each SNP or not?


alphas = 2.**(sp.linspace(-10,10,10)) #list of alphas to test


snp_weights=np.empty((1,snps.shape[1]))
snp_errors=np.empty((1,snps.shape[1]))


for i in range(snps.shape[1]):
	X=snps[:,i]
	X=np.reshape(X,(-1,1))
	MSE_train,MSE_test,W_nonzero, rsquared = lmmlasso.runCrossValidation(lasso,X,y,alphas,n_splits=10,K=K,verbose=True) #X needs to have dimension (5623, 1)
	MSE_train_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_train.mean(axis=0))).derivative(n=2)#Interpolating the values for alphas within -2 to 10...
	MSE_test_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_test.mean(axis=0))).derivative(n=2)
	alphas_inter = 2.**(sp.linspace(-10,10,100))
	idx_train = sp.argmin(MSE_train_inter(alphas_inter)) 
	idx_test = sp.argmin(MSE_test_inter(alphas_inter))
	alpha_cv = (float(alphas_inter[idx_train])+float(alphas_inter[idx_test]))/2
	os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
	kf12 = KFold(n_splits,shuffle=True,random_state=12) #12 because the random state is 12
	lasso.set_params(alpha=alpha_cv) 
	lasso = lasso.fit(X,y,K=K)
	weights = lasso.coef_
	snp_weights[1,i]=weights[0]
	predictions=lasso.predict(X,K)
	residuals=y-predictions
	residuals=residuals**2
	df=y.shape[0]-2 #degrees of freedom
	xbar=np.mean(X)
	se_square=(
	(np.sum(residuals))
	/
	((df-2)*(np.sum((X-xbar)**2)))
	)
	se=np.sqrt(se_square)
	snp_errors[i]=se






# MSE_train,MSE_test,W_nonzero, rsquared = lmmlasso.runCrossValidation(lasso,X,y,alphas,n_splits=10,K=K,verbose=True) #X needs to have dimension (5623, 1)
# MSE_train_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_train.mean(axis=0))).derivative(n=2)#Interpolating the values for alphas within -2 to 10...
# MSE_test_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_test.mean(axis=0))).derivative(n=2)
# alphas_inter = 2.**(sp.linspace(-10,10,100))
# idx_train = sp.argmin(MSE_train_inter(alphas_inter)) 
# idx_test = sp.argmin(MSE_test_inter(alphas_inter))
# alpha_cv = (float(alphas_inter[idx_train])+float(alphas_inter[idx_test]))/2
# os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
# kf12 = KFold(n_splits,shuffle=True,random_state=12) #12 because the random state is 12
# lasso.set_params(alpha=alpha_cv) 
# lasso = lasso.fit(X,y,K=K)
# weights = lasso.coef_

# predictions=lasso.predict(X,K)
# residuals=y-predictions
# residuals=residuals**2
# df=y.shape[0]-2 #degrees of freedom
# xbar=np.mean(X)

# se_square=(
# 	(np.sum(residuals))
# 	/
# 	((df-2)*(np.sum((X-xbar)**2)))
# 	)





''' standard error of a coefficient
is sqrt[(sum of (residuals^2))/(n-2)(sum of (xi - xmean)^2)


X=np.concatenate((snps,environment),axis=1)
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")

import scipy as sp 
alphas = 2.**(sp.linspace(-10,10,50)) #list of alphas to test
from sklearn.linear_model import Lasso 
from sklearn.linear_model import LassoCV
from sklearn.model_selection import KFold




betas=np.zeros((10,X.shape[1]),dtype=float)

kf=KFold(n_splits=10)

for train_index, test_index in kf.split(X):
	X_train,X_test= X[train_index], X[test_index]
	y_train,y_test= y[train_index], y[test_index]
	modelfit=Lasso(alpha=alpha_cv,tol=0.5,random_state=None).fit(X_train,y_train)
	coefficients=modelfit.coef_
	betas[i,:]=coefficients

'''need to find mean and standard error of betas 


















os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/PLINK")
frequencies=pd.read_csv("plink.frq",sep="\t")


logical=frequencies.loc[:,'SNP'].isin(snp_ids.iloc[:,0])

logical=bim.loc[:,'snp'].isin(pd_ids.iloc[:,0]) 
