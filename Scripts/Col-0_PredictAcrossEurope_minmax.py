#using rasterio to load tif files


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
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")
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
alphas_inter = 2.**(sp.linspace(2,10,100))
idx_train = sp.argmin(MSE_train_inter(alphas_inter)) 
idx_test = sp.argmin(MSE_test_inter(alphas_inter))
alpha_cv = (float(alphas_inter[idx_train])+float(alphas_inter[idx_test]))/2
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
kf12 = KFold(n_splits,shuffle=True,random_state=12) #12 because the random state is 12
kf12.get_n_splits(X)

#Model fitting 

lasso.set_params(alpha=alpha_cv) #alpha_cv=4.0
lasso=lasso.fit(X,y,K=K)


#Making prediction of DTB for Col-0 across Europe
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
ids=np.load('Ecotype_ids_NEW_ORDERED.npy')
logical=ids==6909 #Col-0's ecotype id
index=np.where(logical)
col0_snp=snps[index[0][0],]
col0_k=K[index[0][0],]

''' the temperature records start from 1950-01-01, is daily
	need to calculate the number of days since then for the following key dates:
	HalleFall2006: 277th day of 2006; 2006-10-03 ; 20729 days
	NorwichFall2006: 249th day of 2006; 2006-09-06 ; 20702 days
	NorwichSpring2007: 58th day of 2007; 2007-02-27 ; 20876 days
	NorwichSummer2006: 145th day of 2006; 2006-05-25 ; 20598 days
	NorwichSummer2007: 142nd day of 2007: 2007-05-22 ; 20960 days
	OuluFall2007: 255th day of 2007: 2007-09-12 ; 21073 days
	ValenciaFall2006: 312th day of 2006: 2006-11-8 ; 20765 days

	calculate the difference in days using R.
		DATE1=DEPENDS
		DATE2=01-JAN-1950
		as.Date(DATE1,"%d-%B-%Y") - as.Date(DATE2,"%d-%B-%Y")

'''

#Getting microclim data across Europe for a sowing date of Oct 3 2006
import rasterio
os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/E-OBSv19.0HOM")
tx=rasterio.open('tx_gdal.tif')
tn=rasterio.open('tn_gdal.tif')
arraytx=tx.read() #numpy array of values with dimensions (time, lat, long)
arraytn=tn.read()


'''
env data consists of 203 days of observations, with minimum temp first and then maximum temp
'''

day=20597 #days since Jan 1 1950 (maybe do a -1 since it's python?)
day2=day+203
maxtemp=arraytx[day:day2,:,:] #203 x 465 x 705
mintemp=arraytn[day:day2,:,:]


'''manually extracting temperature points instead of reordering'''



new_maxtemp=maxtemp.transpose((1,2,0)) # 465 x 705 x 203
new_mintemp=mintemp.transpose((1,2,0))


#flatten to 2d

new_maxtemp=new_maxtemp.reshape(-1,203) # 327825 x 203
new_mintemp=new_mintemp.reshape(-1,203) # 327825 x 203


maxtemp_bool=[]

for i in (range(new_maxtemp.shape[0])):
	boolmax=-9999 in new_maxtemp[i,:]
	maxtemp_bool.append(boolmax)

maxtemp_bool=np.array(maxtemp_bool)

#


mintemp_bool=[]

for i in (range(new_mintemp.shape[0])):
	boolmin=-9999 in new_mintemp[i,:]
	mintemp_bool.append(boolmin)

mintemp_bool=np.array(mintemp_bool)


new_maxtemp=new_maxtemp/100
new_mintemp=new_mintemp/100

#col0_snp=np.tile(col0_snp,(new_maxtemp.shape[0],1)) #memory error... might have to do the model without snps
col0_k=np.tile(col0_k,(new_maxtemp.shape[0],1))

predictors=np.concatenate((new_mintemp,new_maxtemp),axis=1)

predictions=lasso.predict(predictors,col0_k)

maxtemp_idx=np.where(maxtemp_bool)
mintemp_idx=np.where(mintemp_bool)

predictions[maxtemp_idx,]=-999
predictions[mintemp_idx,]=-999

#something going on with the reshaping

predictions_reshaped=predictions.reshape(465,-1)
predictions_reshaped_transposed=predictions.transpose().reshape(465,-1)
#toobig=predictions_reshaped>1000
#toobig=np.where(toobig)
#predictions_reshaped[toobig]=0

#toosmall=predictions_reshaped<0
#predictions_reshaped[toosmall]=0

#converting the numpy array into a raster
os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/E-OBSv19.0HOM")

naip_data_ras = arraytx
naip_meta = tx.profile

naip_meta
naip_transform = naip_meta["transform"]
naip_crs = naip_meta["crs"]


naip_meta['count'] = 1
naip_meta['dtype'] = "float64"

with rasterio.open('col0_predictionsmap_norwichsummer2006_sowdatev.tif', 'w', **naip_meta) as dst: #version 2: if any missing data, set to -9999
    dst.write(predictions_reshaped, 1)



with rasterio.open('predictionstest.tif', 'w', **naip_meta) as dst: #version 2: if any missing data, set to -9999
    dst.write(predictions_reshaped, 1)


with rasterio.open('predictionstesttransposed.tif', 'w', **naip_meta) as dst: #version 2: if any missing data, set to -9999
    dst.write(predictions_reshaped_transposed, 1)
