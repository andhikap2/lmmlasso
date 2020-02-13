'''
Making predictions using independent G and E components for the genotype Col-0

'''


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

#Fast model fitting 
alpha_cv=4.0
import lmmlasso
lasso = lmmlasso.LmmLasso(warm_start=True,fit_intercept=False,tol=0.05) #note the tolerance value
lasso.set_params(alpha=alpha_cv) #alpha_cv=4.0
lasso=lasso.fit(X,y,K=K)
weights=lasso.coef_
weights=np.reshape(weights,(-1,1))
w_ridge=lasso.w_ridge



#Getting microclim data across Europe for a sowing date of Oct 3 2006
#env data consists of 203 days of observations, with minimum temp first and then maximum temp

import rasterio
os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/E-OBSv19.0HOM")
tx=rasterio.open('tx_gdal.tif')
tn=rasterio.open('tn_gdal.tif')
arraytx=tx.read() #numpy array of values with dimensions (time, lat, long)
arraytn=tn.read()


day=20597 #days since Jan 1 1950 -1
day2=day+203
maxtemp=arraytx[day:day2,:,:] #203 x 465 x 705 (time, lat, long)
mintemp=arraytn[day:day2,:,:]

'''Reconforming'''


new_maxtemp=maxtemp.transpose((1,2,0)) # 465 x 705 x 203
new_mintemp=mintemp.transpose((1,2,0))


new_maxtemp=new_maxtemp.reshape(-1,203) # 327825 x 203
new_mintemp=new_mintemp.reshape(-1,203) # 327825 x 203
new_maxtemp=new_maxtemp/100 #double check this boolean stuff
new_mintemp=new_mintemp/100




''' Environmental Component '''

predictors=np.concatenate((new_mintemp,new_maxtemp),axis=1)
predictions_e=np.matmul(predictors,weights)
predictions_e_reshaped=predictions_e.reshape(465,-1) # (lat, long)


#converting the numpy array into a raster
os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/E-OBSv19.0HOM")

naip_meta = tx.profile

naip_transform = naip_meta["transform"]
naip_crs = naip_meta["crs"]

naip_meta['count'] = 1
naip_meta['dtype'] = "float64"

with rasterio.open('PredMap_EComponentCol-0.tif', 'w', **naip_meta) as dst: 
    dst.write(predictions_e_reshaped, 1)




''' Genetic Component'''

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
ids=np.load('Ecotype_ids_NEW_ORDERED.npy')
logical=ids==6909 #Col-0's ecotype id
index=np.where(logical)
col0_k=K[index[0][0],]
col0_k=np.tile(col0_k,(new_maxtemp.shape[0],1))
predictors_null=np.zeros((new_maxtemp.shape[0],406))

predictions_g=lasso.predict(predictors_null,col0_k)
predictions_g_reshaped=predictions_g.reshape(465,-1)
os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/E-OBSv19.0HOM")

with rasterio.open('PredMap_GComponentCol-0.tif', 'w', **naip_meta) as dst: #version 2: if any missing data, set to -9999
    dst.write(predictions_g_reshaped, 1)


#G+E component

pred_combined=predictions_g_reshaped+predictions_e_reshaped

with rasterio.open('PredMap_G+EComponentCol-0.tif', 'w', **naip_meta) as dst: #version 2: if any missing data, set to -9999
    dst.write(pred_combined, 1)






'''Compare results with a model that  calculates both components at once'''

predictions_ge=lasso.predict(predictors,col0_k)
predictions_ge_reshaped=predictions_ge.reshape(465,-1)
with rasterio.open('PredMap_Col-0.tif', 'w', **naip_meta) as dst: #version 2: if any missing data, set to -9999
    dst.write(pred_combined, 1)
