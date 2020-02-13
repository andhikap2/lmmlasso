#Testing my predmap function and debugging and stuff

import predmap
from predmap import getIBS
import os
import numpy as np
import pandas as pd 
import scipy as sp
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
y=np.load("DTB_NEW.npy")



os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")
G=np.load('SNPs_0.1.npy')
os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
ids=np.load('Ecotype_ids_NEW_ORDERED.npy')
# a=list(range(G.shape[0]))
# a=np.array_split(a,30)



os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/New_Microclim")
E=pd.read_csv('Microclimate_minmaxDaily_threshold_0_0.csv')
logical=E.columns.str.startswith(('PAR','TT','PTT','daylength') )
env=E.iloc[:,~logical] #get columns that don't start with PAR, TT, PTT, daylength
environment=np.array(env)



thingy=predmap.PredMap(G,environment) #create the predmab object
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
K_emmax=np.load('K_MATRIX_EMMAXIBS.npy')
K=getIBS(G,ids)
thingy=thingy.setKinship(K)
alphas = 2.**(sp.linspace(-10,10,10)) #list of alphas to test
thingy.crossvalidate(y,alphas,n_splits=10)
thingy.fitmodel(y)


import rasterio
os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/E-OBSv19.0HOM")
tx=rasterio.open('tx_gdal.tif')
tn=rasterio.open('tn_gdal.tif')
Kt=K[22,:]
Kt=Kt.reshape(1,-1)
cc=thingy.predict(tn,tx,100,Ktest=Kt)

Xtest=G[33,]
cc2=thingy.predict(tn,tx,100,Ktest=None,X=Xtest)


#Testing to make sure the k-matrices are made properly
np.sum(((K_emmax==1)==(K==1))) #should equal 1 in all the places where it's 1 and not equal to 1 in all the same places where they're not 1...


























