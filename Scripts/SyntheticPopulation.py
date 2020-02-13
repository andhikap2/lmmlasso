#Making a synthetic population of K-matrices

import os 
import pandas as pd 
import numpy as np 
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")


data=pd.read_csv("1001genomes_accessiondata.csv")
K=np.load('K_MATRIX_EMMAXIBS.npy') #Load the preferred K-matrix
ids=np.load('Ecotype_ids_NEW_ORDERED.npy')

'''

The E-OBS dataset covers the following areas:
	latitude: 25N - 71.5N (25.00 -- 71.50) 
	longitude: 25W - 45E (-25.00 -- 45.00)

	resolution is at 0.1 degree grids

	matrix is 465 x 705
'''


array=np.empty((K.shape[0],465,705)) #three-dimensional array where the first dimension is the K-matrix
array[:]=np.nan

'''get the list of 1001 genomes id that are in our emmax dataset'''

a=data.loc[:,'tg_ecotypeid'].isin(ids)
dataf=data.loc[a,:]

kvalues=[]

for i in dataf.loc[:,'tg_ecotypeid']:
	b=i==ids
	ind=np.where(b)
	in1=ind[0][0]
	kval=K[in1,:]
	kvalues.append(kval)

kvalues=np.array(kvalues)

coordinates=dataf.iloc[:,4:6] #get the long and lat values
coordinates=coordinates.dropna()


coordinates_conv=coordinates
coordinates_conv.iloc[:,0]=(coordinates.iloc[:,0]-25)*10 #convert latitude to 1st-dimension-coordinates
coordinates_conv.iloc[:,1]=(coordinates.iloc[:,1]+25)*10 #convert longitude to 2nd-dimension-coordinates
coordinates_conv=np.around(coordinates_conv)
coordinates_conv=coordinates_conv.reset_index() #reset the index
coordinates_conv=coordinates_conv.iloc[:,1:3]
''' Filter out coordinates outside the dataset'''
f1=coordinates_conv.iloc[:,0]<465
f2=coordinates_conv.iloc[:,1]<705

i1=np.where(f1)
i1=np.squeeze(i1)
i2=np.where(f2)
i2=np.squeeze(i2)


coordinates_conv=coordinates_conv.iloc[i1,:]
coordinates_conv=coordinates_conv.iloc[i2,:]
coordinates_conv=coordinates_conv.reset_index() #reset the index

kvalues=kvalues[i1,:]
kvalues=kvalues[i2,:]

#expand the number of cells each accession occupies? like from 0.1x0.1 to 0.3x0.3 ?


for j in coordinates_conv.index:
	lati=(coordinates_conv.iloc[j,0]).astype('int64')
	longi=(coordinates_conv.iloc[j,1]).astype('int64')
	if sum(np.isnan(array[:,lati,longi]))>0:
		array[:,lati,longi]=kvalues[j,:]
	else:
		array[:,lati,longi]=np.sum((array[:,lati,longi],kvalues[j,:]),axis=0)
		array[:,lati,longi]=array[:,lati,longi]/2

array=np.transpose(array,(1,2,0))
array=np.reshape(array,(-1,5623))





#imputing (lots of different possible methods)



from missingpy import KNNImputer
imputer=KNNImputer(n_neighbors=5,weights="distance",col_max_missing=1.0,row_max_missing=1.0)
array_imputed=imputer.fit_transform(array)



from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
imp = IterativeImputer(max_iter=100, random_state=12,n_nearest_features=10,tol=0.05,verbose=1)
imp.fit(array)  

array_imputed=imp.transform(array)
# IterativeImputer(add_indicator=False, estimator=None,
#                  imputation_order='ascending', initial_strategy='mean',
#                  max_iter=10, max_value=None, min_value=None,
#                  missing_values=nan, n_nearest_features=None,
#                  random_state=0, sample_posterior=False, tol=0.001,
#                  verbose=0)


import impyute.imputation as ii 

array_imputed=ii.cs.em(array,loops=10000) #memory error
array_imputed=ii.cs.fast_knn(array,k=3) 





'''Making a G-component "layer" '''

predictors_null=np.zeros((array_imputed.shape[0],406))
predictions_g=lasso.predict(predictors_null,array_imputed)
predictions_g_reshaped=predictions_g.reshape(465,-1)
os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/E-OBSv19.0HOM")


import rasterio
os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/E-OBSv19.0HOM")
tx=rasterio.open('tx_gdal.tif')
arraytx=tx.read() #numpy array of values with dimensions (time, lat, long)

naip_meta = tx.profile

naip_transform = naip_meta["transform"]
naip_crs = naip_meta["crs"]

naip_meta['count'] = 1
naip_meta['dtype'] = "float64"

with rasterio.open('PredMap_GComponentImputed.tif', 'w', **naip_meta) as dst: #version 2: if any missing data, set to -9999
    dst.write(predictions_g_reshaped, 1)














''' Making a synthetic "population" of K-matrices based on k2029 emmax but you only subset the columns that are in the dataset..

so it's a 2029 x 5623 matrix 

'''

import os
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/EMMAX")
import pandas as pd 
ibs=pd.read_csv('k2029_emmax.aIBS.csv',sep='\t',header=None)
tfam=pd.read_csv('k2029_emmax.tfam',sep=' ',header=None)
ecotype_ids=tfam.loc[:,0] 
ibs.columns=ecotype_ids
ibs.index=ecotype_ids

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
import numpy as np 
plantings=np.load('Plantings_NEW.npy')

ids_ordered=np.load('Ecotype_ids_NEW_ORDERED.npy') 
ids_unique=np.unique(ids_ordered)
logical=ibs.columns.isin(ids_unique.flatten()) #just keep those columns whose ecotypes we have in the dataset
ibs=ibs.loc[:,logical] 

k_matrix=pd.DataFrame(data=None,index=ibs.index,columns=ids_ordered)
for i in ibs.columns:
	for j in ibs.index:
		a=ibs.loc[j,i] #get the data value
		k_matrix.loc[j,i]=a