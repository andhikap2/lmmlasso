#Principal Component (can only use unique snps because too big to use all snps apparently)


import os
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")

import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import pandas as pd
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
from rpy2.robjects import default_converter
from rpy2.robjects.conversion import localconverter

ro.r('load("AllPlantings_Corrected_SNPs_unique.RData")') #Load the file

with localconverter(default_converter + pandas2ri.converter) as cv:
	pd_snps=r('snps')

pd_snps.columns=pd_snps_colnames
pd_snps.iloc[0:3,0:3] #indexing by integer
snpss=pd_snps.iloc[:,1:(len(pd_snps.columns)+1)]
X=np.array(snpss)

from sklearn.decomposition import PCA
pca= PCA (n_components=250)
principalComponents=pca.fit_transform(X)
principalDf = pd.DataFrame(data = principalComponents)
print(sum(pca.explained_variance_ratio_))
principalDf=principalDf.set_index(pd_snps.iloc[:,0]) #The index of the dataframe becomes plantings


ro.r('load("Phenotypes_dtbfixed.RData")')

with localconverter(default_converter + pandas2ri.converter) as cv:
	pd_phenotypes=r('phenotypes')

DTB=pd_phenotypes.iloc[:,1]
mean_DTB=DTB.mean()






from sklearn.linear_model import LinearRegression
y=np.array(DTB)

model= LinearRegression()
model_fitted=model.fit(principalComponents,y)

intercept = model_fitted.intercept_