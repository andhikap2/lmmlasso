#Decomposing lasso variance components into microclim and snp 

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

ro.r('load("AllPlantings_Corrected_SNPs_unique_colnames.RData")')

with localconverter(default_converter + pandas2ri.converter) as cv:
	pd_snps_colnames=r('snp_colnames')


pd_snps.columns=pd_snps_colnames
pd_snps.iloc[0:3,0:3] #indexing by integer

ro.r('load("Microclimate_hourly_scaled.RData")')

with localconverter(default_converter + pandas2ri.converter) as cv:
	pd_microclim=r('scaled_microclim') # loads as a numpy array automatically for some reason

import pandas as pd
pd_microclim=pd.DataFrame(pd_microclim)


ro.r('load("Microclimate_hourly_scaled_colnames.RData")')

with localconverter(default_converter + pandas2ri.converter) as cv:
	pd_microclim_colnames=r('hourly_colnames')

pd_microclim.columns=pd_microclim_colnames
pd_microclim.iloc[0:3,0:3] #indexing by integer


ro.r('load("Phenotypes_dtbfixed.RData")')
with localconverter(default_converter + pandas2ri.converter) as cv:
	pd_phenotypes=r('phenotypes')

###
pd_snps.reset_index(drop=True,inplace=True)
pd_microclim.reset_index(drop=True,inplace=True)
X=(pd.concat([pd_snps,pd_microclim],axis=1))
y=pd_phenotypes.iloc[:,0:2] #DTB (and plantings)

###

from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1,random_state=12)

np_Xtrain=np.array(X_train.iloc[:,(-len(X_train.columns)+1):-1])
np_Xtest=np.array(X_test.iloc[:,-len(X_test.columns)+1:-1])
np_ytest=np.array(y_test.iloc[:,1])
np_ytrain=np.array(y_train.iloc[:,1])



#Performing lasso
from sklearn.linear_model import Lasso
alpha=0.001
tolerance=0.5
lasso= Lasso(alpha=alpha,tol=tolerance,max_iter=10000)
model_fitted= lasso.fit(np_Xtrain, np_ytrain)
y2_pred=model_fitted.predict(np_Xtest)
from sklearn.metrics import r2_score
r2_score_model = r2_score(np_ytest, y2_pred)
print(lasso)
print("r^2 on test data : %f" % r2_score_model)

weights=model_fitted.coef_
intercept=model_fitted.intercept_

#Need to convert the weights into pandas then assign them column names then separate them into snps and microclim
betas=pd.DataFrame(weights)
betas=pd.DataFrame.transpose(betas)
betas.columns=X_train.columns[-len(X_train.columns)+1:-1]
snps=betas.loc[:, betas.columns.str.startswith(("1","2","3","4","5"))]
microclim=betas.loc[:,betas.columns.str.startswith(("cum","Grnd"))]
X_test_snps=X_test.loc[:,X_test.columns.str.startswith(("1","2","3","4","5"))]
X_test_snps=X_test_snps.astype(float)
X_test_microclim=X_test.loc[:,betas.columns.str.startswith(("cum","Grnd"))]
X_test_microclim=X_test_microclim.astype(float)

a=np.array(X_test_snps)
b=np.array(snps)
b=np.transpose(b)
c=np.array(X_test_microclim)
d=np.array(microclim)
d=np.transpose(d)

y_pred_snps= intercept + (np.dot(a,b))
y_pred_microclim= intercept + (np.dot(c,d))

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

ranges=range(250)
plt.plot(ranges,ranges)
plt.plot(np_ytest,y_pred_snps,'ro',label="SNP data only")
plt.plot(np_ytest,y_pred_microclim,'bX',label="Microclim data only")
plt.plot(np_ytest,y2_pred,'g*',label="SNP + Microclim")
plt.title('Variance partitioning?')
plt.xlabel('DTB_observed')
plt.ylabel('DTB_predicted')
plt.legend(loc='lower right')
plt.show()

r2_score(y_pred_snps,y_test.iloc[:,1])
r2_score(y_pred_microclim,y_test.iloc[:,1])