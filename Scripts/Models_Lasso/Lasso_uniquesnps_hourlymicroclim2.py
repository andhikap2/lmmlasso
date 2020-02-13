#Doing the lasso but perform train-test splitting using pandas first to retain column names

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

np_Xtrain=np.array(X_train.iloc[:,-len(X_train.columns)+1:-1])
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




pd_ypred=pd.DataFrame(y2_pred)
pd_ypred.reset_index(drop=True,inplace=True)
y_test.reset_index(drop=True,inplace=True)
for_plotting=pd.concat([y_test,pd_ypred],axis=1)



#Subsetting!??

HalleFall2006=for_plotting.loc[for_plotting['Planting']=="HalleFall2006"]
NorwichSummer2006=for_plotting.loc[for_plotting['Planting']=="NorwichSummer2006"]
NorwichSummer2007=for_plotting.loc[for_plotting['Planting']=="NorwichSummer2007"]
NorwichSpring2007=for_plotting.loc[for_plotting['Planting']=="NorwichSpring2007"]
NorwichFall2006=for_plotting.loc[for_plotting['Planting']=="NorwichFall2006"]
OuluFall2007=for_plotting.loc[for_plotting['Planting']=="OuluFall2007"]
ValenciaFall2006=for_plotting.loc[for_plotting['Planting']=="ValenciaFall2006"]


#plotting
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines


ranges=range(int(max(y_test.iloc[:,1])))
plt.plot(ranges,ranges)#maybe make the range 250 
plt.title('Lasso alpha=0.001')
plt.xlabel('DTB_observed')
plt.ylabel('DTB_predicted')
plt.plot(HalleFall2006.iloc[:,1],HalleFall2006.iloc[:,2],'ro',label='HalleFall2006')
plt.plot(NorwichSummer2006.iloc[:,1],NorwichSummer2006.iloc[:,2],'yX',label='NorwichSummer2006')
plt.plot(NorwichSummer2007.iloc[:,1],NorwichSummer2007.iloc[:,2],'gp',label='NorwichSummer2007')
plt.plot(NorwichSpring2007.iloc[:,1],NorwichSpring2007.iloc[:,2],'cv',label='NorwichSpring2007')
plt.plot(NorwichFall2006.iloc[:,1],NorwichFall2006.iloc[:,2],'bd',label='NorwichFall2006')
plt.plot(OuluFall2007.iloc[:,1],OuluFall2007.iloc[:,2],'mP',label='OuluFall2007')
plt.plot(ValenciaFall2006.iloc[:,1],ValenciaFall2006.iloc[:,2],'r*',label='ValenciaFall2006')

plt.legend(loc='lower right')
plt.show()
#withinplanting correlations
HalleFall2006.iloc[:,1].corr(HalleFall2006.iloc[:,2])
NorwichSummer2006.iloc[:,1].corr(NorwichSummer2006.iloc[:,2])
NorwichSummer2007.iloc[:,1].corr(NorwichSummer2007.iloc[:,2])
NorwichSpring2007.iloc[:,1].corr(NorwichSpring2007.iloc[:,2])
NorwichFall2006.iloc[:,1].corr(NorwichFall2006.iloc[:,2])
OuluFall2007.iloc[:,1].corr(OuluFall2007.iloc[:,2])
ValenciaFall2006.iloc[:,1].corr(ValenciaFall2006.iloc[:,2])


#Cross validation
X_cv=np.array(X.iloc[:,-len(X.columns)+1:-1])
y_cv=np.array(y.iloc[:,1])

tolerance=0.05


from sklearn.linear_model import LassoCV
model_cv = LassoCV(cv=10,random_state=12,tol=tolerance).fit(X_cv, y_cv)

#alpha selected by lassocv: 0.012359598213280004



#Importing the betas to R for further analysis
coefficients=pd.DataFrame(model_fitted.coef_)
coefficients=coefficients.transpose()
a=X_train.iloc[:,-len(X_train.columns)+1:-1]
coefficients.columns=a.columns


coefficients.to_csv(r'Lasso_Betas.csv')