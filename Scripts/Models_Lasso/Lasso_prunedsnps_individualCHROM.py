#Running Lasso on SNPs from individual chromosomes


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

#Load the snp data
snp1=np.load("SNPs_pruned_CHROM1.npy")
snp2=np.load("SNPs_pruned_CHROM2.npy")
snp3=np.load("SNPs_pruned_CHROM3.npy")
snp4=np.load("SNPs_pruned_CHROM4.npy")
snp5=np.load("SNPs_pruned_CHROM5.npy")


#Load the hourly microclim data
ro.r('load("Microclimate_hourly_scaled.RData")')

with localconverter(default_converter + pandas2ri.converter) as cv:
	pd_microclim=r('scaled_microclim') 

np_microclim=np.array(pd_microclim)


#Load the phenotypes

ro.r('load("Phenotypes_dtbfixed.RData")')
with localconverter(default_converter + pandas2ri.converter) as cv:
	pd_phenotypes=r('phenotypes')

dtb=pd_phenotypes.iloc[:,0:2] #DTB (and plantings)


#Run Lasso on Chromosome 
X=np.concatenate((snp5,np_microclim),axis=1)

from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, dtb, test_size=0.1,random_state=12)

y_trains=np.array(y_train.iloc[:,1])
y_tests=np.array(y_test.iloc[:,1])
from sklearn.linear_model import Lasso
alpha=0.001
tolerance=0.5
lasso= Lasso(alpha=alpha,tol=tolerance,max_iter=10000)
model_fitted= lasso.fit(X_train, y_trains)
y2_pred=model_fitted.predict(X_test)
from sklearn.metrics import r2_score
r2_score_model = r2_score(y_tests, y2_pred)
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
plt.title('Lasso alpha=0.001, Pruned CHROM5 SNPs')
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


