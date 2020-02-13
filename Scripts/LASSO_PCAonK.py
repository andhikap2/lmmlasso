#Regular lasso using PCs derived from K-matrix
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
import matplotlib


os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
plantings=np.load("Plantings_NEW.npy")
K=np.load("K_MATRIX_EMMAXIBS.npy")
y=np.load("DTB_NEW.npy")

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")
snps=np.load('SNPs_0.1.npy') #If x were given as just snps then they would be fixed effects right



os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/New_Microclim")
env=pd.read_csv('Microclimate_minmaxDaily_threshold_0_0.csv',sep=',')
logical=env.columns.str.startswith(('PAR','TT','PTT','daylength') )
env=env.iloc[:,~logical] #get columns that don't start with PAR
environment=np.array(env)

from sklearn.decomposition import PCA
pca1=PCA(0.95,svd_solver='full',whiten=True)

pca1.fit(K)

comps1=pca1.transform(K) #125


#Just take the first 2 PCs...


X1=np.concatenate((snps,comps1[:,0:2],environment),axis=1)

#model fitting
from sklearn.linear_model import LassoCV
from sklearn.linear_model import Lasso
lasso=LassoCV(tol=0.5,cv=10,n_jobs=-1)
reg1=lasso.fit(X1,y)

alpha_cv=reg1.alpha_
model=Lasso(alpha=alpha_cv,tol=0.5)
from sklearn.model_selection import KFold 
kf=KFold(n_splits=10)

import scipy as sp 
from sklearn.metrics import mean_squared_error

n_splits=10


MSE_train_final = sp.zeros((n_splits,))
MSE_test_final  = sp.zeros((n_splits,))
W_nonzero_final = sp.zeros((n_splits,))
rsquared_final  = sp.zeros((n_splits,)) 
kendall_final   = sp.zeros((n_splits,))

ifold=0
for train_index, test_index in kf.split(X1):
	X1_train=X1[train_index,:]
	X1_test=X1[test_index,:]
	y_train=y[train_index]
	y_test=y[test_index]
	model_fit=model.fit(X1_train,y_train) #Fit the model to the training set
	ytrain_star=model_fit.predict(X1_train)
	ytest_star=model_fit.predict(X1_test)
	MSE_train_final[ifold]=mean_squared_error(ytrain_star,y_train)
	MSE_test_final[ifold]=mean_squared_error(ytest_star,y_test)
	W_nonzero_final[ifold]=sp.sum(model_fit.coef_!=0)
	slope, intercept, r_value, p_value, std_err = sp.stats.linregress(y_test, ytest_star)
	rsquared_final[ifold]=r_value**2
	kendall_final[ifold]=sp.stats.kendalltau(y_test,ytest_star)[0]
	ifold +=1

rmse_test=math.sqrt(np.mean(MSE_test_final))
r2_model=np.mean(rsquared_final)
kendall_model=np.mean(kendall_final)
percent_explained=np.sum(pca1.explained_variance_ratio_[0:2])



kf12=KFold(n_splits=10,shuffle=True,random_state=12)
for train_index,test_index in kf12.split(X1):
	X_train, X_test, y_train, y_test, plantings_train, plantings_test= X1[train_index], X1[test_index], y[train_index], y[test_index], plantings[train_index], plantings[test_index] #Split into training and testing data sets

model_fitted = model.fit(X_train,y_train)
Y_hat = lasso.predict(X_test)
pd_Yhat=pd.DataFrame(Y_hat)
pd_Yhat.reset_index(drop=True,inplace=True)
pd_ytest=pd.DataFrame(y_test)
pd_plantings=pd.DataFrame(plantings_test)
for_plotting=pd.concat([pd_plantings,pd_ytest,pd_Yhat],axis=1)
for_plotting.columns=("Planting","Y_test","Y_hat")

HalleFall2006=for_plotting.loc[for_plotting['Planting']=="HalleFall2006"]
NorwichSummer2006=for_plotting.loc[for_plotting['Planting']=="NorwichSummer2006"]
NorwichSummer2007=for_plotting.loc[for_plotting['Planting']=="NorwichSummer2007"]
NorwichSpring2007=for_plotting.loc[for_plotting['Planting']=="NorwichSpring2007"]
NorwichFall2006=for_plotting.loc[for_plotting['Planting']=="NorwichFall2006"]
OuluFall2007=for_plotting.loc[for_plotting['Planting']=="OuluFall2007"]
ValenciaFall2006=for_plotting.loc[for_plotting['Planting']=="ValenciaFall2006"]

import matplotlib.pyplot as plt
matplotlib.use('agg')

plt.figure()


ranges=range(int(np.amax(y_test)))
plt.plot(ranges,ranges)
v=str(alpha_cv)

#plot_title=['LMMLASSO, w/SNPs, KNEWER_BIGK, Daily0_0'] #change title as needed 
#plot_title="".join(plot_title)
#plt.title(plot_title)
plt.xlabel('Observed DTB')
plt.ylabel('Predicted DTB')
plt.plot(HalleFall2006.iloc[:,1],HalleFall2006.iloc[:,2],color='#8f0000',marker='o',linestyle='None',label='HalleFall2006')
plt.plot(NorwichSummer2006.iloc[:,1],NorwichSummer2006.iloc[:,2],color='#c54259',marker='X',linestyle='None',label='NorwichSummer2006')
plt.plot(NorwichSummer2007.iloc[:,1],NorwichSummer2007.iloc[:,2],color='#e881ab',marker='p',linestyle='None',label='NorwichSummer2007')
plt.plot(NorwichSpring2007.iloc[:,1],NorwichSpring2007.iloc[:,2],color='#ffc2f3',marker='v',linestyle='None',label='NorwichSpring2007')
plt.plot(NorwichFall2006.iloc[:,1],NorwichFall2006.iloc[:,2],color='#d890e7',marker='d',linestyle='None',label='NorwichFall2006')
plt.plot(OuluFall2007.iloc[:,1],OuluFall2007.iloc[:,2],color='#9f66e1',marker='P',linestyle='None',label='OuluFall2007')
plt.plot(ValenciaFall2006.iloc[:,1],ValenciaFall2006.iloc[:,2],color='#3a47de',marker='*',linestyle='None',label='ValenciaFall2006')
#plt.figtext(0,2,rmse_string)
#plt.figtext(0,1.5,r2_string)
plt.legend(loc='upper left',frameon=True,fontsize='small')
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/LASSO/Figures")
plt.savefig('PredvsObs_SNP0.1_minmax_KEMMAXIBS_2pc_whiten.png')















# model.fit(X1_train,y_train)
# predictions=model.predict(X1_test)





# reg2=lasso.fit(X2,y)
# reg3=lasso.fit(X3,y)

