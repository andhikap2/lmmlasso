from sklearn.linear_model import Lasso
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import cross_validate
from sklearn.model_selection import KFold ## For K-fold cross validation
from sklearn.model_selection import RepeatedKFold
import limix.vardec.vardec as va 
import limix.qtl as qtl
import scipy as sp
import time
import scipy.linalg as LA
import scipy.optimize as OPT
import pdb
import os
import numpy as np
import pandas as pd 

class LmmLasso(Lasso):
    """
    Lmm-Lasso classo
    """
    def __init__(self, alpha=1., **lasso_args):
        """
        Extension to the sklearn's LASSO to model population structure

        alpha: l1-regularization parameter
        """
        super(LmmLasso, self).__init__(alpha=alpha, **lasso_args)
        self.msg = 'lmmlasso'

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
K=np.load('K_NEWER.npy') 
K1=np.load('K_CHR1.npy')
K2=np.load('K_CHR2.npy')
K3=np.load('K_CHR3.npy')
K4=np.load('K_CHR4.npy')
K5=np.load('K_CHR5.npy')
plantings=np.load("Plantings_NEW.npy")
y=np.load("DTB_NEW.npy")


os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")
snps=np.load('SNPs_0.1.npy') #If x were given as just snps then they would be fixed effects right
env=pd.read_csv('Microclimate_Daily_threshold_0_0.csv',sep=',')
#logical=env.columns.str.startswith('PAR')
#env=env.iloc[:,~logical] #get columns that don't start with PAR
environment=np.array(env)
X=np.concatenate((snps,environment),axis=1)
assert K.shape[0] == K.shape[1], 'K MATRIX IS NOT SYMMETRICAL'
assert K.shape[0] == env.shape[0], 'NO. OF INDIVIDUALS DOES NOT MATCH K MATRIX DIMENSIONS'

from sklearn.model_selection import KFold 
n_splits=10
kf12 = KFold(n_splits,shuffle=True,random_state=12) #12 because the random state is 12
kf12.get_n_splits(X)

for train_index,test_index in kf12.split(X):
    X_train, Xstar, y_train, y_test, plantings_train, plantings_test= X[train_index], X[test_index], y[train_index], y[test_index], plantings[train_index], plantings[test_index] #Split into training and testing data sets
    K_train = K[train_index][:,train_index]
    Kstar  = K[test_index][:,train_index]

K1_train = K1[train_index][:,train_index]
K2_train = K2[train_index][:,train_index]
K3_train = K3[train_index][:,train_index]
K4_train = K4[train_index][:,train_index]
K5_train = K5[train_index][:,train_index]




if y.ndim == 2:
    assert y.shape[1]==1, 'Only one phenotype can be processed at at time.'
    y = y.flatten()


n_s=X_train.shape[0]
n_f=X_train.shape[1]
assert X.shape[0] == y.shape[0], 'dimensions do not match'
assert K.shape[0] == K.shape[1], 'dimensions do not match'
assert K.shape[0] == X.shape[0], 'dimensions do not match'

""" standardizing genotypes and phenotypes """
#if standardize:
    X -= X.mean(axis=0)
    X /= X.std(axis=0)
    y -= y.mean(axis=0)
    y /= y.std(axis=0)
0
""" training null model """
import limix


#vardec=limix.vardec.VarDec(y_train,"normal")
#vardec.append_iid("noise")
#vardec.append(K1_train, "CHR1")
#vardec.append(K2_train, "CHR2")
#vardec.append(K3_train,"CHR3")
#vardec.append(K4_train,"CHR4")
#vardec.append(K5_train,"CHR5")
#vardec.fit(verbose=True)
#delta0   = 2842.330/(151.746+ 0.002+ 610.676+ 0.055+72.454+2842.330) #based on the limix vardec results
delta0= 0.9410828164448579 #this was the delta0 from a BIGK Variance Decomposition
S,U = LA.eigh(K_train) #S = eigenvalues, U = matrix of eigenvectors




#model fitting steps [USE TRAINING DATA]


""" rotating data """
Sdi = 1. / (S + delta0) #using delta0 from multiple k-matrices 
Sdi_sqrt = sp.sqrt(Sdi)
SUX = sp.dot(U.T, X_train)
SUX = SUX * sp.tile(Sdi_sqrt, (n_f, 1)).T
SUy = sp.dot(U.T, y_train)
SUy = Sdi_sqrt * SUy

alpha=65.8177152


""" fitting lasso """
#super(LmmLasso, self).fit(SUX, SUy, alpha=10.3656406044,warm_start=True,fit_intercept=False,tol=0.025) #alpha value is based on cv results from table
lasso=Lasso(alpha=alpha,warm_start=True,fit_intercept=False,tol=0.0005)
fitting=lasso.fit(SUX,SUy)#fit lasso on rotated data
yhat=lasso.predict(X_train)
#yhat = super(LmmLasso, self).predict(X)

#self.w_ridge = LA.solve(K + delta0 * sp.eye(n_s), y - yhat)

w_ridge1= LA.solve(K1_train + K2_train + K3_train + K4_train + K5_train + delta0 * sp.eye(n_s), y_train-yhat) #get weights based on multiple k matrices

#w_ridge2= LA.solve(K + delta0 * sp.eye(n_s), y-yhat)


#fixed_effect = super(LmmLasso, self).predict(Xstar)
#fixed_effect + sp.dot(Kstar, self.w_ridge)

fixed_effect = lasso.predict(Xstar)
fixed_effect1 = fixed_effect + sp.dot(Kstar, w_ridge1) #fixed effects calculated based on multiple k matrices [these are actually predictions]
#fixed_effect2 = fixed_effect + sp.dot(Kstar, w_ridge2) #fixed effects calculated based on a single k matrix... kinda (the initial delta0 should also be based on a single k matrix tbh)

#Y_hat = lasso.predict(X_test, K_test)
pd_Yhat=pd.DataFrame(fixed_effect1)
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
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
plt.figure()


ranges=range(int(np.amax(y_test)))
plt.plot(ranges,ranges)
v=str(alpha_cv)

plot_title=['LMMLASSO, w/SNPs, KNEWER_PERCHRK, Daily0_0'] #change title as needed 
plot_title="".join(plot_title)
plt.title(plot_title)
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
plt.savefig('PredvsObs_KNEWER_PERCHRK.png')

alpha_cv_np=np.array(alpha_cv)
c=np.array(sum(weights!=0))
W_nonzero=np.array(c)
rmse_test=math.sqrt(np.mean(MSE_test_final))
r2_model=np.mean(rsquared_final)
results=np.array((alpha_cv_np,W_nonzero,rmse_test,r2_model))



















vd=va.VarianceDecomposition(y)
vd.addRandomEffect(is_noise=True)
vd.addRandomEffect(K)
vd.addFixedEffect(X)
vd.optimize(init_method='random',n_times=1000,verbose=True) 
varComps = vd.getVarianceComps()
delta0   = varComps[0,0]/varComps.sum()
self.varComps = varComps

#Too big,, probably do it on spartan



vd=va.VarianceDecomposition(y)
vd.addRandomEffect(is_noise=True)
vd.addRandomEffect(K1)
vd.addRandomEffect(K2)
vd.optimize(init_method='random',n_times=1000,verbose=True,inference='gp2KronSum') 
varComps = vd.getVarianceComps()
delta0   = varComps[0,0]/varComps.sum()
self.varComps = varComps

from limix.covar.kronecker import KronCov