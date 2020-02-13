#LMMLASSO 82k SNPs
import limix.vardec.vardec as va 
import scipy.linalg as LA
from sklearn.linear_model import Lasso
import numpy as np 
from sklearn.metrics import mean_squared_error
import math
from sklearn.model_selection import KFold 
import gc
import scipy as sp

X=np.load("SNPs_82K.npy")
y=np.load("DTB_NEW.npy")
K=np.load("K_PLANTINGS_STANDARDIZED.npy") #K matrix calculated by limix
alphas = 2.**(sp.linspace(2,20,10))
n_splits=10

N = X.shape[0]
kf = KFold(n_splits,shuffle=True,random_state=None)
n_alphas = len(alphas)
MSE_train = sp.zeros((n_splits,n_alphas))
MSE_test  = sp.zeros((n_splits,n_alphas))
W_nonzero = sp.zeros((n_splits,n_alphas))
kf.get_n_splits(X)

for train_index,test_index in kf.split(X):
	X_train, X_test, y_train, y_test= X[train_index], X[test_index], y[train_index], y[test_index] #Split into training and testing data sets
	K_train = K[train_index][:,train_index]
	K_test  = K[test_index][:,train_index] 
###############################################################################################

gc.collect()
###############################################################################################
#Doing this with cross validation
"""Lasso Parameters"""
alpha=0.0001
tolerance=0.05
lasso= Lasso(alpha=alpha,tol=tolerance,max_iter=10000,warm_start=True)

#1. Fitting on the training dataset
[n_s, n_f] = X_train.shape
vd = va.VarianceDecomposition(y_train)
vd.addRandomEffect(is_noise=True)
vd.addRandomEffect(K_train)
vd.optimize(init_method='default',n_times=1000,verbose=True) 
varComps = vd.getVarianceComps()
delta0   = varComps[0,0]/varComps.sum()
S,U = LA.eigh(K_train)


""" rotating data """
Sdi = 1. / (S + delta0)
Sdi_sqrt = sp.sqrt(Sdi)
SUX_train = sp.dot(U.T, X_train)
SUX_train = SUX_train * sp.tile(Sdi_sqrt, (n_f, 1)).T
SUy_train = sp.dot(U.T, y_train)
SUy_train = Sdi_sqrt * SUy_train
training_fit=lasso.fit(SUX_train,SUy_train)
ytrain_star=training_fit.predict(X_train) #This is to calculate w_ridge
w_ridge = LA.solve(K_train + delta0 * sp.eye(n_s), y_train - ytrain_star) #Betas
yhat_train= (training_fit.predict(X_train))+sp.dot(K_train, w_ridge) #This is the predicted DTBs based on the training data


#Making predictions on testng set
ytest_star=training_fit.predict(X_test) #Not sure what the point of this is
yhat_test = training_fit.predict(X_test)+sp.dot(K_test, w_ridge) #This gives you the predictions in terms of the original rather than projected coordinatesfor the testing data


#Calculate MSE
RMSE_test=math.sqrt(mean_squared_error(y_test,yhat_test))
RMSE_train=math.sqrt(mean_squared_error(y_train,yhat_train))



#Compare with a regular lasso

reg=lasso.fit(X_train,y_train)
train_pred=reg.predict(X_train)
test_pred=reg.predict(X_test)
RMSE_test_lasso=math.sqrt(mean_squared_error(y_test,test_pred))
RMSE_train_lasso=math.sqrt(mean_squared_error(y_train,train_pred))








