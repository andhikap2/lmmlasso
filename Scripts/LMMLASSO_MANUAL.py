#Manual LMMLASSO using Kinship matrix and 340k snps

X=np.load("SNPs_340K.npy")
y=np.load("DTB_NEW.npy")
K=np.load("K_NEW.npy") #K matrix calculated by limix
import scipy as sp
alphas = 2.**(sp.linspace(2,20,10))
n_splits=10

from sklearn.model_selection import KFold 
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
	K_test  = K[test_index][:,train_index] #Do the same for the Kinship matrices

#Model Fitting (without standardizing) on the whole set

import limix.vardec.vardec as va 

vd = va.VarianceDecomposition(y)
vd.addRandomEffect(is_noise=True)
vd.addRandomEffect(K)
vd.optimize(init_method='random',n_times=1000,verbose=True)
varComps = vd.getVarianceComps() #The first term is the noise and the second term is the K matrix
delta0   = varComps[0,0]/varComps.sum()

import scipy.linalg as LA
S,U = LA.eigh(K) # S = eigenvalues, W = eigenvector
Sdi = 1. / (S + delta0)
Sdi_sqrt = sp.sqrt(Sdi)
SUX = sp.dot(U.T, X)
SUX = SUX * sp.tile(Sdi_sqrt, (n_f, 1)).T
SUy = sp.dot(U.T, y)
SUy = Sdi_sqrt * SUy

#Fitting to a Lasso

from sklearn.linear_model import Lasso
alpha=1
tolerance=0.5
lasso= Lasso(alpha=alpha,tol=tolerance,max_iter=10000)
model_fitted= lasso.fit(SUX, SUy)

#Prediction step
yhat=model_fitted.predict(X)
n_s=X.shape[0]
w_ridge = LA.solve(K + delta0 * sp.eye(n_s), y - yhat)
fixed_effect = lasso.predict(X)+sp.dot(K, w_ridge) #This gives you the predictions in terms of the original rather than projected coordinates


r2_score_model = r2_score(y,fixed_effect)

###############################################################################################
#Doing this with cros validation
lasso= Lasso(alpha=alpha,tol=tolerance,max_iter=10000)

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
from sklearn.metrics import mean_squared_error
import math
RMSE_test=math.sqrt(mean_squared_error(y_test,yhat_test))
RMSE_train=math.sqrt(mean_squared_error(y_train,yhat_train))



#Compare with a regular lasso

reg=lasso.fit(X_train,y_train)
train_pred=reg.predict(X_train)
test_pred=reg.predict(X_test)
RMSE_test_lasso=math.sqrt(mean_squared_error(y_test,test_pred))
RMSE_train_lasso=math.sqrt(mean_squared_error(y_train,train_pred))
























#How they do it in ccross-validation so convert the above to cross validaiton
else:
                estimator.fit(X_train,y_train,K_train)
                ytrain_star = estimator.predict(X_train,K_train)
                ytest_star  = estimator.predict(X_test, K_test)


assert self.w_ridge.shape[0]==Kstar.shape[1], 'number of training samples is not consistent.'
        assert self.coef_.shape[0]==Xstar.shape[1],   'number of SNPs is not consistent.'
        assert Xstar.shape[0]==Kstar.shape[0],   'number of test samples is not consistent.'


