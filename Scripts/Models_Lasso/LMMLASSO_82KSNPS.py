#LMMLASSO IN A LOOP
import limix.vardec.vardec as va q
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
K_plantings=np.load("K_PLANTINGS.npy") #K matrix calculated by limix


#Running LMM lasso (100k SNPs + covariance matrix)
import lmmlasso
#Cross validation to get the optimal parameters
alphas = 2.**(sp.linspace(2,12,10))
alphas = alphas[::-1]
from lmmlasso import runCrossValidation
lasso=lmmlasso.LmmLasso() #No need to set parameters because these will be decided through cross-validation [may need to set tolerance higher, use tol=0.05 as a baseline]
MSE_train,MSE_test,W_nonzero = lmmlasso.runCrossValidation(lasso,X,y,alphas,n_splits=10,K=K_plantings,verbose=True)



import pylab as pl

MSE_train_inter=sp.interpolate.UnivariateSpline(x=np.flip(alphas,axis=0), y=np.flip(MSE_train.mean(axis=0),axis=0)).derivative(n=2)
MSE_test_inter=sp.interpolate.UnivariateSpline(x=np.flip(alphas,axis=0), y=np.flip(MSE_test.mean(axis=0),axis=0)).derivative(n=2)
alphas_inter = 2.**(sp.linspace(2,12,100))
idx_train = sp.argmin(MSE_train_inter(alphas_inter)) 
idx_test = sp.argmin(MSE_test_inter(alphas_inter)) 
alpha_cv = (float(alphas_inter[idx_train])+float(alphas_inter[idx_test]))/2 #Selects an alpha of 4 even though idx_train and idx_test both equal zero?



from sklearn.model_selection import KFold 
kf = KFold(10,shuffle=True,random_state=None)
for train_index,test_index in kf.split(X):
	X_train, X_test, y_train, y_test= X[train_index], X[test_index], y[train_index], y[test_index] #Split into training and testing data sets
	K_train = K_plantings[train_index][:,train_index]
	K_test  = K_plantings[test_index][:,train_index] 

lasso=lmmlasso.LmmLasso(alpha=alpha_cv)
model_fit = lasso.fit(X_train,y_train,K=K_train) #fit the lmm lasso model to the whole data (shouldn't this be training data? Look into it.)
weights = model_fit.coef_ #the b




import matplotlib.pyplot as plt
# cross_val_predict returns an array of the same size as `y` where each entry
# is a prediction obtained by cross validation:
Y_hat = model_fit.predict(X_test, K_test) #create predictions for the phenotype. Should be predictions based on testing set though, not whole dataset.
fig, ax = plt.subplots()
ax.scatter(y_test, Y_hat, edgecolors=(0, 0, 0))
ax.plot([Pheno_data.min(), Pheno_data.max()], [Y_hat.min(), Y_hat.max()], 'k--', lw=4)
ax.set_xlabel('Measured')
ax.set_ylabel('Predicted')
plt.show()
plt.savefig("PredvsObs_"+PHENO_file+".png")



from sklearn.linear_model import Lasso

ls=Lasso(alpha=alpha_cv)
reg=ls.fit(X_train,y_train)
predictions=reg.predict(X_test)
fig, ax = plt.subplots()
ax.scatter(y_test, predictions, edgecolors=(0, 0, 0))
ax.plot([Pheno_data.min(), Pheno_data.max()], [predictions.min(), predictions.max()], 'k--', lw=4)
ax.set_xlabel('Measured')
ax.set_ylabel('Predicted')
plt.show()