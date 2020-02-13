#Setting the working directory

import os
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")


from sklearn import linear_model

import numpy as np
#Dataset
snps=np.load('AllPlantings_Corrected_SNPs.npy')
phenotypes=np.load('AllPlantings_Corrected_Phenotypes.npy')
X = snps
y = phenotypes [:,3] #for cumPTT at bolt
# Train test split
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=12,)

#Performing lasso
from sklearn.linear_model import Lasso
alpha=0.001
tolerance=0.001
lasso= Lasso(alpha=alpha,tol=tolerance,max_iter=10000)
model_fitted= lasso.fit(X_train, y_train)
y2_pred=model_fitted.predict(X_test)
from sklearn.metrics import r2_score
r2_score_model = r2_score(y_test, y2_pred)
print(model_fitted)
print("r^2 on test data : %f" % r2_score_model)


#Saving the trained model
import pickle
from joblib import dump, load
dump(model_fitted, 'Lasso_AllPlantings_Corrected_DTB_0.05tol.joblib') 

#Loading the trained model
clf = load('Lasso_AllPlantings_Corrected_0.05tol.joblib') 


#Performing lassoCV
from sklearn.linear_model import LassoCV
alpha=0.001
tolerance=0.001
model= Lasso(alpha=alpha,tol=tolerance,max_iter=10000)
model=LassoCV(cv=10,n_jobs=5,tol=tolerance,random_state=12,max_iter=10000)
reg = LassoCV(cv=10,n_jobs=5,tol=tolerance,random_state=12,max_iter=10000).fit(X, y)
reg.score(X, y) #Returns R2 value
reg.predict(X[:,]) #Predicts values of y using the linear model



#Loading SNP names

pickle_in = open("AllPlantings_Corrected_SNPs_colnames.pickle","rb")
snp_names = pickle.load(pickle_in)
