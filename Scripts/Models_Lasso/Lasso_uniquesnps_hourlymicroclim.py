# Ridge regression on model using unique SNPs + hourly microclimate data

import os
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")

from sklearn import linear_model
import numpy as np
#Dataset
snps=np.load('SNPs_unique.npy')
phenotypes=np.load('Phenotypes_dtbfixed.npy')
microclimate=np.load('Microclimate_hourly_scaled.npy')
X = np.concatenate((snps,microclimate),axis=1) #axis = 1 means concatenate along columns
y = phenotypes [:,3] #for cumPTT at bolt
# Train test split
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1,)

#Performing ridge
from sklearn.linear_model import Lasso
alpha=0.001
tolerance=0.05
lasso= Lasso(alpha=alpha,tol=tolerance,max_iter=10000)
model_fitted= lasso.fit(X_train, y_train)
y2_pred=model_fitted.predict(X_test)
from sklearn.metrics import r2_score
r2_score_model = r2_score(y_test, y2_pred)
print(lasso)
print("r^2 on test data : %f" % r2_score_model)


#Plotting results

import matplotlib.pyplot as plt 
plt.plot(y2_pred,y_test,'ro')
plt.show()