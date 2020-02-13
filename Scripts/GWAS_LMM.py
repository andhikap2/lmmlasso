
from datetime import *
import sys
import os
import re
import h5py
import numpy as np
import scipy as sp
import limix.qtl
from scipy import stats
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import pandas as pd
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
from rpy2.robjects import default_converter
from rpy2.robjects.conversion import localconverter

#Get the snp data
ro.r('load("AllPlantings_Corrected_SNPs_unique.RData")') #Load the file

with localconverter(default_converter + pandas2ri.converter) as cv:
	pd_snps=r('snps')
	snps=pd_snps.iloc[:,1:len(pd_snps.columns)+1]
	SNP_data=np.array(snps)

#get the phenotype data

ro.r('load("Phenotypes_dtbfixed.RData")')

with localconverter(default_converter + pandas2ri.converter) as cv:
	pd_phenotypes=r('phenotypes')
	Pheno_data=np.array(pd_phenotypes.iloc[:,1]) #Days-to-bolting phenotypes only

#include the kinship matrix

ro.r('load("Kinship_matrix.RData")')

with localconverter(default_converter + pandas2ri.converter) as cv:
	pd_K=r('K_matrix')	
	K_data=np.array(pd_K)

#or make it using limix
from limix.stats import linear_kinship
K = linear_kinship(SNP_data, verbose=True)
K_data=K

#Another way to make the kinship matrix
from numpy import dot
from limix.qc import normalise_covariance
X=SNP_data
K=dot(X, X.T)
Kn=normalise_covariance(K)

#Missing Value Threshold
Miss_Tol=.3
#Minor Allele Frequency Threshold
MAF_Tol=.05
#estimating the allele frequencies in the data
SNPsum=np.nansum(SNP_data,axis=0)
nInd=np.sum(~np.isnan(SNP_data),axis=0)
freq_hat=np.array(SNPsum,dtype="float")/(2*nInd)
mask=np.ndarray.flatten(np.array(np.all([freq_hat>MAF_Tol,freq_hat<(1-MAF_Tol)],axis=0)).astype("bool"))
SNP_data=SNP_data[:,mask]
SNP_names=SNP_names[mask,:] #Too many indices for array error
MAF=freq_hat[mask]
SNP_in=np.sum(mask.astype("int"))


#Running LMM
import limix.qtl.lmm
SNP_data=SNP_data.astype('float') #The three arrays have to be the same data type
mm_results=limix.qtl.lmm.LMM(SNP_data,Pheno_data,K_data,test='f') #run a linear mixed model. Changing the test from maximum likelihood to F-test prevents getting NaN for p-values.
p_val=mm_results.getPv() #These are NaNs..
betas=mm_results.getBetaSNP()
logp=-np.log10(p_val)
import matplotlib.pyplot as plt 
plt.plot(logp[0,:])
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
plt.ylabel('-log10(p-value)')
plt.xlabel('SNPs')
plt.title('GWAS using LMM; F-test; Unique SNPs only; kinship from popkin')
plt.show()


mm_results2=limix.qtl.lmm.LMM(SNP_data,Pheno_data,Kn,test='f') #run a linear mixed model. Changing the test from maximum likelihood to F-test prevents getting NaN for p-values.
plt.plot(logp2[0,:])
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
plt.ylabel('-log10(p-value)')
plt.xlabel('SNPs')
plt.title('GWAS using LMM; F-test; Unique SNPs only')
plt.show()


#Running LMM lasso (100k SNPs + covariance matrix)
import lmmlasso
#Cross validation to get the optimal parameters
alphas = 2.**(sp.linspace(2,12,10))
alphas = alphas[::-1]
from lmmlasso import runCrossValidation
lasso=lmmlasso.LmmLasso() #No need to set parameters because these will be decided through cross-validation [may need to set tolerance higher, use tol=0.05 as a baseline]
MSE_train,MSE_test,W_nonzero = lmmlasso.runCrossValidation(lasso,SNP_data,Pheno_data,alphas,n_splits=10,K=K_data,verbose=True)

#Then from Alex's code..
import pylab as pl

MSE_train_inter=sp.interpolate.UnivariateSpline(x=np.flip(alphas,axis=0), y=np.flip(MSE_train.mean(axis=0),axis=0)).derivative(n=2)
MSE_test_inter=sp.interpolate.UnivariateSpline(x=np.flip(alphas,axis=0), y=np.flip(MSE_test.mean(axis=0),axis=0)).derivative(n=2)
alphas_inter = 2.**(sp.linspace(2,12,100))
idx_train = sp.argmin(MSE_train_inter(alphas_inter)) #Why would you want this to be max when the goal should be to minimise both errors... maybe because it's impossible?
idx_test = sp.argmin(MSE_test_inter(alphas_inter)) 
alpha_cv = (float(alphas_inter[idx_train])+float(alphas_inter[idx_test]))/2 #The halfway point between the alpha value giving max training error and minimum testing error

pl.figure(figsize=[20,4])
plt = pl.subplot(1,3,1)
plt.plot(sp.log2(alphas),MSE_train.mean(axis=0),linewidth=2)
pl.axvline(sp.log2(alpha_cv),color='r')
pl.xlabel('log alpha')
pl.ylabel('training error')
pl.grid(True)

plt = pl.subplot(1,3,2)
plt.plot(sp.log2(alphas),MSE_test.mean(axis=0),linewidth=2)
pl.axvline(sp.log2(alpha_cv),color='r')
pl.xlabel('log alpha')
pl.ylabel('test error')
pl.grid(True)

plt = pl.subplot(1,3,3)
plt.plot(sp.log2(alphas),W_nonzero.mean(axis=0),linewidth=2)
pl.axvline(sp.log2(alpha_cv),color='r')
pl.xlabel('log alpha')
pl.ylabel('number of nonzero coefficients')
pl.grid(True)

pl.show()


lasso.setparams(alpha=alpha_cv) #Set the alpha as the ideal based on cross-validation
model_fit = lasso.fit(SNP_data,Pheno_data,K=K_data) #fit the lmm lasso model to the whole data (shouldn't this be training data? Look into it.)
weights = model_fit.coef_ #the betas

import matplotlib.pyplot as plt
# cross_val_predict returns an array of the same size as `y` where each entry
# is a prediction obtained by cross validation:
Y_hat = lasso.predict(SNP_data, K_data) #create predictions for the phenotype. Should be predictions based on testing set though, not whole dataset.
fig, ax = plt.subplots()
ax.scatter(Pheno_data, Y_hat, edgecolors=(0, 0, 0))
ax.plot([Pheno_data.min(), Pheno_data.max()], [Y_hat.min(), Y_hat.max()], 'k--', lw=4)
ax.set_xlabel('Measured')
ax.set_ylabel('Predicted')
plt.savefig("PredvsObs_"+PHENO_file+".png")


#GWAS using LMM Lasso
p_val=np.repeat(-999,weights.shape[0])
test=np.vstack((SNP_names[:,0].astype("int"),SNP_names[:,1].astype("int"),weights,p_val,MAF)) #Need SNP names hmmmhmmgmrrgmgrmgrmgrmrgmgrm #Need to calculate MAF
test=np.transpose(test)
header='Chromosome,Position,beta,p-val,MAF'
filename="GWAS_LASSO_for_"+PHENO_file+"_out.csv"
np.savetxt(filename, test, delimiter=",",header=header,fmt=["%i","%i","%s","%s","%s"])
END_LMM=datetime.now()
print('LMM-LASSO completed in '+str(END_LMM-START))


#It runs but the predictions are absolutely terrible!!!_!_!_!_!!_!_!_!_!_! 
#But I'm probably doing it properly because if I try doing regular lasso with a kinship matrix it fails!
#Try making a lmmlasso vs lasso for a subset of the population only
