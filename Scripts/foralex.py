#LMMLASSO for alex

import os
import numpy as np 
import pandas as pd 
from sklearn.model_selection import KFold 


os.chdir("/data/GWAS/At/")

pc1=pd.read_csv("FIBR_PTUtB_PC1.csv")
pc2=pd.read_csv("FIBR_PTUtB_PC2.csv")
eco_ids=pc1.iloc[:,0]
eco_ids=np.array(eco_ids)

kmaf01=np.load("k2029_IBSMatrix_MAF0.1_R20.5.npy")
ids_kmaf01=pd.read_csv("k2029_IBS_MAF0.1_R20.5.prune.in",header=None)
ecoids_kmaf01=pd.read_csv("k2029_IBS_MAF0.1_R20.5.mibs.id",header=None,sep='\t')
kmaf01_pd=pd.DataFrame(kmaf01)
kmaf01_pd.index=ecoids_kmaf01.iloc[:,0]
kmaf01_pd.columns=ecoids_kmaf01.iloc[:,0]




kmaf005=np.load("k2029_IBSMatrix_MAF0.05_R20.9.npy")
ids_kmaf005=pd.read_csv("k2029_IBS_MAF0.05_R20.9.prune.in",header=None)
ecoids_kmaf005=pd.read_csv("k2029_IBS_MAF0.05_R20.9.mibs.id",header=None,sep='\t')
kmaf005_pd=pd.DataFrame(kmaf005)
kmaf005_pd.index=ecoids_kmaf005.iloc[:,0]
kmaf005_pd.columns=ecoids_kmaf005.iloc[:,0]


from pandas_plink import example_file_prefix
import limix
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/PLINK")


(bim, fam, bed) = limix.io.plink.read("k2029", verbose=False)
fid=fam.loc[:,'fid'].astype('int')
# logical=bim.loc[:,'snp'].isin(ids_kmaf01.iloc[:,0]) 
# logic=np.array(logical)
# logical_eco=fid.isin(eco_ids) #211
# logice=np.array(logical_eco)
# indices_e = np.flatnonzero(logice) #Ecotype indices

# indices = np.flatnonzero(logic) #SNP id indices
# b=bed[indices,:] #Slice by SNPs
# b=b[:,indices_e] #Slice by ecotype id

# b=b.T
# beds=b.compute()
# beds[beds==2]=1


# np.save("k2029_SNPs_MAF0.1_R20.5.npy",beds)


v=pc1.iloc[:,0].astype('int').isin(fid)
pc1=pc1.loc[v,:]
pc1=pc1.sort_values(by=['ecotype_id'])
pc2=pc2.loc[v,:]
pc2=pc2.sort_values(by=['ecotype_id'])

kfilter=kmaf01_pd.index.isin(pc1.iloc[:,0])
k01=kmaf01_pd.loc[kfilter,:]
k01=k01.loc[:,kfilter]

kfilter=kmaf005_pd.index.isin(pc1.iloc[:,0])
k005=kmaf005_pd.loc[kfilter,:]
k005=k005.loc[:,kfilter]



snps_kmaf01=np.load("k2029_SNPs_MAF0.1_R20.5.npy")
snps_kmaf005=np.load("k2029_SNPs_MAF0.05_R20.9.npy")


import scipy as sp 
alphas = 2.**(sp.linspace(-10,10,10)) #list of alphas to test
n_alphas = len(alphas)
n_splits=10
kf = KFold(n_splits,shuffle=True,random_state=None)
MSE_train = sp.zeros((n_splits,n_alphas))
MSE_test  = sp.zeros((n_splits,n_alphas))
W_nonzero = sp.zeros((n_splits,n_alphas))
rsquared  = sp.zeros((n_splits,n_alphas))
os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
import lmmlasso
lasso = lmmlasso.LmmLasso(warm_start=True,fit_intercept=False,tol=0.5) #note the tolerance value

K=np.array(k01)
y=np.array(pc1.iloc[:,1])
X=snps_kmaf01
MSE_train,MSE_test,W_nonzero, rsquared = lmmlasso.runCrossValidation(lasso,X,y,alphas,n_splits=10,K=K,verbose=True)
MSE_train_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_train.mean(axis=0))).derivative(n=2)#Interpolating the values for alphas within -2 to 10...
MSE_test_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_test.mean(axis=0))).derivative(n=2)
alphas_inter = 2.**(sp.linspace(-10,10,100))
idx_train = sp.argmin(MSE_train_inter(alphas_inter)) 
idx_test = sp.argmin(MSE_test_inter(alphas_inter))
alpha_cv = (float(alphas_inter[idx_train])+float(alphas_inter[idx_test]))/2
lasso.set_params(alpha=alpha_cv)

lasso = lasso.fit(X,y,K=K)

weights_01pc1=lasso.coef_ #all of these are zero
varcomps_01pc1=lasso.varComps
kprop_01pc1=lasso.varComps[0,1]/(lasso.varComps[0,1]+lasso.varComps[0,0]) #proportion of variance explained by the K matrix

######

MSE_train = sp.zeros((n_splits,n_alphas))
MSE_test  = sp.zeros((n_splits,n_alphas))
W_nonzero = sp.zeros((n_splits,n_alphas))
rsquared  = sp.zeros((n_splits,n_alphas))
os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
lasso = lmmlasso.LmmLasso(warm_start=True,fit_intercept=False,tol=0.5) #note the tolerance value

K=np.array(k005)
y=np.array(pc1.iloc[:,1])
X=snps_kmaf005
MSE_train,MSE_test,W_nonzero, rsquared = lmmlasso.runCrossValidation(lasso,X,y,alphas,n_splits=10,K=K,verbose=True)
MSE_train_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_train.mean(axis=0))).derivative(n=2)#Interpolating the values for alphas within -2 to 10...
MSE_test_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_test.mean(axis=0))).derivative(n=2)
alphas_inter = 2.**(sp.linspace(-10,10,100))
idx_train = sp.argmin(MSE_train_inter(alphas_inter)) 
idx_test = sp.argmin(MSE_test_inter(alphas_inter))
alpha_cv = (float(alphas_inter[idx_train])+float(alphas_inter[idx_test]))/2
lasso.set_params(alpha=alpha_cv)

lasso = lasso.fit(X,y,K=K)

weights_005pc1=lasso.coef_ #all of these are zero
varcomps_005pc1=lasso.varComps
kprop_005pc1=lasso.varComps[0,1]/(lasso.varComps[0,1]+lasso.varComps[0,0]) #proportion of variance explained by the K matrix


########


MSE_train = sp.zeros((n_splits,n_alphas))
MSE_test  = sp.zeros((n_splits,n_alphas))
W_nonzero = sp.zeros((n_splits,n_alphas))
rsquared  = sp.zeros((n_splits,n_alphas))
os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
lasso = lmmlasso.LmmLasso(warm_start=True,fit_intercept=False,tol=0.5) #note the tolerance value

K=np.array(k005)
y=np.array(pc2.iloc[:,1])
X=snps_kmaf005
MSE_train,MSE_test,W_nonzero, rsquared = lmmlasso.runCrossValidation(lasso,X,y,alphas,n_splits=10,K=K,verbose=True)
MSE_train_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_train.mean(axis=0))).derivative(n=2)#Interpolating the values for alphas within -2 to 10...
MSE_test_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_test.mean(axis=0))).derivative(n=2)
alphas_inter = 2.**(sp.linspace(-10,10,100))
idx_train = sp.argmin(MSE_train_inter(alphas_inter)) 
idx_test = sp.argmin(MSE_test_inter(alphas_inter))
alpha_cv = (float(alphas_inter[idx_train])+float(alphas_inter[idx_test]))/2
lasso.set_params(alpha=alpha_cv)

lasso = lasso.fit(X,y,K=K)

weights_005pc2=lasso.coef_ #all of these are zero
varcomps_005pc2=lasso.varComps
kprop_005pc2=lasso.varComps[0,1]/(lasso.varComps[0,1]+lasso.varComps[0,0]) #proportion of variance explained by the K matrix



#########

MSE_train = sp.zeros((n_splits,n_alphas))
MSE_test  = sp.zeros((n_splits,n_alphas))
W_nonzero = sp.zeros((n_splits,n_alphas))
rsquared  = sp.zeros((n_splits,n_alphas))
os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
lasso = lmmlasso.LmmLasso(warm_start=True,fit_intercept=False,tol=0.5) #note the tolerance value

K=np.array(k01)
y=np.array(pc2.iloc[:,1])
X=snps_kmaf01
MSE_train,MSE_test,W_nonzero, rsquared = lmmlasso.runCrossValidation(lasso,X,y,alphas,n_splits=10,K=K,verbose=True)
MSE_train_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_train.mean(axis=0))).derivative(n=2)#Interpolating the values for alphas within -2 to 10...
MSE_test_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_test.mean(axis=0))).derivative(n=2)
alphas_inter = 2.**(sp.linspace(-10,10,100))
idx_train = sp.argmin(MSE_train_inter(alphas_inter)) 
idx_test = sp.argmin(MSE_test_inter(alphas_inter))
alpha_cv = (float(alphas_inter[idx_train])+float(alphas_inter[idx_test]))/2
lasso.set_params(alpha=alpha_cv)

lasso = lasso.fit(X,y,K=K)

weights_01pc2=lasso.coef_ #all of these are zero
varcomps_01pc2=lasso.varComps
kprop_01pc2=lasso.varComps[0,1]/(lasso.varComps[0,1]+lasso.varComps[0,0]) #proportion of variance explained by the K matrix
