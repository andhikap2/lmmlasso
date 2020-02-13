import sys
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
#from limix_legacy.deprecated.modules import varianceDecomposition as va
import limix.vardec.vardec as va 
import limix
import limix.vardec.vardec

#Load the required data

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
plantings=np.load("Plantings_NEW.npy")
y=np.load("DTB_NEW.npy")

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")

#env=pd.read_csv('Microclimate_Daily_threshold_0.csv')
#snps=np.load("SNPs_0.1.npy")
#environment=np.array(env)

#X=np.concatenate((snps,environment),axis=1)
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")

K=np.load('K_NEWER.npy')
K1=np.load('K_CHR1.npy')
K2=np.load('K_CHR2.npy')
K3=np.load('K_CHR3.npy')
K4=np.load('K_CHR4.npy')
K5=np.load('K_CHR5.npy')






#VarDec for a single K-matrix
vardecs=[]
vardec=limix.vardec.VarDec(y,"normal")
vardec.append(K,"K")
vardec.append_iid("noise") #add noise
vardec.fit(verbose=True)
vardecs.append(vardec)

vardecs[0].plot() #Generate the variance decomposition plot






#Variance Decomposition multiple K Matrices
vardecs = []
vardec = limix.vardec.VarDec(y,"normal")

vardec.append(K1, "CHR1")
vardec.append(K2, "CHR2")
vardec.append(K3,"CHR3")
vardec.append(K4,"CHR4")
vardec.append(K5,"CHR5")
vardec.append_iid("noise")
vardec.fit(verbose=True)
vardecs.append(vardec)

vardecs[0].plot()


#Variance Decomposition for multiple K-matrices + chr 3 is split
from limix.stats import linear_kinship

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")
snps=np.load('SNPs_0.1.npy')
snp_names=pd.read_csv('SNPs_0.1.prune.in',header=None)
snp_names=pd.DataFrame.transpose(snp_names)


pd_snps=pd.DataFrame(snps)
pd_snps.columns=snp_names.iloc[0,:]

chr3=pd_snps.loc[:,pd_snps.columns.str.startswith(pat='3_')]
chr3=np.array(chr3)

np.split(chr3,3,axis=1) #split chromosome 3 into three arrays
chr3_1=np.split(chr3,3,axis=1)[0]
chr3_2=np.split(chr3,3,axis=1)[1]
chr3_3=np.split(chr3,3,axis=1)[2]


chr3_1_std=np.std(chr3_1,axis=0)
logical=chr3_1_std!=0
chr3_1_trimmed=chr3_1[:,logical] #get rid of invariant columns
K3_1=linear_kinship(chr3_1_trimmed,verbose=True)


chr3_2_std=np.std(chr3_2,axis=0)
logical=chr3_2_std!=0
chr3_2_trimmed=chr3_2[:,logical] #get rid of invariant columns
K3_2=linear_kinship(chr3_2_trimmed,verbose=True)


chr3_3_std=np.std(chr3_3,axis=0)
logical=chr3_3_std!=0
chr3_3_trimmed=chr3_3[:,logical] #get rid of invariant columns
K3_3=linear_kinship(chr3_3_trimmed,verbose=True)

vardecs=[]
vardec=limix.vardec.vardec.VarDec(y,"normal")
vardec.append(K1, "CHR1")
vardec.append(K2, "CHR2")
vardec.append(K3_1,"CHR3.1")
vardec.append(K3_2,"CHR3.2")
vardec.append(K3_3,"CHR3.3")
vardec.append(K4,"CHR4")
vardec.append(K5,"CHR5")
vardec.append_iid("noise") #add noise
vardec.fit(verbose=True)
vardecs.append(vardec)

vardecs[0].plot() #Generate the variance decomposition plot

#Other method


vd=va.VarianceDecomposition(y)
vd.addRandomEffect(is_noise=True)
vd.addRandomEffect(K1,trait_covar_type='freeform')
vd.addRandomEffect(K2,trait_covar_type='freeform')
vd.addRandomEffect(K3,trait_covar_type='freeform')
vd.addRandomEffect(K4,trait_covar_type='freeform')
vd.addRandomEffect(K5,trait_covar_type='freeform')
#vd.addFixedEffect(X)
vd.optimize()

