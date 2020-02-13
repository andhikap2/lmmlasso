#Cutting up chromosome 3 
import numpy as np
import pandas as pd
from limix.stats import linear_kinship
from numpy.random import RandomState


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
chr3_1_normalised=((chr3_1_trimmed)-(np.mean(chr3_1_trimmed,axis=0)))/np.std(chr3_1_trimmed,axis=0)
K3_1=linear_kinship(chr3_1_normalised,verbose=True)


chr3_2_std=np.std(chr3_2,axis=0)
logical=chr3_2_std!=0
chr3_2_trimmed=chr3_2[:,logical] #get rid of invariant columns
chr3_2_normalised=((chr3_2_trimmed)-(np.mean(chr3_2_trimmed,axis=0)))/np.std(chr3_2_trimmed,axis=0)
K3_2=linear_kinship(chr3_2_normalised,verbose=True)


chr3_3_std=np.std(chr3_3,axis=0)
logical=chr3_3_std!=0
chr3_3_trimmed=chr3_3[:,logical] #get rid of invariant columns
chr3_3_normalised=((chr3_3_trimmed)-(np.mean(chr3_3_trimmed,axis=0)))/np.std(chr3_3_trimmed,axis=0)
K3_3=linear_kinship(chr3_3_normalised,verbose=True)

import os
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
y=np.load("DTB_NEW.npy")
import limix.vardec

#VarDec for a single K-matrix
vardecs=[]
vardec=limix.vardec.VarDec(y,"normal")
vardec.append(K3_1,"K3.1")
vardec.append(K3_2,"K3.2")
vardec.append(K3_3,"K3.3")
vardec.append_iid("noise") #add noise
vardec.fit(verbose=True)
vardecs.append(vardec)

vardecs[0].plot() #Generate the variance decomposition plot

#Splitting the 3rd split of the 3rd

chr3_3_1=np.split(chr3_3_trimmed,3,axis=1)[0] #split chromosome 3 into three arrays
chr3_3_2=np.split(chr3_3_trimmed,3,axis=1)[1]