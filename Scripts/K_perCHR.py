#Making K-matrix for each chromosome

import numpy as np
import pandas as pd

snps=np.load('SNPs_0.1.npy')
snp_names=pd.read_csv('SNPs_0.1.prune.in',header=None)
snp_names=pd.DataFrame.transpose(snp_names)


pd_snps=pd.DataFrame(snps)
pd_snps.columns=snp_names.iloc[0,:]


#Separating the snps by chromosome
chr1=pd_snps.loc[:,pd_snps.columns.str.startswith(pat='1_')]
chr2=pd_snps.loc[:,pd_snps.columns.str.startswith(pat='2_')]
chr3=pd_snps.loc[:,pd_snps.columns.str.startswith(pat='3_')]
chr4=pd_snps.loc[:,pd_snps.columns.str.startswith(pat='4_')]
chr5=pd_snps.loc[:,pd_snps.columns.str.startswith(pat='5_')]

from limix.stats import linear_kinship
from numpy.random import RandomState

chr1=np.array(chr1)
chr2=np.array(chr2)
chr3=np.array(chr3)
chr4=np.array(chr4)
chr5=np.array(chr5)

#Need to standardise them


chr1_std=np.std(chr1,axis=1)
logical=chr1_std!=0
chr1_trimmed=chr1[:,logical] #get rid of invariant columns
chr1_normalised=((chr1_trimmed)-(np.mean(chr1_trimmed,axis=0)))/np.std(chr1_trimmed,axis=0)
K1=linear_kinship(chr1_normalised,verbose=True)
np.save('K_CHR1.npy',K1)


chr2_std=np.std(chr2,axis=0)
logical=chr2_std!=0
chr2_trimmed=chr2[:,logical] #get rid of invariant columns
chr2_normalised=((chr2_trimmed)-(np.mean(chr2_trimmed,axis=0)))/np.std(chr2_trimmed,axis=0)
K2=linear_kinship(chr2_normalised,verbose=True)
np.save('K_CHR2.npy',K2)


chr3_std=np.std(chr3,axis=0)
logical=chr3_std!=0
chr3_trimmed=chr3[:,logical] #get rid of invariant columns
chr3_normalised=((chr3_trimmed)-(np.mean(chr3_trimmed,axis=0)))/np.std(chr3_trimmed,axis=0)
K3=linear_kinship(chr3_normalised,verbose=True)
np.save('K_CHR3.npy',K3)


chr4_std=np.std(chr4,axis=0)
logical=chr4_std!=0
chr4_trimmed=chr4[:,logical] #get rid of invariant columns
chr4_normalised=((chr4_trimmed)-(np.mean(chr4_trimmed,axis=0)))/np.std(chr4_trimmed,axis=0)
K4=linear_kinship(chr4_normalised,verbose=True)
np.save('K_CHR4.npy',K4)

chr5_std=np.std(chr5,axis=0)
logical=chr5_std!=0
chr5_trimmed=chr5[:,logical] #get rid of invariant columns
chr5_normalised=((chr5_trimmed)-(np.mean(chr5_trimmed,axis=0)))/np.std(chr5_trimmed,axis=0)
K5=linear_kinship(chr5_normalised,verbose=True)
np.save('K_CHR5.npy',K5)




snps_std=np.std(snps,axis=0)
snp_logical=snps_std!=0
snps_trimmed=snps[:,snp_logical]
snps_normalised=((snps_trimmed)-(np.mean(snps_trimmed,axis=0)))/np.std(snps_trimmed,axis=0)
K_snps=linear_kinship(snps_normalised,verbose=True)



import os
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")

np.save('K_NEWER.npy',K_snps)

chr1means=np.mean(chr1,axis=0)





random = RandomState(1)
X = random.randn(4, 100)
K = linear_kinship(X, verbose=False)



K1= linear_kinship(chr1,verbose=True)
K2= linear_kinship(np.array(chr2),verbose=True)
K3= linear_kinship(np.array(chr3),verbose=True)
K4= linear_kinship(np.array(chr4),verbose=True)
K5= linear_kinship(np.array(chr5),verbose=True)

linear_kinship(chr1,K1,verbose=False)


np.isfinite(chr1).all()
np.isfinite(chr1).all()