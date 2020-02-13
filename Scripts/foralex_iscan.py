from limix.qtl import iscan

import os
import numpy as np 
import pandas as pd 


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
logical=bim.loc[:,'snp'].isin(ids_kmaf01.iloc[:,0]) 
logic=np.array(logical)
logical_eco=fid.isin(eco_ids) #211
logice=np.array(logical_eco)
indices_e = np.flatnonzero(logice) #Ecotype indices

indices = np.flatnonzero(logic) #SNP id indices
b=bed[indices,:] #Slice by SNPs
b=b[:,indices_e] #Slice by ecotype id

b=b.T
beds=b.compute()
beds[beds==2]=1


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

K=np.array(k01)

iscan_model=iscan(snps_kmaf01,pc1.iloc[:,1],K=K,M=None)


from limix._data import asarray as _asarray, conform_dataset, normalize_likelihood
from limix._display import session_block
from limix.qtl._assert import assert_finite
from limix.qtl._result import IScanResultFactory
from numpy_sugar.linalg import economic_qs
from xarray import concat
from numpy import asarray, empty, ones
lik = normalize_likelihood("normal")
lik_name = lik[0]
y=pc1.iloc[:,1]
M=np.ones((211,2))
G=snps_kmaf01

k01=np.array(k01)
K=k01iscan_model=iscan(G,y,K=k01,M=none)
