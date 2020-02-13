#Timetest k2029 in emmax versus in python

#in python:
import limix
import os 
import numpy as np 
import time
import itertools

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/PLINK")

(bim, fam, bed) = limix.io.plink.read("k2029", verbose=False)
G=bed
G=G.T
G[G==2]=1
ids=fam.iloc[:,0]
ids=np.array(ids,dtype='int')

n=G.shape[0]
p=G.shape[0]
K=np.zeros((n,p))
nsnps=G.shape[1]

start=time.time()

for i in np.unique(ids):
	for j in np.unique(ids):
		if i==j:
			prop=1
			v=i==ids
			vo=np.where(v)[0][0]
			idx_i=np.where(v)
			idx_i=np.array(idx_i).flatten()
			s=j==ids
			so=np.where(s)[0][0]
			idx_j=np.where(s)
			idx_j=np.array(idx_j).flatten()
			for t in itertools.product(idx_i,idx_j):
				print(t)
				K[t]=prop
		else:
			v=i==ids #logical for i
			vo=np.where(v)[0][0]
			idx_i=np.where(v)
			idx_i=np.array(idx_i).flatten()
			isnp=G[vo,:].compute()
			s=j==ids
			so=np.where(s)[0][0]
			idx_j=np.where(s)
			idx_j=np.array(idx_j).flatten()
			jsnp=G[so,:].compute()
			total=sum(isnp==jsnp)
			prop=total/nsnps
			for t in itertools.product(idx_i,idx_j):
				print(t)
				K[t]=prop


stop=time.time()
duration=stop-start

start=time.time()
G[vo,:].compute()
stop=time.time()
duration=stop-start