#GetIBS could probably be parallelized if executed from bash
'''
The argument name would be the file name, the number of processes you wanna split it over, and the ith process


So something like


getIBS.py genotype_file snps_file 30 i
Some sort of script to replace the i with range (0-(30-1))
the output of each would be a partial numpy array you could concanetae

within python 

'''
import os 
import numpy as np 
import sys 
filename=str(sys.argv[1])
n_cores=sys.argv[2] #no. of cores to split process over
nth_split=sys.argv[3] #the nth split?

G=np.load(filename)
print(G)

"""
calculate identity-by-state matrix
G: SNP data
ids: accession / individual ids
"""
n=G.shape[0]
K=np.zeros((n,n))
nsnps=G.shape[1]
for i in np.unique(ids):
	for j in np.unique(ids):
		if i==j:
			prop=1
			v=i==ids
			vo=np.where(v)[0][0]
			idx_i=np.where(v)
			idx_i=np.array(idx_i).flatten()
			isnp=G[vo,:]
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
			isnp=G[vo,:]
			s=j==ids
			so=np.where(s)[0][0]
			idx_j=np.where(s)
			idx_j=np.array(idx_j).flatten()
			jsnp=G[so,:]
			total=sum(isnp==jsnp)
			prop=total/nsnps
			for t in itertools.product(idx_i,idx_j):
				print(t)
				K[t]=prop
return K 