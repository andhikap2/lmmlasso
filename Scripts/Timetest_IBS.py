#Version 1

import time
start=time.time()

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


stop=time.time()
duration=stop-start
#9392 seconds -> a little less than 3 hours

#can you convert an iterator object into a list of tuples?
list(itertools.product(idx_id,idx_j))

#do a timetest of emmax ibs versus python using the same exact k2029 matrix

#Version 2

start=time.time()

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
			prop=total/nsnps
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


stop=time.time()
duration=stop-start