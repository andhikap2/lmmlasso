#K Emmax ibs 2029 * 5623

import os
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/EMMAX")
import pandas as pd 
ibs=pd.read_csv('k2029_emmax.aIBS.csv',sep='\t',header=None)
tfam=pd.read_csv('k2029_emmax.tfam',sep=' ',header=None)
ecotype_ids=tfam.loc[:,0] 
ibs.columns=ecotype_ids
ibs.index=ecotype_ids



os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
import numpy as np 
ids_ordered=np.load('Ecotype_ids_NEW_ORDERED.npy') 
ids_unique=np.unique(ids_ordered)

k_matrix=pd.DataFrame(data=None,index=ibs.index,columns=ids_ordered)

for i in ibs.columns:
	for j in ibs.index:
		a=ibs.loc[j,i] #get the data value
		k_matrix.loc[j,i]=a