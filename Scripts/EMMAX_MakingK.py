#Conforming the EMMAX IBS matrix 

import os
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/EMMAX")
import pandas as pd 
ibs=pd.read_csv('k2030_emmax.aIBS.csv',sep='\t')
tfam=pd.read_csv('k2030_emmax.tfam',sep=' ',header=None)
ecotype_ids=tfam.loc[:,0] 
ibs.columns=ecotype_ids
ibs.index=ecotype_ids

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
import numpy as np 
plantings=np.load('Plantings_NEW.npy')

ids_ordered=np.load('Ecotype_ids_NEW_ORDERED.npy') 
ids_unique=np.unique(ids_ordered)
logical=ibs.columns.isin(ids_unique.flatten()) #just keep those columns whose ecotypes we have in the dataset
ibs=ibs.loc[logical,logical] 

k_matrix=pd.DataFrame(data=None,index=ids_ordered,columns=ids_ordered)

for i in ibs.columns:
	for j in ibs.index:
		a=ibs.loc[j,i] #get the data value
		k_matrix.loc[j,i]=a

pd.write_csv('k2030_emmaxibsmatrix.csv',k_matrix)
k_np=np.array(k_matrix)
np.save('k2030_emmaxibsmatrix.npy',k_matrix)

#Conforming the EMMAX aBN (Balding-Nichols) matrix 
import os
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/EMMAX")
import pandas as pd 
ibs=pd.read_csv('k2029_emmax.aBN.csv',sep='\t',header=None)
tfam=pd.read_csv('k2029_emmax.tfam',sep=' ',header=None)
ecotype_ids=tfam.loc[:,0] 
ibs.columns=ecotype_ids
ibs.index=ecotype_ids

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
import numpy as np 
plantings=np.load('Plantings_NEW.npy')

ids_ordered=np.load('Ecotype_ids_NEW_ORDERED.npy') 
ids_unique=np.unique(ids_ordered)
logical=ibs.columns.isin(ids_unique.flatten()) #just keep those columns whose ecotypes we have in the dataset
ibs=ibs.loc[logical,logical] 

k_matrix=pd.DataFrame(data=None,index=ids_ordered,columns=ids_ordered)

for i in ibs.columns:
	for j in ibs.index:
		a=ibs.loc[j,i] #get the data value
		k_matrix.loc[j,i]=a

'''
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
from rpy2.robjects import default_converter
from rpy2.robjects.conversion import localconverter

ro.r('load("Multitrait_NEW.RData")') #Load the file

with localconverter(default_converter + pandas2ri.converter) as cv:
    phenotypes=r('Multitrait_NEW')

Ecotype_ids_NEW_ORDERED=phenotypes.iloc[:,1]
Ecotype_ids_NEW_ORDERED_np=np.array(Ecotype_ids_NEW_ORDERED,dtype=int)
np.save('Ecotype_ids_NEW_ORDERED.npy',Ecotype_ids_NEW_ORDERED_np)
'''