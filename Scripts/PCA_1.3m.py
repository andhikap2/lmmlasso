#Incremental PCA on 500k trimmed SNPs


import os
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")

import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import pandas as pd
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
from rpy2.robjects import default_converter
from rpy2.robjects.conversion import localconverter

ro.r('load("SNPs_pruned.RData")') #Load the file

with localconverter(default_converter + pandas2ri.converter) as cv:
	pd_snps=r('snps_in')

snps=pd_snps.iloc[:,1:len(pd_snps.columns)+1]

a=np.array(snps)

from sklearn.decomposition import IncrementalPCA
ipca= IncrementalPCA(n_components=125,batch_size=125) #Keep n_components and batch_size the same

pcas=ipca.fit_transform(a)
print(sum(ipca.explained_variance_ratio_))

###########################