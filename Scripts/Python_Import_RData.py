#Setting the working directory

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
ro.r('load("AllPlantings_Corrected_SNPs.RData")') #Load the file

###


with localconverter(default_converter + pandas2ri.converter) as cv:
	pd_scaled_microclimate=r('scaled_microclimate')
#Saving SNP data

import numpy as np
np_AllPlantings_Corrected_SNPs=np.array(pd_AllPlantings_Corrected_SNPs)
np.save('AllPlantings_Corrected_SNPs.npy',np_AllPlantings_Corrected_SNPs)

#Saving SNP column names
import pickle
AllPlantings_Corrected_SNPs_colnames = list(pd_AllPlantings_Corrected_SNPs.columns)
pickle_out=open("AllPlantings_Corrected_SNPs_colnames.pickle","wb")
pickle.dump(AllPlantings_Corrected_SNPs_colnames, pickle_out)
pickle_out.close()

#Loading the column names
pickle_in = open("AllPlantings_Corrected_SNPs_colnames.pickle","rb")
example_dict = pickle.load(pickle_in)



#Saving Phenotype data

ro.r('load("AllPlantings_Corrected_Phenotypes.RData")') #Load the file

with localconverter(default_converter + pandas2ri.converter) as cv:
    pd_AllPlantings_Corrected_Phenotypes = r('AllPlantings_Corrected_Phenotypes')


np_AllPlantings_Corrected_Phenotypes=np.array(pd_AllPlantings_Corrected_Phenotypes)
np.save('AllPlantings_Corrected_Phenotypes.npy',np_AllPlantings_Corrected_Phenotypes)


#Saving Phenotype column names
import pickle
AllPlantings_Corrected_Phenotypes_colnames = list(pd_AllPlantings_Corrected_Phenotypes.columns)
pickle_out=open("AllPlantings_Corrected_Phenotypes_colnames.pickle","wb")
pickle.dump(AllPlantings_Corrected_Phenotypes_colnames, pickle_out)
pickle_out.close()



#Saving accession IDs
ro.r('load("AllPlantings_Corrected_IDs.RData")') #Load the file


with localconverter(default_converter + pandas2ri.converter) as cv:
    pd_AllPlantings_Corrected_IDs = r('AllPlantings_Corrected_IDs')


np_AllPlantings_Corrected_IDs=np.array(pd_AllPlantings_Corrected_IDs)
np.save('AllPlantings_Corrected_IDs.npy',np_AllPlantings_Corrected_IDs)

