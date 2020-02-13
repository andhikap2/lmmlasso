#!/usr/bin/env python3

#SNP PROCESSING
import sys
from pandas_plink import example_file_prefix
import limix
import os
import numpy as np
import csv
import pandas as pd
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
from rpy2.robjects import default_converter
from rpy2.robjects.conversion import localconverter
import time

tic=time.clock()
print("PROCESSING SNPS...")

'''arguments'''
r2=str(sys.argv[1]) #The argument

''''''
if r2 =="K":
    sys.exit()

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/PLINK")

(bim, fam, bed) = limix.io.plink.read("k2029", verbose=False)

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")

#Loading the snp ids into python

filename=["SNPs_",r2,".prune.in"]
filename="".join(filename)
filename=str(filename)

with open(filename) as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    ids=[]
    for row in readCSV:
    	yeah=row
    	ids.append(yeah)

pd_ids=pd.DataFrame(ids)

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")

with open('Ecotype_ids_UNIQUE.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    ecotypes=list()
    for row in readCSV:
    	yeah=row[1]
    	ecotypes.append(yeah)


os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/PLINK")

#Convert this list to a dataframe
pd_ecotypes=pd.DataFrame(ecotypes) #Unique ecotype IDs

logical=bim.loc[:,'snp'].isin(pd_ids.iloc[:,0]) 

if (sum(logical)!=(pd_ids.shape[0])):
    print("Incorrect SNP IDs")

logical_eco=fam.loc[:,'fid'].isin(pd_ecotypes.iloc[:,0]) 


if (sum(logical_eco)!=222):
    print("Incorrect Ecotype IDs")


logic=np.array(logical)
logice=np.array(logical_eco)

#Subset the bed file
indices = np.flatnonzero(logic) #SNP id indices
indices_e = np.flatnonzero(logice) #Ecotype indices

b=bed[indices,:] #Slice by SNPs
b=b[:,indices_e] #Slice by ecotype id
b=b.T


beds=b.compute()
pd_beds=pd.DataFrame(beds)


pd_logice=pd.DataFrame(logice)
pd_logic=pd.DataFrame(logic)


#Need to reset indices before concatenating
fam.reset_index(drop=True, inplace=True)
pd_logice.reset_index(drop=True, inplace=True)
ye=pd.concat((fam,pd_logice),axis=1,ignore_index=True)
ye.columns=['fid','iid','father','mother','gender','trait','i','logical']
ye.head()

yes=ye[ye['logical']==True]

yes_ecotypes=yes['fid'] #Get the ecotype ids (lengths don't match??)
pd_ecotypes=pd.DataFrame(yes_ecotypes)
pd_ecotypes.reset_index(drop=True,inplace=True)



#assign column names
pd_beds=pd.concat((pd_ecotypes,pd_beds),axis=1,ignore_index=True)


bim.reset_index(drop=True,inplace=True)
pd_logic.reset_index(drop=True,inplace=True)
ah=pd.concat((bim,pd_logic),axis=1,ignore_index=True)
ah.columns=['chrom','snp','cm','pos','a0','a1','i','logical']
ah.head()

ahs=ah[ah['logical']==True]
ahs.head()

ahs_ids=np.array(ahs['snp'])
ahs_ids=ahs_ids.astype(str)
dd=np.array('ecotype_id')
ccc=np.hstack((dd,ahs_ids))
pd_beds.columns=[ccc] #Name the columns


#replacing the 2's with one's
pd_beds=pd_beds.replace(to_replace=2.0,value=1.0)



os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")

#Conforming in Python


ro.r('load("Multitrait_NEW.RData")') #Load the file

with localconverter(default_converter + pandas2ri.converter) as cv:
    phenotypes=r('Multitrait_NEW')


snp_array=np.empty((phenotypes.shape[0],(pd_beds.shape[1]-1)))

for i in range(snp_array.shape[0]):
    id=phenotypes.iloc[i,1]
    boolean=pd_beds['ecotype_id']==id
    boo=np.array(boolean)
    index=np.nonzero(boo)
    snp_row=pd_beds.iloc[index[0],1:]
    snp_array[i]=snp_row

print(snp_array[0:5,0:5])


os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")






savefilename=["SNPs_",r2,".npy"]
savefilename="".join(savefilename)
np.save(savefilename,snp_array)



toc=time.clock()
timediff=toc-tic
print("SNP PROCESSING COMPLETED IN",timediff,"SECONDS")

