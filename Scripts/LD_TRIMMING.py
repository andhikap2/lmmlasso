#Bed files stuff


#Step 1: Get the list of SNP_ids that we want to keep (THIS IS IN PLINK)


sudo ./plink --bfile k2029 --geno 0.05 --maf 0.01 --indep 300 60 1.3 --out SNPs_whogem #WHOGEM parameters 222,969 SNPs
sudo ./plink --bfile k2029 --geno 0.05 -maf 0.05 --indep-pairwise 300 60 0.2 --out SNPs_0.2 # 82,231 SNPs
sudo ./plink --bfile k2029 --geno 0.05 -maf 0.05 --indep-pairwise 300 60 0.1 --out SNPs_0.1 # 29,171 SNPs
sudo ./plink --bfile k2029 --geno 0.05 --maf 0.05 --indep-pairwise 500 100 0.6 --out SNPs_0.6 # 341,615 SNPs
sudo ./plink --bfile k2029 --geno 0.05 --maf 0.01 --indep-pairwise 500 100 0.9 --out SNPs_many #Really loose criteria... 1,416,222 SNPs (too loose to run computations)
sudo ./plink --bfile k2029 --geno 0.05 --maf 0.01 --indep-pairwise 500 100 0.85 --out SNPs_085 #1,229,867 SNPs

#Import the PLINK file reader
from pandas_plink import example_file_prefix
import limix
import os
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/PLINK")

(bim, fam, bed) = limix.io.plink.read("k2029", verbose=False)

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")

#To compute the bed file:
#Example

#Loading the snp ids into python
import csv

with open('SNPs_1.2M_pruneIN.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    ids=list()
    for row in readCSV:
    	yeah=row[1]
    	ids.append(yeah)

with open('Ecotype_ids_UNIQUE.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    ecotypes=list()
    for row in readCSV:
    	yeah=row[1]
    	ecotypes.append(yeah)


os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/PLINK")

#Convert this list to a dataframe
import pandas as pd 
import numpy as np
pd_ids=pd.DataFrame(ids) #Ids of runed snps
pd_ecotypes=pd.DataFrame(ecotypes) #Unique ecotype IDs

logical=bim.loc[:,'snp'].isin(pd_ids.iloc[:,0]) #Create a logical to match snp ids
sum(logical) #Should equal the dimensions of pd_ids - 1

logical_eco=fam.loc[:,'fid'].isin(pd_ecotypes.iloc[:,0]) #Create a logical to match ecotype ids
sum(logical_eco) #Should equal 222

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

import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import pandas as pd
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
from rpy2.robjects import default_converter
from rpy2.robjects.conversion import localconverter

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

np.save("SNPs_1.2M.npy",snp_array)

from limix.stats import linear_kinship
K=linear_kinship(snp_array)




#Sanity check: Ecotype id 6909 (Col-0) should have all SNPs equal zero because it is the reference.
#None of the row sums equal the number of SNPs, but Col-0 should be reference.. There's probably something wrong with the indexing..

#Saving as a csv file
pd_beds.to_csv("SNPs_0.2_pruned.csv",sep='\t') #tab-delimited file
