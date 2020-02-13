#Conforming the Korves Plantings
import os
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
import pandas as pd 
korves=pd.read_csv('korves_rhodeisland_fixed.csv',sep='\\')

import numpy as np 
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/PLINK")


'''get ecotype id indices'''

import limix
(bim, fam, bed) = limix.io.plink.read("k2029", verbose=False)
logical_eco=fam.loc[:,'fid'].astype('int64').isin(korves.iloc[:,0]) #need the astype because fam fid was an Object!
logice=np.array(logical_eco)
indices_e = np.flatnonzero(logice) #ecotype id indices


'''get snp indices'''

import csv

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")
filename="SNPs_0.1.prune.in"
filename=str(filename)

with open(filename) as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    ids=[]
    for row in readCSV:
    	yeah=row
    	ids.append(yeah)

pd_ids=pd.DataFrame(ids)
logical=bim.loc[:,'snp'].isin(pd_ids.iloc[:,0]) 
logic=np.array(logical)
indices = np.flatnonzero(logic) #SNP id indices


b=bed[indices,:] #Slice by SNPs
b=b[:,indices_e] #Slice by ecotype id
b=b.T

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/PLINK")


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
pd_beds.to_csv('Korves_SNPs.csv')

arrays=np.array(pd_beds.iloc[:,1:],dtype=int)
np.save('Korves_SNPs.npy',arrays)