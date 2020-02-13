#Genotype data k2029 whole matrix

from pandas_plink import example_file_prefix
import pandas as pd
import limix
import os
import numpy as np
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/PLINK")
(bim, fam, bed) = limix.io.plink.read("k2029", verbose=False)

accession_ids=fam.iloc[:,0]
accessions=np.array(accession_ids)

np.save('k2029_accessions.npy',accessions)

#10709466 / 12 = 892455.5
''' 

Let's make it 850,000 SNPs per file


'''
for i in range(12):
	print((i*892455),((i*892455)+892454))
	

part01=bed[0:850000,:]
part02=bed[850000:1700000,:]
part03=bed[1700000:2550000,:]
part04=bed[2550000:3400000,:]
part05=bed[3400000:4250000,:]
part06=bed[4250000:5100000,:]
part07=bed[5100000:5950000,:]
part08=bed[5950000:6800000,:]
part09=bed[6800000:7650000,:]
part10=bed[7650000:8500000,:]
part11=bed[8500000:9350000,:]
part12=bed[9350000:10200000,:]

part13=bed[10200000:,:]

part01.shape==part02.shape==part03.shape==part04.shape==part05.shape==part06.shape==part07.shape==part08.shape==part09.shape==part10.shape==part11.shape==part12.shape #should be true

bed01=part01.compute()
bed02=part02.compute()
bed03=part03.compute()
bed04=part04.compute()
bed05=part05.compute()
bed06=part06.compute()
bed07=part07.compute()
bed08=part08.compute()
bed09=part09.compute()
bed10=part10.compute()
bed11=part11.compute()
bed12=part12.compute()

bed=part01.compute()
pdb=pd.DataFrame(bed)
pdb=pdb.replace(to_replace=2.0,value=1)
bed=np.array(pdb)
bed.shape
bed
np.save('k2029_SNPdata_part01.npy',bed)


snp01=bim.iloc[0:850000,1]
snp02=bim.iloc[850000:1700000,1]
snp03=bim.iloc[1700000:2550000,1]
snp04=bim.iloc[2550000:3400000,1]
snp05=bim.iloc[3400000:4250000,1]
snp06=bim.iloc[4250000:5100000,1]
snp07=bim.iloc[5100000:5950000,1]
snp08=bim.iloc[5950000:6800000,1]
snp09=bim.iloc[6800000:7650000,1]
snp10=bim.iloc[7650000:8500000,1]
snp11=bim.iloc[8500000:9350000,1]
snp12=bim.iloc[9350000:10200000,1]
snp13=bim.iloc[10200000:,1]

snps01=np.array(snp01)
snps01
np.save('k2029_SNPnames_part01.npy',snps01)


