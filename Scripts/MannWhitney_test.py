#Mann-Whitney Rank Test

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

with localconverter(default_converter + pandas2ri.converter) as cv:
	pd_snps=r('AllPlantings_Corrected_SNPs')

snps=pd_snps.iloc[:,1:len(pd_snps.columns)+1]

ro.r('load("Phenotypes_dtbfixed.RData")')

with localconverter(default_converter + pandas2ri.converter) as cv:
	pd_phenotypes=r('phenotypes')


DTB=pd_phenotypes.iloc[:,1]

X=np.array(snps)
y=np.array(DTB)
ylog=np.log(y) #make it more normally distributed


teststat=np.zeros(1)
pvalues=np.zeros(1)

from scipy.stats import mannwhitneyu

for i in range(len(snps.columns)):
	c=X[:,i]
	DTB_snp=np.vstack((ylog,c))
	DTB_snp=np.transpose(DTB_snp)
	is0=DTB_snp[:,1]==0
	is1=DTB_snp[:,1]==1
	a = DTB_snp[is0]
	dtb0 = a[:,0]
	b = DTB_snp[is1]
	dtb1 = b[:,0]
	ts = mannwhitneyu(dtb0,dtb1,alternative='two-sided')
	stat = ts.statistic
	pvalue = ts.pvalue
	teststat=np.append(teststat,stat)
	pvalues=np.append(pvalues,pvalue)
	print(i)

pvaluess=pvalues[1:(len(pvalues)+1)]
logp=-np.log10(pvaluess)

import matplotlib.pyplot as plt 

plt.plot(logp)
plt.xlabel('SNPs')
plt.ylabel('-log10 p values')
plt.title('Mann-Whitney Rank Test')
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
plt.show()