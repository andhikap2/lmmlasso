#iscan with covariates?
from limix.qtl import iscan 

import sys
#import limix.vardec.vardec as va 
import scipy.linalg as LA
from sklearn.linear_model import Lasso
import numpy as np 
from sklearn.metrics import mean_squared_error
import math
from sklearn.model_selection import KFold 
import gc
import scipy as sp
import os
import csv
import time
import pandas as pd
import scipy.stats

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
plantings=np.load("Plantings_NEW.npy")

import matplotlib
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
#K=np.load('K_NEWER.npy') #Load the preferred K-matrix
K=np.load("K_MATRIX_EMMAXIBS.npy")

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")
snps=np.load('SNPs_0.1.npy') #If x were given as just snps then they would be fixed effects right
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/New_Microclim")
env=pd.read_csv('Microclimate_minmaxDaily_threshold_0_0.csv',sep=',')
logical=env.columns.str.startswith(('PAR','TT','PTT','daylength') )
env=env.iloc[:,~logical] #get columns that don't start with PAR
environment=np.array(env)
assert K.shape[0] == K.shape[1], 'K MATRIX IS NOT SYMMETRICAL'
assert K.shape[0] == env.shape[0], 'NO. OF INDIVIDUALS DOES NOT MATCH K MATRIX DIMENSIONS'
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
y=np.load("DTB_NEW.npy")

iscan_model=iscan(snps,y,K=K,M=env)

#𝐲 ~ 𝓝(𝙼𝜶, 4777.187⋅𝙺 + 394.918⋅𝙸) --prop explained by K: 0.9236446282509733


h1_effsizes=iscan_model.effsizes['h1'] #test every combination between snp and covariate (dim= 86760 * 406)
h1_effsizes.to_csv('h1_effsizes_minmaxEMMAXIBS.csv')
h2_effsizes=iscan_model.effsizes['h2']
h2_effsizes.to_csv('h2_effsizes_minmaxEMMAXIBS.csv')
stats=iscan_model.stats
stats.to_csv('stats_minmaxEMMAXIBS.csv')