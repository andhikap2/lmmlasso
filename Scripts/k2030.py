'''
Making a k2030 to figure out how adding a new line changes the scaling of the original k2029 matrix


'''

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
import dask.array as da
import dask.dataframe as dd


os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/PLINK")

(bim, fam, bed) = limix.io.plink.read("k2029", verbose=False)




'''
read a ped file
modify the ped file
save and merge


./plink --file k2029_pedfile --keep k2029_idfile.csv --make-bed --out test #keeps only a single row of data
./plink --bfile test --recode tab --out test 
./plink --bfile k2029 --merge test --make-bed --allow-no-sex --out k2030_test #to get the merge to work, I changed a few of the bases and the fid and iid

'''
(bim, fam, bed) = limix.io.plink.read("k2030_test", verbose=False)





























#Adding a new line to bed
newline=np.zeros((10709466,1),dtype='float64')
np.save('k2030_newline.bed',newline)
newline=pd.DataFrame(newline)
newline.to_csv('k2030_newline.bed',header=False)

template=fam.iloc[2028,:] #get the last line, replace the values, append this new line to make 2030
template[0]='999999'
template[1]='999999'
template[6]=2029
template.to_csv('k2030_newline.fam',header=False)


# newline_pd=pd.DataFrame(newline)
# newline_pd.to_csv('k2030_newline.txt')


'''
then in on plink run the following commands:

./plink --bfile k2029 --merge k2030_newline --make-bed --out k2030




# bed_df=dd.from_dask_array(bed_new) #too big

#Adding a new line to fam [a pandas DataFrame]
template=fam.iloc[2028,:] #get the last line, replace the values, append this new line to make 2030
template[0]='999999'
template[1]='999999'
template[6]=2029

newfam=fam.append(template,ignore_index=False)
