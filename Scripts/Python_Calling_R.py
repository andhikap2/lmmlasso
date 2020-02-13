#Calling R in python?

#Setting the working directory

import os
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")


#Following the tutorial.. (slightly modified)

from rpy2.robjects.packages import importr
import rpy2.robjects as ro

#Load dataset

ro.r('load("AllPlantings_Corrected.RData")') #Load the file
ro.r('library(glmnet)')
ro.r('library(caret)')

ro.r('rownames(AllPlantings_Corrected)=NULL') 
ro.r('lambdas_to_try <- seq(0, 1, length.out = 100)')
ro.r('lasso_cv <- cv.glmnet(x=as.matrix(AllPlantings_Corrected[,-c(1:14)]),y=as.matrix(AllPlantings_Corrected$DTB), alpha = 0.001,lambda=lambdas_to_try, nfolds=10)')
ro.r('lasso_lambda<- lasso_cv$lambda.min')
print(ro.r(lasso_lambda))




pydf = pd.DataFrame(Data)