setwd("~/Clim_GWAS/Clim_GWAS_2")
library(caret)
library(glmnet)
library(tidyverse)
library(broom)
load("AllPlantings_Corrected.RData")

for (i in unique(AllPlantings_Corrected$Planting)) {

	assign(paste0("",i),AllPlantings_Corrected[AllPlantings_Corrected$Planting==i,])
}



for ( i in unique(AllPlantings_Corrected$Planting)){

	lasso_cv <- cv.glmnet(x=data.matrix(eval(as.name(i)))[,-c(1:14)],y=data.matrix(eval(as.name(i))[,4]), alpha = 0.001, nfolds=10)
	opt_lambda<- lasso_cv$lambda.min #optimum lambda
	model_fit<- lasso_cv$glmnet.fit
	glm_betas <- tidy(model_fit) %>%
    filter(term != "(Intercept)", lambda == opt_lambda)  
    assign(paste0("Beta_",i),glm_betas)
}