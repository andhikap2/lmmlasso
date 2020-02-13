#!/usr/bin/env python3

#LMMLASSO PIPELINE
import sys
import limix.vardec.vardec as va 
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

tic=time.clock()
print("RUNNING LMMLASSO...")
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
plantings=np.load("Plantings_NEW.npy")

r2=str(sys.argv[1])
Th=str(sys.argv[2])
N=int(sys.argv[3])

if(N==720):
	b="Monthly"
elif (N==168):
	b="Weekly"
elif (N==24):
	b="Daily"
elif (N==1):
	b="Hourly"
else:
	print("RECONSIDER PERIOR CHOICE")




if (r2=="K"):
	import matplotlib
	matplotlib.use('Agg')

	os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")

	K=np.load('K_NEW.npy')

	os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")

	xfilename=["Microclimate_",b,"_threshold_",Th,".csv"]
	xfilename="".join(xfilename)
	env=pd.read_csv(xfilename,sep=',')
	X=np.array(env)

	assert K.shape[0] == K.shape[1], 'K MATRIX IS NOT SYMMETRICAL'
	assert K.shape[0] == env.shape[0], 'NO. OF INDIVIDUALS DOES NOT MATCH K MATRIX DIMENSIONS'

	os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")

	y=np.load("DTB_NEW.npy")

	assert y.shape[0] == X.shape[0], 'NO. OF OBSERVATIONS DOES NOT MATCH NO. OF INDIVIDUALS'
	
	#Running LMMLASSO
	alphas = 2.**(sp.linspace(-2,10,10)) #list of alphas to test
	n_splits=10
	N = X.shape[0]
	kf = KFold(n_splits,shuffle=True,random_state=None)
	n_alphas = len(alphas)
	MSE_train = sp.zeros((n_splits,n_alphas))
	MSE_test  = sp.zeros((n_splits,n_alphas))
	W_nonzero = sp.zeros((n_splits,n_alphas))
	kf.get_n_splits(X)

	os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")

	import lmmlasso

	lasso = lmmlasso.LmmLasso(warm_start=True,fit_intercept=False,tol=0.5) #note the tolerance value

	MSE_train,MSE_test,W_nonzero = lmmlasso.runCrossValidation(lasso,X,y,alphas,n_splits=10,K=K,verbose=True)
	MSE_train_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_train.mean(axis=0))).derivative(n=2) #something about the rotation here is different from the original script
	MSE_test_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_test.mean(axis=0))).derivative(n=2)

	alphas_inter = 2.**(sp.linspace(2,12,100))
	idx_train = sp.argmin(MSE_train_inter(alphas_inter)) 
	idx_test = sp.argmin(MSE_test_inter(alphas_inter))
	alpha_cv = (float(alphas_inter[idx_train])+float(alphas_inter[idx_test]))/2

	import pylab as pl
	pl.figure(figsize=[20,4])
	pls = pl.subplot(1,3,1)
	pls.plot(sp.log2(alphas),MSE_train.mean(axis=0),linewidth=2)
	pl.axvline(sp.log2(alpha_cv),color='r')
	pl.xlabel('log alpha')
	pl.ylabel('training error')
	pl.grid(True)

	pls = pl.subplot(1,3,2)
	pls.plot(sp.log2(alphas),MSE_test.mean(axis=0),linewidth=2)
	pl.axvline(sp.log2(alpha_cv),color='r')
	pl.xlabel('log alpha')
	pl.ylabel('test error')
	pl.grid(True)

	pls = pl.subplot(1,3,3)
	pls.plot(sp.log2(alphas),W_nonzero.mean(axis=0),linewidth=2)
	pl.axvline(sp.log2(alpha_cv),color='r')
	pl.xlabel('log alpha')
	pl.ylabel('number of nonzero coefficients')
	pl.grid(True)

	os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Figures")
	pl.savefig("Alphas_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+".png") #save the figure. Example file name: Alphas_SNPs_0.95_Microclimate_Monthly.png

	
	#Now we make the predictions
	os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
	kf12 = KFold(n_splits,shuffle=True,random_state=12) #12 because the random state is 12
	kf12.get_n_splits(X)

	for train_index,test_index in kf12.split(X):
		X_train, X_test, y_train, y_test, plantings_train, plantings_test= X[train_index], X[test_index], y[train_index], y[test_index], plantings[train_index], plantings[test_index] #Split into training and testing data sets
		K_train = K[train_index][:,train_index]
		K_test  = K[test_index][:,train_index]

	assert X_train.shape[0] == y_train.shape[0] == plantings_train.shape[0] == K_train.shape[0], 'FAILED TO PROPERLY SPLIT DATA'
	lasso.set_params(alpha=alpha_cv)
	lasso = lasso.fit(X_train,y_train,K=K_train)
	weights = lasso.coef_
	Y_hat = lasso.predict(X_test, K_test)
	pd_Yhat=pd.DataFrame(Y_hat)
	pd_Yhat.reset_index(drop=True,inplace=True)
	pd_ytest=pd.DataFrame(y_test)
	pd_plantings=pd.DataFrame(plantings_test)
	for_plotting=pd.concat([pd_plantings,pd_ytest,pd_Yhat],axis=1)
	for_plotting.columns=("Planting","Y_test","Y_hat")

	mse_test=mean_squared_error(y_test,Y_hat)
	rmse_test=np.sqrt(mse_test)
	yeah=str(rmse_test)
	rmse_string=["RMSE=",yeah]
	rmse_string="".join(rmse_string)

	from sklearn.metrics import r2_score

	r2_model=r2_score(y_test,Y_hat)
	yeahs=str(yeah)
	r2_string=["R-SQUARED = ",yeahs]
	r2_string="".join(r2_string)

	HalleFall2006=for_plotting.loc[for_plotting['Planting']=="HalleFall2006"]
	NorwichSummer2006=for_plotting.loc[for_plotting['Planting']=="NorwichSummer2006"]
	NorwichSummer2007=for_plotting.loc[for_plotting['Planting']=="NorwichSummer2007"]
	NorwichSpring2007=for_plotting.loc[for_plotting['Planting']=="NorwichSpring2007"]
	NorwichFall2006=for_plotting.loc[for_plotting['Planting']=="NorwichFall2006"]
	OuluFall2007=for_plotting.loc[for_plotting['Planting']=="OuluFall2007"]
	ValenciaFall2006=for_plotting.loc[for_plotting['Planting']=="ValenciaFall2006"]

	import matplotlib.pyplot as plt

	plt.figure()


	ranges=range(int(np.amax(y_test)))
	plt.plot(ranges,ranges)
	v=str(alpha_cv)

	plot_title=['LMMLASSO_SNPs_',r2,'_Microclimate_',b,'_alpha=',v]
	plot_title="".join(plot_title)
	plt.title(plot_title)
	plt.xlabel('DTB_observed')
	plt.ylabel('DTB_predicted')
	plt.plot(HalleFall2006.iloc[:,1],HalleFall2006.iloc[:,2],color='#8f0000',marker='o',linestyle='None',label='HalleFall2006')
	plt.plot(NorwichSummer2006.iloc[:,1],NorwichSummer2006.iloc[:,2],color='#c54259',marker='X',linestyle='None',label='NorwichSummer2006')
	plt.plot(NorwichSummer2007.iloc[:,1],NorwichSummer2007.iloc[:,2],color='#e881ab',marker='p',linestyle='None',label='NorwichSummer2007')
	plt.plot(NorwichSpring2007.iloc[:,1],NorwichSpring2007.iloc[:,2],color='#ffc2f3',marker='v',linestyle='None',label='NorwichSpring2007')
	plt.plot(NorwichFall2006.iloc[:,1],NorwichFall2006.iloc[:,2],color='#d890e7',marker='d',linestyle='None',label='NorwichFall2006')
	plt.plot(OuluFall2007.iloc[:,1],OuluFall2007.iloc[:,2],color='#9f66e1',marker='P',linestyle='None',label='OuluFall2007')
	plt.plot(ValenciaFall2006.iloc[:,1],ValenciaFall2006.iloc[:,2],color='#3a47de',marker='*',linestyle='None',label='ValenciaFall2006')
	plt.figtext(0,2,rmse_string)
	plt.figtext(0,1.5,r2_string)
	plt.legend(loc='lower right')
	
	os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Figures")
	plt.savefig("PredvsObs_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+".png")

	os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/LMMLASSO_RESULTS")
	np.save("W_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+".npy",weights) #save the coefficients

	alpha_cv_np=np.array(alpha_cv)
	c=np.array(sum(weights!=0))
	W_nonzero=np.array(c)

	

	


	results=np.array((alpha_cv_np,W_nonzero,rmse_test,r2_model))
	np.savetxt("Parameters_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+".csv",results)













else:
	import matplotlib
	matplotlib.use('Agg')

	os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")

	K=np.load('K_NEW.npy')

	os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")

	snpfilename=["SNPs_",r2,".npy"]
	snpfilename="".join(snpfilename)
	snps=np.load(snpfilename) #If x were given as just snps then they would be fixed effects right

	envfilename=["Microclimate_",b,"_threshold_",Th,".csv"]
	envfilename="".join(envfilename)
	env=pd.read_csv(envfilename,sep=',')
	environment=np.array(env)

	X=np.concatenate((snps,environment),axis=1)

	assert K.shape[0] == K.shape[1], 'K MATRIX IS NOT SYMMETRICAL'
	assert K.shape[0] == env.shape[0], 'NO. OF INDIVIDUALS DOES NOT MATCH K MATRIX DIMENSIONS'

	os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")

	y=np.load("DTB_NEW.npy")

	assert y.shape[0] == X.shape[0], 'NO. OF OBSERVATIONS DOES NOT MATCH NO. OF INDIVIDUALS'
	
	#Running LMMLASSO
	alphas = 2.**(sp.linspace(-2,10,10)) #list of alphas to test
	n_splits=10
	N = X.shape[0]
	kf = KFold(n_splits,shuffle=True,random_state=None)
	n_alphas = len(alphas)
	MSE_train = sp.zeros((n_splits,n_alphas))
	MSE_test  = sp.zeros((n_splits,n_alphas))
	W_nonzero = sp.zeros((n_splits,n_alphas))
	kf.get_n_splits(X)

	os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")

	import lmmlasso

	lasso = lmmlasso.LmmLasso(warm_start=True,fit_intercept=False,tol=0.5) #note the tolerance value

	MSE_train,MSE_test,W_nonzero = lmmlasso.runCrossValidation(lasso,X,y,alphas,n_splits=10,K=K,verbose=True)
	MSE_train_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_train.mean(axis=0))).derivative(n=2) #something about the rotation here is different from the original script
	MSE_test_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_test.mean(axis=0))).derivative(n=2)

	alphas_inter = 2.**(sp.linspace(2,12,100))
	idx_train = sp.argmin(MSE_train_inter(alphas_inter)) 
	idx_test = sp.argmin(MSE_test_inter(alphas_inter))
	alpha_cv = (float(alphas_inter[idx_train])+float(alphas_inter[idx_test]))/2

	import pylab as pl
	pl.figure(figsize=[20,4])
	pls = pl.subplot(1,3,1)
	pls.plot(sp.log2(alphas),MSE_train.mean(axis=0),linewidth=2)
	pl.axvline(sp.log2(alpha_cv),color='r')
	pl.xlabel('log alpha')
	pl.ylabel('training error')
	pl.grid(True)

	pls = pl.subplot(1,3,2)
	pls.plot(sp.log2(alphas),MSE_test.mean(axis=0),linewidth=2)
	pl.axvline(sp.log2(alpha_cv),color='r')
	pl.xlabel('log alpha')
	pl.ylabel('test error')
	pl.grid(True)

	pls = pl.subplot(1,3,3)
	pls.plot(sp.log2(alphas),W_nonzero.mean(axis=0),linewidth=2)
	pl.axvline(sp.log2(alpha_cv),color='r')
	pl.xlabel('log alpha')
	pl.ylabel('number of nonzero coefficients')
	pl.grid(True)

	os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Figures")
	pl.savefig("Alphas_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+".png") #save the figure. Example file name: Alphas_SNPs_0.95_Microclimate_Monthly.png

	
	#Now we make the predictions
	os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
	kf12 = KFold(n_splits,shuffle=True,random_state=12) #12 because the random state is 12
	kf12.get_n_splits(X)

	for train_index,test_index in kf12.split(X):
		X_train, X_test, y_train, y_test, plantings_train, plantings_test= X[train_index], X[test_index], y[train_index], y[test_index], plantings[train_index], plantings[test_index] #Split into training and testing data sets
		K_train = K[train_index][:,train_index]
		K_test  = K[test_index][:,train_index]

	assert X_train.shape[0] == y_train.shape[0] == plantings_train.shape[0] == K_train.shape[0], 'FAILED TO PROPERLY SPLIT DATA'
	lasso.set_params(alpha=alpha_cv)
	lasso = lasso.fit(X_train,y_train,K=K_train)
	weights = lasso.coef_
	Y_hat = lasso.predict(X_test, K_test)
	pd_Yhat=pd.DataFrame(Y_hat)
	pd_Yhat.reset_index(drop=True,inplace=True)
	pd_ytest=pd.DataFrame(y_test)
	pd_plantings=pd.DataFrame(plantings_test)
	for_plotting=pd.concat([pd_plantings,pd_ytest,pd_Yhat],axis=1)
	for_plotting.columns=("Planting","Y_test","Y_hat")

	mse_test=mean_squared_error(y_test,Y_hat)
	rmse_test=np.sqrt(mse_test)
	yeah=str(rmse_test)
	rmse_string=["RMSE=",yeah]
	rmse_string="".join(rmse_string)

	from sklearn.metrics import r2_score

	r2_model=r2_score(y_test,Y_hat)
	yeahs=str(yeah)
	r2_string=["R-SQUARED = ",yeahs]
	r2_string="".join(r2_string)

	HalleFall2006=for_plotting.loc[for_plotting['Planting']=="HalleFall2006"]
	NorwichSummer2006=for_plotting.loc[for_plotting['Planting']=="NorwichSummer2006"]
	NorwichSummer2007=for_plotting.loc[for_plotting['Planting']=="NorwichSummer2007"]
	NorwichSpring2007=for_plotting.loc[for_plotting['Planting']=="NorwichSpring2007"]
	NorwichFall2006=for_plotting.loc[for_plotting['Planting']=="NorwichFall2006"]
	OuluFall2007=for_plotting.loc[for_plotting['Planting']=="OuluFall2007"]
	ValenciaFall2006=for_plotting.loc[for_plotting['Planting']=="ValenciaFall2006"]

	import matplotlib.pyplot as plt

	plt.figure()


	ranges=range(int(np.amax(y_test)))
	plt.plot(ranges,ranges)
	v=str(alpha_cv)

	plot_title=['LMMLASSO_SNPs_',r2,'_Microclimate_',b,'_alpha=',v]
	plot_title="".join(plot_title)
	plt.title(plot_title)
	plt.xlabel('DTB_observed')
	plt.ylabel('DTB_predicted')
	plt.plot(HalleFall2006.iloc[:,1],HalleFall2006.iloc[:,2],color='#8f0000',marker='o',linestyle='None',label='HalleFall2006')
	plt.plot(NorwichSummer2006.iloc[:,1],NorwichSummer2006.iloc[:,2],color='#c54259',marker='X',linestyle='None',label='NorwichSummer2006')
	plt.plot(NorwichSummer2007.iloc[:,1],NorwichSummer2007.iloc[:,2],color='#e881ab',marker='p',linestyle='None',label='NorwichSummer2007')
	plt.plot(NorwichSpring2007.iloc[:,1],NorwichSpring2007.iloc[:,2],color='#ffc2f3',marker='v',linestyle='None',label='NorwichSpring2007')
	plt.plot(NorwichFall2006.iloc[:,1],NorwichFall2006.iloc[:,2],color='#d890e7',marker='d',linestyle='None',label='NorwichFall2006')
	plt.plot(OuluFall2007.iloc[:,1],OuluFall2007.iloc[:,2],color='#9f66e1',marker='P',linestyle='None',label='OuluFall2007')
	plt.plot(ValenciaFall2006.iloc[:,1],ValenciaFall2006.iloc[:,2],color='#3a47de',marker='*',linestyle='None',label='ValenciaFall2006')
	plt.figtext(0,2,rmse_string)
	plt.figtext(0,1.5,r2_string)
	plt.legend(loc='lower right')
	
	os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Figures")
	plt.savefig("PredvsObs_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+".png")

	os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/LMMLASSO_RESULTS")
	np.save("W_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+".npy",weights) #save the coefficients

	alpha_cv_np=np.array(alpha_cv)
	c=np.array(sum(weights!=0))
	W_nonzero=np.array(c)

	

	


	results=np.array((alpha_cv_np,W_nonzero,rmse_test,r2_model))
	np.savetxt("Parameters_SNPs_"+r2+"_Microclimate_"+b+"_Temp_"+Th+".csv",results)


	









toc=time.clock()
timediff=toc-tic
print("LMMLASSO COMPLETED IN",timediff,"SECONDS")
