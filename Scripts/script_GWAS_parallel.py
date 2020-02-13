#!/usr/bin/env python

from datetime import *
import sys
import os
import re
import h5py
import numpy as np
import scipy as sp
import limix.qtl
from scipy import stats
if os.path.isdir("mixmogam"):
    sys.path.append("mixmogam")
elif os.path.isdir("/data/GWAS/At/mixmogam"):
    sys.path.append("/data/GWAS/At/mixmogam")
else: sys.exit("Can't find a path to mixmogam libraries") 
import linear_models as lm

if not os.path.exists("1001genomes_accessions.npy"):
    prefix="/data/GWAS/At/"
else: prefix="" #sys.exit("Can't find the genotype file(s)")

START=datetime.now()

#Missing Value Threshold
Miss_Tol=.3
#Minor Allele Frequency Threshold
MAF_Tol=.05
#Which Models should be computed?

#sys.argv=['','DTB_HalleFall2006.csv','1001genomes_SNPdata_part01.npy','LASSO']

PHENO_file=str(sys.argv[1]).replace('.csv','')
if re.search('/',PHENO_file):
	PHENO_file=PHENO_file.rsplit('/', 1)[1].split(".")[0]

if "MT-LMM" not in sys.argv and "LASSO" not in sys.argv:

	if "GLM" not in sys.argv:
		GLM=False
	else: GLM=True

	if "EMMAX" in sys.argv:
		EMMAX=True
	else: EMMAX=False

	if "KW" in sys.argv:
		KW=True
	else: KW=False

	if "LMM" in sys.argv:
		LMM=True
	else: LMM=False

	#get the SNP data
	SNP_data=np.load(sys.argv[2])
	SNP_names=np.load(sys.argv[2].replace('SNPdata','SNPnames'))
	n_snps_start=SNP_data.shape[0]
	Lines_geno=np.load(prefix+"1001genomes_accessions.npy")

	#get the phenotype data
	Pheno_data = np.genfromtxt(str(sys.argv[1]), delimiter=',')
	Pheno_handle=open(sys.argv[1],'r')
	Pheno_data = Pheno_handle.readlines()
	Pheno_data = np.array([x.replace('\n','').split(',') for x in Pheno_data])
	Lines_pheno=Pheno_data[1:,0]
	Pheno_data=Pheno_data[1:,1].astype("float")

	#remove lines with missing phenotype
	Pheno_data=Pheno_data[~np.isnan(Pheno_data)]
	Lines_pheno=Lines_pheno[~np.isnan(Pheno_data)]
	#sort the lines in Phenotype
	Pheno_sort=np.argsort(Lines_pheno)
	Pheno_data=Pheno_data[Pheno_sort]
	#remove the lines with no genotypes
	Lines_in_pheno=np.in1d(Lines_pheno,Lines_geno)
	Pheno_data=Pheno_data[Lines_in_pheno]
	Lines_pheno=Lines_pheno[Lines_in_pheno]

	#remove lines with missing phenotype
	Lines_in_geno=np.in1d(Lines_geno,Lines_pheno)
	SNP_data=SNP_data[:,Lines_in_geno]
	print('=============================================================')
	print('FROM PHENOTYPE FILE # '+str(sys.argv[1]))
	print('Line/Ecotype ID included:')
	print(Lines_pheno)

	#sort the lines in SNP
	SNP_sort=np.argsort(Lines_geno[Lines_in_geno])
	SNP_data=np.transpose(SNP_data[:,SNP_sort])
	SNP_data=SNP_data.astype('float')

	#include the kinship matrix
	with h5py.File(prefix+'kinship_ibs_mac5.hdf5','r') as hf:
		K_data = hf.get('kinship')
		K_data = np.array(K_data)
		Lines_K = hf.get('accessions')
		Lines_K = np.array(Lines_K)

	Lines_in_K=np.in1d(Lines_K,Lines_pheno)
	K_data=K_data[Lines_in_K,:][:,Lines_in_K]
	K_sort=np.argsort(Lines_K[Lines_in_K])
	K_data=K_data[K_sort,:][:,K_sort]
	K_data=K_data/np.min(K_data)

	#estimating the allele frequencies in the data
	SNPsum=np.nansum(SNP_data,axis=0)
	nInd=np.sum(~np.isnan(SNP_data),axis=0)
	freq_hat=np.array(SNPsum,dtype="float")/(2*nInd)
	mask=np.ndarray.flatten(np.array(np.all([freq_hat>MAF_Tol,freq_hat<(1-MAF_Tol)],axis=0)).astype("bool"))
	SNP_data=SNP_data[:,mask]
	SNP_names=SNP_names[mask,:]
	MAF=freq_hat[mask]
	SNP_in=np.sum(mask.astype("int"))

	print('=============================================================')
	print('FROM GENOTYPE FILE # '+str(sys.argv[2]))
	if SNP_in==0:
		sys.exit('ERROR: No SNP satisfied Missing Value and/or Minor Allele Frequency Threshold, please edit')
	else:
		print(str(n_snps_start-SNP_in)+'/'+str(n_snps_start)+' SNPs did not meet the MAF threshold')
		print(str(SNP_data.shape[1])+'/'+str(SNP_in)+' SNPs included in this GWAS')
		END_load=datetime.now()
		print('Files were loaded in '+str(END_load-START))

	if GLM:
		START=datetime.now()
		test=[]
		print('_____________________________________________________________')
		print('Starting the computation of GLM')
		for i in range(SNP_data.shape[1]):
			test_data=np.transpose((SNP_data[:,i],Pheno_data))
			test_data=test_data[~np.isnan(test_data).any(axis=1)]
			try:
				slope, intercept, r_value, p_value, std_err = stats.linregress(test_data[:,0],test_data[:,1])
			except ValueError:
				slope, p_value = float('nan'), float('nan')
			test.append([SNP_names[i,0],SNP_names[i,1],slope,p_value,MAF[i]])
			sys.stdout.write('\r')
			sys.stdout.write("[%-20s] %d%%" % ('='*(100*i/SNP_in), 100*i/SNP_in))
			sys.stdout.flush()

		test=np.array(test)
		header='Chromosome,Position,beta,p-val,MAF'
		filename="GWAS_GLM_for_"+PHENO_file+"_"+str(sys.argv[2]).split("_")[2].replace(".npy","")+".csv"
		np.savetxt(filename, test, delimiter=",",header=header,fmt="%s")
		END_GLM=datetime.now()
		print('\rGLM completed in '+str(END_GLM-START))

	if KW:
		START=datetime.now()
		test_H=[]
		test_p=[]
		print('_____________________________________________________________')
		print('Starting the computation of KW')
		for i in range(SNP_data.shape[1]):
			H, p_value = stats.kruskal(Pheno_data[SNP_data[:,i]==0],Pheno_data[SNP_data[:,i]==1])
			test_H.append(H)
			test_p.append(p_value)
			sys.stdout.write('\r')
			sys.stdout.write("[%-20s] %d%%" % ('='*(100*i/SNP_in), 100*i/SNP_in))
			sys.stdout.flush()

		test=np.vstack((SNP_names[:,0].astype("int"),SNP_names[:,1].astype("int"),np.array(test_H),np.array(test_p),MAF))
		test=np.transpose(test)
		header='Chromosome,Position,H,p-val,MAF'
		filename="GWAS_KW_for_"+PHENO_file+"_"+str(sys.argv[2]).split("_")[2].replace(".npy","")+".csv"
		np.savetxt(filename, test, delimiter=",",header=header,fmt="%s")
		END_KW=datetime.now()
		print('\rKW completed in '+str(END_KW-START))


	if EMMAX:
		START=datetime.now()
		print('_____________________________________________________________')
		print('Starting the computation of EMMAX using mixmogam')
		mm_results = lm.emmax(np.transpose(SNP_data), Pheno_data, K_data,with_betas=True)
		betas=mm_results['betas']
		betas=[[0,0] if i is None else i for i in betas]
		betas=[[0,0] if len(i)!=2 else i for i in betas]
		betas=np.hstack(betas).reshape(len(betas),2)
		test=np.vstack((SNP_names[:,0].astype("int"),SNP_names[:,1].astype("int"),betas[:,1],mm_results['ps'],MAF))
		test=np.transpose(test)
		header='Chromosome,Position,beta,p-val,MAF'
		filename="GWAS_EMMAX_for_"+PHENO_file+"_"+str(sys.argv[2]).split("_")[2].replace(".npy","")+".csv"
		np.savetxt(filename, test, delimiter=",",header=header,fmt=["%i","%i","%s","%s","%s"])
		END_EMMAX=datetime.now()
		print('\rEMMAX completed in '+str(END_EMMAX-START))

	if LMM:
		START=datetime.now()
		print('_____________________________________________________________')
		print('Starting the computation of LMM using limix')
		mm_results = limix.qtl.qtl_test_lmm(SNP_data, Pheno_data, K=K_data)
		p_val=mm_results.getPv()
		betas=mm_results.getBetaSNP()
		test=np.vstack((SNP_names[:,0].astype("int"),SNP_names[:,1].astype("int"),betas,p_val,MAF))
		test=np.transpose(test)
		header='Chromosome,Position,beta,p-val,MAF'
		filename="GWAS_LMM_for_"+PHENO_file+"_"+str(sys.argv[2]).split("_")[2].replace(".npy","")+".csv"
		np.savetxt(filename, test, delimiter=",",header=header,fmt=["%i","%i","%s","%s","%s"])
		END_LMM=datetime.now()
		print('\LMM completed in '+str(END_LMM-START))

elif "MT-LMM" in sys.argv:

	import pandas as pd

	SNP_data=np.load(sys.argv[2])
	SNP_names=np.load(sys.argv[2].replace('SNPdata','SNPnames'))
	n_snps_start=SNP_data.shape[0]
	Lines_geno=np.load("/data/GWAS/At/1001genomes_accessions.npy")

	Data=pd.read_csv(str(sys.argv[1]),header=0)
	Data=Data.values

	#Extract the Phenotype and convert it to numerical
	Pheno_data=Data[:,1:]
	Lines_pheno=Data[:,0].astype('int').astype('str')
	#remove lines with missing phenotype
	Pheno_data=Pheno_data[~np.all(pd.isna(Pheno_data),axis=1),:]
	Lines_pheno=Lines_pheno[~np.all(pd.isna(Pheno_data),axis=1)]
	#sort the lines in Phenotype
	Pheno_sort=np.argsort(Lines_pheno)
	Pheno_data=Pheno_data[Pheno_sort,:]
	#remove the lines with no genotypes
	Lines_in_pheno=np.in1d(Lines_pheno,Lines_geno)
	Pheno_data=Pheno_data[Lines_in_pheno]
	Lines_pheno=Lines_pheno[Lines_in_pheno]

	#remove lines with missing phenotype
	Lines_in_geno=np.in1d(Lines_geno,Lines_pheno)
	SNP_data=SNP_data[:,Lines_in_geno]
	print('=============================================================')
	print('FROM PHENOTYPE FILE # '+str(sys.argv[1]))
	print(str(Lines_pheno.shape[0])+' Line/Ecotype included with ID:')
	print(Lines_pheno)

	#sort the lines in SNP
	SNP_sort=np.argsort(Lines_geno[Lines_in_geno])
	SNP_data=np.transpose(SNP_data[:,SNP_sort])
	SNP_data=SNP_data.astype('float')

	#include the kinship matrix
	with h5py.File('/data/GWAS/At/kinship_ibs_mac5.hdf5','r') as hf:
		K_data = hf.get('kinship')
		K_data = np.array(K_data)
		Lines_K = hf.get('accessions')
		Lines_K = np.array(Lines_K)

	Lines_in_K=np.in1d(Lines_K,Lines_pheno)
	K_data=K_data[Lines_in_K,:][:,Lines_in_K]
	K_sort=np.argsort(Lines_K[Lines_in_K])
	K_data=K_data[K_sort,:][:,K_sort]
	K_data=K_data/np.min(K_data)

	#estimating the allele frequencies in the data
	SNPsum=np.nansum(SNP_data,axis=0)
	nInd=np.sum(~np.isnan(SNP_data),axis=0)
	freq_hat=np.array(SNPsum,dtype="float")/(nInd)
	mask=np.ndarray.flatten(np.array(np.all([freq_hat>MAF_Tol,freq_hat<(1-MAF_Tol)],axis=0)).astype("bool"))
	SNP_data=SNP_data[:,mask]
	SNP_names=SNP_names[mask,:]
	MAF=freq_hat[mask]
	SNP_in=np.sum(mask.astype("int"))

	print('=============================================================')
	print('FROM GENOTYPE FILE # '+str(sys.argv[2]))
	if SNP_in==0:
		sys.exit('ERROR: No SNP satisfied Missing Value and/or Minor Allele Frequency Threshold, please edit')
	else:
		print(str(n_snps_start-SNP_in)+'/'+str(n_snps_start)+' SNPs did not meet the MAF threshold')
		print(str(SNP_data.shape[1])+'/'+str(SNP_in)+' SNPs included in this GWAS')
		END_load=datetime.now()
		print('Files were loaded in '+str(END_load-START))


	N, P = Pheno_data.shape

	covs = sp.ones((N, 1))      #covariates
	Acovs = sp.eye(P)           #the design matrix for the covariates   
	#Asnps = sp.eye(P)           #the design matrix for the SNPs having an effect on any trait
	# or
	Asnps = sp.ones((1,P))           #the design matrix for the SNPs having an effect on every trait

	from limix.qtl import qtl_test_lmm_kronecker

	lmm, p_val = qtl_test_lmm_kronecker(snps=SNP_data,
	                                      phenos=Pheno_data,
	                                      covs=covs,
	                                      Acovs=Acovs,
	                                      Asnps=Asnps,
	                                      K1r=K_data)

	betas=np.repeat(-999,p_val.shape[1])
	test=np.vstack((SNP_names[:,0].astype("int"),SNP_names[:,1].astype("int"),betas,p_val,MAF))
	test=np.transpose(test)
	header='Chromosome,Position,beta,p-val,MAF'
	filename="GWAS_MT-LMM_for_"+PHENO_file+"_"+str(sys.argv[2]).split("_")[2].replace(".npy","")+".csv"
	np.savetxt(filename, test, delimiter=",",header=header,fmt=["%i","%i","%s","%s","%s"])
	END_LMM=datetime.now()
	print('MT-LMM completed in '+str(END_LMM-START))

elif "LASSO" in sys.argv:

	import lmmlasso
	import pylab as pl

		#get the SNP data
	with h5py.File('/data/GWAS/At/imputed_snps_binary.hdf5','r') as hf:
		SNP_data = hf.get('snps')
		SNP_data = np.array(SNP_data)
		Lines_geno = hf.get('accessions')
		Lines_geno = np.array(Lines_geno)
		positions = hf.get('positions')
		chr_regions = positions.attrs['chr_regions']
		Chr=np.array([])
		k=1
		for i in chr_regions[:,1]-chr_regions[:,0]:
			Chr=np.concatenate((Chr,np.repeat(k,i)))
			k+=1
		SNP_names = np.transpose(np.vstack((Chr,positions))).astype("int")
	
	n_snps_start=SNP_data.shape[0]
	hf=[]

	#get the phenotype data
	Pheno_data = np.genfromtxt(str(sys.argv[1]), delimiter=',')
	Pheno_handle=open(sys.argv[1],'r')
	Pheno_data = Pheno_handle.readlines()
	Pheno_data = np.array([x.replace('\n','').split(',') for x in Pheno_data])
	Lines_pheno=Pheno_data[1:,0]
	Pheno_data=Pheno_data[1:,1].astype("float")

	#remove lines with missing phenotype
	Pheno_data=Pheno_data[~np.isnan(Pheno_data)]
	Lines_pheno=Lines_pheno[~np.isnan(Pheno_data)]
	#sort the lines in Phenotype
	Pheno_sort=np.argsort(Lines_pheno)
	Pheno_data=Pheno_data[Pheno_sort]
	#remove the lines with no genotypes
	Lines_in_pheno=np.in1d(Lines_pheno,Lines_geno)
	Pheno_data=Pheno_data[Lines_in_pheno]
	Lines_pheno=Lines_pheno[Lines_in_pheno]

	#remove lines with missing phenotype
	Lines_in_geno=np.in1d(Lines_geno,Lines_pheno)
	SNP_data=SNP_data[:,Lines_in_geno]
	print('=============================================================')
	print('FROM PHENOTYPE FILE # '+str(sys.argv[1]))
	print('Line/Ecotype ID included:')
	print(Lines_pheno)

	#sort the lines in SNP
	SNP_sort=np.argsort(Lines_geno[Lines_in_geno])
	SNP_data=np.transpose(SNP_data[:,SNP_sort])
	SNP_data=SNP_data.astype('float')

	#include the kinship matrix
	with h5py.File('/data/GWAS/At/kinship_ibs_mac5.hdf5','r') as hf:
		K_data = hf.get('kinship')
		K_data = np.array(K_data)
		Lines_K = hf.get('accessions')
		Lines_K = np.array(Lines_K)

	Lines_in_K=np.in1d(Lines_K,Lines_pheno)
	K_data=K_data[Lines_in_K,:][:,Lines_in_K]
	K_sort=np.argsort(Lines_K[Lines_in_K])
	K_data=K_data[K_sort,:][:,K_sort]
	K_data=K_data/np.min(K_data)

	#estimating the allele frequencies in the data
	SNPsum=np.nansum(SNP_data,axis=0)
	nInd=np.sum(~np.isnan(SNP_data),axis=0)
	freq_hat=np.array(SNPsum,dtype="float")/(2*nInd)
	mask=np.ndarray.flatten(np.array(np.all([freq_hat>MAF_Tol,freq_hat<(1-MAF_Tol)],axis=0)).astype("bool"))
	SNP_data=SNP_data[:,mask]
	SNP_names=SNP_names[mask,:]
	MAF=freq_hat[mask]
	SNP_in=np.sum(mask.astype("int"))

	print('=============================================================')
	print('FROM GENOTYPE FILE # imputed_snps_binary.hdf5')
	if SNP_in==0:
		sys.exit('ERROR: No SNP satisfied Missing Value and/or Minor Allele Frequency Threshold, please edit')
	else:
		print(str(n_snps_start-SNP_in)+'/'+str(n_snps_start)+' SNPs did not meet the MAF threshold')
		print(str(SNP_data.shape[1])+'/'+str(SNP_in)+' SNPs included in this GWAS')
		END_load=datetime.now()
		print('Files were loaded in '+str(END_load-START))


	# running cross-validation
	alphas = 2.**(sp.linspace(2,12,10))
	alphas = alphas[::-1]

	# running LMM-Lasso
	lasso = lmmlasso.LmmLasso(warm_start=True,fit_intercept=False)
	MSE_train,MSE_test,W_nonzero = lmmlasso.runCrossValidation(lasso,SNP_data,Pheno_data,alphas,n_folds=10,K=K_data,verbose=True)

	# the dirty version
	#idx = sp.argmin(MSE_test.mean(axis=0))
	#alpha_cv = alphas[idx]
	#OR
	# the verion with secondary derivatives
	MSE_train_inter=sp.interpolate.UnivariateSpline(x=np.flip(alphas,axis=0), y=np.flip(MSE_train.mean(axis=0),axis=0)).derivative(n=2)
	MSE_test_inter=sp.interpolate.UnivariateSpline(x=np.flip(alphas,axis=0), y=np.flip(MSE_test.mean(axis=0),axis=0)).derivative(n=2)
	alphas_inter = 2.**(sp.linspace(2,12,100))
	idx_train = sp.argmax(MSE_train_inter(alphas_inter))
	idx_test = sp.argmin(MSE_test_inter(alphas_inter))
	alpha_cv = (float(alphas_inter[idx_train])+float(alphas_inter[idx_test]))/2

	pl.figure(figsize=[20,4])
	plt = pl.subplot(1,3,1)
	plt.plot(sp.log2(alphas),MSE_train.mean(axis=0),linewidth=2)
	pl.axvline(sp.log2(alpha_cv),color='r')
	pl.xlabel('log alpha')
	pl.ylabel('training error')
	pl.grid(True)

	plt = pl.subplot(1,3,2)
	plt.plot(sp.log2(alphas),MSE_test.mean(axis=0),linewidth=2)
	pl.axvline(sp.log2(alpha_cv),color='r')
	pl.xlabel('log alpha')
	pl.ylabel('test error')
	pl.grid(True)

	plt = pl.subplot(1,3,3)
	plt.plot(sp.log2(alphas),W_nonzero.mean(axis=0),linewidth=2)
	pl.axvline(sp.log2(alpha_cv),color='r')
	pl.xlabel('log alpha')
	pl.ylabel('number of nonzero coefficients')
	pl.grid(True)

	pl.savefig("Alphas_"+PHENO_file+".png")

	lasso.set_params(alpha=alpha_cv)
	lasso = lasso.fit(SNP_data,Pheno_data,K=K_data)
	weights = lasso.coef_

	import matplotlib.pyplot as plt
	# cross_val_predict returns an array of the same size as `y` where each entry
	# is a prediction obtained by cross validation:
	Y_hat = lasso.predict(SNP_data, K_data)

	fig, ax = plt.subplots()
	ax.scatter(Pheno_data, Y_hat, edgecolors=(0, 0, 0))
	ax.plot([Pheno_data.min(), Pheno_data.max()], [Y_hat.min(), Y_hat.max()], 'k--', lw=4)
	ax.set_xlabel('Measured')
	ax.set_ylabel('Predicted')
	plt.savefig("PredvsObs_"+PHENO_file+".png")

	p_val=np.repeat(-999,weights.shape[0])
	test=np.vstack((SNP_names[:,0].astype("int"),SNP_names[:,1].astype("int"),weights,p_val,MAF))
	test=np.transpose(test)
	header='Chromosome,Position,beta,p-val,MAF'
	filename="GWAS_LASSO_for_"+PHENO_file+"_out.csv"
	np.savetxt(filename, test, delimiter=",",header=header,fmt=["%i","%i","%s","%s","%s"])
	END_LMM=datetime.now()
	print('LMM-LASSO completed in '+str(END_LMM-START))

