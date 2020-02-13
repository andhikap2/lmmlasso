'''

Function for making predictions across a raster given some environment and SNP data (but only uses the IBS matrix to make actual predictions??)

'''
import limix.vardec.vardec as va 
import scipy.linalg as LA
from sklearn.linear_model import Lasso
import numpy as np 
# from sklearn.metrics import mean_squared_error
import math
# from sklearn.model_selection import KFold 
# import gc
import scipy as sp
# import os
# import csv
# import time
import pandas as pd
import lmmlasso
import operator
import itertools
# from tqdm import tqdm


class PredMap(Lasso):
	"""
	prediction map object class
	"""
	def __init__(self, G, E):
		"""
		G = SNPs
		E = environmental data of daily min followed by daily max temp.

		"""
		self.G = G
		self.E = E
		self.msg = 'PredMap'


	def setKinship(self,K):
		"""
		manually set a k-matrix of choice
		"""
		assert K.shape[0] == K.shape[1], 'Kinship matrix is not square'
		assert K.shape[0] == self.G.shape[0], 'K-matrix dimensions not consistent with no. of observations'
		self.K=K
		return self


	def crossvalidate(self,y,alphas,n_splits=10):
		"""
		lmmlasso cross-validation to get optimal alpha
		alphas = list of alphas to perform cross-validation over
		y = phenotype
		"""
		lasso=lmmlasso.LmmLasso(warm_start=True,fit_intercept=False,tol=0.5)
		X=self.E
		K=self.K

		assert K is not None, 'no kinship matrix defined'
		MSE_train,MSE_test,W_nonzero, rsquared = lmmlasso.runCrossValidation(lasso,self.E,y,alphas,n_splits=n_splits,K=K,verbose=True)
		train_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_train.mean(axis=0))).derivative(n=2)#Interpolating the values for alphas within the range
		test_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_test.mean(axis=0))).derivative(n=2)
		alphas_inter = (sp.linspace(min(alphas),max(alphas),100))
		idx_train = sp.argmin(train_inter(alphas_inter))  # :/
		idx_test = sp.argmin(test_inter(alphas_inter))   # :/
		alpha_cv = (float(alphas_inter[idx_train])+float(alphas_inter[idx_test]))/2
		self.alpha=alpha_cv
		return self.alpha

	def fitmodel(self, y, alpha = None,tol=0.05):
		"""
		fit the model
		"""
		if alpha is None:
			assert self.alpha is not None, 'Set an alpha value'
		else:
			self.alpha = alpha
		assert y.shape[0]==self.G.shape[0], 'No. of observations does not match'
		lasso=lmmlasso.LmmLasso(warm_start=True,fit_intercept=False,tol=tol)
		X=self.E
		lasso.set_params(alpha=self.alpha)
		assert self.K is not None, 'Include a Kinship matrix'
		lasso=lasso.fit(X,y,K=self.K)
		self.feweights=lasso.coef_ #fixed effect weights
		self.feweights=np.reshape(self.feweights,(-1,1))
		self.lasso=lasso

	def predict(self, tmin,tmax,day,Ktest=None,X=None,div=100):
		"""
		predict phenotype across raster given the snp data or the k-matrix for one accession
		you only need to provide either the snps or the K-matrix

		tmin = tif object of minimum daily temperature (time, lat, long)
		tmax = tif object of maximum daily temperature (time, lat, long)
		X = SNP data of the accession you want to predict
		Kstar = kinship matrix of the ecotype [relative to the training dataset] of length (y x 1)
		day = planting date (the first day = 1)
		div = value to divide raster observations by. default is 100
			then you could just divide the length of the weights which are just days of environment x 2 by 2 to get the end date

		"""
		if Ktest is None:
			assert X is not None, 'Need to provide either SNPs or K-matrix'
		if X is not None:
			if Ktest is None:
				X=X.reshape(-1,self.G.shape[1])
				Ktest=np.zeros((1,self.G.shape[0]))
				for i in range(X.shape[0]):
					for j in range(self.G.shape[0]):
						total = sum (operator.eq(self.G[i,:],self.G[j,:]))
						prop = operator.truediv(total,self.G.shape[1])
						Ktest[i,j]=prop
		le=len(self.feweights) #number of days of data
		ndays=le//2
		day2=day+ndays
		day2=int(day2)
		a=range(day,(day2))
		arraytmax=tmax.read((a)) #numpy array of values with dimensions (time, lat, long)
		lati=arraytmax.shape[1]
		longi=arraytmax.shape[2]
		arraytmin=tmin.read((a))
		maxtemp=arraytmax.transpose((1,2,0))
		mintemp=arraytmin.transpose((1,2,0))
		maxtemp=maxtemp.reshape((-1,ndays))
		mintemp=mintemp.reshape((-1,ndays))
		maxtemp=maxtemp/div
		mintemp=mintemp/div
		predictors=np.concatenate((mintemp,maxtemp),axis=1)
		pe=np.matmul(predictors,self.feweights)


		pnull=np.zeros((maxtemp.shape[0],(ndays*2)))
		Kte=np.tile(Ktest,((lati*longi),1))
		assert Kte.shape[0]==(lati*longi),'HMM'

		pg=self.lasso.predict(pnull,Kte)
		pg=pg.reshape(-1,1)

		pred=pg+pe
		pred=pred.reshape(lati,longi)
		self.predictions=pred
		return pred

	def tifsave(arr,ref,filename):
		"""
		save numpy array into a tif file with coordinate system(?) based on a reference file
		ref = reference tif file
		"""
		meta = tx.profile
		naip_transform = meta["transform"]
		naip_crs = meta["crs"]
		meta['count'] = 1
		meta['dtype'] = "float64"
		with rasterio.open(filename, 'w', **meta) as dst: #version 2: if any missing data, set to -9999
			dst.write(arr, 1)

def getIBS(G,ids):
	"""
	calculate identity-by-state matrix
	G: SNP data
	ids: accession / individual ids
	"""
	n=G.shape[0]
	K=np.zeros((n,n))
	nsnps=G.shape[1]
	for i in np.unique(ids):
		for j in np.unique(ids):
			if i==j:
				prop=1
				v=i==ids
				vo=np.where(v)[0][0]
				idx_i=np.where(v)
				idx_i=np.array(idx_i).flatten()
				isnp=G[vo,:]
				s=j==ids
				so=np.where(s)[0][0]
				idx_j=np.where(s)
				idx_j=np.array(idx_j).flatten()
				for t in itertools.product(idx_i,idx_j):
					print(t)
					K[t]=prop
			else:
				v=i==ids #logical for i
				vo=np.where(v)[0][0]
				idx_i=np.where(v)
				idx_i=np.array(idx_i).flatten()
				isnp=G[vo,:]
				s=j==ids
				so=np.where(s)[0][0]
				idx_j=np.where(s)
				idx_j=np.array(idx_j).flatten()
				jsnp=G[so,:]
				total=sum(isnp==jsnp)
				prop=total/nsnps
				for t in itertools.product(idx_i,idx_j):
					print(t)
					K[t]=prop
	return K 