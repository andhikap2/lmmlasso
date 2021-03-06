#!/usr/bin/env python3
from numpy import loadtxt
from keras.models import Sequential
from keras.layers import Dense
import numpy as np 
import os
import pandas as pd 
import sys
from sklearn.preprocessing import StandardScaler

hl1=int(sys.argv[1])
hl2=int(sys.argv[2])
hl3=int(sys.argv[3])
hl4=int(sys.argv[4])
ep=int(sys.argv[5])

# load the dataset
# split into input (X) and output (y) variables
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
K=np.load("K_MATRIX_EMMAXIBS.npy")
dtb=np.load('DTB_NEW.npy')
dtb=dtb.reshape(-1,1)
seednum=np.load('SeedNum.npy')
seednum=seednum.reshape(-1,1)
lfl=np.load('LeafLength.npy')
lfl=lfl.reshape(-1,1)




os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/New_Microclim")
env=pd.read_csv('Microclimate_minmaxDaily_threshold_0_0.csv',sep=',')
logical=env.columns.str.startswith(('PAR','TT','PTT','daylength') )
env=env.iloc[:,~logical] #get columns that don't start with PAR
environment=np.array(env)
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")
snps=np.load('SNPs_0.1.npy')

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
pca=PCA(n_components=0.95)
pca.fit(K)
K=pca.transform(K)

X = np.concatenate((snps,environment,K),axis=1)
y = np.concatenate((dtb,lfl,seednum),axis=1)
scaler=StandardScaler(copy=False)
scaler.fit(y)
y = scaler.transform(y) #scale so it has unit variance and unit sd

# define the keras model
model = Sequential()
model.add(Dense(hl1, input_dim=X.shape[1], activation='relu'))
model.add(Dense(hl2, activation='relu'))
model.add(Dense(hl3,activation='relu'))
model.add(Dense(hl4,activation='relu'))
model.add(Dense(3, activation='relu')) #corresponds to the number of traits to be predicted

# compile the keras model
model.compile(loss='mean_squared_error', optimizer='adam',metrics=['mse']) #mse is mean squared error
# fit the keras model on the dataset
model.fit(X, y, epochs=ep, batch_size=10)
# evaluate the keras model (on the same dataset)
predictions=model.predict(X)
y_back=scaler.inverse_transform(y) #bring it back
y_predictions=scaler.inverse_transform(predictions)
print(predictions)


from sklearn.metrics import mean_squared_error
rmse_dtb=np.sqrt(mean_squared_error(y_back[:,0],y_predictions[:,0]))
rmse_lfl=np.sqrt(mean_squared_error(y_back[:,1],y_predictions[:,1]))
rmse_sn=np.sqrt(mean_squared_error(y_back[:,2],y_predictions[:,2]))

keras_results=np.array((hl1,hl2,hl3,hl4,ep,rmse_dtb,rmse_lfl,rmse_sn))
kpd=pd.DataFrame(keras_results)
kpd=kpd.transpose()
kpd.columns=(("HL1","HL2","HL3","HL4","Epochs","DTBrmse","LfLrmse","SNrmse"))
kpd.to_csv("Keras_"+str(hl1)+"."+str(hl2)+"."+str(hl3)+"."+str(hl4)+"."+str(ep)+".csv")


#Montesinos-Lopez used either one, two, or three hidden layers and optimized via grid search; could probably do this in SLURM and choose based on minimizing the loss function