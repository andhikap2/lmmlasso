#Plotting the weights

import os
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/New_Microclim")

#Load the dataframe to get the column names
import pandas as pd 

microclim=pd.read_csv('Microclimate_Hourly_threshold_7_0.csv')

colnames=microclim.columns

#Load the weights
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/LMMLASSO/Weights_Corrected")

import numpy as np 
weights=np.load("Weights_SNPs_K_Microclimate_Hourly_threshold_7_0.npy")

#Convert the numpy array into a dataframe
weights_pd=pd.DataFrame(weights)
weights_pd=pd.DataFrame.transpose(weights_pd)
weights_pd.columns=colnames




#Separating the data into the weights for each microclimatic variable
groundtemp=weights_pd.loc[:,weights_pd.columns.str.startswith(pat = 'groundtemp')]
daylength=weights_pd.loc[:,weights_pd.columns.str.startswith(pat = 'daylength')]
PAR=weights_pd.loc[:,weights_pd.columns.str.startswith(pat = 'PAR')]
TT=weights_pd.loc[:,weights_pd.columns.str.startswith(pat = 'TT_')]
PTT=weights_pd.loc[:,weights_pd.columns.str.startswith(pat = 'PTT_')]



#Standardizing
groundtemp=groundtemp.sub(groundtemp.mean(1), axis=0).div(groundtemp.std(1), axis=0)
daylength=daylength.sub(daylength.mean(1), axis=0).div(daylength.std(1), axis=0)
PAR=PAR.sub(PAR.mean(1), axis=0).div(PAR.std(1), axis=0)
TT=TT.sub(TT.mean(1), axis=0).div(TT.std(1), axis=0)
PTT=PTT.sub(PTT.mean(1), axis=0).div(PTT.std(1), axis=0)


os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Figures")

#Plotting them

import matplotlib.pyplot as plt 
matplotlib.use('TkAgg')
ranges=range(groundtemp.shape[1])
fig, axs = plt.subplots(3, 2, gridspec_kw={'hspace': 0.5})
axs[0,0].plot(ranges,pd.DataFrame.transpose(groundtemp))
axs[0,0].title.set_text('Ground Temp.')
axs[0,0].set_ylabel('β',rotation=0)
axs[0,0].set_xticks([])
axs[0,1].plot(ranges,pd.DataFrame.transpose(daylength))
axs[0,1].title.set_text('Daylength')
axs[0,1].set_ylabel('β',rotation=0)
axs[0,1].set_xticks([])
axs[1,0].plot(ranges,pd.DataFrame.transpose(PAR))
axs[1,0].title.set_text('PAR')
axs[1,0].set_ylabel('β',rotation=0)
axs[1,0].set_xticks([])
axs[1,1].plot(ranges,pd.DataFrame.transpose(TT))
axs[1,1].title.set_text('Thermal Time')
axs[1,1].set_ylabel('β',rotation=0)
axs[1,1].set_xlabel('Hours since sowing')
axs[2,0].plot(ranges,pd.DataFrame.transpose(PTT))
axs[2,0].title.set_text('Photothermal Time')
axs[2,0].set_xlabel('Hours since sowing')
axs[2,0].set_ylabel('β',rotation=0)
fig.delaxes(axs[2,1])

axs[0,0].yaxis.set_label_coords(-0.1,0.5)
axs[0,1].yaxis.set_label_coords(-0.125,0.5)
axs[1,0].yaxis.set_label_coords(-0.1,0.5)
axs[1,1].yaxis.set_label_coords(-0.125,0.5)
axs[2,0].yaxis.set_label_coords(-0.1,0.5)

plt.savefig('Betas_K_Hourly_Temp3.png',bbox_inches='tight')

















#Plotting the weights but for datasets containing SNPs

import os

#Load the dataframe to get the column names
import pandas as pd 
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/New_Microclim")

microclim=pd.read_csv('Microclimate_Hourly_threshold_7_0.csv')
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")

snp_names=pd.read_csv('SNPs_0.1.prune.in',header=None)

snp_names=pd.DataFrame.transpose(snp_names)
np_snps=np.array(snp_names)

colnames_microclim=microclim.columns
# groundtemp_logical=microclim.columns.str.startswith(pat='groundtemp')
# groundtemp=microclim.loc[:,groundtemp_logical]
# groundtemp_colnames=groundtemp.columns

# np_groundtemp_colnames=np.array(groundtemp_colnames)

colnames=np.concatenate([np_snps[0,:],colnames_microclim])

#Load the weights
os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/LMMLASSO/Weights_Corrected")

import numpy as np 
Hourly_0=np.load("Weights_SNPs_0.1_Microclimate_Hourly_threshold_7_0.npy")

#Convert the numpy array into a dataframe
weights_pd=pd.DataFrame(Hourly_0)
weights_pd=pd.DataFrame.transpose(weights_pd)
weights_pd.columns=colnames



#Separating the snps by chromosome
chr1=weights_pd.loc[:,weights_pd.columns.str.startswith(pat='1_')]
chr2=weights_pd.loc[:,weights_pd.columns.str.startswith(pat='2_')]
chr3=weights_pd.loc[:,weights_pd.columns.str.startswith(pat='3_')]
chr4=weights_pd.loc[:,weights_pd.columns.str.startswith(pat='4_')]
chr5=weights_pd.loc[:,weights_pd.columns.str.startswith(pat='5')]

#Separating the data into the weights for each microclimatic variable
groundtemp=weights_pd.loc[:,weights_pd.columns.str.startswith(pat = 'groundtemp')]
daylength=weights_pd.loc[:,weights_pd.columns.str.startswith(pat = 'daylength')]
PAR=weights_pd.loc[:,weights_pd.columns.str.startswith(pat = 'PAR')]
TT=weights_pd.loc[:,weights_pd.columns.str.startswith(pat = 'TT_')]
PTT=weights_pd.loc[:,weights_pd.columns.str.startswith(pat = 'PTT_')]

os.chdir("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/LMMLASSO/Weights_Corrected")


#Plotting the SNPs
import matplotlib.pyplot as plt 
matplotlib.use('TkAgg')
range1=range(chr1.shape[1])
range2=range(chr2.shape[1])
range3=range(chr3.shape[1])
range4=range(chr4.shape[1])
range5=range(chr5.shape[1])


fig, axs = plt.subplots(3, 2, gridspec_kw={'hspace': 0.5})
axs[0,0].plot(range1,pd.DataFrame.transpose(chr1))
axs[0,0].title.set_text('CHR. 1')
axs[0,0].set_ylabel('β',rotation=0)
axs[0,0].set_xticks([])
axs[0,1].plot(range2,pd.DataFrame.transpose(chr2))
axs[0,1].title.set_text('CHR. 2')
axs[0,1].set_ylabel('β',rotation=0)
axs[0,1].set_xticks([])
axs[1,0].plot(range3,pd.DataFrame.transpose(chr3))
axs[1,0].title.set_text('CHR. 3')
axs[1,0].set_ylabel('β',rotation=0)
axs[1,0].set_xticks([])
axs[1,1].plot(range4,pd.DataFrame.transpose(chr4))
axs[1,1].title.set_text('CHR. 4')
axs[1,1].set_xticks([])
axs[1,1].set_ylabel('β',rotation=0)
axs[2,0].plot(range5,pd.DataFrame.transpose(chr5))
axs[2,0].title.set_text('CHR. 5')
axs[2,0].set_xticks([])
axs[2,0].set_ylabel('β',rotation=0)
fig.delaxes(axs[2,1])

axs[0,0].yaxis.set_label_coords(-0.1,0.5)
axs[0,1].yaxis.set_label_coords(-0.125,0.5)
axs[1,0].yaxis.set_label_coords(-0.1,0.5)
axs[1,1].yaxis.set_label_coords(-0.125,0.5)
axs[2,0].yaxis.set_label_coords(-0.1,0.5)

plt.savefig('BestModel_CHRBetas.png')
#Plotting the microclim
import matplotlib.pyplot as plt 
import matplotlib
matplotlib.use('TkAgg')
ranges=range(groundtemp.shape[1])
fig, axs = plt.subplots(3, 2, gridspec_kw={'hspace': 0.5})
axs[0,0].plot(ranges,pd.DataFrame.transpose(groundtemp))
axs[0,0].title.set_text('Ground Temp.')
axs[0,0].set_ylabel('β',rotation=0)
axs[0,0].set_xticks([])
axs[0,1].plot(ranges,pd.DataFrame.transpose(daylength))
axs[0,1].title.set_text('Daylength')
axs[0,1].set_ylabel('β',rotation=0)
axs[0,1].set_xticks([])
axs[1,0].plot(ranges,pd.DataFrame.transpose(PAR))
axs[1,0].title.set_text('PAR')
axs[1,0].set_ylabel('β',rotation=0)
axs[1,0].set_xticks([])
axs[1,1].plot(ranges,pd.DataFrame.transpose(TT))
axs[1,1].title.set_text('Thermal Time')
axs[1,1].set_ylabel('β',rotation=0)
axs[1,1].set_xlabel('Hours since sowing')
axs[2,0].plot(ranges,pd.DataFrame.transpose(PTT))
axs[2,0].title.set_text('Photothermal Time')
axs[2,0].set_xlabel('Hours since sowing')
axs[2,0].set_ylabel('β',rotation=0)
fig.delaxes(axs[2,1])

axs[0,0].yaxis.set_label_coords(-0.1,0.5)
axs[0,1].yaxis.set_label_coords(-0.125,0.5)
axs[1,0].yaxis.set_label_coords(-0.1,0.5)
axs[1,1].yaxis.set_label_coords(-0.125,0.5)
axs[2,0].yaxis.set_label_coords(-0.1,0.5)
plt.savefig('BestModel_CLIMBetas.png')
