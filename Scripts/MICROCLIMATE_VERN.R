#!/usr/bin/env Rscript

#Microclimate with vernalisation parameter included
##STEP ONE: CREATING THE RAW HOURLY MICROCLIMATIC VARIABLES, SEPARATED BY PLANTING

setwd("~/Clim_GWAS/Clim_GWAS_2")
load('full_environ_data_list.Robj') #load environmental data
load('Multitrait_NEW.RData') #Load the phenotype data
library(tictoc)
tic("PROCESSING MICROCLIMATIC DATA")
print("Starting...")

args = commandArgs(trailingOnly=TRUE)
#Parameters
Th=as.numeric(args[1]) #induction temperature threshold
N=as.numeric(args[2]) #period in hours

message(sprintf("Induction temp = %s", args[1]))#Prints a message to bash
message(sprintf("Period %s", args[2]))#Prints a message to bash



seasons=c("HalleFall2006","NorwichSpring2007","NorwichSummer2006","NorwichSummer2007","OuluFall2007","ValenciaFall2006","NorwichFall2006")
daylengths=c("short","long","long","long","long","short","short","short") #I'm assuming summer & spring = long, fall = short
envs=c(3,5,6,7,8,10,4)



site_ptus = do.call(rbind, lapply(1:length(envs), function(i){
  
  #i=1
  env = environ_data[[envs[i]]] #load corresponding environmental data file
  
  env$TT = with(env,pmax(Grnd.Tmp- Th,0)) #thermal time
  env$PTT = with(env,pmax(Grnd.Tmp - Th,0)*Hrs.light) #photothermal time
  env$cumPTT = cumsum(env$PTT) #accumulated thermal time
  
  #env$daylength=daylengths[i]
  env$season=seasons[i]
  env$hour = seq(1:nrow(env))
  env$day=floor(env$hour/24) #get days since sowing
  
  envhrs = do.call(rbind, lapply(1:nrow(env), function(j){
    #j=1
    if(env$Daylength[j]>=12){"long"} else {"short"}
  }))
  env$daylength = as.vector(envhrs)
  env
}))

for (i in unique(site_ptus$season)){
	assign(paste0(i,"_site_ptus"),site_ptus[site_ptus$season==i,])
}


##








##STEP TWO: PROCESSING THIS HOURLY MICROCLIM DATA TO GET THE LEVEL OF RESOLUTION WE DESIRE

#Parameter
	if(N==720){
		b="Monthly"
	} else if (N==168){
		b="Weekly"
		}else if (N==24){
			b="Daily"
			}else if (N==1){
				b="Hourly"
				}else(print("Reconsider period choice"))
				


df_list = mget(ls(pattern = "\\_site_ptus"))

#Ground Temp

	for (ENV in df_list){
		df=c()
		for (i in 1:floor((nrow(ENV)/N))) {
		ENV.sub=ENV[(i-1)*N+1:(i*N),]
		ENV.sub=na.omit(ENV.sub)
		muT=mean(ENV.sub$Grnd.Tmp)
		df=c(df,muT)
		}
		assign(paste0(ENV$season[1],"_groundtemp_",b),df)}

	#Making the matrix
	lengths=lapply(mget(ls(pattern=paste0("_groundtemp_",b))),length)
	lengths=unlist(lengths)


	assign(paste0("groundtemp_",b),matrix(NA,nrow=nrow(Multitrait_NEW),ncol=(min(lengths)+1)))
	assign(paste0("groundtemp_",b),as.data.frame(eval(as.name(paste0("groundtemp_",b)))))
	groundtemp_colnames=c()

	for (i in 1:(ncol(eval(as.name(paste0("groundtemp_",b))))-1)){
		names=paste0("groundtemp_",b,"_",i)
		groundtemp_colnames=c(groundtemp_colnames,names)
	}
	a=eval(as.name(paste0("groundtemp_",b)))
	colnames(a)[-1]=groundtemp_colnames
	colnames(a)[1]="Planting"
	a[,1]=Multitrait_NEW$Planting



	for (i in 1:nrow(a)){
		a[i,-c(1)]=eval(as.name(paste0(a$Planting[i],"_groundtemp_",b)))[1:(ncol(a)-1)]
	}
	assign(paste0("groundtemp_",b),a)


#Daylength

	for (ENV in df_list){
			df=c()
			for (i in 1:floor((nrow(ENV)/N))) {
			ENV.sub=ENV[(i-1)*N+1:(i*N),]
			ENV.sub=na.omit(ENV.sub)
			muT=mean(ENV.sub$Daylength)
			df=c(df,muT)
			}
			assign(paste0(ENV$season[1],"_daylength_",b),df)}

	#Making the matrix
	lengths=lapply(mget(ls(pattern=paste0("_daylength_",b))),length)
		lengths=unlist(lengths)


	assign(paste0("daylength_",b),matrix(NA,nrow=nrow(Multitrait_NEW),ncol=(min(lengths)+1)))
	assign(paste0("daylength_",b),as.data.frame(eval(as.name(paste0("daylength_",b)))))
	daylength_colnames=c()

	for (i in 1:(ncol(eval(as.name(paste0("daylength_",b))))-1)){
		names=paste0("daylength_",b,"_",i)
		daylength_colnames=c(daylength_colnames,names)
	}
	a=eval(as.name(paste0("daylength_",b)))
	colnames(a)[-1]=daylength_colnames
	colnames(a)[1]="Planting"
	a[,1]=Multitrait_NEW$Planting



	for (i in 1:nrow(a)){
		a[i,-c(1)]=eval(as.name(paste0(a$Planting[i],"_daylength_",b)))[1:(ncol(a)-1)]
	}
	assign(paste0("daylength_",b),a)


#PAR

for (ENV in df_list){
			df=c()
			for (i in 1:floor((nrow(ENV)/N))) {
			ENV.sub=ENV[(i-1)*N+1:(i*N),]
			ENV.sub=na.omit(ENV.sub)
			muT=mean(ENV.sub$PAR)
			df=c(df,muT)
			}
			assign(paste0(ENV$season[1],"_PAR_",b),df)}

	#Making the matrix
	lengths=lapply(mget(ls(pattern=paste0("_PAR_",b))),length)
		lengths=unlist(lengths)


	assign(paste0("PAR_",b),matrix(NA,nrow=nrow(Multitrait_NEW),ncol=(min(lengths)+1)))
	assign(paste0("PAR_",b),as.data.frame(eval(as.name(paste0("PAR_",b)))))
	PAR_colnames=c()

	for (i in 1:(ncol(eval(as.name(paste0("PAR_",b))))-1)){
		names=paste0("PAR_",b,"_",i)
		PAR_colnames=c(PAR_colnames,names)
	}
	a=eval(as.name(paste0("PAR_",b)))
	colnames(a)[-1]=PAR_colnames
	colnames(a)[1]="Planting"
	a[,1]=Multitrait_NEW$Planting



	for (i in 1:nrow(a)){
		a[i,-c(1)]=eval(as.name(paste0(a$Planting[i],"_PAR_",b)))[1:(ncol(a)-1)]
	}
	assign(paste0("PAR_",b),a)




#TT

for (ENV in df_list){
			df=c()
			for (i in 1:floor((nrow(ENV)/N))) {
			ENV.sub=ENV[(i-1)*N+1:(i*N),]
			ENV.sub=na.omit(ENV.sub)
			muT=mean(ENV.sub$TT)
			df=c(df,muT)
			}
			assign(paste0(ENV$season[1],"_TT_",b),df)}

	#Making the matrix
	lengths=lapply(mget(ls(pattern=paste0("_TT_",b))),length)
		lengths=unlist(lengths)


	assign(paste0("TT_",b),matrix(NA,nrow=nrow(Multitrait_NEW),ncol=(min(lengths)+1)))
	assign(paste0("TT_",b),as.data.frame(eval(as.name(paste0("TT_",b)))))
	TT_colnames=c()

	for (i in 1:(ncol(eval(as.name(paste0("TT_",b))))-1)){
		names=paste0("TT_",b,"_",i)
		TT_colnames=c(TT_colnames,names)
	}
	a=eval(as.name(paste0("TT_",b)))
	colnames(a)[-1]=TT_colnames
	colnames(a)[1]="Planting"
	a[,1]=Multitrait_NEW$Planting



	for (i in 1:nrow(a)){
		a[i,-c(1)]=eval(as.name(paste0(a$Planting[i],"_TT_",b)))[1:(ncol(a)-1)]
	}
	assign(paste0("TT_",b),a)

#PTT


for (ENV in df_list){
			df=c()
			for (i in 1:floor((nrow(ENV)/N))) {
			ENV.sub=ENV[(i-1)*N+1:(i*N),]
			ENV.sub=na.omit(ENV.sub)
			muT=mean(ENV.sub$PTT)
			df=c(df,muT)
			}
			assign(paste0(ENV$season[1],"_PTT_",b),df)}

	#Making the matrix
	lengths=lapply(mget(ls(pattern=paste0("_PTT_",b))),length)
		lengths=unlist(lengths)


	assign(paste0("PTT_",b),matrix(NA,nrow=nrow(Multitrait_NEW),ncol=(min(lengths)+1)))
	assign(paste0("PTT_",b),as.data.frame(eval(as.name(paste0("PTT_",b)))))
	PTT_colnames=c()

	for (i in 1:(ncol(eval(as.name(paste0("PTT_",b))))-1)){
		names=paste0("PTT_",b,"_",i)
		PTT_colnames=c(PTT_colnames,names)
	}
	a=eval(as.name(paste0("PTT_",b)))
	colnames(a)[-1]=PTT_colnames
	colnames(a)[1]="Planting"
	a[,1]=Multitrait_NEW$Planting



	for (i in 1:nrow(a)){
		a[i,-c(1)]=eval(as.name(paste0(a$Planting[i],"_PTT_",b)))[1:(ncol(a)-1)]
	}
	assign(paste0("PTT_",b),a)

#Vernalisation
is.nan.data.frame <- function(x)
do.call(cbind, lapply(x, is.nan))

Tv_min=-3.5 #Minimum vernalisation temperature [this is not the same thing as the threshold temp]
Tv_max=6 #Maximum vernalisation temperature
Vsat=960 #Saturation point where FLC is permanently repressed. 960 h = 40 days; based on literature.
kappa=-5.17 #some sort of parameter; value from Chew (2012) Table S2
omega=2.23 #some sort of parameter; value from Chew (2012) Table S2
xi=1.00 #some sort of parameter; value from Chew (2012) Table S2

vernalisation_effectiveness <- function(T) { #T as in Grnd.Tmp
	ve= exp(kappa)* ((T - Tv_min)^(omega))*((Tv_max - T)^xi)
	return(ve)
}


for (ENV in df_list){
			vernalisation=c()
			for (i in 1:floor((nrow(ENV)/N))) {
			ENV.sub=ENV[(i-1)*N+1:(i*N),]
			ENV.sub=na.omit(ENV.sub)
			vern.sub=vernalisation_effectiveness(ENV.sub$Grnd.Tmp[i])								
			vernalisation=c(vernalisation,vern.sub)
			vernalisation[is.nan(vernalisation)] <- 0
			cumsumv=cumsum(vernalisation)
			}
			assign(paste0(ENV$season[1],"_Vernalisation_",b),vernalisation)
			assign(paste0(ENV$season[1],"_cumsumVernalisation_",b),cumsumv)}

#Making the matrix
	lengths=lapply(mget(ls(pattern="_Vernalisation_")),length)
		lengths=unlist(lengths)

	assign(paste0("Vernalisation_",b),matrix(NA,nrow=nrow(Multitrait_NEW),ncol=(min(lengths)+1)))
	assign(paste0("Vernalisation_",b),as.data.frame(eval(as.name(paste0("Vernalisation_",b)))))
	Vernalisation_colnames=c()


	for (i in 1:(ncol(eval(as.name(paste0("Vernalisation_",b))))-1)){
		names=paste0("Vernalisation_",b,"_",i)
		Vernalisation_colnames=c(Vernalisation_colnames,names)
	}
	a=eval(as.name(paste0("Vernalisation_",b)))
	colnames(a)[-1]=Vernalisation_colnames
	colnames(a)[1]="Planting"
	a[,1]=Multitrait_NEW$Planting

	for (i in 1:nrow(a)){
		a[i,-c(1)]=eval(as.name(paste0(a$Planting[i],"_Vernalisation_",b)))[1:(ncol(a)-1)]
	}
	assign(paste0("Vernalisation_",b),a)

#Making the cumsum matrix

	lengths=lapply(mget(ls(pattern="_cumsumVernalisation_")),length)
		lengths=unlist(lengths)

	assign(paste0("cumsumVernalisation_",b),matrix(NA,nrow=nrow(Multitrait_NEW),ncol=(min(lengths)+1)))
	assign(paste0("cumsumVernalisation_",b),as.data.frame(eval(as.name(paste0("cumsumVernalisation_",b)))))
	cumsumVernalisation_colnames=c()


	for (i in 1:(ncol(eval(as.name(paste0("cumsumVernalisation_",b))))-1)){
		names=paste0("cumsumVernalisation_",b,"_",i)
		cumsumVernalisation_colnames=c(cumsumVernalisation_colnames,names)
	}
	a=eval(as.name(paste0("cumsumVernalisation_",b)))
	colnames(a)[-1]=cumsumVernalisation_colnames
	colnames(a)[1]="Planting"
	a[,1]=Multitrait_NEW$Planting

	for (i in 1:nrow(a)){
		a[i,-c(1)]=eval(as.name(paste0(a$Planting[i],"_cumsumVernalisation_",b)))[1:(ncol(a)-1)]
	}
	assign(paste0("cumsumVernalisation_",b),a)






	combined = cbind(
	(eval(as.name(paste0("groundtemp_",b)))),
	(eval(as.name(paste0("daylength_",b))))[,-1],
	(eval(as.name(paste0("PAR_",b))))[,-1],
	(eval(as.name(paste0("TT_",b))))[,-1],
	(eval(as.name(paste0("PTT_",b))))[,-1],
	(eval(as.name(paste0("Vernalisation_",b))))
	)
combined=combined[,-1]
print(combined[1:5,1:5])

#Save the output

setwd("~/Clim_GWAS/Clim_GWAS_2/Temp_Files")

filename=paste0("Microclimate_vern_",b,"_threshold_",Th)

write.csv(combined,file=paste0(filename,".csv"),row.names=FALSE)

print("Finished...")
toc()

#Vernalisation




























































#Determining vernalisation


#Calculating hourly vernalisation

#Should add the vernalisation values to each planting's _site_ptus.

