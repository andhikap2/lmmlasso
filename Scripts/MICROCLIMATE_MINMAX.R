#!/usr/bin/env Rscript

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
PAR_Th=as.numeric(args[3]) #PAR threshold

message(sprintf("Induction temp = %s", args[1]))#Prints a message to bash
message(sprintf("Period %s", args[2]))#Prints a message to bash
message(sprintf("PAR threshold %s",args[3]))


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


# Setting the PAR threshold

df_list = mget(ls(pattern = "\\_site_ptus"))

for (ENV in df_list){
	ENV$PAR[ENV$PAR<=PAR_Th]<-0
	assign(paste0(ENV$season[1],"_site_ptus"),ENV)
}

df_list = mget(ls(pattern = "\\_site_ptus"))




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
				




#Max Ground Temp

	for (ENV in df_list){
		df=c()
		for (i in 1:floor((nrow(ENV)/N))) {
		ENV.sub=ENV[((i-1)*N+1):(i*N),]
		ENV.sub=na.omit(ENV.sub)
		ENV.sub=na.omit(ENV.sub)
		maxT=max(ENV.sub$Grnd.Tmp)
		df=c(df,maxT)
		}
		assign(paste0(ENV$season[1],"_maxgroundtemp_",b),df)}

	#Making the matrix
	lengths=lapply(mget(ls(pattern=paste0("_maxgroundtemp_",b))),length)
	lengths=unlist(lengths)


	assign(paste0("maxgroundtemp_",b),matrix(NA,nrow=nrow(Multitrait_NEW),ncol=(min(lengths)+1)))
	assign(paste0("maxgroundtemp_",b),as.data.frame(eval(as.name(paste0("maxgroundtemp_",b)))))
	groundtemp_colnames=c()

	for (i in 1:(ncol(eval(as.name(paste0("maxgroundtemp_",b))))-1)){
		names=paste0("maxgroundtemp_",b,"_",i)
		groundtemp_colnames=c(groundtemp_colnames,names)
	}
	a=eval(as.name(paste0("maxgroundtemp_",b)))
	colnames(a)[-1]=groundtemp_colnames
	colnames(a)[1]="Planting"
	a[,1]=Multitrait_NEW$Planting



	for (i in 1:nrow(a)){
		a[i,-c(1)]=eval(as.name(paste0(a$Planting[i],"_maxgroundtemp_",b)))[1:(ncol(a)-1)]
	}
	assign(paste0("maxgroundtemp_",b),a)




#Min Ground Temp

	for (ENV in df_list){
		df=c()
		for (i in 1:floor((nrow(ENV)/N))) {
		ENV.sub=ENV[((i-1)*N+1):(i*N),]
		ENV.sub=na.omit(ENV.sub)
		minT=min(ENV.sub$Grnd.Tmp)
		df=c(df,minT)
		}
		assign(paste0(ENV$season[1],"_mingroundtemp_",b),df)}

	#Making the matrix
	lengths=lapply(mget(ls(pattern=paste0("_mingroundtemp_",b))),length)
	lengths=unlist(lengths)


	assign(paste0("mingroundtemp_",b),matrix(NA,nrow=nrow(Multitrait_NEW),ncol=(min(lengths)+1)))
	assign(paste0("mingroundtemp_",b),as.data.frame(eval(as.name(paste0("maxgroundtemp_",b)))))
	groundtemp_colnames=c()

	for (i in 1:(ncol(eval(as.name(paste0("mingroundtemp_",b))))-1)){
		names=paste0("mingroundtemp_",b,"_",i)
		groundtemp_colnames=c(groundtemp_colnames,names)
	}
	a=eval(as.name(paste0("mingroundtemp_",b)))
	colnames(a)[-1]=groundtemp_colnames
	colnames(a)[1]="Planting"
	a[,1]=Multitrait_NEW$Planting



	for (i in 1:nrow(a)){
		a[i,-c(1)]=eval(as.name(paste0(a$Planting[i],"_mingroundtemp_",b)))[1:(ncol(a)-1)]
	}
	assign(paste0("mingroundtemp_",b),a)


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


	combined = cbind(
	(eval(as.name(paste0("mingroundtemp_",b)))),
	(eval(as.name(paste0("maxgroundtemp_",b))))[,-1],
	(eval(as.name(paste0("daylength_",b))))[,-1],
	(eval(as.name(paste0("PAR_",b))))[,-1],
	(eval(as.name(paste0("TT_",b))))[,-1],
	(eval(as.name(paste0("PTT_",b))))[,-1]
	)
combined=combined[,-1]
print(combined[1:5,1:5])

#Save the output

setwd("~/Clim_GWAS/New_Microclim")

filename=paste0("Microclimate_minmax",b,"_threshold_",Th,"_",PAR_Th)

write.csv(combined,file=paste0(filename,".csv"),row.names=FALSE)

print("Finished...")
toc()



#Save as something
















