#Making monthly microclimate variables

setwd("~/Clim_GWAS/Clim_GWAS_2")
load('full_environ_data_list.Robj') #load environmental data
load("AllPlantings_Corrected.RData")
seasons=c("HalleFall2006","NorwichSpring2007","NorwichSummer2006","NorwichSummer2007","OuluFall2007","ValenciaFall2006","NorwichFall2006")
daylengths=c("short","long","long","long","long","short","short","short") #I'm assuming summer & spring = long, fall = short
envs=c(3,5,6,7,8,10,4)
names(environ_data)[c(3,5,6,7,8,10,4)]



site_ptus = do.call(rbind, lapply(1:length(envs), function(i){
  
  #i=1
  env = environ_data[[envs[i]]] #load corresponding environmental data file
  
  env$TT = with(env,pmax(Grnd.Tmp-3,0)) #thermal time
  env$PTT = with(env,pmax(Grnd.Tmp - 3,0)*Hrs.light) #photothermal time
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

#Separating microclimate data into each planting
for (i in unique(site_ptus$season)){
	assign(paste0(i,"_site_ptus"),site_ptus[site_ptus$season==i,])
}
####################################################################################################################################################################################################################



######## muT monthly 

N=720 #  hourly period
Th=3 # Temp threshold
hnu=0

df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	muT=mean(ENV.sub$Grnd.Tmp)
	df=c(df,muT)
	}
	assign(paste0(ENV$season[1],"_muT_monthly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_muT_monthly")),length)
lengths=unlist(lengths)

muT_monthly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
muT_monthly=data.frame(muT_monthly)
muT_monthly_colnames=c()

for (i in 1:(ncol(muT_monthly)-1)){
	names=paste0("muT_Month",i)
	muT_monthly_colnames=c(muT_monthly_colnames,names)
}
colnames(muT_monthly)[-c(1)]=muT_monthly_colnames
colnames(muT_monthly)[1]="Planting"
muT_monthly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(muT_monthly)){
	muT_monthly[i,-c(1)]=eval(as.name(paste0(muT_monthly$Planting[i],"_muT_monthly")))[1:(ncol(muT_monthly)-1)]
}


#Mean monthly daylength


df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	muDL=mean(ENV.sub$Daylength)
	df=c(df,muDL)
	}
	assign(paste0(ENV$season[1],"_muDL_monthly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_muDL_monthly")),length)
lengths=unlist(lengths)

muDL_monthly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
muDL_monthly=data.frame(muDL_monthly)
muDL_monthly_colnames=c()

for (i in 1:(ncol(muDL_monthly)-1)){
	names=paste0("muDL_Month",i)
	muDL_monthly_colnames=c(muDL_monthly_colnames,names)
}
colnames(muDL_monthly)[-c(1)]=muDL_monthly_colnames
colnames(muDL_monthly)[1]="Planting"
muDL_monthly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(muDL_monthly)){
	muDL_monthly[i,-c(1)]=eval(as.name(paste0(muDL_monthly$Planting[i],"_muDL_monthly")))[1:(ncol(muDL_monthly)-1)]
}

#Mean monthly thermal time


df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	muTT=mean(ENV.sub$TT)
	df=c(df,muTT)
	}
	assign(paste0(ENV$season[1],"_muTT_monthly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_muTT_monthly")),length)
lengths=unlist(lengths)

muTT_monthly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
muTT_monthly=data.frame(muTT_monthly)
muTT_monthly_colnames=c()

for (i in 1:(ncol(muTT_monthly)-1)){
	names=paste0("muTT_Month",i)
	muTT_monthly_colnames=c(muTT_monthly_colnames,names)
}
colnames(muTT_monthly)[-c(1)]=muTT_monthly_colnames
colnames(muTT_monthly)[1]="Planting"
muTT_monthly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(muTT_monthly)){
	muTT_monthly[i,-c(1)]=eval(as.name(paste0(muTT_monthly$Planting[i],"_muTT_monthly")))[1:(ncol(muTT_monthly)-1)]
}
#Mean monthly [PAR (Photosynthetically active radiation) * Hrs.Light]; a bit different from PTT which is (temp >3 * hours daylight)
#muPART (PAR time)





df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	muPART=mean((ENV.sub$PAR)*(ENV.sub$Hrs.light))
	df=c(df,muPART)
	}
	assign(paste0(ENV$season[1],"_muPART_monthly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_muPART_monthly")),length)
lengths=unlist(lengths)

muPART_monthly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
muPART_monthly=data.frame(muPART_monthly)
muPART_monthly_colnames=c()

for (i in 1:(ncol(muPART_monthly)-1)){
	names=paste0("muPART_Month",i)
	muPART_monthly_colnames=c(muPART_monthly_colnames,names)
}
colnames(muPART_monthly)[-c(1)]=muPART_monthly_colnames
colnames(muPART_monthly)[1]="Planting"
muPART_monthly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(muPART_monthly)){
	muPART_monthly[i,-c(1)]=eval(as.name(paste0(muPART_monthly$Planting[i],"_muPART_monthly")))[1:(ncol(muPART_monthly)-1)]
}



#Mean monthly PAR


df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	muPAR=mean((ENV.sub$PAR))
	df=c(df,muPAR)
	}
	assign(paste0(ENV$season[1],"_muPAR_monthly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_muPAR_monthly")),length)
lengths=unlist(lengths)

muPAR_monthly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
muPAR_monthly=data.frame(muPAR_monthly)
muPAR_monthly_colnames=c()

for (i in 1:(ncol(muPAR_monthly)-1)){
	names=paste0("muPAR_Month",i)
	muPAR_monthly_colnames=c(muPAR_monthly_colnames,names)
}
colnames(muPAR_monthly)[-c(1)]=muPAR_monthly_colnames
colnames(muPAR_monthly)[1]="Planting"
muPAR_monthly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(muPAR_monthly)){
	muPAR_monthly[i,-c(1)]=eval(as.name(paste0(muPAR_monthly$Planting[i],"_muPAR_monthly")))[1:(ncol(muPAR_monthly)-1)]
}


#Mean Monthly Thermal Time
df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	muTT=mean((ENV.sub$TT))
	df=c(df,muTT)
	}
	assign(paste0(ENV$season[1],"_muTT_monthly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_muTT_monthly")),length)
lengths=unlist(lengths)

muTT_monthly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
muTT_monthly=data.frame(muTT_monthly)
muTT_monthly_colnames=c()

for (i in 1:(ncol(muTT_monthly)-1)){
	names=paste0("muTT_Month",i)
	muTT_monthly_colnames=c(muTT_monthly_colnames,names)
}
colnames(muTT_monthly)[-c(1)]=muTT_monthly_colnames
colnames(muTT_monthly)[1]="Planting"
muTT_monthly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(muTT_monthly)){
	muTT_monthly[i,-c(1)]=eval(as.name(paste0(muTT_monthly$Planting[i],"_muTT_monthly")))[1:(ncol(muTT_monthly)-1)]
}


#Mean Monthly Phototermal Time
df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	muPTT=mean((ENV.sub$TT))
	df=c(df,muPTT)
	}
	assign(paste0(ENV$season[1],"_muPTT_monthly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_muPTT_monthly")),length)
lengths=unlist(lengths)

muPTT_monthly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
muPTT_monthly=data.frame(muPTT_monthly)
muPTT_monthly_colnames=c()

for (i in 1:(ncol(muPTT_monthly)-1)){
	names=paste0("muPTT_Month",i)
	muPTT_monthly_colnames=c(muPTT_monthly_colnames,names)
}
colnames(muPTT_monthly)[-c(1)]=muPTT_monthly_colnames
colnames(muPTT_monthly)[1]="Planting"
muPTT_monthly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(muPTT_monthly)){
	muPTT_monthly[i,-c(1)]=eval(as.name(paste0(muPTT_monthly$Planting[i],"_muPTT_monthly")))[1:(ncol(muPTT_monthly)-1)]
}



microclim_monthly=cbind(muDL_monthly,muT_monthly,muPAR_monthly,muPART_monthly,muTT_monthly,muPTT_monthly)
microclim_monthly=microclim_monthly[,colnames(microclim_monthly)!="Planting"]
save(microclim_monthly,file="Microclimate_monthly.RData")