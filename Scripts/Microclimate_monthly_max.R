#Making max microclimate variables

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



######## maxT monthly 

N=720 #  hourly period
Th=3 # Temp threshold
hnu=0

df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	maxT=max(ENV.sub$Grnd.Tmp)
	df=c(df,maxT)
	}
	assign(paste0(ENV$season[1],"_maxT_monthly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_maxT_monthly")),length)
lengths=unlist(lengths)

maxT_monthly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
maxT_monthly=data.frame(maxT_monthly)
maxT_monthly_colnames=c()

for (i in 1:(ncol(maxT_monthly)-1)){
	names=paste0("maxT_Month",i)
	maxT_monthly_colnames=c(maxT_monthly_colnames,names)
}
colnames(maxT_monthly)[-c(1)]=maxT_monthly_colnames
colnames(maxT_monthly)[1]="Planting"
maxT_monthly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(maxT_monthly)){
	maxT_monthly[i,-c(1)]=eval(as.name(paste0(maxT_monthly$Planting[i],"_maxT_monthly")))[1:(ncol(maxT_monthly)-1)]
}


#max monthly daylength


df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	maxDL=max(ENV.sub$Daylength)
	df=c(df,maxDL)
	}
	assign(paste0(ENV$season[1],"_maxDL_monthly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_maxDL_monthly")),length)
lengths=unlist(lengths)

maxDL_monthly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
maxDL_monthly=data.frame(maxDL_monthly)
maxDL_monthly_colnames=c()

for (i in 1:(ncol(maxDL_monthly)-1)){
	names=paste0("maxDL_Month",i)
	maxDL_monthly_colnames=c(maxDL_monthly_colnames,names)
}
colnames(maxDL_monthly)[-c(1)]=maxDL_monthly_colnames
colnames(maxDL_monthly)[1]="Planting"
maxDL_monthly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(maxDL_monthly)){
	maxDL_monthly[i,-c(1)]=eval(as.name(paste0(maxDL_monthly$Planting[i],"_maxDL_monthly")))[1:(ncol(maxDL_monthly)-1)]
}

#max monthly thermal time


df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	maxTT=max(ENV.sub$TT)
	df=c(df,maxTT)
	}
	assign(paste0(ENV$season[1],"_maxTT_monthly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_maxTT_monthly")),length)
lengths=unlist(lengths)

maxTT_monthly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
maxTT_monthly=data.frame(maxTT_monthly)
maxTT_monthly_colnames=c()

for (i in 1:(ncol(maxTT_monthly)-1)){
	names=paste0("maxTT_Month",i)
	maxTT_monthly_colnames=c(maxTT_monthly_colnames,names)
}
colnames(maxTT_monthly)[-c(1)]=maxTT_monthly_colnames
colnames(maxTT_monthly)[1]="Planting"
maxTT_monthly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(maxTT_monthly)){
	maxTT_monthly[i,-c(1)]=eval(as.name(paste0(maxTT_monthly$Planting[i],"_maxTT_monthly")))[1:(ncol(maxTT_monthly)-1)]
}
#max monthly [PAR (Photosynthetically active radiation) * Hrs.Light]; a bit different from PTT which is (temp >3 * hours daylight)
#maxPART (PAR time)





df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	maxPART=max((ENV.sub$PAR)*(ENV.sub$Hrs.light))
	df=c(df,maxPART)
	}
	assign(paste0(ENV$season[1],"_maxPART_monthly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_maxPART_monthly")),length)
lengths=unlist(lengths)

maxPART_monthly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
maxPART_monthly=data.frame(maxPART_monthly)
maxPART_monthly_colnames=c()

for (i in 1:(ncol(maxPART_monthly)-1)){
	names=paste0("maxPART_Month",i)
	maxPART_monthly_colnames=c(maxPART_monthly_colnames,names)
}
colnames(maxPART_monthly)[-c(1)]=maxPART_monthly_colnames
colnames(maxPART_monthly)[1]="Planting"
maxPART_monthly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(maxPART_monthly)){
	maxPART_monthly[i,-c(1)]=eval(as.name(paste0(maxPART_monthly$Planting[i],"_maxPART_monthly")))[1:(ncol(maxPART_monthly)-1)]
}



#max monthly PAR


df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	maxPAR=max((ENV.sub$PAR))
	df=c(df,maxPAR)
	}
	assign(paste0(ENV$season[1],"_maxPAR_monthly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_maxPAR_monthly")),length)
lengths=unlist(lengths)

maxPAR_monthly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
maxPAR_monthly=data.frame(maxPAR_monthly)
maxPAR_monthly_colnames=c()

for (i in 1:(ncol(maxPAR_monthly)-1)){
	names=paste0("maxPAR_Month",i)
	maxPAR_monthly_colnames=c(maxPAR_monthly_colnames,names)
}
colnames(maxPAR_monthly)[-c(1)]=maxPAR_monthly_colnames
colnames(maxPAR_monthly)[1]="Planting"
maxPAR_monthly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(maxPAR_monthly)){
	maxPAR_monthly[i,-c(1)]=eval(as.name(paste0(maxPAR_monthly$Planting[i],"_maxPAR_monthly")))[1:(ncol(maxPAR_monthly)-1)]
}


#max Monthly Thermal Time
df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	maxTT=max((ENV.sub$TT))
	df=c(df,maxTT)
	}
	assign(paste0(ENV$season[1],"_maxTT_monthly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_maxTT_monthly")),length)
lengths=unlist(lengths)

maxTT_monthly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
maxTT_monthly=data.frame(maxTT_monthly)
maxTT_monthly_colnames=c()

for (i in 1:(ncol(maxTT_monthly)-1)){
	names=paste0("maxTT_Month",i)
	maxTT_monthly_colnames=c(maxTT_monthly_colnames,names)
}
colnames(maxTT_monthly)[-c(1)]=maxTT_monthly_colnames
colnames(maxTT_monthly)[1]="Planting"
maxTT_monthly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(maxTT_monthly)){
	maxTT_monthly[i,-c(1)]=eval(as.name(paste0(maxTT_monthly$Planting[i],"_maxTT_monthly")))[1:(ncol(maxTT_monthly)-1)]
}


#max Monthly Phototermal Time
df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	maxPTT=max((ENV.sub$TT))
	df=c(df,maxPTT)
	}
	assign(paste0(ENV$season[1],"_maxPTT_monthly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_maxPTT_monthly")),length)
lengths=unlist(lengths)

maxPTT_monthly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
maxPTT_monthly=data.frame(maxPTT_monthly)
maxPTT_monthly_colnames=c()

for (i in 1:(ncol(maxPTT_monthly)-1)){
	names=paste0("maxPTT_Month",i)
	maxPTT_monthly_colnames=c(maxPTT_monthly_colnames,names)
}
colnames(maxPTT_monthly)[-c(1)]=maxPTT_monthly_colnames
colnames(maxPTT_monthly)[1]="Planting"
maxPTT_monthly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(maxPTT_monthly)){
	maxPTT_monthly[i,-c(1)]=eval(as.name(paste0(maxPTT_monthly$Planting[i],"_maxPTT_monthly")))[1:(ncol(maxPTT_monthly)-1)]
}



microclim_monthly_max=cbind(maxDL_monthly,maxT_monthly,maxPAR_monthly,maxPART_monthly,maxTT_monthly,maxPTT_monthly)
microclim_monthly_max=microclim_monthly_max[,colnames(microclim_monthly_max)!="Planting"]
save(microclim_monthly_max,file="Microclimate_monthly_max.RData")