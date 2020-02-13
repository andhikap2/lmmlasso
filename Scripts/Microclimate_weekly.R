#Microclimate variables weekly

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
  
  env$TT = with(env,pmu(Grnd.Tmp-3,0)) #thermal time
  env$PTT = with(env,pmu(Grnd.Tmp - 3,0)*Hrs.light) #photothermal time
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



N=168 #  hourly period
Th=3 # Temp threshold
hnu=0


######## muT weekly 

df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	muT=mean(ENV.sub$Grnd.Tmp)
	df=c(df,muT)
	}
	assign(paste0(ENV$season[1],"_muT_weekly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_muT_weekly")),length)
lengths=unlist(lengths)

muT_weekly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or mu
muT_weekly=data.frame(muT_weekly)
muT_weekly_colnames=c()

for (i in 1:(ncol(muT_weekly)-1)){
	names=paste0("muT_Week",i)
	muT_weekly_colnames=c(muT_weekly_colnames,names)
}
colnames(muT_weekly)[-c(1)]=muT_weekly_colnames
colnames(muT_weekly)[1]="Planting"
muT_weekly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(muT_weekly)){
	muT_weekly[i,-c(1)]=eval(as.name(paste0(muT_weekly$Planting[i],"_muT_weekly")))[1:(ncol(muT_weekly)-1)]
}


#Mean weekly daylength

N=168 #  hourly period
Th=NA # Just looking at daylength ~
hnu=0

df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	muDL=mean(ENV.sub$Daylength)
	df=c(df,muDL)
	}
	assign(paste0(ENV$season[1],"_muDL_weekly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_muDL_weekly")),length)
lengths=unlist(lengths)

muDL_weekly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or mu
muDL_weekly=data.frame(muDL_weekly)
muDL_weekly_colnames=c()

for (i in 1:(ncol(muDL_weekly)-1)){
	names=paste0("muDL_Week",i)
	muDL_weekly_colnames=c(muDL_weekly_colnames,names)
}
colnames(muDL_weekly)[-c(1)]=muDL_weekly_colnames
colnames(muDL_weekly)[1]="Planting"
muDL_weekly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(muDL_weekly)){
	muDL_weekly[i,-c(1)]=eval(as.name(paste0(muDL_weekly$Planting[i],"_muDL_weekly")))[1:(ncol(muDL_weekly)-1)]
}


#Mean weekly thermal time


N=168 #  hourly period
Th=3 
hnu=0

df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	muTT=mean(ENV.sub$TT)
	df=c(df,muTT)
	}
	assign(paste0(ENV$season[1],"_muTT_weekly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_muTT_weekly")),length)
lengths=unlist(lengths)

muTT_weekly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or mu
muTT_weekly=data.frame(muTT_weekly)
muTT_weekly_colnames=c()

for (i in 1:(ncol(muTT_weekly)-1)){
	names=paste0("muTT_Week",i)
	muTT_weekly_colnames=c(muTT_weekly_colnames,names)
}
colnames(muTT_weekly)[-c(1)]=muTT_weekly_colnames
colnames(muTT_weekly)[1]="Planting"
muTT_weekly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(muTT_weekly)){
	muTT_weekly[i,-c(1)]=eval(as.name(paste0(muTT_weekly$Planting[i],"_muTT_weekly")))[1:(ncol(muTT_weekly)-1)]
}


#Mean weekly [PAR (Photosynthetically active radiation) * Hrs.Light]; a bit different from PTT which is (temp >3 * hours daylight)
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
	assign(paste0(ENV$season[1],"_muPART_weekly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_muPART_weekly")),length)
lengths=unlist(lengths)

muPART_weekly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or mu
muPART_weekly=data.frame(muPART_weekly)
muPART_weekly_colnames=c()

for (i in 1:(ncol(muPART_weekly)-1)){
	names=paste0("muPART_Week",i)
	muPART_weekly_colnames=c(muPART_weekly_colnames,names)
}
colnames(muPART_weekly)[-c(1)]=muPART_weekly_colnames
colnames(muPART_weekly)[1]="Planting"
muPART_weekly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(muPART_weekly)){
	muPART_weekly[i,-c(1)]=eval(as.name(paste0(muPART_weekly$Planting[i],"_muPART_weekly")))[1:(ncol(muPART_weekly)-1)]
}


#mu weekly Phototermal Time
df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
df=c()
for (i in 1:floor((nrow(ENV)/N))) {
ENV.sub=ENV[(i-1)*N+1:(i*N),]
ENV.sub=na.omit(ENV.sub)
muPTT=mean((ENV.sub$TT))
df=c(df,muPTT)
}
assign(paste0(ENV$season[1],"_muPTT_weekly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_muPTT_weekly")),length)
lengths=unlist(lengths)

muPTT_weekly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or mu
muPTT_weekly=data.frame(muPTT_weekly)
muPTT_weekly_colnames=c()

for (i in 1:(ncol(muPTT_weekly)-1)){
names=paste0("muPTT_week",i)
muPTT_weekly_colnames=c(muPTT_weekly_colnames,names)
}
colnames(muPTT_weekly)[-c(1)]=muPTT_weekly_colnames
colnames(muPTT_weekly)[1]="Planting"
muPTT_weekly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(muPTT_weekly)){
muPTT_weekly[i,-c(1)]=eval(as.name(paste0(muPTT_weekly$Planting[i],"_muPTT_weekly")))[1:(ncol(muPTT_weekly)-1)]
}



#mu weekly PAR


df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	muPAR=mean((ENV.sub$PAR))
	df=c(df,muPAR)
	}
	assign(paste0(ENV$season[1],"_muPAR_weekly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_muPAR_weekly")),length)
lengths=unlist(lengths)

muPAR_weekly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or mu
muPAR_weekly=data.frame(muPAR_weekly)
muPAR_weekly_colnames=c()

for (i in 1:(ncol(muPAR_weekly)-1)){
	names=paste0("muPAR_week",i)
	muPAR_weekly_colnames=c(muPAR_weekly_colnames,names)
}
colnames(muPAR_weekly)[-c(1)]=muPAR_weekly_colnames
colnames(muPAR_weekly)[1]="Planting"
muPAR_weekly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(muPAR_weekly)){
	muPAR_weekly[i,-c(1)]=eval(as.name(paste0(muPAR_weekly$Planting[i],"_muPAR_weekly")))[1:(ncol(muPAR_weekly)-1)]
}

microclim_weekly=cbind(muDL_weekly,muT_weekly,muPAR_weekly,muPART_weekly,muTT_weekly,muPTT_weekly)
microclim_weekly=microclim_weekly[,colnames(microclim_weekly)!="Planting"]
save(microclim_weekly,file="Microclimate_weekly.RData")