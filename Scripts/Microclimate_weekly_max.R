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



######## maxT weekly 

N=168 #  hourly period
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
	assign(paste0(ENV$season[1],"_maxT_weekly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_maxT_weekly")),length)
lengths=unlist(lengths)

maxT_weekly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
maxT_weekly=data.frame(maxT_weekly)
maxT_weekly_colnames=c()

for (i in 1:(ncol(maxT_weekly)-1)){
	names=paste0("maxT_week",i)
	maxT_weekly_colnames=c(maxT_weekly_colnames,names)
}
colnames(maxT_weekly)[-c(1)]=maxT_weekly_colnames
colnames(maxT_weekly)[1]="Planting"
maxT_weekly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(maxT_weekly)){
	maxT_weekly[i,-c(1)]=eval(as.name(paste0(maxT_weekly$Planting[i],"_maxT_weekly")))[1:(ncol(maxT_weekly)-1)]
}


#max weekly daylength


df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	maxDL=max(ENV.sub$Daylength)
	df=c(df,maxDL)
	}
	assign(paste0(ENV$season[1],"_maxDL_weekly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_maxDL_weekly")),length)
lengths=unlist(lengths)

maxDL_weekly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
maxDL_weekly=data.frame(maxDL_weekly)
maxDL_weekly_colnames=c()

for (i in 1:(ncol(maxDL_weekly)-1)){
	names=paste0("maxDL_week",i)
	maxDL_weekly_colnames=c(maxDL_weekly_colnames,names)
}
colnames(maxDL_weekly)[-c(1)]=maxDL_weekly_colnames
colnames(maxDL_weekly)[1]="Planting"
maxDL_weekly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(maxDL_weekly)){
	maxDL_weekly[i,-c(1)]=eval(as.name(paste0(maxDL_weekly$Planting[i],"_maxDL_weekly")))[1:(ncol(maxDL_weekly)-1)]
}

#max weekly thermal time


df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	maxTT=max(ENV.sub$TT)
	df=c(df,maxTT)
	}
	assign(paste0(ENV$season[1],"_maxTT_weekly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_maxTT_weekly")),length)
lengths=unlist(lengths)

maxTT_weekly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
maxTT_weekly=data.frame(maxTT_weekly)
maxTT_weekly_colnames=c()

for (i in 1:(ncol(maxTT_weekly)-1)){
	names=paste0("maxTT_week",i)
	maxTT_weekly_colnames=c(maxTT_weekly_colnames,names)
}
colnames(maxTT_weekly)[-c(1)]=maxTT_weekly_colnames
colnames(maxTT_weekly)[1]="Planting"
maxTT_weekly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(maxTT_weekly)){
	maxTT_weekly[i,-c(1)]=eval(as.name(paste0(maxTT_weekly$Planting[i],"_maxTT_weekly")))[1:(ncol(maxTT_weekly)-1)]
}
#max weekly [PAR (Photosynthetically active radiation) * Hrs.Light]; a bit different from PTT which is (temp >3 * hours daylight)
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
	assign(paste0(ENV$season[1],"_maxPART_weekly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_maxPART_weekly")),length)
lengths=unlist(lengths)

maxPART_weekly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
maxPART_weekly=data.frame(maxPART_weekly)
maxPART_weekly_colnames=c()

for (i in 1:(ncol(maxPART_weekly)-1)){
	names=paste0("maxPART_week",i)
	maxPART_weekly_colnames=c(maxPART_weekly_colnames,names)
}
colnames(maxPART_weekly)[-c(1)]=maxPART_weekly_colnames
colnames(maxPART_weekly)[1]="Planting"
maxPART_weekly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(maxPART_weekly)){
	maxPART_weekly[i,-c(1)]=eval(as.name(paste0(maxPART_weekly$Planting[i],"_maxPART_weekly")))[1:(ncol(maxPART_weekly)-1)]
}



#max weekly PAR


df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	maxPAR=max((ENV.sub$PAR))
	df=c(df,maxPAR)
	}
	assign(paste0(ENV$season[1],"_maxPAR_weekly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_maxPAR_weekly")),length)
lengths=unlist(lengths)

maxPAR_weekly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
maxPAR_weekly=data.frame(maxPAR_weekly)
maxPAR_weekly_colnames=c()

for (i in 1:(ncol(maxPAR_weekly)-1)){
	names=paste0("maxPAR_week",i)
	maxPAR_weekly_colnames=c(maxPAR_weekly_colnames,names)
}
colnames(maxPAR_weekly)[-c(1)]=maxPAR_weekly_colnames
colnames(maxPAR_weekly)[1]="Planting"
maxPAR_weekly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(maxPAR_weekly)){
	maxPAR_weekly[i,-c(1)]=eval(as.name(paste0(maxPAR_weekly$Planting[i],"_maxPAR_weekly")))[1:(ncol(maxPAR_weekly)-1)]
}


#max weekly Thermal Time
df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	maxTT=max((ENV.sub$TT))
	df=c(df,maxTT)
	}
	assign(paste0(ENV$season[1],"_maxTT_weekly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_maxTT_weekly")),length)
lengths=unlist(lengths)

maxTT_weekly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
maxTT_weekly=data.frame(maxTT_weekly)
maxTT_weekly_colnames=c()

for (i in 1:(ncol(maxTT_weekly)-1)){
	names=paste0("maxTT_week",i)
	maxTT_weekly_colnames=c(maxTT_weekly_colnames,names)
}
colnames(maxTT_weekly)[-c(1)]=maxTT_weekly_colnames
colnames(maxTT_weekly)[1]="Planting"
maxTT_weekly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(maxTT_weekly)){
	maxTT_weekly[i,-c(1)]=eval(as.name(paste0(maxTT_weekly$Planting[i],"_maxTT_weekly")))[1:(ncol(maxTT_weekly)-1)]
}


#max weekly Phototermal Time
df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	df=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	ENV.sub=na.omit(ENV.sub)
	maxPTT=max((ENV.sub$TT))
	df=c(df,maxPTT)
	}
	assign(paste0(ENV$season[1],"_maxPTT_weekly"),df)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_maxPTT_weekly")),length)
lengths=unlist(lengths)

maxPTT_weekly=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
maxPTT_weekly=data.frame(maxPTT_weekly)
maxPTT_weekly_colnames=c()

for (i in 1:(ncol(maxPTT_weekly)-1)){
	names=paste0("maxPTT_week",i)
	maxPTT_weekly_colnames=c(maxPTT_weekly_colnames,names)
}
colnames(maxPTT_weekly)[-c(1)]=maxPTT_weekly_colnames
colnames(maxPTT_weekly)[1]="Planting"
maxPTT_weekly[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(maxPTT_weekly)){
	maxPTT_weekly[i,-c(1)]=eval(as.name(paste0(maxPTT_weekly$Planting[i],"_maxPTT_weekly")))[1:(ncol(maxPTT_weekly)-1)]
}



microclim_weekly_max=cbind(maxDL_weekly,maxT_weekly,maxPAR_weekly,maxPART_weekly,maxTT_weekly,maxPTT_weekly)
microclim_weekly_max=microclim_weekly_max[,colnames(microclim_weekly_max)!="Planting"]
save(microclim_weekly_max,file="Microclimate_weekly_max.RData")