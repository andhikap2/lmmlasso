#### Making microenvironmental variables


setwd("~/Clim_GWAS/Clim_GWAS_2")
load('full_environ_data_list.Robj') #load environmental data

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


####################################################################################################################################################################################################################











#Separating microclimate data into each planting
for (i in unique(site_ptus$season)){
	assign(paste0(i,"_site_ptus"),site_ptus[site_ptus$season==i,])
}

#############################################


##Getting hourly microclim values for all plantings at once
for (a in unique(site_ptus$season)){

assign(paste0(a,"_microclim"),matrix(NA,nrow=4,ncol=nrow(eval(as.name(paste0(a,"_site_ptus"))))))
y=data.frame(eval(as.name(paste0(a,"_microclim"))))
rownames(y)=c("Grnd.Tmp","cumPTT","Hrs.Light","PAR")

ndays=eval(as.name(paste0(a,"_site_ptus")))$Date[nrow(eval(as.name(paste0(a,"_site_ptus"))))]-eval(as.name(paste0(a,"_site_ptus")))$Date[1]+1

colnameslist=c()


for (i in 1:ndays)
for(j in 1:24){{
column_name=paste0("Day",i,"Hour",j)
colnameslist=c(colnameslist,column_name)

}}


colnames(y)=colnameslist

#Filling up ground temp
ground_temp=eval(as.name(paste0(a,"_site_ptus")))$Grnd.Tmp
y[1,]=ground_temp

#Filling up cumPTT
cumulative_PTT=eval(as.name(paste0(a,"_site_ptus")))$cumPTT
y[2,]=cumulative_PTT

#filling up hours light
hours_light=eval(as.name(paste0(a,"_site_ptus")))$Hrs.light
y[3,]=hours_light

#filling up PAR
par=eval(as.name(paste0(a,"_par")))$PAR
y[4,]=par


assign(paste0(a,"_microclim"),y)

}



## Creating a matrix of individual values

########################## 

##Making a function

#Creating matrices

create_matrix=function(x) {
	a=matrix(NA,nrow=1,ncol=(nrow(eval(as.name(paste0(x,"_microclim"))))*ncol(eval(as.name(paste0(x,"_microclim"))))))
	a=data.frame(a)
	assign(paste0(x,"_microclimvalues"),a,envir=.GlobalEnv)
}

season=data.frame(seasons)
for (i in unique(season)){
results=apply(unique(season),1,create_matrix)
}

#Assigning column names
for (i in unique(seasons)){
	a_colnameslist=c()
for (j in rownames(eval(as.name(paste0(i,"_microclim")))))
for (k in colnames(eval(as.name(paste0(i,"_microclim"))))){{


	a_colname=paste0(j,"_",k)
	a_colnameslist=c(a_colnameslist,a_colname)


	}
	d=get(paste0(i,"_microclimvalues"))
	colnames(d)=a_colnameslist
	assign(paste0(i,"_microclimvalues"),d)
	
	}}




####################

for (i in unique(seasons)) {
	values=c(eval(as.name(paste0(i,"_microclim")))["Grnd.Tmp",],eval(as.name(paste0(i,"_microclim")))["cumPTT",],eval(as.name(paste0(i,"_microclim")))["daylength",])
	b=eval(as.name(paste0(i,"_microclimvalues")))
	b[]=values
	assign(paste0(i,"_microclimvalues"),b)	
}

######################################################################################################################################################################################################################


#### Combining with genetic data


load("Multitrait_NEW.RData")
Multitrait_NEW=data.frame(Multitrait_NEW)

###################################

for (i in unique(seasons)){

empty_matrix=matrix(NA,nrow=nrow(eval(as.name(paste0(i)))),ncol=ncol(eval(as.name(paste0(i,"_microclimvalues")))))
empty_matrix=data.frame(empty_matrix)
colnames(empty_matrix)=colnames(eval(as.name(paste0(i,"_microclimvalues"))))
for (j in 1:ncol(empty_matrix)){
	empty_matrix[,j]=eval(as.name(pacdste0(i,"_microclimvalues")))[,j]

}
assign(paste0(i,"_emptymatrix"),empty_matrix)
}

####################################
for (i in unique(seasons)){

	g=cbind(eval(as.name(paste0(i))),eval(as.name(paste0(i,"_emptymatrix"))))
	assign(paste0(i,"_hourly_env"),g)
}

#Creating matrix of purely microclimatic variables (no SNPs)
library(dplyr)
microclimvalues=bind_rows(mget(ls(pattern="_microclimvalues")))
rownames(microclimvalues)=c("HalleFall2006","NorwichFall2006","NorwichSpring2007","NorwichSummer2006","NorwichSummer2007","OuluFall2007","ValenciaFall2006")
microclimvalues=data.frame(microclimvalues)
microclimvalues=microclimvalues %>%
	mutate_all(as.numeric)
microclimvalues=microclimvalues[,colSums(is.na(microclimvalues))==0] #Remove columns with NA...
#####


microclimate=matrix(NA,nrow=nrow(Multitrait_NEW),ncol=(ncol(microclimvalues)+1))
microclimate=data.frame(microclimate)
microclimate[,1]=(Multitrait_NEW$Planting)
colnames(microclimate)[1]="Planting"

#####

Geno=matrix(NA,nrow=nrow(Multitrait_NEW),ncol=(ncol(microclimvalues)+1))
rownames(Geno)=(Multitrait_NEW$Planting)
colnames(Geno)[2:ncol(Geno)]=colnames(microclimvalues)
colnames(Geno)[1]="Planting"
Geno[,1]=rownames(Geno)
Geno=data.frame(Geno,stringsAsFactors=FALSE)
rownames(Geno)=NULL

for (i in 1:nrow(Geno)){
	planting_site=Geno[i,1]
	Geno[i,-1]=unlist(microclimvalues[(rownames(microclimvalues)==planting_site),])

}

microclimate=Geno
save(microclimate,file="Multitrait_NEW_microclimate.RData")






####################################

N=24 #  hourly period
Th=3 # Temp threshold
hnu=0

cumT=c()
ENV=HalleFall2006_site_ptus
for (i in 1:(nrow(ENV)/N)) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	muT=mean((ENV.sub$Grnd.Tmp*ENV.sub$Hrs.light)[ENV.sub$Grnd.Tmp>Th])
	cumT=c(cumT,muT+sum(cumT[length(cumT)],na.omit=T))
}

plot(ENV$cumPTT)
lines(cumT,col=2)

###### cumT 6-hourly


N=6 #  hourly period
Th=3 # Temp threshold
hnu=0



df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	cumT=c()
	for (i in 1:floor(nrow(ENV)/N)) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	muT=mean((ENV.sub$Grnd.Tmp*ENV.sub$Hrs.light)[ENV.sub$Grnd.Tmp>Th])
	cumT=c(cumT,muT+sum(cumT[length(cumT)],na.omit=T))
}
	assign(paste0(ENV$season[1],"_cumT_6hours"),cumT)
}
#Making the matrix

lengths_cumT_6hours=lapply(mget(ls(pattern="_cumT_6hours")),length)
lengths_cumT_6hours=unlist(lengths_cumT_6hours)

cumT_6hours=matrix(NA,nrow=nrow(Multitrait_NEW),ncol=(min(lengths_cumT_6hours)+1))
cumT_6hours=data.frame(cumT_6hours)
cumT_6hours_colnames=c()

for (i in 1:(ncol(cumT_6hours)-1)){
	names=paste0("cumT_6hours",i)
	cumT_6hours_colnames=c(cumT_6hours_colnames,names)
}
colnames(cumT_6hours)[-c(1)]=cumT_6hours_colnames
colnames(cumT_6hours)[1]="Planting"
cumT_6hours[,1]=Multitrait_NEW$Planting


for (i in 1:nrow(cumT_6hours)){
	cumT_6hours[i,-c(1)]=eval(as.name(paste0(cumT_6hours$Planting[i],"_cumT_6hours")))[1:(ncol(cumT_6hours)-1)]
}

###### muT 6-hourly


N=6 #  hourly period
Th=3 # Temp threshold
hnu=0



df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	muT=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	muT=mean((ENV.sub$Grnd.Tmp*ENV.sub$Hrs.light)[ENV.sub$Grnd.Tmp>Th])
	muT=c(muT,muT+sum(muT[length(muT)],na.omit=T))
}
	assign(paste0(ENV$season[1],"_muT_6hours"),muT)
}

#Making the matrix

lengths_muT_6hours=lapply(mget(ls(pattern="_muT_6hours")),length)
lengths_muT_6hours=unlist(lengths_muT_6hours)

muT_6hours=matrix(NA,nrow=nrow(Multitrait_NEW),ncol=(min(lengths_muT_6hours)+1))
muT_6hours=data.frame(muT_6hours)
muT_6hours_colnames=c()

for (i in 1:(ncol(muT_6hours)-1)){
	names=paste0("muT_6hours",i)
	muT_6hours_colnames=c(muT_6hours_colnames,names)
}
colnames(muT_6hours)[-c(1)]=muT_6hours_colnames
colnames(muT_6hours)[1]="Planting"
muT_6hours[,1]=Multitrait_NEW$Planting


for (i in 1:nrow(muT_6hours)){
	muT_6hours[i,-c(1)]=eval(as.name(paste0(muT_6hours$Planting[i],"_muT_6hours")))[1:(ncol(muT_6hours)-1)]
}

###### cumT daily


N=24 #  hourly period
Th=3 # Temp threshold
hnu=0



df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	cumT=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	muT=mean((ENV.sub$Grnd.Tmp*ENV.sub$Hrs.light)[ENV.sub$Grnd.Tmp>Th])
	cumT=c(cumT,muT+sum(cumT[length(cumT)],na.omit=T))
}
	assign(paste0(ENV$season[1],"_cumT_daily"),cumT)
}
#Making the matrix

lengths_cumT_daily=lapply(mget(ls(pattern="_cumT_daily")),length)
lengths_cumT_daily=unlist(lengths_cumT_daily)

cumT_daily=matrix(NA,nrow=nrow(Multitrait_NEW),ncol=(min(lengths_cumT_daily)+1))
cumT_daily=data.frame(cumT_daily)
cumT_daily_colnames=c()

for (i in 1:(ncol(cumT_daily)-1)){
	names=paste0("cumT_Day",i)
	cumT_daily_colnames=c(cumT_daily_colnames,names)
}
colnames(cumT_daily)[-c(1)]=cumT_daily_colnames
colnames(cumT_daily)[1]="Planting"
cumT_daily[,1]=Multitrait_NEW$Planting


for (i in 1:nrow(cumT_daily)){
	cumT_daily[i,-c(1)]=eval(as.name(paste0(cumT_daily$Planting[i],"_cumT_daily")))[1:(ncol(cumT_daily)-1)]
}







######## cumT weekly 

N=168 #  hourly period
Th=3 # Temp threshold
hnu=0


df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	cumT=c()
	for (i in 1:floor((nrow(ENV)/N))) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	muT=mean((ENV.sub$Grnd.Tmp*ENV.sub$Hrs.light)[ENV.sub$Grnd.Tmp>Th])
	cumT=c(cumT,muT+sum(cumT[length(cumT)],na.omit=T))
}
	assign(paste0(ENV$season[1],"_cumT_weekly"),cumT)}



#Making the matrix
lengths=lapply(mget(ls(pattern="_cumT_weekly")),length)
lengths=unlist(lengths)

cumT_weekly=matrix(NA,nrow=nrow(Multitrait_NEW),ncol=(min(lengths)+1)) ##Should I do min or max
cumT_weekly=data.frame(cumT_weekly)
cumT_weekly_colnames=c()

for (i in 1:(ncol(cumT_weekly)-1)){
	names=paste0("cumT_Week",i)
	cumT_weekly_colnames=c(cumT_weekly_colnames,names)
}
colnames(cumT_weekly)[-c(1)]=cumT_weekly_colnames
colnames(cumT_weekly)[1]="Planting"
cumT_weekly[,1]=Multitrait_NEW$Planting


for (i in 1:nrow(cumT_weekly)){
	cumT_weekly[i,-c(1)]=eval(as.name(paste0(cumT_weekly$Planting[i],"_cumT_weekly")))[1:(ncol(cumT_weekly)-1)]
}





#### Meean number of daylight hours per week?

