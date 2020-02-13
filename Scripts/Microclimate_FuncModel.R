#Functional modelling of microclimate data based on Chew (2012) paper

#Get the conditions for each location
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

for (i in unique(site_ptus$season)){
	assign(paste0(i,"_site_ptus"),site_ptus[site_ptus$season==i,])
}

##################

N=1 #  hourly period
Th=3 # Temperature threshold
hnu=0
CSDL=10 #Critical short daylength, based on Chew(2012) Table S2
CLDL=14 #Critical long daylength, based on Chew(2012) Table S2
DLD=1 #Maximum rate to bolting
DSD=0.6626 #Minimum rate to bolting; based on Chew(2012) Table S2- Wilczek et al model




#Defining the Photoperiod
df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){

	photoperiod=c()
	for (i in 1:floor(nrow(ENV)/N)) {

		if (ENV$Daylength[i]<=CSDL){

			pp=DSD
			} else if (ENV$Daylength[i]>=CLDL){

				pp=DLD
			} else {

				pp=  DSD + ( ((ENV$Daylength[i]-CSDL)*(DLD-DSD)) / (CLDL-CSDL) )
			}	

			photoperiod=c(photoperiod,pp)
			if (length(photoperiod)==nrow(ENV)) {
				ENV=cbind(ENV,photoperiod)
			}

	} 
}

#Adding the photoperiods to each planting
photoperiod_list= mget(ls(pattern="\\_photoperiod")) #List of photoperiod objects


#Making the matrix
lengths=lapply(mget(ls(pattern="_photoperiod")),length)
lengths=unlist(lengths)

photoperiod=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1))
photoperiod=data.frame(photoperiod)
photoperiod_colnames=c()

for (i in 1:(ncol(photoperiod)-1)){
	names=paste0("photoperiod",i)
	photoperiod_colnames=c(photoperiod_colnames,names)
}
colnames(photoperiod)[-c(1)]=photoperiod_colnames
colnames(photoperiod)[1]="Planting"
photoperiod[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(photoperiod)){
	photoperiod[i,-c(1)]=eval(as.name(paste0(photoperiod$Planting[i],"_photoperiod")))[1:(ncol(photoperiod)-1)]
}



#Defining Thermal time
df_list=mget(ls(pattern="\\_site_ptus"))

for (ENV in df_list){
	thermaltime=c()
	for (i in 1:floor(nrow(ENV)/N)) {
			if (ENV$Grnd.Tmp[i]>=Th) {
				tt= (ENV$Grnd.Tmp[i]-Th)*(ENV$Hrs.light[i]) #Multiply by hours light
					} else {
						tt=0
					}
					thermaltime=c(thermaltime,tt)
					assign(paste0(ENV$season[1],"_thermaltime"),thermaltime)
					}}

#Making the matrix
lengths=lapply(mget(ls(pattern="_thermaltime")),length)
lengths=unlist(lengths)

thermaltime=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1))
thermaltime=data.frame(thermaltime)
thermaltime_colnames=c()

for (i in 1:(ncol(thermaltime)-1)){
	names=paste0("thermaltime",i)
	thermaltime_colnames=c(thermaltime_colnames,names)
}
colnames(thermaltime)[-c(1)]=thermaltime_colnames
colnames(thermaltime)[1]="Planting"
thermaltime[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(thermaltime)){
	thermaltime[i,-c(1)]=eval(as.name(paste0(thermaltime$Planting[i],"_thermaltime")))[1:(ncol(thermaltime)-1)]
}

#Determining vernalisation

Tv_min=-3.5 #Minimum vernalisation temperature 
Tv_max=6 #Maximum vernalisation temperature
kappa=-5.17 #some sort of parameter; value from Chew (2012) Table S2
omega=2.23 #some sort of parameter; value from Chew (2012) Table S2
xi=1 #some sort of parameter; value from Chew (2012) Table S2

vernalisation_effectiveness <- function(T) { #T as in Grnd.Tmp
	ve= exp(kappa)*((T - Tv_min)^(omega))*((Tv_max - T)^xi)
	return(ve)
}

#Calculating hourly vernalisation
df_list=mget(ls(pattern="\\_site_ptus"))
vernalisation=c()

for (ENV in df_list){
	for (i in 1:floor(nrow(ENV)/N)) {
			v=vernalisation_effectiveness(ENV$Grnd.Tmp[i])
								
					vernalisation=c(vernalisation,v)
				}
					assign(paste0(ENV$season[1],"_vernalisation"),vernalisation)
					}







#Making the matrix
lengths=lapply(mget(ls(pattern="_vernalisation")),length)
lengths=unlist(lengths)

vernalisation=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1))
vernalisation=data.frame(vernalisation)
vernalisation_colnames=c()

for (i in 1:(ncol(vernalisation)-1)){
	names=paste0("vernalisation",i)
	vernalisation_colnames=c(vernalisation_colnames,names)
}
colnames(vernalisation)[-c(1)]=vernalisation_colnames
colnames(vernalisation)[1]="Planting"
vernalisation[,1]=AllPlantings_Corrected$Planting


for (i in 1:nrow(vernalisation)){
	vernalisation[i,-c(1)]=eval(as.name(paste0(vernalisation$Planting[i],"_vernalisation")))[1:(ncol(vernalisation)-1)]
}

#Should add the vernalisation values to each planting's _site_ptus.




#Calculating cumulative vernalisation
Vh=c()
ENV=NorwichSummer2007_site_ptus 
for (i in 1:(nrow(ENV)/N)) {
	ENV.sub=ENV[(i-1)*N+1:(i*N),]
	muT=mean((ENV.sub$Grnd.Tmp*ENV.sub$Hrs.light)[ENV.sub$Grnd.Tmp>Th])
	cumT=c(cumT,muT+sum(cumT[length(cumT)],na.omit=T))
}
