#Making interactions for the microclim data

setwd("~/Clim_GWAS/Clim_GWAS_2")
load('full_environ_data_list.Robj') #load environmental data
load('Multitrait_NEW.RData') #Load the phenotype data
Th=as.numeric(args[1]) #induction temperature threshold
N=as.numeric(args[2]) #period in hours
Th=0
N=1

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


groundtemp=c(HalleFall2006_groundtemp_Daily,NorwichFall2006_groundtemp_Daily,NorwichSpring2007_groundtemp_Daily,NorwichSummer2006_groundtemp_Daily,NorwichSummer2007_groundtemp_Daily,
	ValenciaFall2006_groundtemp_Daily)
daylength=c(HalleFall2006_daylength_Daily,NorwichFall2006_daylength_Daily,NorwichSpring2007_daylength_Daily,NorwichSummer2006_daylength_Daily,NorwichSummer2007_daylength_Daily,
	ValenciaFall2006_daylength_Daily)
PAR=c(HalleFall2006_PAR_Daily,NorwichFall2006_PAR_Daily,NorwichSpring2007_PAR_Daily,NorwichSummer2006_PAR_Daily,NorwichSummer2007_PAR_Daily,
	ValenciaFall2006_PAR_Daily)
TT=c(HalleFall2006_TT_Daily,NorwichFall2006_TT_Daily,NorwichSpring2007_TT_Daily,NorwichSummer2006_TT_Daily,NorwichSummer2007_TT_Daily,
	ValenciaFall2006_TT_Daily)
PTT=c(HalleFall2006_PTT_Daily,NorwichFall2006_PTT_Daily,NorwichSpring2007_PTT_Daily,NorwichSummer2006_PTT_Daily,NorwichSummer2007_PTT_Daily,
	ValenciaFall2006_PTT_Daily)

dailymicroclim=rbind(groundtemp,daylength,PAR,TT,PTT)
dailymicroclims=t(dailymicroclim)
corrr=cor(dailymicroclims)

for (ENV in df_list){


assign((paste0(ENV$season[1],"_",b))
	,

	c(
	(eval(as.name((paste0(ENV$season[1],"_groundtemp_",b))))),
	(eval(as.name((paste0(ENV$season[1],"_daylength_",b))))),
	(eval(as.name(paste0(ENV$season[1],"_PAR_",b)))),
	(eval(as.name(paste0(ENV$season[1],"_TT_",b)))),
	(eval(as.name(paste0(ENV$season[1],"_PTT_",b))))
	))

}