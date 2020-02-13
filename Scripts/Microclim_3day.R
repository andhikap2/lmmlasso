#Making microclim in 3-day intervals

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


#Separating microclimate data into each planting
for (i in unique(site_ptus$season)){
  assign(paste0(i,"_site_ptus"),site_ptus[site_ptus$season==i,])
}


temperatures=NorwichSpring2007_site_ptus$Grnd.Tmp
temp_array<-matrix(temperatures,nrow=251,byrow=T)


norwichspring2007_temps=data.frame(matrix(NA,nrow=83,ncol=2))

for (i in 1:nrow(norwichspring2007_temps)){

  b=((i-1)*72)+1
  mintemp=min(temperatures[(b):(b+71)]) #the minimum temperature in a 3-day period
  maxtemp=max(temperatures[(b):(b+71)]) #the maximum temperature in a 3-day period
  norwichspring2007_temps[i,1]<-mintemp
  norwichspring2007_temps[i,2]<-maxtemp

}

norwichspring2007_when=data.frame(matrix(NA,nrow=83,ncol=2))

for (i in 1:nrow(norwichspring2007_when)){
  b=((i-1)*72)+1
  when_min=(which.min(temperatures[(b):(b+71)]))%%24
  when_max=(which.max(temperatures[(b):(b+71)]))%%24
  norwichspring2007_when[i,1]<-when_min
  norwichspring2007_when[i,2]<-when_max
}

norwichspring2007_when[norwichspring2007_when==0]<-24





daylengths=NorwichSpring2007_site_ptus$Daylength
daylength_array<-matrix(daylengths,nrow=251,byrow=T)

norwichspring2007_dls=data.frame(matrix(NA,nrow=83,ncol=1))

for (i in 1:nrow(norwichspring2007_dls)){
  b=((i-1)*72)+1
  meandayl=mean(daylengths[(b):(b+71)])
  norwichspring2007_dls[i,1]<-meandayl
}

##Interpolating the other hours

Hn=7 #time of dawn

intertemps=data.frame(matrix(NA,nrow=83,ncol=24))

for (i in 1:nrow(intertemps)){
  mintemp=norwichspring2007_temps[i,1]
  maxtemp=norwichspring2007_temps[i,2]
  when_min=norwichspring2007_when[i,1]
  when_max=norwichspring2007_when[i,2]
  intertemps[i,when_min]<-mintemp
  intertemps[i,when_max]<-maxtemp
}

for (i in 1:nrow(intertemps)){
  for (t in 1:ncol(intertemps)){
    Hs=Hn + norwichspring2007_dls[i,1] #time of dusk
    Hm = norwichspring2007_when[i,2] #time of maximum temperature
    Tn = norwichspring2007_temps[i,1] #the day's minimum temperature
    Tm = norwichspring2007_temps[i,2] #the day's maximum temperature
    Tp = norwichspring2007_temps[(i+1),1] #the next day's minimum temperature
    Hp = 31 #7+24
    s= 0.227538
    Ts = Tm - s*(Tm - Tp) #temperature at sunset (dusk)
    c= (Tm + Tn) * 0.5 #mean of Tm and Tn
    a = Tm - Tn 
    k = Tm - Ts 
    j = 1 + Hs - Hm 
    L = j - (t - Hm)
    b= (Tp - Ts) / (sqrt(Hp - Hs))
    if (!is.na(intertemps[i,t])){
      next
    }
    if (Hn <= t & t <= Hm){
      T= c +(a/2)*(cos(pi +(((t-Hn)/(Hm-Hn)*pi))))
      intertemps[i,t]<-T
    }
    if (Hm <= t & t <= Hs){
      T= Ts + k*log(L,base=j)
      intertemps[i,t]<-T
    }
    if (Hs <= (t +24) & (t+24) <= Hp){
      T=Ts + b * sqrt((t+24)-Hs)
      intertemps[i,t]<-T
    }
    if (Hs<= t & t <= Hp){
      T = Ts + b * sqrt((t)-Hs)
      intertemps[i,t]<-T
    }
  }
}

for (i in 1:nrow(intertemps)){
  mintemp=norwichspring2007_temps[i,1]
  maxtemp=norwichspring2007_temps[i,2]
  when_min=norwichspring2007_when[i,1]
  when_max=norwichspring2007_when[i,2]
  intertemps[i,when_min]<-mintemp
  intertemps[i,when_max]<-maxtemp
}


############################################



temperatures=ValenciaFall2006_site_ptus$Grnd.Tmp
temp_array<-matrix(temperatures,nrow=203,byrow=TRUE)


valenciafall2006_temps=data.frame(matrix(NA,nrow=68,ncol=2))

for (i in 1:nrow(valenciafall2006_temps)){

  b=((i-1)*72)+1
  mintemp=min(temperatures[(b):(b+71)]) #the minimum temperature in a 3-day period
  maxtemp=max(temperatures[(b):(b+71)]) #the maximum temperature in a 3-day period
  valenciafall2006_temps[i,1]<-mintemp
  valenciafall2006_temps[i,2]<-maxtemp

}

valenciafall2006_when=data.frame(matrix(NA,nrow=68,ncol=2))

for (i in 1:nrow(valenciafall2006_when)){
  b=((i-1)*72)+1
  when_min=(which.min(temperatures[(b):(b+71)]))%%24
  when_max=(which.max(temperatures[(b):(b+71)]))%%24
  valenciafall2006_when[i,1]<-when_min
  valenciafall2006_when[i,2]<-when_max
}

valenciafall2006_when[valenciafall2006_when==0]<-24





daylengths=ValenciaFall2006_site_ptus$Daylength
daylength_array<-matrix(daylengths,nrow=251,byrow=T)

valenciafall2006_dls=data.frame(matrix(NA,nrow=68,ncol=1))

for (i in 1:nrow(valenciafall2006_dls)){
  b=((i-1)*72)+1
  meandayl=mean(daylengths[(b):(b+71)])
  valenciafall2006_dls[i,1]<-meandayl
}

##Interpolating the other hours

Hn=7 #time of dawn

intertemps=data.frame(matrix(NA,nrow=68,ncol=24))

for (i in 1:nrow(intertemps)){
  mintemp=valenciafall2006_temps[i,1]
  maxtemp=valenciafall2006_temps[i,2]
  when_min=valenciafall2006_when[i,1]
  when_max=valenciafall2006_when[i,2]
  intertemps[i,when_min]<-mintemp
  intertemps[i,when_max]<-maxtemp
}

for (i in 1:nrow(intertemps)){
  for (t in 1:ncol(intertemps)){
    Hs=Hn + valenciafall2006_dls[i,1] #time of dusk
    Hm = valenciafall2006_when[i,2] #time of maximum temperature
    Tn = valenciafall2006_temps[i,1] #the day's minimum temperature
    Tm = valenciafall2006_temps[i,2] #the day's maximum temperature
    Tp = valenciafall2006_temps[(i+1),1] #the next day's minimum temperature
    Hp = 31 #7+24
    s= 0.227538
    Ts = Tm - s*(Tm - Tp) #temperature at sunset (dusk)
    c= (Tm + Tn) * 0.5 #mean of Tm and Tn
    a = Tm - Tn 
    k = Tm - Ts 
    j = 1 + Hs - Hm 
    L = j - (t - Hm)
    b= (Tp - Ts) / (sqrt(Hp - Hs))
    if (!is.na(intertemps[i,t])){
      next
    }
    if (Hn <= t & t <= Hm){
      T= c +(a/2)*(cos(pi +(((t-Hn)/(Hm-Hn)*pi))))
      intertemps[i,t]<-T
    }
    if (Hm <= t & t <= Hs){
      T= Ts + k*log(L,base=j)
      intertemps[i,t]<-T
    }
    if (Hs <= (t +24) & (t+24) <= Hp){
      T=Ts + b * sqrt((t+24)-Hs)
      intertemps[i,t]<-T
    }
    if (Hs<= t & t <= Hp){
      T = Ts + b * sqrt((t)-Hs)
      intertemps[i,t]<-T
    }
  }
}

for (i in 1:nrow(intertemps)){
  mintemp=valenciafall2006_temps[i,1]
  maxtemp=valenciafall2006_temps[i,2]
  when_min=valenciafall2006_when[i,1]
  when_max=valenciafall2006_when[i,2]
  intertemps[i,when_min]<-mintemp
  intertemps[i,when_max]<-maxtemp
}