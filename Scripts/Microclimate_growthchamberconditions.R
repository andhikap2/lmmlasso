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


norwichspring2007_temps=data.frame(matrix(NA,nrow=29,ncol=24))

for (i in 1:ncol(norwichspring2007_temps)){
  for (j in 1:nrow(norwichspring2007_temps)){
    b=((j-1)*168)+i
tempvalues=temperatures[c(b,b+24,b+48,b+72,b+96,b+120,b+144)] # 7 days of values
meantemp=mean(tempvalues)
norwichspring2007_temps[j,i]<-meantemp
}}

colnameslist=c()
rownamelist=c()

for (i in 1:nrow(norwichspring2007_temps)){
row_name=paste0("Week",i)
rownamelist=c(rownamelist,row_name)
}

for(j in 1:24){
column_name=paste0("Hour",j)
colnameslist=c(colnameslist,column_name)
}

rownames(norwichspring2007_temps)=rownamelist
colnames(norwichspring2007_temps)=colnameslist
write.csv(norwichspring2007_temps,file='NorwichSpring2007_grndtemp.csv',row.names=TRUE)


daylengths=NorwichSpring2007_site_ptus$Daylength


norwichspring2007_daylengths=data.frame(matrix(NA,nrow=29,ncol=24))

for (i in 1:ncol(norwichspring2007_daylengths)){
  for (j in 1:nrow(norwichspring2007_daylengths)){
    b=((j-1)*168)+i
daylvalues=daylengths[c(b,b+24,b+48,b+72,b+96,b+120,b+144)] # 7 days of values
meandayl=mean(daylvalues)
norwichspring2007_daylengths[j,i]<-meandayl
}}
rownames(norwichspring2007_daylengths)=rownamelist
colnames(norwichspring2007_daylengths)=colnameslist
write.csv(norwichspring2007_daylengths,file='NorwichSpring2007_daylength.csv',row.names=TRUE)








temperatures=ValenciaFall2006_site_ptus$Grnd.Tmp


valenciafall2006_temps=data.frame(matrix(NA,nrow=29,ncol=24))

for (i in 1:ncol(valenciafall2006_temps)){
  for (j in 1:nrow(valenciafall2006_temps)){
    b=((j-1)*168)+i
tempvalues=temperatures[c(b,b+24,b+48,b+72,b+96,b+120,b+144)] # 7 days of values
meantemp=mean(tempvalues)
valenciafall2006_temps[j,i]<-meantemp
}}

colnameslist=c()
rownamelist=c()

for (i in 1:nrow(valenciafall2006_temps)){
row_name=paste0("Week",i)
rownamelist=c(rownamelist,row_name)
}

for(j in 1:24){
column_name=paste0("Hour",j)
colnameslist=c(colnameslist,column_name)
}

rownames(valenciafall2006_temps)=rownamelist
colnames(valenciafall2006_temps)=colnameslist
write.csv(valenciafall2006_temps,file='ValenciaFall2006_grndtemp.csv',row.names=TRUE)


daylengths=ValenciaFall2006_site_ptus$Daylength


valenciafall2006_daylengths=data.frame(matrix(NA,nrow=29,ncol=24))

for (i in 1:ncol(valenciafall2006_daylengths)){
  for (j in 1:nrow(valenciafall2006_daylengths)){
    b=((j-1)*168)+i
daylvalues=daylengths[c(b,b+24,b+48,b+72,b+96,b+120,b+144)] # 7 days of values
meandayl=mean(daylvalues)
valenciafall2006_daylengths[j,i]<-meandayl
}}
rownames(valenciafall2006_daylengths)=rownamelist
colnames(valenciafall2006_daylengths)=colnameslist
write.csv(valenciafall2006_daylengths,file='ValenciaFall2006_daylength.csv',row.names=TRUE)

