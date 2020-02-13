#Plotting predictions in R
setwd("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/data/Clim_GWAS/Clim_GWAS_2")
for_plotting=read.csv("DTBMODEL1_wholedataset.csv")
multitrait=read.csv("Multitrait_NEW.csv",sep='\t')
multitrait_filtered=multitrait #no need to filter if using whole dataset

for_plotting$Residual<-for_plotting$Y_hat - for_plotting$Y_test
for_plotting=cbind(multitrait_filtered[,1:2],for_plotting)
for_plotting=for_plotting[,c(1,2,4,5,6,7)]
accessiondata=read.csv("1001genomes_accessions.csv")


latitude=c()
longitude=c()

for (j in (1:nrow(for_plotting))){
	b=for_plotting[j,'Ecotype_id']
	if (b %in% accessiondata$tg_ecotypeid){
		latvalue=accessiondata[accessiondata$tg_ecotypeid==b,'latitude']
	} else {
		latvalue=NA
	}
	latitude=c(latitude,latvalue)

	if (b %in% accessiondata$tg_ecotypeid){
		longvalue=accessiondata[accessiondata$tg_ecotypeid==b,'longitude']
	} else {
		longvalue=NA
	}
	longitude=c(longitude,longvalue)
}

for_plotting=cbind(for_plotting[,c(1,2)],latitude,longitude,for_plotting[,c(3,4,5,6)])


sum(multitrait_filtered$DTB!=for_plotting$Y_test)==0 #should be true



#manually adding coordinates (if coordinates were given as a range, e.g. N49-N50, we took the mid-point value (N49.5))
for_plotting[for_plotting$Ecotype=='Ema-1','latitude']<-51.3
for_plotting[for_plotting$Ecotype=='Ema-1','longitude']<-0.5

for_plotting[for_plotting$Ecotype=='Ksk-1','latitude']<-54.6
for_plotting[for_plotting$Ecotype=='Ksk-1','longitude']<--3.1


for_plotting[for_plotting$Ecotype=='CIBC-5','latitude']<-51.4083
for_plotting[for_plotting$Ecotype=='CIBC-5','longitude']<--0.6383


for_plotting[for_plotting$Ecotype=='Bay-0','latitude']<-49.5
for_plotting[for_plotting$Ecotype=='Bay-0','longitude']<-11



for_plotting[for_plotting$Ecotype=='Bur-0','latitude']<-53.5
for_plotting[for_plotting$Ecotype=='Bur-0','longitude']<-8


for_plotting[for_plotting$Ecotype=='C24','latitude']<-41.25
for_plotting[for_plotting$Ecotype=='C24','longitude']<--8.45


for_plotting[for_plotting$Ecotype=='Ct-1','latitude']<-37.5
for_plotting[for_plotting$Ecotype=='Ct-1','longitude']<-15


for_plotting[for_plotting$Ecotype=='Edi-0','latitude']<-56
for_plotting[for_plotting$Ecotype=='Edi-0','longitude']<-3


for_plotting[for_plotting$Ecotype=='Est-1','latitude']<-58.5
for_plotting[for_plotting$Ecotype=='Est-1','longitude']<-25.5


for_plotting[for_plotting$Ecotype=='Goettingen-7','latitude']<-51.5
for_plotting[for_plotting$Ecotype=='Goettingen-7','longitude']<-10


for_plotting[for_plotting$Ecotype=='KZ1','latitude']<-NA
for_plotting[for_plotting$Ecotype=='KZ1','longitude']<-NA


for_plotting[for_plotting$Ecotype=='Lz-0','latitude']<-46
for_plotting[for_plotting$Ecotype=='Lz-0','longitude']<-3.5


for_plotting[for_plotting$Ecotype=='Mrk-0','latitude']<-49
for_plotting[for_plotting$Ecotype=='Mrk-0','longitude']<-9.5


for_plotting[for_plotting$Ecotype=='Nd-1','latitude']<-50
for_plotting[for_plotting$Ecotype=='Nd-1','longitude']<-10


for_plotting[for_plotting$Ecotype=='Sha','latitude']<-38.35
for_plotting[for_plotting$Ecotype=='Sha','longitude']<-68.48


for_plotting[for_plotting$Ecotype=='Spr1-2','latitude']<-56.31635
for_plotting[for_plotting$Ecotype=='Spr1-2','longitude']<-16.03529


for_plotting[for_plotting$Ecotype=='Van-0','latitude']<-49.5
for_plotting[for_plotting$Ecotype=='Van-0','longitude']<-123


for_plotting[for_plotting$Ecotype=='Wa-1','latitude']<-52.5
for_plotting[for_plotting$Ecotype=='Wa-1','longitude']<-21


for_plotting[for_plotting$Ecotype=='Ws-0','latitude']<-52.5
for_plotting[for_plotting$Ecotype=='Ws-0','longitude']<-30


for_plotting[for_plotting$Ecotype=='Yo-0','latitude']<-37.45
for_plotting[for_plotting$Ecotype=='Yo-0','longitude']<--119.35


for_plotting[for_plotting$Ecotype=='Zdr-6','latitude']<-49.3853
for_plotting[for_plotting$Ecotype=='Zdr-6','longitude']<-16.2544


for_plotting[for_plotting$Ecotype=='Alc-0','latitude']<-40
for_plotting[for_plotting$Ecotype=='Alc-0','longitude']<-3


for_plotting[for_plotting$Ecotype=='Bla-1','latitude']<-41.5
for_plotting[for_plotting$Ecotype=='Bla-1','longitude']<-3


for_plotting[for_plotting$Ecotype=='Blh-2','latitude']<-48
for_plotting[for_plotting$Ecotype=='Blh-2','longitude']<-19


for_plotting[for_plotting$Ecotype=='Cha-0','latitude']<-46
for_plotting[for_plotting$Ecotype=='Cha-0','longitude']<-7


for_plotting[for_plotting$Ecotype=='Cit-0','latitude']<-43.5
for_plotting[for_plotting$Ecotype=='Cit-0','longitude']<-2.5


for_plotting[for_plotting$Ecotype=='Db-0','latitude']<-50.5
for_plotting[for_plotting$Ecotype=='Db-0','longitude']<-8.5 


for_plotting[for_plotting$Ecotype=='Ei-4','latitude']<-50.5
for_plotting[for_plotting$Ecotype=='Ei-4','longitude']<-6.5


for_plotting[for_plotting$Ecotype=='Eil-0','latitude']<-51.5
for_plotting[for_plotting$Ecotype=='Eil-0','longitude']<-12.5


for_plotting[for_plotting$Ecotype=='Ep-0','latitude']<-50.5
for_plotting[for_plotting$Ecotype=='Ep-0','longitude']<-8.5


for_plotting[for_plotting$Ecotype=='Ge-1','latitude']<-46.5
for_plotting[for_plotting$Ecotype=='Ge-1','longitude']<-6


for_plotting[for_plotting$Ecotype=='Go-2','latitude']<-51.5
for_plotting[for_plotting$Ecotype=='Go-2','longitude']<-10


for_plotting[for_plotting$Ecotype=='HOG','latitude']<-NA
for_plotting[for_plotting$Ecotype=='HOG','longitude']<-NA


for_plotting[for_plotting$Ecotype=='Hl-3','latitude']<-51.5
for_plotting[for_plotting$Ecotype=='Hl-3','longitude']<-9.5


for_plotting[for_plotting$Ecotype=='Kl-0','latitude']<-51
for_plotting[for_plotting$Ecotype=='Kl-0','longitude']<-7


for_plotting[for_plotting$Ecotype=='La-1','latitude']<-52.5
for_plotting[for_plotting$Ecotype=='La-1','longitude']<-15.5


for_plotting[for_plotting$Ecotype=='Li-6:1','latitude']<-50.5
for_plotting[for_plotting$Ecotype=='Li-6:1','longitude']<-8


for_plotting[for_plotting$Ecotype=='Litva','latitude']<-NA
for_plotting[for_plotting$Ecotype=='Litva','longitude']<-NA


for_plotting[for_plotting$Ecotype=='Ll-2','latitude']<-42
for_plotting[for_plotting$Ecotype=='Ll-2','longitude']<-3


for_plotting[for_plotting$Ecotype=='Ma-0','latitude']<-50.5
for_plotting[for_plotting$Ecotype=='Ma-0','longitude']<-8.5


for_plotting[for_plotting$Ecotype=='Mc-0','latitude']<-57
for_plotting[for_plotting$Ecotype=='Mc-0','longitude']<-54


for_plotting[for_plotting$Ecotype=='Mh-1','latitude']<-53.5
for_plotting[for_plotting$Ecotype=='Mh-1','longitude']<-20.5


for_plotting[for_plotting$Ecotype=='No-0','latitude']<-51
for_plotting[for_plotting$Ecotype=='No-0','longitude']<-13.5


for_plotting[for_plotting$Ecotype=='Pa-2','latitude']<-38
for_plotting[for_plotting$Ecotype=='Pa-2','longitude']<-13.5


for_plotting[for_plotting$Ecotype=='Pf-0','latitude']<-48.5
for_plotting[for_plotting$Ecotype=='Pf-0','longitude']<-9


for_plotting[for_plotting$Ecotype=='Pla-0','latitude']<-41.5
for_plotting[for_plotting$Ecotype=='Pla-0','longitude']<-2.5 


for_plotting[for_plotting$Ecotype=='Sh-0','latitude']<-51.5
for_plotting[for_plotting$Ecotype=='Sh-0','longitude']<-10.5


for_plotting[for_plotting$Ecotype=='Ty-0','latitude']<-56.4278
for_plotting[for_plotting$Ecotype=='Ty-0','longitude']<--5.23439


for_plotting[for_plotting$Ecotype=='Te-0','latitude']<-63
for_plotting[for_plotting$Ecotype=='Te-0','longitude']<-25.5


for_plotting[for_plotting$Ecotype=='Uk-4','latitude']<-48
for_plotting[for_plotting$Ecotype=='Uk-4','longitude']<-7.5


for_plotting[for_plotting$Ecotype=='Vi-0','latitude']<-50.5
for_plotting[for_plotting$Ecotype=='Vi-0','longitude']<-8.5


for_plotting[for_plotting$Ecotype=='Wei-1','latitude']<-47.25
for_plotting[for_plotting$Ecotype=='Wei-1','longitude']<-8.26


for_plotting[for_plotting$Ecotype=='Wc-2','latitude']<-53
for_plotting[for_plotting$Ecotype=='Wc-2','longitude']<-10


for_plotting[for_plotting$Ecotype=='Ko-2','latitude']<-55.5
for_plotting[for_plotting$Ecotype=='Ko-2','longitude']<-11.5 


for_plotting[for_plotting$Ecotype=='N13','latitude']<-61.36
for_plotting[for_plotting$Ecotype=='N13','longitude']<-34.15


for_plotting[for_plotting$Ecotype=='N7','latitude']<-61.5
for_plotting[for_plotting$Ecotype=='N7','longitude']<-34


for_plotting[for_plotting$Ecotype=='Rld-2','latitude']<-56.25
for_plotting[for_plotting$Ecotype=='Rld-2','longitude']<-34.3167


for_plotting[for_plotting$Ecotype=='PHW-26','latitude']<-50.6728
for_plotting[for_plotting$Ecotype=='PHW-26','longitude']<--3.8404


for_plotting[for_plotting$Ecotype=='PHW-31','latitude']<-51.4666
for_plotting[for_plotting$Ecotype=='PHW-31','longitude']<--3.2


for_plotting[for_plotting$Ecotype=='Mr-0','latitude']<-44.45
for_plotting[for_plotting$Ecotype=='Mr-0','longitude']<-9.5


for_plotting[for_plotting$Ecotype=='Rmx-A02','latitude']<-42.036
for_plotting[for_plotting$Ecotype=='Rmx-A02','longitude']<--86.511


for_plotting[for_plotting$Ecotype=='Pna-10','latitude']<-42.0945
for_plotting[for_plotting$Ecotype=='Pna-10','longitude']<--86.3253


for_plotting[for_plotting$Ecotype=='Pro-0','latitude']<-43.25
for_plotting[for_plotting$Ecotype=='Pro-0','longitude']<--6


for_plotting[for_plotting$Ecotype=='Fei-0','latitude']<-40.92
for_plotting[for_plotting$Ecotype=='Fei-0','longitude']<--8.54 


for_plotting[for_plotting$Ecotype=='Ang-0','latitude']<-50.5
for_plotting[for_plotting$Ecotype=='Ang-0','longitude']<-5.5


for_plotting[for_plotting$Ecotype=='Can-0','latitude']<-28
for_plotting[for_plotting$Ecotype=='Can-0','longitude']<-15.5


for_plotting[for_plotting$Ecotype=='Cen-0','latitude']<-49
for_plotting[for_plotting$Ecotype=='Cen-0','longitude']<-0.5


for_plotting[for_plotting$Ecotype=='Hi-0','latitude']<-54
for_plotting[for_plotting$Ecotype=='Hi-0','longitude']<-34


for_plotting[for_plotting$Ecotype=='Hs-0','latitude']<-52.5
for_plotting[for_plotting$Ecotype=='Hs-0','longitude']<-9.5


for_plotting[for_plotting$Ecotype=='Ka-0','latitude']<-46.5
for_plotting[for_plotting$Ecotype=='Ka-0','longitude']<-14.5


for_plotting[for_plotting$Ecotype=='Lc-0','latitude']<-57.5
for_plotting[for_plotting$Ecotype=='Lc-0','longitude']<-4.5


for_plotting[for_plotting$Ecotype=='Lip-0','latitude']<-50
for_plotting[for_plotting$Ecotype=='Lip-0','longitude']<-19.5 


for_plotting[for_plotting$Ecotype=='Pa-1','latitude']<-38
for_plotting[for_plotting$Ecotype=='Pa-1','longitude']<-13.5


for_plotting[for_plotting$Ecotype=='Stw-0','latitude']<-52.5
for_plotting[for_plotting$Ecotype=='Stw-0','longitude']<-36.5


for_plotting[for_plotting$Ecotype=='Ta-0','latitude']<-49.5
for_plotting[for_plotting$Ecotype=='Ta-0','longitude']<-14.5


for_plotting[for_plotting$Ecotype=='Tu-0','latitude']<-45
for_plotting[for_plotting$Ecotype=='Tu-0','longitude']<-7.5


for_plotting[for_plotting$Ecotype=='Sap-0','latitude']<-49.49
for_plotting[for_plotting$Ecotype=='Sap-0','longitude']<-14.24


for_plotting[for_plotting$Ecotype=='Co-2','latitude']<-40.5
for_plotting[for_plotting$Ecotype=='Co-2','longitude']<--8.5


for_plotting[for_plotting$Ecotype=='Arby-1','latitude']<-56.7
for_plotting[for_plotting$Ecotype=='Arby-1','longitude']<-16.1


for_plotting[for_plotting$Ecotype=='Sap-0','latitude']<-49.49
for_plotting[for_plotting$Ecotype=='Sap-0','longitude']<-14.24


for_plotting[for_plotting$Ecotype=='Spr1-6','latitude']<-56.31635
for_plotting[for_plotting$Ecotype=='Spr1-6','longitude']<-16.03529






for (i in unique(multitrait_filtered$Planting)){
	assign(paste0(i),(for_plotting[for_plotting$Planting==i,]))
}


HalleFall2006$Residual<-scale(HalleFall2006$Residual)
NorwichFall2006$Residual<-scale(NorwichFall2006$Residual)
NorwichSpring2007$Residual<-scale(NorwichSpring2007$Residual)
NorwichSummer2006$Residual<-scale(NorwichSummer2006$Residual)
NorwichSummer2007$Residual<-scale(NorwichSummer2007$Residual)
OuluFall2007$Residual<-scale(OuluFall2007$Residual)
ValenciaFall2006$Residual<-scale(ValenciaFall2006$Residual)


library(raster)
library(maps)

#Extract max temp. data
tmax <- raster(paste(getwd(), "/Bioclim/bio6.bil", sep = ""))
tmax[tmax>1000]=NA
#Divide by 10 because data has been recorded in decimals
tmax=tmax/10


dev.new(width=10,height=10)
par(mfrow=c(1,1))
#Plot temperatures on map
plot(tmax,xlim=c(-25,90),ylim=c(-35,75),main="Halle in Fall 2006",xaxt='n',yaxt='n',legend.args=list(text="Min. temperature coldest quarter",side=4,line=2.5,cex=0.8))
map('world',add=TRUE)
#Plot location of known data points
points(HalleFall2006$latitude~HalleFall2006$longitude,col='black',bg='#db70ff',pch=21,cex=(2/(1+abs(HalleFall2006$Residual)))) #bg is the fill color
points(x=11.9688,y=51.4970,col='black',bg='#ffec70',pch=23,cex=2)#Halle coordinates
legend("bottomright",legend=c("Planting site"),pch=23,col="black",pt.bg="#ffec70")
dev.copy(png,"HalleFall2006DTBmap.png")
dev.off()
dev.off()


#Instead of using the bioclim, should use the data from the netCDF4




dev.new(width=10,height=10)
par(mfrow=c(1,1))
#Plot temperatures on map
plot(tmax,xlim=c(-25,90),ylim=c(-35,75),main="Norwich in Summer 2006",xaxt='n',yaxt='n',legend.args=list(text="Min. temperature coldest quarter",side=4,line=2.5,cex=0.8))
map('world',add=TRUE)
#Plot location of known data points
points(NorwichSummer2006$latitude~NorwichSummer2006$longitude,col='black',bg='#db70ff',pch=21,cex=(2/(1+abs(NorwichSummer2006$Residual)))) #bg is the fill color
points(x=1.2974,y=52.6309,col='black',bg='#ffec70',pch=23,cex=2)#Norwich coordinates
legend("bottomright",legend=c("Planting site"),pch=23,col="black",pt.bg="#ffec70")
dev.copy(png,"NorwichSummer2006DTBmap.png")
dev.off()
dev.off()



dev.new(width=10,height=10)
par(mfrow=c(1,1))
#Plot temperatures on map
plot(tmax,xlim=c(-25,90),ylim=c(-35,75),main="Norwich in Summer 2007",xaxt='n',yaxt='n',legend.args=list(text="Min. temperature coldest quarter",side=4,line=2.5,cex=0.8))
map('world',add=TRUE)
#Plot location of known data points
points(NorwichSummer2007$latitude~NorwichSummer2007$longitude,col='black',bg='#db70ff',pch=21,cex=(2/(1+abs(NorwichSummer2007$Residual)))) #bg is the fill color
points(x=1.2974,y=52.6309,col='black',bg='#ffec70',pch=23,cex=2)#Norwich coordinates
legend("bottomright",legend=c("Planting site"),pch=23,col="black",pt.bg="#ffec70")
dev.copy(png,"NorwichSummer2007DTBmap.png")
dev.off()
dev.off()




dev.new(width=10,height=10)
par(mfrow=c(1,1))
#Plot temperatures on map
plot(tmax,xlim=c(-25,90),ylim=c(-35,75),main="Norwich in Fall 2006",xaxt='n',yaxt='n',legend.args=list(text="Min. temperature coldest quarter",side=4,line=2.5,cex=0.8))
map('world',add=TRUE)
#Plot location of known data points
points(NorwichFall2006$latitude~NorwichFall2006$longitude,col='black',bg='#db70ff',pch=21,cex=(2/(1+abs(NorwichFall2006$Residual)))) #bg is the fill color
points(x=1.2974,y=52.6309,col='black',bg='#ffec70',pch=23,cex=2)#Norwich coordinates
legend("bottomright",legend=c("Planting site"),pch=23,col="black",pt.bg="#ffec70")
dev.copy(png,"NorwichFall2006DTBmap.png")
dev.off()
dev.off()



dev.new(width=10,height=10)
par(mfrow=c(1,1))
#Plot temperatures on map
plot(tmax,xlim=c(-25,90),ylim=c(-35,75),main="Norwich in Spring 2007",xaxt='n',yaxt='n',legend.args=list(text="Min. temperature coldest quarter",side=4,line=2.5,cex=0.8))
map('world',add=TRUE)
#Plot location of known data points
points(NorwichSpring2007$latitude~NorwichSpring2007$longitude,col='black',bg='#db70ff',pch=21,cex=(2/(1+abs(NorwichSpring2007$Residual)))) #bg is the fill color
points(x=1.2974,y=52.6309,col='black',bg='#ffec70',pch=23,cex=2)#Norwich coordinates
legend("bottomright",legend=c("Planting site"),pch=23,col="black",pt.bg="#ffec70")
dev.copy(png,"NorwichSpring2007DTBmap.png")
dev.off()
dev.off()




dev.new(width=10,height=10)
par(mfrow=c(1,1))
#Plot temperatures on map
plot(tmax,xlim=c(-25,90),ylim=c(-35,75),main="Valencia in Fall 2006",xaxt='n',yaxt='n',legend.args=list(text="Min. temperature coldest quarter",side=4,line=2.5,cex=0.8))
map('world',add=TRUE)
#Plot location of known data points
points(ValenciaFall2006$latitude~ValenciaFall2006$longitude,col='black',bg='#db70ff',pch=21,cex=(2/(1+abs(ValenciaFall2006$Residual)))) #bg is the fill color
points(x=-0.3763,y=39.4699,col='black',bg='#ffec70',pch=23,cex=2)#Valencia coordinates
legend("bottomright",legend=c("Planting site"),pch=23,col="black",pt.bg="#ffec70")
dev.copy(png,"ValenciaFall2006DTBmap.png")
dev.off()
dev.off()



dev.new(width=10,height=10)
par(mfrow=c(1,1))
#Plot temperatures on map
plot(tmax,xlim=c(-25,90),ylim=c(-35,75),main="Oulu in Fall 2007",xaxt='n',yaxt='n',legend.args=list(text="Min. temperature coldest quarter",side=4,line=2.5,cex=0.8))
map('world',add=TRUE)
#Plot location of known data points
points(OuluFall2007$latitude~OuluFall2007$longitude,col='black',bg='#db70ff',pch=21,cex=(2/(1+abs(OuluFall2007$Residual)))) #bg is the fill color
points(x=25.4651,y=65.0121,col='black',bg='#ffec70',pch=23,cex=2)#Norwich coordinates
legend("bottomright",legend=c("Planting site"),pch=23,col="black",pt.bg="#ffec70")
dev.copy(png,"OuluFall2007DTBmap.png")
dev.off()
dev.off()



#Calculating physical distance between accession origin site and planting site
library(geosphere)


distancesHalleFall2006=c()
for (i in 1:nrow(HalleFall2006)){
	dist=distHaversine(c(11.9688,51.4970),c(HalleFall2006[i,'longitude'],HalleFall2006[i,'latitude']))
	dist=dist/1000
	distancesHalleFall2006=c(distancesHalleFall2006,dist)
}
HalleFall2006=cbind(HalleFall2006,distancesHalleFall2006)



distancesNorwichSummer2007=c()
for (i in 1:nrow(NorwichSummer2007)){
	dist=distHaversine(c(1.2974,52.6309),c(NorwichSummer2007[i,'longitude'],NorwichSummer2007[i,'latitude']))
	dist=dist/1000
	distancesNorwichSummer2007=c(distancesNorwichSummer2007,dist)
}
NorwichSummer2007=cbind(NorwichSummer2007,distancesNorwichSummer2007)



distancesNorwichFall2006=c()
for (i in 1:nrow(NorwichFall2006)){
	dist=distHaversine(c(1.2974,52.6309),c(NorwichFall2006[i,'longitude'],NorwichFall2006[i,'latitude']))
	dist=dist/1000
	distancesNorwichFall2006=c(distancesNorwichFall2006,dist)
}
NorwichFall2006=cbind(NorwichFall2006,distancesNorwichFall2006)


distancesNorwichSpring2007=c()
for (i in 1:nrow(NorwichSpring2007)){
	dist=distHaversine(c(1.2974,52.6309),c(NorwichSpring2007[i,'longitude'],NorwichSpring2007[i,'latitude']))
	dist=dist/1000
	distancesNorwichSpring2007=c(distancesNorwichSpring2007,dist)
}
NorwichSpring2007=cbind(NorwichSpring2007,distancesNorwichSpring2007)


distancesOuluFall2007=c()
for (i in 1:nrow(OuluFall2007)){
	dist=distHaversine(c(25.4651,65.0121),c(OuluFall2007[i,'longitude'],OuluFall2007[i,'latitude']))
	dist=dist/1000
	distancesOuluFall2007=c(distancesOuluFall2007,dist)
}
OuluFall2007=cbind(OuluFall2007,distancesOuluFall2007)


distancesNorwichSummer2006=c()
for (i in 1:nrow(NorwichSummer2006)){
	dist=distHaversine(c(1.2974,52.6309),c(NorwichSummer2006[i,'longitude'],NorwichSummer2006[i,'latitude']))
	dist=dist/1000
	distancesNorwichSummer2006=c(distancesNorwichSummer2006,dist)
}
NorwichSummer2006=cbind(NorwichSummer2006,distancesNorwichSummer2006)

distancesValenciaFall2006=c()
for (i in 1:nrow(ValenciaFall2006)){
	dist=distHaversine(c(0.3763,39.4699),c(ValenciaFall2006[i,'longitude'],ValenciaFall2006[i,'latitude']))
	dist=dist/1000
	distancesValenciaFall2006=c(distancesValenciaFall2006,dist)
}
ValenciaFall2006=cbind(ValenciaFall2006,distancesValenciaFall2006)




#Calculating difference in min temp. coldest quarter in bioclim at planting site and original accession sites (because the other model suggested these were enough)


tempdiffHalleFall2006=c()
for (i in 1:nrow(HalleFall2006)){
	tempsite=extract(tmax,cbind(11.9688,51.4970)) #temperature at planting site
	tempacc=extract(tmax,cbind(HalleFall2006[i,'longitude'],HalleFall2006[i,'latitude']))
	diff=tempacc-tempsite
	tempdiffHalleFall2006=c(tempdiffHalleFall2006,diff)
}
HalleFall2006=cbind(HalleFall2006,tempdiffHalleFall2006)




HalleFall2006=HalleFall2006[order(HalleFall2006$tempdiffHalleFall2006),]
plot(HalleFall2006$tempdiffHalleFall2006,abs(HalleFall2006$Residual),type='l',pch=19,xlab="Temperature difference (deg C)", ylab="|Standardized residual|",main="Halle in Fall 2006")
points(HalleFall2006$tempdiffHalleFall2006,abs(HalleFall2006$Residual),type='p',pch=19,col='black')






tempdiffNorwichSummer2007=c()
for (i in 1:nrow(NorwichSummer2007)){
	tempsite=extract(tmax,cbind(1.2974,52.6309))
	tempacc=extract(tmax,cbind(NorwichSummer2007[i,'longitude'],NorwichSummer2007[i,'latitude']))
	diff=tempacc-tempsite
	tempdiffNorwichSummer2007=c(tempdiffNorwichSummer2007,diff)
}
NorwichSummer2007=cbind(NorwichSummer2007,tempdiffNorwichSummer2007)


NorwichSummer2007=NorwichSummer2007[order(NorwichSummer2007$tempdiffNorwichSummer2007),]
plot(NorwichSummer2007$tempdiffNorwichSummer2007,abs(NorwichSummer2007$Residual),type='l',pch=19,xlab="Temperature difference (deg C)", ylab="|Standardized residual|",main="Norwich in Summer 2007")
points(NorwichSummer2007$tempdiffNorwichSummer2007,abs(NorwichSummer2007$Residual),type='p',pch=19,col='black')







tempdiffNorwichFall2006=c()
for (i in 1:nrow(NorwichFall2006)){
	dist=distHaversine(c(1.2974,52.6309),c(NorwichFall2006[i,'longitude'],NorwichFall2006[i,'latitude']))
	dist=dist/1000
	tempdiffNorwichFall2006=c(tempdiffNorwichFall2006,dist)
}
NorwichFall2006=cbind(NorwichFall2006,tempdiffNorwichFall2006)


tempdiffNorwichSpring2007=c()
for (i in 1:nrow(NorwichSpring2007)){
	dist=distHaversine(c(1.2974,52.6309),c(NorwichSpring2007[i,'longitude'],NorwichSpring2007[i,'latitude']))
	dist=dist/1000
	tempdiffNorwichSpring2007=c(tempdiffNorwichSpring2007,dist)
}
NorwichSpring2007=cbind(NorwichSpring2007,tempdiffNorwichSpring2007)


tempdiffOuluFall2007=c()
for (i in 1:nrow(OuluFall2007)){
	dist=distHaversine(c(25.4651,65.0121),c(OuluFall2007[i,'longitude'],OuluFall2007[i,'latitude']))
	dist=dist/1000
	tempdiffOuluFall2007=c(tempdiffOuluFall2007,dist)
}
OuluFall2007=cbind(OuluFall2007,tempdiffOuluFall2007)


tempdiffNorwichSummer2006=c()
for (i in 1:nrow(NorwichSummer2006)){
	dist=distHaversine(c(1.2974,52.6309),c(NorwichSummer2006[i,'longitude'],NorwichSummer2006[i,'latitude']))
	dist=dist/1000
	tempdiffNorwichSummer2006=c(tempdiffNorwichSummer2006,dist)
}
NorwichSummer2006=cbind(NorwichSummer2006,tempdiffNorwichSummer2006)

tempdiffValenciaFall2006=c()
for (i in 1:nrow(ValenciaFall2006)){
	dist=distHaversine(c(0.3763,39.4699),c(ValenciaFall2006[i,'longitude'],ValenciaFall2006[i,'latitude']))
	dist=dist/1000
	tempdiffValenciaFall2006=c(tempdiffValenciaFall2006,dist)
}
ValenciaFall2006=cbind(ValenciaFall2006,tempdiffValenciaFall2006)







#Plotting standardized residuals against distance

HalleFall2006=HalleFall2006[order(HalleFall2006$distancesHalleFall2006),]
plot(HalleFall2006$distancesHalleFall2006,abs(HalleFall2006$Residual),type='p',pch=19,xlab="Distance of accession origin from planting site (km)", ylab="|Standardized residual|",main="Halle in Fall 2006")
#points(HalleFall2006$distancesHalleFall2006,abs(HalleFall2006$Residual),type='p',pch=19,col='black')

dev.copy(png,"HalleFall2006MODEL1_residvdist.png")
dev.off()
dev.off()


OuluFall2007=OuluFall2007[order(OuluFall2007$distancesOuluFall2007),]
plot(OuluFall2007$distancesOuluFall2007,abs(OuluFall2007$Residual),type='p',pch=19,xlab="Distance of accession origin from planting site (km)", ylab="|Standardized residual|",main="Oulu in Fall 2007")
#points(OuluFall2007$distancesOuluFall2007,abs(OuluFall2007$Residual),type='p',pch=19,col='black')

dev.copy(png,"OuluFall2007MODEL1_residvdist.png")
dev.off()
dev.off()


ValenciaFall2006=ValenciaFall2006[order(ValenciaFall2006$distancesValenciaFall2006),]
plot(ValenciaFall2006$distancesValenciaFall2006,abs(ValenciaFall2006$Residual),type='p',pch=19,xlab="Distance of accession origin from planting site (km)", ylab="|Standardized residual|",main="Valencia in Fall 2006")
#points(ValenciaFall2006$distancesValenciaFall2006,abs(ValenciaFall2006$Residual),type='p',pch=19,col='black')

dev.copy(png,"ValenciaFall2006MODEL1_residvdist.png")
dev.off()
dev.off()


NorwichFall2006=NorwichFall2006[order(NorwichFall2006$distancesNorwichFall2006),]
plot(NorwichFall2006$distancesNorwichFall2006,abs(NorwichFall2006$Residual),type='p',pch=19,xlab="Distance of accession origin from planting site (km)", ylab="|Standardized residual|",main="Norwich in Fall 2006")
#points(NorwichFall2006$distancesNorwichFall2006,abs(NorwichFall2006$Residual),type='p',pch=19,col='black')

dev.copy(png,"NorwichFall2006MODEL1_residvdist.png")
dev.off()
dev.off()



NorwichSpring2007=NorwichSpring2007[order(NorwichSpring2007$distancesNorwichSpring2007),]
plot(NorwichSpring2007$distancesNorwichSpring2007,abs(NorwichSpring2007$Residual),type='p',pch=19,xlab="Distance of accession origin from planting site (km)", ylab="|Standardized residual|",main="Norwich in Spring 2007")
#points(NorwichSpring2007$distancesNorwichSpring2007,abs(NorwichSpring2007$Residual),type='p',pch=19,col='black')

dev.copy(png,"NorwichSpring2007MODEL1_residvdist.png")
dev.off()
dev.off()


NorwichSummer2006=NorwichSummer2006[order(NorwichSummer2006$distancesNorwichSummer2006),]
plot(NorwichSummer2006$distancesNorwichSummer2006,abs(NorwichSummer2006$Residual),type='p',pch=19,xlab="Distance of accession origin from planting site (km)", ylab="|Standardized residual|",main="Norwich in Summer 2006")
#points(NorwichSummer2006$distancesNorwichSummer2006,abs(NorwichSummer2006$Residual),type='p',pch=19,col='black')

dev.copy(png,"NorwichSummer2006MODEL1_residvdist.png")
dev.off()
dev.off()



NorwichSummer2007=NorwichSummer2007[order(NorwichSummer2007$distancesNorwichSummer2007),]
plot(NorwichSummer2007$distancesNorwichSummer2007,abs(NorwichSummer2007$Residual),type='p',pch=19,xlab="Distance of accession origin from planting site (km)", ylab="|Standardized residual|",main="Norwich in Summer 2007")
#points(NorwichSummer2007$distancesNorwichSummer2007,abs(NorwichSummer2007$Residual),type='p',pch=19,col='black')

dev.copy(png,"NorwichSummer2007MODEL1_residvdist.png")
dev.off()
dev.off()
