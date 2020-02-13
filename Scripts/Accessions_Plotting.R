#Plotting where those accessions came from

a=read.csv('1001genomes_accessions.csv',sep=',')




#Top 10

logical=a$tg_ecotypeid==6937 | a$tg_ecotypeid==8254 | a$tg_ecotypeid==6916 | a$tg_ecotypeid==6930 | a$tg_ecotypeid==8323 | a$tg_ecotypeid==7062 | a$tg_ecotypeid==6958 | a$tg_ecotypeid==6943 |a$tg_ecotypeid==7218 | a$tg_ecotypeid==7003
b=a[logical,]


#Bottom 10

logical=a$tg_ecotypeid== 7438 | a$tg_ecotypeid==7328 | a$tg_ecotypeid==6976 | a$tg_ecotypeid==7217 | a$tg_ecotypeid==7319 | a$tg_ecotypeid==6929 | a$tg_ecotypeid== 6976 | a$tg_ecotypeid==8213 | a$tg_ecotypeid== 7239 | a$tg_ecotypeid==6897

c=a[logical,]

#Plotting them

library(raster)
library(maps)
tmax <- raster(paste(getwd(), "/Bioclim/bio5.bil", sep = "")) #max temp of warmest quarter
tmax[tmax>1000]=NA
#Divide by 10 because data has been recorded in decimals
tmax=tmax/10
png("NorwichSpring2007_Accessions_BestWorst.png")
plot(tmax,xlim=c(-30,80),ylim=c(-30,75),legend.args=list(text='Max temp warmest quarter',side=4,font=2,line=2.5,cex=0.8))
map('world',add=TRUE)
points(b$latitude~b$longitude,col="#d81c45",pch=19,cex=1.0) #top 10 points
points(c$latitude~c$longitude,col="#6059c9",pch=19,cex=1.0) #bottom 10
points(1.2974,52.6309,col="#42a650",cex=1.4,pch=15)



##Based on percentage difference


#Top 10

logical=a$tg_ecotypeid==6937 | a$tg_ecotypeid==8254 | a$tg_ecotypeid==6916 | a$tg_ecotypeid==7062 | a$tg_ecotypeid==8323 | a$tg_ecotypeid==6930 | a$tg_ecotypeid==6958 | a$tg_ecotypeid==6943 |a$tg_ecotypeid==7218 | a$tg_ecotypeid==7003
b=a[logical,]


#Bottom 10

logical=a$tg_ecotypeid== 7328 | a$tg_ecotypeid==7217 | a$tg_ecotypeid==6976 | a$tg_ecotypeid==6929 | a$tg_ecotypeid==7319 | a$tg_ecotypeid==8213 | a$tg_ecotypeid== 7239 | a$tg_ecotypeid==6897 | a$tg_ecotypeid== 7239 

c=a[logical,]

tmax <- raster(paste(getwd(), "/Bioclim/bio5.bil", sep = "")) #max temp of warmest quarter
tmax[tmax>1000]=NA
#Divide by 10 because data has been recorded in decimals
tmax=tmax/10
png("NorwichSpring2007_Accessions_BestWorst_percentdiff.png")
plot(tmax,xlim=c(-30,80),ylim=c(-30,75),legend.args=list(text='Max temp warmest quarter',side=4,font=2,line=2.5,cex=0.8))
map('world',add=TRUE)
points(b$latitude~b$longitude,col="#d81c45",pch=19,cex=1.0) #top 10 points
points(c$latitude~c$longitude,col="#6059c9",pch=19,cex=1.0) #bottom 10
points(1.2974,52.6309,col="#42a650",cex=1.4,pch=15)


png("NorwichSpring2007_Accessions_BestWorst.png")
par(mar=c(1,1,1,5))
plot(tmax,xlim=c(-30,80),ylim=c(-40,75),legend.args=list(text='Max temp warmest quarter',side=4,font=2,line=2.5,cex=0.8),xaxt='n',yaxt='n')
map('world',add=TRUE)
points(top10$LAT~top10$LON,col="#d81c45",pch=19,cex=1.0) #top 10 points
points(bottom10$LAT~bottom10$LON,col="#6059c9",pch=19,cex=1.0) #bottom 10
legend("bottomright", legend=c("Top10", "Bottom10", "Norwich"),
       col=c("#d81c45", "#6059c9","#42a650"), pch=c(19,19,15), cex=0.8,
       title="Accessions", text.font=4, bg='lightblue')
points(1.2974,52.6309,col="#42a650",cex=1.4,pch=15)
dev.off()