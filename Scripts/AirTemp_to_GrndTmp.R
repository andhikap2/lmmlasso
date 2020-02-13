#Converting air to ground temperature

# N=1 # calendar day of the year
# delta= 0.409*sin((2*pi)*N/365 - 1.39) #declination angle [37]
# phi = 52.6309 *pi/180  #latitude in radians (Norwich in this example)
# hs= acos(-tan(phi) * tan(delta)) #hour angle at sunrise and sunset [47]
# dr= 1 + 0.033*cos((2*pi)*N/365) #distance factor [53]
# Gsc=0.082 #solar constant 

# Ra=(24*60/pi) *Gsc*dr*((hs*sin(delta)*sin(phi)) + (cos(phi)*cos(delta)*sin(hs)))
# alt= 2 #altitude or elevation in m (2 in this case because dataset measures at 2m)
# Rso= (0.75 + (2 * (10^-5)) * alt)*Ra #daily clear-sky global irradiance

# #Converting air temp to ground temp
# a=0.004099
# c=0.920493
# d=22.466179
# e=-1.861643
# f=1.549941




# Ta= 20+(273.15) #air temp in Kelvin
# t = 0 #time since midnight january 1 in days
# W=Rso

# Tg = (a*W) + (c*Ta) + e*(sin((2*pi*t/365)+f)) + d 
# gtemp=Tg-273.15 #ground temp in Kelvin

setwd("/srv/uom-data1-q.unimelb.edu.au/6300-afournier/data/Clim_GWAS/Clim_GWAS_2/CMIP5")
library(raster)
min2071=stack("tasmin_day_CMCC-CM_rcp85_r1i1p1_20710101-20711231.nc")
min2071points=rasterToPoints(min2071)
dayseq=seq(1:365)
colnames(min2071points)[-2:-1]=dayseq

## Constants


Gsc=0.082
ac=0.004099
c=0.920493
d=22.466179
e=-1.861643
f=1.549941
alt=2


groundtemps=min2071points


# for (i in 1:365){
# 	for (j in 1:nrow(min2071points)){
# 		N=i
# 		delta= 0.409*sin((2*pi)*N/365 - 1.39) #declination angle [37]
# 		phi= (min2071points[j,2]) *pi/180
# 		hs = acos(-tan(phi) * tan(delta))
# 		dr = 1 + 0.033*cos((2*pi)*N/365)
# 		Ra=(24*60/pi) *Gsc*dr*((hs*sin(delta)*sin(phi)) + (cos(phi)*cos(delta)*sin(hs)))
# 		Rso= (0.75 + (2 * (10^-5)) * alt)*Ra #daily clear-sky global irradiance
# 		Ta= min2071points[j,colnames(min2071points)==i]
# 		t=i
# 		W=Rso
# 		Tg = (a*W) + (c*Ta) + e*(sin((2*pi*t/365)+f)) + d 
# 		gtemp = Tg-273.15
# 		groundtemps[j,colnames(groundtemps==i)]<-gtemp
# 	}
# }


delta<-function(day){
	N=day
	delta= 0.409*sin((2*pi)*N/365 - 1.39) #declination angle [37]
	return(delta)
}

deltas=lapply(c(1:365),delta)
deltas=unlist(deltas)

latitudes=min2071points[,2]
phi<-function(latitude){
	phi=latitude*pi/180
	return(phi)
}
phis=lapply(latitudes,phi)
phis=unlist(phis)

# hs_matrix=matrix(NA,nrow=length(phis),ncol=length(deltas))

# for (i in 1:length(phis)){
# 	for (j in 1:length(deltas)){
# 		v=acos(-tan(phis[i])*tan(deltas[j]))
# 		hs_matrix[i,j]=v
# 	}
# }

aa=-tan(phis)
bb=tan(deltas)

aa=matrix(aa,ncol=1)
bb=matrix(bb,nrow=1)

hs_matrix=aa%*%bb
hs_acos=acos(hs_matrix)
# hs1=hs_acos[,1]
# hs1=matrix(hs1,nrow=1)

# N=1
# dr = 1 + 0.033*cos((2*pi)*N/365)
# Ra=(24*60/pi) *Gsc*dr*((hs1*sin(deltas[1])*sin(phis)) + (cos(phis)*cos(deltas[1])*sin(hs1)))
# Rso= (0.75 + (2 * (10^-5)) * alt)*Ra #daily clear-sky global irradiance


newraster=min2071
# Ta= values(newraster[[1]])
# t=N
# W=Rso
# Tg = (ac*W) + (c*Ta) + e*(sin((2*pi*t/365)+f)) + d 
# gtemp = Tg-273.15


# values(newraster[[1]])<-gtemp
# plot(newraster[[i]],zlim=c(0,50))





for (day in 1:365){
	hs_day=hs_acos[,day]
	hs_day=matrix(hs_day,nrow=1)
	N=day
	dr = 1 + 0.033*cos((2*pi)*N/365)
	Ra=(24*60/pi) *Gsc*dr*((hs_day*sin(deltas[day])*sin(phis)) + (cos(phis)*cos(deltas[day])*sin(hs_day)))
	Rso= (0.75 + (2 * (10^-5)) * alt)*Ra #daily clear-sky global irradiance
	Ta= values(min2071[[day]])
	t=N
	W=Rso
	Tg = (ac*W) + (c*Ta) + e*(sin((2*pi*t/365)+f)) + d 
	gtemp = Tg-273.15
	values(newraster[[day]])<-gtemp



}




#Efficiently do this with apply?
#Switch to Python..