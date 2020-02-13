#Korves Air temp to ground temp

#Constants
Gsc=0.082
ac=0.004099
c=0.920493
d=22.466179
e=-1.861643
f=1.549941
alt=2

setwd("~/Clim_GWAS/Daymet/")
korves_env=read.csv("korves_env_2002-2003.csv")

#Fall planting date: October 28
#Spring planting date: April 1

falldays=korves_env[301:(301+202),]
springdays=korves_env[(365+91):(365+91+202),]

rhlat=41.6771 #rhode island latitude
rhlong=71.2662 #rhode island longitude


tmin_fall=falldays[,4]
tmax_fall=falldays[,3]

tmin_spring=springdays[,4]
tmax_spring=springdays[,3]


delta<-function(day){
	N=day
	delta= 0.409*sin((2*pi)*N/365 - 1.39) #declination angle [37]
	return(delta)
}

fall_deltas=lapply(falldays[,2],delta)
fall_deltas=unlist(fall_deltas)

spring_deltas=lapply(springdays[,2],delta)
spring_deltas=unlist(spring_deltas)


phis=rhlat*pi/180


#Spring min airtemp to gtemp

aa=-tan(phis)
bb=tan(spring_deltas)

hs_matrix=aa*bb
hs_acos=acos(hs_matrix)

tmin_spring_gtemp=c()

for (day in 1:length(tmin_spring)){
	hs_day=hs_acos[day]
	hs_day=matrix(hs_day,nrow=1)
	N=springdays[day,2]
	dr = 1 + 0.033*cos((2*pi)*N/365)
	Ra=(24*60/pi) *Gsc*dr*((hs_day*sin(spring_deltas[day])*sin(phis)) + (cos(phis)*cos(spring_deltas[day])*sin(hs_day)))
	Rso= (0.75 + (2 * (10^-5)) * alt)*Ra #daily clear-sky global irradiance
	Ta= tmin_spring[day]+273.15 #air temp
	t=N
	W=Rso
	Tg = (ac*W) + (c*Ta) + e*(sin((2*pi*t/365)+f)) + d #ground temp (K)
	gtemp = Tg-273.15 #ground temp (c)
	tmin_spring_gtemp=c(tmin_spring_gtemp,gtemp)

}

#Spring max airtemp to gtemp

aa=-tan(phis)
bb=tan(spring_deltas)

hs_matrix=aa*bb
hs_acos=acos(hs_matrix)

tmax_spring_gtemp=c()

for (day in 1:length(tmax_spring)){
	hs_day=hs_acos[day]
	hs_day=matrix(hs_day,nrow=1)
	N=springdays[day,2]
	dr = 1 + 0.033*cos((2*pi)*N/365)
	Ra=(24*60/pi) *Gsc*dr*((hs_day*sin(spring_deltas[day])*sin(phis)) + (cos(phis)*cos(spring_deltas[day])*sin(hs_day)))
	Rso= (0.75 + (2 * (10^-5)) * alt)*Ra #daily clear-sky global irradiance
	Ta= tmax_spring[day]+273.15 #air temp
	t=N
	W=Rso
	Tg = (ac*W) + (c*Ta) + e*(sin((2*pi*t/365)+f)) + d #ground temp (K)
	gtemp = Tg-273.15 #ground temp (c)
	tmax_spring_gtemp=c(tmax_spring_gtemp,gtemp)

}


#Fall min airtemp to gtemp

aa=-tan(phis)
bb=tan(fall_deltas)

hs_matrix=aa*bb
hs_acos=acos(hs_matrix)

tmin_fall_gtemp=c()

for (day in 1:length(tmin_fall)){
	hs_day=hs_acos[day]
	hs_day=matrix(hs_day,nrow=1)
	N=falldays[day,2]
	dr = 1 + 0.033*cos((2*pi)*N/365)
	Ra=(24*60/pi) *Gsc*dr*((hs_day*sin(fall_deltas[day])*sin(phis)) + (cos(phis)*cos(fall_deltas[day])*sin(hs_day)))
	Rso= (0.75 + (2 * (10^-5)) * alt)*Ra #daily clear-sky global irradiance
	Ta= tmin_fall[day]+273.15 #air temp
	t=N
	W=Rso
	Tg = (ac*W) + (c*Ta) + e*(sin((2*pi*t/365)+f)) + d #ground temp (K)
	gtemp = Tg-273.15 #ground temp (c)
	tmin_fall_gtemp=c(tmin_fall_gtemp,gtemp)

}


#Fall max airtemp to gtemp

aa=-tan(phis)
bb=tan(fall_deltas)

hs_matrix=aa*bb
hs_acos=acos(hs_matrix)

tmax_fall_gtemp=c()

for (day in 1:length(tmax_fall)){
	hs_day=hs_acos[day]
	hs_day=matrix(hs_day,nrow=1)
	N=falldays[day,2]
	dr = 1 + 0.033*cos((2*pi)*N/365)
	Ra=(24*60/pi) *Gsc*dr*((hs_day*sin(fall_deltas[day])*sin(phis)) + (cos(phis)*cos(fall_deltas[day])*sin(hs_day)))
	Rso= (0.75 + (2 * (10^-5)) * alt)*Ra #daily clear-sky global irradiance
	Ta= tmax_fall[day]+273.15 #air temp
	t=N
	W=Rso
	Tg = (ac*W) + (c*Ta) + e*(sin((2*pi*t/365)+f)) + d #ground temp (K)
	gtemp = Tg-273.15 #ground temp (c)
	tmax_fall_gtemp=c(tmax_fall_gtemp,gtemp)

}


fall_gtemp=c(tmin_fall_gtemp,tmax_fall_gtemp)
spring_gtemp=c(tmin_spring_gtemp,tmax_spring_gtemp)





