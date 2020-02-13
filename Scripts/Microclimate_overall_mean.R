#Creating single arithmetic average of each microclim variable 

#Mean of Grnd Tmp, Hrs. Light, PAR, TT, etc


df_list = mget(ls(pattern = "\\_site_ptus"))
for (ENV in df_list){
	muT=mean(ENV[,"Grnd.Tmp"])
	muLight=mean(ENV[,"Hrs.light"])
	muDaylength=mean(ENV[,"Daylength"])
	muPAR=mean(ENV[,"PAR"])
	muTT=mean(ENV[,"TT"])

	results=cbind(muT,muLight,muDaylength,muPAR,muTT)
	assign(paste0(ENV$season[1],"grand_mean"),results)
}

#Making the matrix
lengths=lapply(mget(ls(pattern="grand_mean")),length)
lengths=unlist(lengths)

grand_mean=matrix(NA,nrow=nrow(AllPlantings_Corrected),ncol=(min(lengths)+1)) ##Should I do min or max
grand_mean=data.frame(grand_mean)
grand_mean_colnames=c("muT","muLight","muDayLength","muPAR","muTT")
colnames(grand_mean)[-c(1)]=grand_mean_colnames
colnames(grand_mean)[1]="Planting"
grand_mean[,1]=AllPlantings_Corrected$Planting



for (i in 1:nrow(grand_mean)){
	grand_mean[i,-c(1)]=eval(as.name(paste0(grand_mean$Planting[i],"grand_mean")))[1:(ncol(grand_mean)-1)]
}
overall_mean=grand_mean
save(overall_mean,file="Microclimate_overall_mean")
