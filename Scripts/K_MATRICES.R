#Making kinship matrices

K_Env=matrix(NA,nrow=nrow(Multitrait_NEW),ncol=nrow(Multitrait_NEW))

for (i in (1:nrow(K_Env))) {
	for (j in (1:ncol(K_Env))){
		if(Multitrait_NEW[i,"Planting"]==Multitrait_NEW[j,"Planting"]){
			K_Env[i,j]=1}
			else {
				K_Env[i,j]=0
			}
	}
}

#Standardizing kinship matrix

K_Env_st=standardize(K_Env)


K_Id=matrix(NA,nrow=nrow(Multitrait_NEW),ncol=nrow(Multitrait_NEW))

for (i in (1:nrow(K_Id))) {
	for (j in (1:ncol(K_Id))){
		if(Multitrait_NEW[i,"Ecotype_id"]==Multitrait_NEW[j,"Ecotype_id"]){
			K_Id[i,j]=1}
			else {
				K_Id[i,j]=0
			}
	}
}

K_Id_st=standardize(K_Id)