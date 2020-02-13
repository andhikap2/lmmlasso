#Lasso beta shenanigans

library(dplyr)

a=read.csv("Lasso_Betas.csv")
a=data.frame(a)
colnames(a) = gsub("X", "", colnames(a))

grnd_temp=a %>% select(starts_with("Grnd."))
snps1 =a %>% select(starts_with("1_"))
snps2 =a %>% select(starts_with("2_"))
snps3 =a %>% select(starts_with("3_"))
snps4 =a %>% select(starts_with("4_"))
snps5 =a %>% select(starts_with("5_"))
snps=cbind(snps1,snps2,snps3,snps4,snps5)
cum_ptt=a %>% select(starts_with("cum"))


cum_ptt=cum_ptt[,colSums(cum_ptt)!=0] #Remove non-zero betas
nz_snps=snps[,colSums(snps)!=0]

#Rank snp effects by absolute values
d=sort(abs(unlist(nz_snps[1,])))

snp_ids=colnames(nz_snps)
write.table(snp_ids,file="SNP_ids.txt",sep="\t")