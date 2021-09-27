# Individual SNP data 

#BST1 SNP data extraction from plink
plink --bfile mergedageatonset_pooled --snps 4:15737348:G:A --recode --out BSTA1SNP

bsta1_rs4698412 <- fread("BSTA1SNP.ped")
bsta1_rs4698412$SNP <- paste(bsta1_rs4698412$V7, bsta1_rs4698412$V8, sep="")
bsta1_rs4698412 <-  bsta1_rs4698412 %>% rename(FID=V1) %>% select(-V2,-V3, -V4, -V7, -V8, -V5, -V6)
bsta1_aao_rs4698412 <- left_join(clinical_cases_sporadic, bsta1_rs4698412) 

bsta1_aao_rs4698412 <- bsta1_aao_rs4698412 %>% filter(SNP != "00") 
bsta1_aao_rs4698412 %>% group_by(SNP) %>% summarise(n=n(), mean = mean(age_onset), sd = sd(age_onset))


# Graph
p<-ggplot(bsta1_aao_rs4698412, aes(x=SNP, y=age_onset))+ geom_density(adjust=1.5, alpha=.4) +
  theme_ipsum()

ggplot(bsta1_aao_rs4698412, aes(x = age_onset)) +
  geom_density(aes(color = SNP))