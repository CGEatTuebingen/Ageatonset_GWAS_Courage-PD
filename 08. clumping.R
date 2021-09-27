# Clumping

# Create input files for running clumping in Plink
uniq <- unique(unlist(sumstat_snpid$CHR))
dataset <- "courage_all"
for (i in 1:length(uniq)){
  data_1 <- subset(sumstat_snpid, CHR == uniq[i])
  myfile <- paste0("clump_input19df_mychr","_",uniq[i],"_", dataset,".txt" )
  write.table(data_1, myfile,row.names=FALSE,sep="\t", quote = FALSE) 
}
plink --bfile ./HRC_ld/chr1 --clump clump_input19df_mychr_1_courage_all.txt --clump-field P --clump-kb 250 --clump-p1 1 --clump-r2 0.1 --clump-snp-field SNP --memory 120000 --out chr1_clump_courage19df_all 
plink --bfile ./HRC_ld/chr2 --clump clump_input19df_mychr_2_courage_all.txt --clump-field P --clump-kb 250 --clump-p1 1 --clump-r2 0.1 --clump-snp-field SNP --memory 120000 --out chr2_clump_courage19df_all 
plink --bfile ./HRC_ld/chr3 --clump clump_input19df_mychr_3_courage_all.txt --clump-field P --clump-kb 250 --clump-p1 1 --clump-r2 0.1 --clump-snp-field SNP --memory 120000 --out chr3_clump_courage19df_all 
plink --bfile ./HRC_ld/chr4 --clump clump_input19df_mychr_4_courage_all.txt --clump-field P --clump-kb 250 --clump-p1 1 --clump-r2 0.1 --clump-snp-field SNP --memory 120000 --out chr4_clump_courage19df_all 
plink --bfile ./HRC_ld/chr5 --clump clump_input19df_mychr_5_courage_all.txt --clump-field P --clump-kb 250 --clump-p1 1 --clump-r2 0.1 --clump-snp-field SNP --memory 120000 --out chr5_clump_courage19df_all 
plink --bfile ./HRC_ld/chr7 --clump clump_input19df_mychr_7_courage_all.txt --clump-field P --clump-kb 250 --clump-p1 1 --clump-r2 0.1 --clump-snp-field SNP --memory 120000 --out chr7_clump_courage19df_all 
plink --bfile ./HRC_ld/chr8 --clump clump_input19df_mychr_8_courage_all.txt --clump-field P --clump-kb 250 --clump-p1 1 --clump-r2 0.1 --clump-snp-field SNP --memory 120000 --out chr8_clump_courage19df_all
plink --bfile ./HRC_ld/chr10 --clump clump_input19df_mychr_10_courage_all.txt --clump-field P --clump-kb 250 --clump-p1 1 --clump-r2 0.1 --clump-snp-field SNP --memory 120000 --out chr10_clump_courage19df_all
plink --bfile ./HRC_ld/chr14 --clump clump_input19df_mychr_14_courage_all.txt --clump-field P --clump-kb 250 --clump-p1 1 --clump-r2 0.1 --clump-snp-field SNP --memory 120000 --out chr14_clump_courage19df_all
plink --bfile ./HRC_ld/chr16 --clump clump_input19df_mychr_16_courage_all.txt --clump-field P --clump-kb 250 --clump-p1 1 --clump-r2 0.1 --clump-snp-field SNP --memory 120000 --out chr16_clump_courage19df_all
plink --bfile ./HRC_ld/chr21 --clump clump_input19df_mychr_21_courage_all.txt --clump-field P --clump-kb 250 --clump-p1 1 --clump-r2 0.1 --clump-snp-field SNP --memory 120000 --out chr21_clump_courage19df_all



masterlist =list.files(pattern = "19df_all\\.clumped$")
l <- lapply(masterlist, fread)
dt <- rbindlist(l)
dt <- dt %>% select(CHR,SNP, BP, P, SP2) %>%  mutate(SP2 = gsub("\\(1\\)", "", SP2))
write.csv(dt, "courage_clumping19df_all_results.csv")

female_sumstat_snpid  <- female_sumstat_snpid %>% rename(P = `P-value`)
female_sumstat_snpid  <- female_sumstat_snpid  %>% mutate(Allele1 = toupper(Allele1))
female_sumstat_snpid  <- female_sumstat_snpid  %>% mutate(Allele2 = toupper(Allele2))
female_sumstat_snpid  <- female_sumstat_snpid %>% separate("MarkerName", c("chr", "bp"), ":")
female_sumstat_snpid$chr <- as.numeric(as.character(female_sumstat_snpid$chr))
female_sumstat_snpid$bp <- as.numeric(as.character(female_sumstat_snpid$bp))
female_sumstat_snpid <- female_sumstat_snpid %>% filter(P < 0.00001)
#write.csv(female_significant_snps, "female_significant_snps.csv")
#write.csv(female_significant_snps, "female_significant_snps_050521.csv")
female_19df_sumstat_snpid <- female_sumstat_snpid %>% filter(HetDf>18) %>% select(chr, SNP, bp, Allele1, StdErr, P) %>% rename(CHR=chr,BP=bp, A1=Allele1, SE=StdErr)



# SNP annotation tools used for functional annotation and  nearby genes
# http://db.systemsbiology.net/kaviar/
# https://snipa.helmholtz-muenchen.de/snipa3/
# https://www.snp-nexus.org/v4/