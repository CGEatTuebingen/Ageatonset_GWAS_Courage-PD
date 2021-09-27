# PRS SNP extraction from case control analysis
## A. Load relevant libraries
library(tidyverse)
library(data.table)
library(readxl)

## Load top hits from case control analysis
df <- read_excel("courage+ipdgc_010721.xlsx")
mergedageatonset_bim <- fread("mergedageatonset_pooled.bim")

write.table(prs_bim, "snpforprs.bim", row.names = F, quote = F, sep = "\t")

patterns <- paste0("^", df$`CHR:BP`)
result <- filter(mergedageatonset_bim, grepl(paste(patterns, collapse="|"), mergedageatonset_bim$V2))
write.table(result$V2, "snplist.txt", col.names = F, row.names = F, quote = F)


filter(mergedageatonset_bim, grepl("^1:88362", mergedageatonset_bim$V2))

plink --bfile mergedageatonset_pooled --extract snplist.txt  --recode --out snpforprs
prs_ped <- fread("snpforprs.ped")
prs_map <- fread("snpforprs.map")
## Merge alleles into genotypes for each SNP
prs_ped <- prs_ped %>% unite("1:152192927",  V7:V8, sep = "",  remove = T)
write.csv(prs_ped, "prs_input_aao.csv")

prs_ped <- prs_ped %>% rename(FID = V1)
prs_input_cov <- left_join(clinical_cases_sporadic_pcs, prs_ped)
prs_input_cov <- prs_input_cov %>% rename(SEX = gender) %>% select(FID, IID, SEX, PC1:PC5) %>% write.table("prs_cov_aao.txt", row.names = F, quote = F, sep = "\t")
prs_input_cov <- left_join(clinical_cases_sporadic_pcs, prs_ped)

prs_input_cov <-  prs_input_cov %>% rename(PHENOTYPE = age_onset) %>% select(FID, IID, PHENOTYPE) %>% write.table("prs_pheno_aao.txt", row.names = F, quote = F, sep = "\t")
df <- read_excel("courage+ipdgc_010721.xlsx")

df <- df %>% rename(snp = MarkerName, A1 = Allele1, A2= Allele2, beta = Effect, pvalue = 'P-value of Courage + IPDGC meta-analysis') %>% separate('CHR:BP', into = c("chr", "pos"), sep = ":")
df <- df %>% select(snp, chr, pos, A1, A2, beta, pvalue)
write.table(df, "prs_base_ipdgc_courage.assoc", row.names = F, quote = F, sep = "\t")
prs_ped <- fread("snpforprs.ped")
prs_ped %>% 
  prs_bim <- fread("snpforprs.bim")
prs_bim$V2 <- str_sub(prs_bim$V2, end=-5)
write.table(prs_bim, "snpforprs.bim", row.names = F, quote = F, sep = "\t")


# extract SNPdata
plink --bfile mergedageatonset_pooled --extract snplist.txt --make-bed --out snpforprs
#update rsid
plink --bfile snpforprs --update-name updatersid.txt --make-bed --out snpforprs


setwd("D:/courage_data/data/metanalysis_05_05_2021")
Rscript PRSice.R --dir . \
--prsice ./PRSice_linux \
--base prs_base_ipdgc_courage.assoc \
--ignore-fid \
--snp snp --chr chr --bp bp --A1 A1 --A2 A2 --stat beta --pvalue pvalue \
--target snpforprs \
--pheno prs_pheno_aao.txt \
--pheno-col PHENOTYPE \
--cov prs_cov_aao.txt \
--cov-col SEX,PC1,PC2,PC3,PC4,PC5 \
--score sum \
--stat beta \
--beta \
--binary-target F \
--perm 10000 \
--cov-factor SEX \


prsice_data <- fread("PRSice.summary")
