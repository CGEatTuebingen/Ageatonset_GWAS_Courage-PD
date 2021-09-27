#Heitability and correlation between GWAS datasets using LDSC

## A.Load relevant libraries
library(data.table)
library(tidyverse)

## B. Load the input meta-analysis file and format it for LDSC; On terminal/ Windows using R
dataset_A <-
  fread("cosortium_sumstat_A.tbl")
hrc_data <-
  fread("hrc_annotations.txt")
consortiumA_sumstat_annotated <-
  left_join(dataset_A, hrc_data) %>% mutate(Allele1 = toupper(cosortiumA_sumstat_annotated$Allele1)) %>% mutate(Allele2 = toupper(cosortiumA_sumstat_annotated$Allele2)) %>%
  rename(P.value = `P-value`) %>% filter(HetISq < 50) %>% filter(HetDf > 16) %>% select(SNP, Allele1, Allele2, P.value, Effect, StdErr) %>%
  write.table("consortium_sumstat_ldsc_input.txt",
              row.names = F,
              quote = FALSE)

## C. Generate h2 for different datasets; On terminal; Tool: LDSC; N: Sample size
munge_sumstats.py --sumstats consortiumA_sumstat_ldsc_input.txt --N XXXX --merge -alleles eur_w_ld_chr/ w_hm3.snplist --chunksize 500000 --out consortiumA_sumstat_ldsc_input
ldsc.py --h2 consortiumA_sumstat_ldsc_input.sumstats.gz --ref -ld -chr eur_w_ld_chr/ --w -ld -chr eur_w_ld_chr/ --out consortiumA_sumstat_h2

munge_sumstats.py --sumstats consortiumB_sumstat_ldsc_input.txt --N YYYY --merge -alleles eur_w_ld_chr/ w_hm3.snplist --chunksize 500000 --out consortiumB_sumstat_ldsc_input
ldsc.py --h2 consortiumB_sumstat_ldsc_input.sumstats.gz --ref -ld -chr eur_w_ld_chr/ --w -ld -chr eur_w_ld_chr/ --out consortiumB_sumstat_h2

## D. Generate correlation between different datasets
ldsc.py --rg consortiumA_sumstat_ldsc_input.sumstats.gz,consortiumB_sumstat_ldsc_input.sumstats.gz  --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out  correaltion_bip 