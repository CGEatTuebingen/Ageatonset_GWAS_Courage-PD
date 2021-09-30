#Heritability and correlation between GWAS datasets using LDSC

## A.Load relevant libraries
library(data.table)
library(tidyverse)

## B. Load the input meta-analysis file and format it for LDSC; On terminal/ Windows using R
dataset_A <-
  fread("cosortium_sumstat_A.tbl")
### e.g of meta-analysis file - "meta_analysis_courage_05052021.tbl" (See 04. Linear regression.R)
head(dataset_A)
# Classes 'data.table' and 'data.frame':  20473650 obs. of  15 variables:
#   $ MarkerName: chr  "5:85928892" "14:36082010" "11:107819621" "21:41565796" ...
# $ Allele1   : chr  "t" "t" "a" "t" ...
# $ Allele2   : chr  "c" "c" "c" "c" ...
# $ Freq1     : num  0.066 0.0019 0.0103 0.0023 0.989 ...
# $ FreqSE    : num  0.0146 0.002 0.0034 0 0.0029 0.0264 0.004 0.0006 0.0013 0.0025 ...
# $ MinFreq   : num  0.0362 0.0007 0.0023 0.0023 0.9841 ...
# $ MaxFreq   : num  0.0939 0.0054 0.012 0.0023 0.9972 ...
# $ Effect    : num  0.1941 0.0961 3.4662 -21.2955 -1.7573 ...
# $ StdErr    : num  0.344 5.64 2.843 10.041 1.008 ...
# $ P-value   : num  0.5727 0.9864 0.2227 0.0339 0.0812 ...
# $ Direction : chr  "--+-++++++---+-+++-+-+-+?++-++" "?????????-?????????????+???+??" "??????????????+-???+??????????" "???????????????-??????????????" ...
# $ HetISq    : num  22.6 37.8 0 0 0 0 0 0 66.8 0 ...
# $ HetChiSq  : num  36.19 3.22 1.3 0 16.37 ...
# $ HetDf     : int  28 2 2 0 19 29 3 1 2 22 ...
# $ HetPVal   : num  0.138 0.2 0.523 1 0.632 ...
# - attr(*, ".internal.selfref")=<externalptr>

hrc_data <-
  fread("hrc_annotations.txt")

consortiumA_sumstat_annotated <-
  left_join(dataset_A, hrc_data) %>% mutate(Allele1 = toupper(cosortiumA_sumstat_annotated$Allele1)) %>% mutate(Allele2 = toupper(cosortiumA_sumstat_annotated$Allele2)) %>%
  rename(P.value = `P-value`) %>% filter(HetISq < 50) %>% filter(HetDf > 16) %>% select(SNP, Allele1, Allele2, P.value, Effect, StdErr) %>%
  write.table("consortiumA_sumstat_ldsc_input.txt",
              row.names = F,
              quote = FALSE)

### head of ldsc input file
head(fread("consortiumA_sumstat_ldsc_input.txt"))
# SNP Allele1 Allele2 P.value  Effect StdErr
# 1: rs113534962       T       C  0.5727  0.1941 0.3442
# 2: rs559397866       T       C  0.0812 -1.7573 1.0077
# 3:   rs2366866       T       C  0.1377  0.2504 0.1687
# 4:  rs79253331       T       C  0.2688  1.4123 1.2772
# 5:  rs62099898       T       C  0.4961 -0.1264 0.1857
# 6:  rs11725240       T       C  0.1402  0.3328 0.2256

### Similarly generate ldsc input for other consortium using metal ouput from COURAGE PD GWAS or IPDGC AAO GWAS

## C. Generate h2 for different datasets; On terminal; Tool: LDSC; N: Sample size
munge_sumstats.py--sumstats consortiumA_sumstat_ldsc_input.txt--N 16443--merge -
  alleles eur_w_ld_chr / w_hm3.snplist--chunksize 500000--out consortiumA_sumstat_ldsc_input
ldsc.py--h2 consortiumA_sumstat_ldsc_input.sumstats.gz--ref - ld - chr eur_w_ld_chr / --w -
  ld - chr eur_w_ld_chr / --out consortiumA_sumstat_h2


munge_sumstats.py--sumstats consortiumB_sumstat_ldsc_input.txt--N 14919 --merge -
  alleles eur_w_ld_chr / w_hm3.snplist--chunksize 500000--out consortiumB_sumstat_ldsc_input
ldsc.py--h2 consortiumB_sumstat_ldsc_input.sumstats.gz--ref - ld - chr eur_w_ld_chr / --w -
  ld - chr eur_w_ld_chr / --out consortiumB_sumstat_h2

### Sample output from consortiumA_sumstat_h2.log file
# Reading summary statistics from consortiumA_sumstat_ldsc_input.sumstats.gz ...
# Read summary statistics for 1061219 SNPs.
# Reading reference panel LD Score from eur_w_ld_chr/[1-22] ... (ldscore_fromlist)
# Read reference panel LD Scores for 1290028 SNPs.
# Removing partitioned LD Scores with zero variance.
# Reading regression weight LD Score from eur_w_ld_chr/[1-22] ... (ldscore_fromlist)
# Read regression weight LD Scores for 1290028 SNPs.
# After merging with reference panel LD, 1048041 SNPs remain.
# After merging with regression SNP LD, 1048041 SNPs remain.
# Using two-step estimator with cutoff at 30.
# Total Observed scale h2: 0.0835 (0.0568)
# Lambda GC: 0.9927
# Mean Chi^2: 0.9885
# Intercept: 0.9736 (0.0072)
# Ratio: NA (mean chi^2 < 1)

## D. Generate correlation between different datasets
ldsc.py--rg consortiumA_sumstat_ldsc_input.sumstats.gz, consortiumB_sumstat_ldsc_input.sumstats.gz --ref -
  ld - chr eur_w_ld_chr / --w - ld - chr eur_w_ld_chr / --out  correaltion_bip

### Sample output of correlation between different datasets
  # Reading summary statistics from consortiumA_sumstat_ldsc_input.sumstats.gz ...
 # Read summary statistics for 1061219 SNPs.
 # Reading reference panel LD Score from eur_w_ld_chr/[1-22] ... (ldscore_fromlist)
 # Read reference panel LD Scores for 1290028 SNPs.
 # Removing partitioned LD Scores with zero variance.
 # Reading regression weight LD Score from eur_w_ld_chr/[1-22] ... (ldscore_fromlist)
 # Read regression weight LD Scores for 1290028 SNPs.
 # After merging with reference panel LD, 1048041 SNPs remain.
 # After merging with regression SNP LD, 1048041 SNPs remain.
 # Computing rg for phenotype 2/2
 # Reading summary statistics from consortiumB_sumstat_ldsc_input.sumstats.gz ...
 # Read summary statistics for 1217311 SNPs.
 # After merging with summary statistics, 1048041 SNPs remain.
 # 1047234 SNPs with valid alleles.
 # 
 # Heritability of phenotype 1
 # ---------------------------
 #   Total Observed scale h2: 0.0833 (0.0577)
 # Lambda GC: 0.9927
 # Mean Chi^2: 0.9885
 # Intercept: 0.9737 (0.0072)
 # Ratio: NA (mean chi^2 < 1)
 # 
 # Heritability of phenotype 2/2
 # -----------------------------
 #   Total Observed scale h2: 0.2914 (0.041)
 # Lambda GC: 1.0315
 # Mean Chi^2: 1.0489
 # Intercept: 0.9485 (0.0082)
 # Ratio < 0 (usually indicates GC correction).
 # 
 # Genetic Covariance
 # ------------------
 #   Total Observed scale gencov: -0.0462 (0.0304)
 # Mean z1*z2: -0.011
 # Intercept: 0.0008 (0.0052)
 # 
 # Genetic Correlation
 # -------------------
 #   Genetic Correlation: -0.2964 (0.2245)
 # Z-score: -1.3204
 # P: 0.1867
 # 
 # Summary of Genetic Correlation Results
 # p1                                    p2      rg      se       z       p  h2_obs  h2_obs_se  h2_int  h2_int_se  gcov_int  gcov_int_se
 # consortiumA_sumstat_ldsc_input.sumstats.gz  consortiumB_sumstat_ldsc_input.sumstats.gz -0.2964  0.2245 -1.3204  0.1867  0.2914      0.041  0.9485     0.0082    0.0008       0.0052
