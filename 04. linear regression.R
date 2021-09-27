# Continued from 03.add PCs_phenotype.R
## A. Create cohort specific sample list from finalized phenotypic sheet for filtering samples in GWAS files; Tool: Command-line on linux terminal
study_list<-unique(clinical_cases_sporadic_pcs$study_recoded)
head(study_list)
# [1] kruger   kim      stefanis tolosa   mellick  annesi  
# 30 Levels: aasly annesi bardien brice brighina carmine_belin chartier_harlin chung deutschlander elbaz farrer ferreira gasser ... zimprich

for(study in study_list)
{
  file=paste0(study, "_sample_list.txt")
  clinical_cases_sporadic_pcs %>% filter(study_recoded == study) %>% select("FID", "IID") %>%  write.table(file, row.names = F, quote = FALSE)
}



## B. Create GWAS data for age at onset GWAS (bim, fam and bed) from in-house PD GWAS with finalized for age at onset; Tool: Plink on linux terminal
### Repeat the process for each study cohort
plink --bfile study_gwasdata --keep study_sample_list.txt --keep-allele-order --make-bed --out study_ageatonset

## C. Convert plink files to vcf files; Tool: Plink on linux terminal
### Repeat the process for each study cohort
plink --bfile study_ageatonset --recode vcf --keep-allele-order --out study_ageatonset --threads 3 --memory 30000 

## D. Run Linear regression (Cohort specific); Tool: rvtest on Linux terminal 
### Repeat the process for each study cohort and for each subgroup (e.g.males as shown below in gender specific regression command)

### Example header from Cohort specific input file
head aasly_pheno.ped
# fid iid fatid matid sex y1 y2 y3 y4
# TRO_174580550_TRO_174580550_TRO_174580550_TRO_174580550 TRO_174580550_TRO_174580550_TRO_174580550_TRO_174580550 0 0 1 59 0 0 0
# TRO_177275912_TRO_177275912_TRO_177275912_TRO_177275912 TRO_177275912_TRO_177275912_TRO_177275912_TRO_177275912 0 0 2 65 0 0 0
# TRO_177275913_TRO_177275913_TRO_177275913_TRO_177275913 TRO_177275913_TRO_177275913_TRO_177275913_TRO_177275913 0 0 2 42 0 0 0
# TRO_177275914_TRO_177275914_TRO_177275914_TRO_177275914 TRO_177275914_TRO_177275914_TRO_177275914_TRO_177275914 0 0 1 55 0 0 0
# TRO_177275920_TRO_177275920_TRO_177275920_TRO_177275920 TRO_177275920_TRO_177275920_TRO_177275920_TRO_177275920 0 0 1 59 0 0 0
# TRO_177275934_TRO_177275934_TRO_177275934_TRO_177275934 TRO_177275934_TRO_177275934_TRO_177275934_TRO_177275934 0 0 1 48 0 0 0

### Example header from Cohort specific input file
head aasly_covar.covar
# fid iid sex pc1 pc2 pc3 pc4 pc5
# TRO_174580550_TRO_174580550_TRO_174580550_TRO_174580550 TRO_174580550_TRO_174580550_TRO_174580550_TRO_174580550 1 -0.000204557 -0.00407448 0.000559757 -0.0125841 0.0150817
# TRO_177275912_TRO_177275912_TRO_177275912_TRO_177275912 TRO_177275912_TRO_177275912_TRO_177275912_TRO_177275912 2 -0.014812 0.00362692 -0.00469811 -0.0148678 -0.00525114
# TRO_177275913_TRO_177275913_TRO_177275913_TRO_177275913 TRO_177275913_TRO_177275913_TRO_177275913_TRO_177275913 2 -0.00189356 -0.0070984 -0.00639011 0.00384111 0.00424962
# TRO_177275914_TRO_177275914_TRO_177275914_TRO_177275914 TRO_177275914_TRO_177275914_TRO_177275914_TRO_177275914 1 0.000965011 0.00375134 -0.00787583 -0.00737997 -0.00302513
# TRO_177275920_TRO_177275920_TRO_177275920_TRO_177275920 TRO_177275920_TRO_177275920_TRO_177275920_TRO_177275920 1 -0.00371203 -0.0169196 -0.006055 0.00474924 -0.00565091
# TRO_177275934_TRO_177275934_TRO_177275934_TRO_177275934 TRO_177275934_TRO_177275934_TRO_177275934_TRO_177275934 1 -0.00671023 -0.0127435 0.00804586 0.00582786 0.00148317

### Example header from cohort specific input file
header aasly_sample_list.txt
# FID IID
# TRO_174580550_TRO_174580550 TRO_174580550_TRO_174580550
# TRO_177275912_TRO_177275912 TRO_177275912_TRO_177275912
# TRO_177275913_TRO_177275913 TRO_177275913_TRO_177275913
# TRO_177275914_TRO_177275914 TRO_177275914_TRO_177275914
# TRO_177275920_TRO_177275920 TRO_177275920_TRO_177275920
# TRO_177275934_TRO_177275934 TRO_177275934_TRO_177275934

### Cohort specific linear regression; Tool: rvtest on Linux terminal
rvtest --inVcf study_ageatonset.vcf --pheno study_pheno.ped --pheno-name y1 --covar study_covar.covar --covar-name sex,pc1,pc2,pc3,pc4,pc5 --out study --single wald,score

### Example header from Cohort specific output file
head aasly.SingleWald.assoc
# CHROM   POS     REF     ALT     N_INFORMATIVE   Test    Beta    SE      Pvalue
# 1       751580  T       C       500     1:751580        20.9432 10.7101 0.0505291
# 1       751580  T       C       500     sex     2.69895 0.981511        0.00596325
# 1       751580  T       C       500     pc1     -14.3461        46.6215 0.758301
# 1       751580  T       C       500     pc2     87.5936 49.6013 0.0774036
# 1       751580  T       C       500     pc3     5.80774 56.4363 0.918036
# 1       751580  T       C       500     pc4     12.2833 46.6568 0.792344

### Example header from Cohort specific output file
head aasly.SingleScore.assoc
# CHROM   POS     REF     ALT     N_INFORMATIVE   AF      U       V       STAT    DIRECTION       EFFECT  SE      PVALUE
# 1       751580  T       C       500     0.999   20.8555 114.618 3.79479 +       20.9432 10.751  0.0514123
# 1       752566  A       G       500     0.177   134.317 16860.2 1.07005 +       0.916952        0.886431        0.300935
# 1       752721  G       A       500     0.178   113.462 17119.5 0.751986        +       0.762844        0.879693        0.385848
# 1       753541  A       G       500     0.853   -145.484        15916.9 1.32976 -       -1.05204        0.912318        0.248848
# 1       755274  T       C       500     0.999   20.8555 114.618 3.79479 +       20.9432 10.751  0.0514123
# 1       755868  T       C       500     0.999   20.8555 114.618 3.79479 +       20.9432 10.751  0.0514123
# 1       757120  G       A       500     0.999   20.8555 114.618 3.79479 +       20.9432 10.751  0.0514123


### Gender specific linear regression (Cohort specific); Tool: rvtest on Linux terminal 
rvtest --inVcf study_ageatonset.vcf --pheno study_pheno.ped --pheno-name y1 --covar study_covar.covar --covar-name pc1,pc2,pc3,pc4,pc5 --peopleIncludeFile study_males_sample_list.txt --out study_male --single wald,score

## E. Format the regression output (Cohort specific); Tool: Standard R package on Windows/Terminal  
### Repeat the process for each study cohort
library(data.table)
library(tidyverse)
fread("study.SingleScore.assoc") %>% mutate(SNP = paste(CHROM, POS, sep = ":")) %>% rename(N = N_INFORMATIVE) %>% fwrite("study_ageatonset_summstat.txt")

### Structure of an example regression output file
str(fread("aasly_ageatonset_summstat.txt"))
# Classes 'data.table' and 'data.frame':  10891340 obs. of  14 variables:
#   $ CHROM    : int  1 1 1 1 1 1 1 1 1 1 ...
# $ POS      : int  751580 752566 752721 753541 755274 755868 757120 762856 766007 769223 ...
# $ REF      : chr  "T" "A" "G" "A" ...
# $ ALT      : chr  "C" "G" "A" "G" ...
# $ N        : int  500 500 500 500 500 500 500 500 500 500 ...
# $ AF       : num  0.999 0.177 0.178 0.853 0.999 0.999 0.999 0.999 0.876 0.854 ...
# $ U        : num  20.9 134.3 113.5 -145.5 20.9 ...
# $ V        : num  115 16860 17120 15917 115 ...
# $ STAT     : num  3.795 1.07 0.752 1.33 3.795 ...
# $ DIRECTION: chr  "+" "+" "+" "-" ...
# $ EFFECT   : num  20.943 0.917 0.763 -1.052 20.943 ...
# $ SE       : num  10.751 0.886 0.88 0.912 10.751 ...
# $ PVALUE   : num  0.0514 0.3009 0.3858 0.2488 0.0514 ...
# $ SNP      : chr  "1:751580" "1:752566" "1:752721" "1:753541" ...

## F. Meta-analysis; Tool: Metal on Terminal
### Sample metal script (metal_script_ageatonset_consortium_date.txt)
# SCHEME STDERR
# GENOMICCONTROL ON
# MINMAXFREQ ON
# MARKER   SNP
# ALLELE   ALT REF
# PVALUELABEL     PVALUE
# EFFECTLABEL   EFFECT
# STDERR   SE
# ANALYZE HETEROGENEITY
# COLUMNCOUNTING LENIENT
# SEPARATOR  COMMA
# FREQ AF 
# PROCESS study1_ageatonset_summstat.txt
# PROCESS study2_ageatonset_summstat.txt
# PROCESS study3_ageatonset_summstat.txt
# ...
# PROCESS study30_ageatonset_summstat.txt
# OUTFILE meta_analysis_conosrtium_date.tbl
# ANALYZE HETEROGENEITY
# QUIT

### Run the meta-analysis; Tool: Metal on Linux terminal
metal metal_script_ageatonset_consortium_date.txt

### Read the output file in R
sumstat <- fread("meta_analysis_courage_05052021.tbl")
head(sumstat)
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

### Annotate the output file 
hrc_data <- fread("hrc_annotations.txt")
# hrc_annotations.txt was obtained from HRC which contains SNP rsids and corresponding positions
hrc_data <- hrc_data %>% select(MarkerName, SNP)
sumstat_snpid <- left_join(sumstat, hrc_data); rm(sumstat)
sumstat_snpid   <- sumstat_snpid %>% rename(P = `P-value`)
sumstat_snpid   <- sumstat_snpid  %>% separate("MarkerName", c("CHR", "BP"), ":")

### Remove bad quality SNP
sumstat_snpid <- sumstat_snpid %>% filter(HetDf>18 & HetISq<=50)

# In a similar fashion, stratified meta-analysis was conducted on the basis of gender and ethnicity 

# All the summary statitistics files for courage cohorts will be available on request 
# Grover2022_courage_aao_sumstat.txt
# Grover2022_courage_aao_caucasian_sumstat.txt
# Grover2022_courage_aao_asian_sumstat.txt
# Grover2022_courage_aao_males_sumstat.txt
# Grover2022_courage_aao_females_sumstat.txt





