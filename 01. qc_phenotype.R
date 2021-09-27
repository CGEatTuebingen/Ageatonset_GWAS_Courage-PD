# Inclusion and exclusion criteria 

## A. Load relevant libraries
library(tidyverse)
library(data.table)


## B. Loading complete phenotypic data from clinical collaborator (Includes cohort/s with age at onset data available only in cases)
data.dir <- "parent_directory_project"
clinical.fn <- sprintf("%s/phenotype/all_studies.txt", data.dir)
clinical <- fread(clinical.fn)
str(clinical)
# #Classes 'data.table' and 'data.frame':	24622 obs. of  15 variables:
# $ FID             : chr  "TRO_174580550_TRO_174580550" "TRO_174581017_TRO_174581017" "TRO_174581062_TRO_174581062" "TRO_177274656_TRO_177274656" ...
# $ IID             : chr  "TRO_174580550_TRO_174580550" "TRO_174581017_TRO_174581017" "TRO_174581062_TRO_174581062" "TRO_177274656_TRO_177274656" ...
# $ nstudy          : int  1 1 1 1 1 1 1 1 1 1 ...
# $ study_recoded   : chr  "aasly" "aasly" "aasly" "aasly" ...
# $ age_study       : num  77 NA NA NA NA NA NA NA NA NA ...
# $ age_onset       : num  59 NA NA NA NA NA NA NA NA NA ...
# $ age_diag        : int  NA NA NA NA NA NA NA NA NA NA ...
# $ pd_duration     : int  18 NA NA NA NA NA NA NA NA NA ...
# $ atc_fam         : int  NA 0 0 0 0 0 0 0 0 0 ...
# $ ethnicity       : int  1 NA NA NA NA NA NA NA NA NA ...
# $ GBA             : int  0 0 0 0 0 0 0 0 0 0 ...
# $ GBA_mutation    : chr  NA NA NA NA ...
# $ LRKK2           : int  0 0 0 0 0 0 0 0 0 0 ...
# $ LRKK2_mutation  : chr  NA NA NA NA ...
# $ atc_fam_mutation: int  NA 0 0 0 0 0 0 0 0 0 ...
# - attr(*, ".internal.selfref")=<externalptr>  
cat ("Total number of samples = ",(nrow(clinical)))
count_sites <- as.data.frame(count(clinical, study_recoded)) %>% arrange(-n)

## C.Loading in-house case-control phenotype data on age and gender (Includes cohort/s with a mixture of different ethnicites)
fam.fn <- sprintf("%s/phenotype/phenotype_fam.csv", data.dir)
clinical_fam <-  fread(fam.fn)
str(clinical_fam)
# Classes 'data.table' and 'data.frame':	24859 obs. of  3 variables:
#   $ sample_id: chr  "TRO_174580550" "TRO_174581017" "TRO_174581062" "TRO_177274656" ...
# $ gender   : int  1 1 2 2 1 2 1 1 2 2 ...
# $ disease  : int  2 1 1 1 1 1 1 1 1 1 ...
# - attr(*, ".internal.selfref")=<externalptr>
clinical_fam$sample_id_dupl <- clinical_fam$sample_id
clinical_fam$sample_id <- paste(clinical_fam$sample_id, clinical_fam$sample_id_dupl, sep="_")
clinical_fam <-  clinical_fam[,-c(4)]
cat ("Total number of samples = ",(nrow(clinical_fam)))

## D. Synchronized both the phenotype datasheets (Sites included with age at onset data in cases, sites excluded with different ethnicities)
clinical <- merge(clinical, clinical_fam, by.x = "FID", by.y = "sample_id", all.x=TRUE) %>% filter(study_recoded != c("XXX")) %>% mutate(age_diag = as.numeric(age_diag))
cat ("Total number of samples = ",(nrow(clinical)))
count_clinical <- as.data.frame(count(clinical, study_recoded)) %>% arrange(-n)

## E. Filter cases only and merge age at onset with age at diagnosis columns 
clinical %>% group_by(disease) %>% summarise(mean_age_onset = mean(age_onset, na.rm = T), number_samples = n(), missing = sum(is.na(age_onset)))

## F. Removal of samples with missing data on age at onset
clinical_cases <- clinical %>% filter(disease == 2) %>% mutate(age_onset_diag = coalesce(age_onset, age_diag)) %>% drop_na(age_onset_diag)
sum(is.na(clinical_cases$age_onset_diag))
clinical %>% group_by(disease) %>% summarise(mean_age_diag = mean(age_diag, na.rm = T), number_samples = n(), missing = sum(is.na(age_diag)))
clinical_cases %>% group_by(disease) %>% summarise(mean_age_onset = mean(age_onset_diag, na.rm = T), number_samples = n(), missing = sum(is.na(age_onset_diag)))
sum(is.na(clinical_cases$age_onset))
sum(is.na(clinical_cases$age_diag))
sum(is.na(clinical_cases$age_onset_diag))
cat ("Total number of samples = ",(nrow(clinical_cases)))
summary(clinical_cases$age_onset)
count_sites <- as.data.frame(count(clinical_cases, study_recoded)) %>% arrange(-n)
view(count_sites)

## F. Removal of cohorts with low number of samples (n<25)
### courageA represents specifc courage cohort with a low number of cases
clinical_cases <- clinical_cases %>% filter((study_recoded != c("courageA"))) 
cat ("Total number of samples = ",(nrow(clinical_cases)))

## G. Removal of completely overlapping IPDGC and COURGAGE cohorts
### ipdgcA and ipdgcB represents specific COURAGE cohorts that were completelely overlapping with IPDGC cohorts
clinical_cases <- clinical_cases %>% filter((study_recoded != c("ipdgcA")))  %>% filter((study_recoded != c("ipdgcB"))) 
cat ("Total number of samples = ",(nrow(clinical_cases)))

## H. Removal of overlapping IPDGC samples in the some of the COURAGE cohort/s
### ipdgcC_samples represents list of overlapping samples from a specific COURAGE cohort that were overlapping with IPDGC cohorts
ipdgcC <- as.vector(fread("ipdgcC_samples.txt"))
ipdgcC <- as.vector(ipdgcD$FID_ipdgcC)
clinical_cases <- filter(clinical_cases, !(FID %in% ipdgcC))
cat ("Total number of samples = ",(nrow(clinical_cases)))

## I. Remove of samples with family history of PD
clinical_cases$atc_fam_mutation[is.na(clinical_cases$atc_fam_mutation)] <- 0
clinical_cases_sporadic <- clinical_cases %>% filter(atc_fam_mutation == 0)
count_sites <- as.data.frame(count(clinical_cases_sporadic, study_recoded)) %>% arrange(-n)
view(count_sites)

## J. Added a column for ethnicity and summarize AAO by different genders and ethnic groups
clinical_cases_sporadic <- clinical_cases_sporadic %>% mutate(ethnicity = if_else(study_recoded == "hattori"|study_recoded == "kim"|study_recoded == "matsuo"|study_recoded == "chung", "Asian", "Caucasian"))
table(clinical_cases_sporadic$ethnicity)
clinical_cases_sporadic %>% group_by(ethnicity) %>% summarise(mean = mean(age_onset), sd = sd(age_onset))
clinical_cases_sporadic %>% group_by(gender) %>% summarise(mean = mean(age_onset), sd = sd(age_onset))