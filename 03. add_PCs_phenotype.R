# Continued from 02. descriptive_phenotype.R
# Merge Principal compnents in the phenotypic data file
# Prinicpal components for each cohort were derived from QC ran during PD GWAS, and stored in cohort specific text files.
# A folder (name enclosed within path below) was created with text files for each cohort names where 1st and 2nd columns correspond 
# to sample id (columns named V1 and V2) and remaining columns correspond to PC1 to PC5 (Columns named V3 to V7)
all.files <- list.files(path = ".", pattern = ".txt")
l <- lapply(all.files, data.table::fread)
dt <- data.table::rbindlist(l)
str(dt)
# Classes 'data.table' and 'data.frame':	24859 obs. of  7 variables:
#   $ V1: chr  "study1_174580550_study1_174580550" "study1__174581017_study1__174581017" "study1__174581062_study1__174581062" "study1__177274656_study1__177274656" ...
# $ V2: chr  "study1__174580550_study1__174580550" "study1__174581017_study1__174581017" "study1__174581062_study1__174581062" "study1__177274656_study1__177274656" ...
# $ V3: num  -0.000205 0.000817 0.01412 0.005044 0.001199 ...
# $ V4: num  -0.004074 -0.002514 0.018957 -0.006775 0.000713 ...
# $ V5: num  0.00056 0.00708 0.00336 -0.00566 0.00605 ...
# $ V6: num  -0.012584 0.000476 0.001673 -0.008008 0.004574 ...
# $ V7: num  0.015082 0.008378 0.007267 0.001258 0.000677 ...
# - attr(*, ".internal.selfref")=<externalptr> 
dt <-
  dt %>% rename(
    FID = V1,
    IID = V2,
    PC1 = V3,
    PC2 = V4,
    PC3 = V5,
    PC4 = V6,
    PC5 = V7
  )
dt <-  dt %>% select(-IID)
clinical_cases_sporadic_pcs <-
  clinical_cases_sporadic %>% left_join(dt)
clinical_cases_sporadic_pcs %>% select(
  FID,
  IID,
  study_recoded,
  gender,
  age_study,
  age_diag,
  age_onset,
  pd_duration,
  PC1,
  PC2,
  PC3,
  PC4,
  PC5
) %>% data.table::fwrite("phenotype_data_ageatonset.txt")
