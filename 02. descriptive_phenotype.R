# Create descriptive table for paper
# eTable1
# Continued from 01.qc_phenotype.R
clinical_cases_sporadic$gender <-
  factor(
    clinical_cases_sporadic$gender,
    levels = c(1, 2),
    labels =  c("M", "F")
  )
clinical_cases_sporadic$study_recoded <-
  factor(clinical_cases_sporadic$study_recoded)
clinical_cases_sporadic$age_study <-
  as.numeric(clinical_cases_sporadic$age_study)
clinical_cases_sporadic$age_onset <-
  as.numeric(clinical_cases_sporadic$age_onset)
clinical_cases_sporadic$pd_duration <-
  as.numeric(clinical_cases_sporadic$pd_duration)
clinical_cases_sporadic$age_diag <-
  as.numeric(clinical_cases_sporadic$age_diag)
table1::table1( ~ gender + age_study + age_onset + age_diag + pd_duration |
                  study_recoded,
                data = clinical_cases_sporadic)