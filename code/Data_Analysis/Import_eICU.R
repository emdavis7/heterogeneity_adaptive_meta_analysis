##############
## Import eICU data and deduplicate.
##############
library(readr)
library(dplyr)
library(janitor)

# Import and reduce vital:
source("./code/data_analysis/import_vital.R")

# Import data
vital <- readRDS("./vital_reduced.rds")
apachePatRes <- read_csv("./apachePatientResult.csv")
hospital <- read_csv("./hospital.csv")
patient <- read_csv("./patient.csv")
apachePredVar <- read_csv("./apachePredVar.csv")



##### Keep the resultId that's highest for each patient
patRes_dedup <- apachePatRes %>%
  group_by(patientunitstayid) %>%
  arrange(patientunitstayid, -apachepatientresultsid) %>%
  ungroup()
patRes_dedup <-
  patRes_dedup[!duplicated(patRes_dedup$patientunitstayid),]

## Deduplicate patient table to keep the first instance of the
# patient in the hospital, as unit transfers
## create a new row.
patient_dedup <- patient %>%
  group_by(uniquepid) %>%
  arrange(uniquepid, unitvisitnumber) %>%
  ungroup()
patient_dedup <- patient_dedup[!duplicated(patient_dedup$uniquepid),]

## There are different numbers in each; keep the minimum set.
full <- patient_dedup %>%
  left_join(apachePredVar,
            by= ("patientunitstayid" = "patientunitstayid"),
            keep = F) %>%
  left_join(patRes_dedup,
            by = ("patientunitstayid" = "patientunitstayid"),
            keep = F) %>%
  left_join(vital,
            by = ("patientunitstayid" = "patientunitstayid"),
            keep = F) %>%
  left_join(hospital, by = ("hospitalid" = "hospitalid"), keep = F)


## Use only observations with APACHE version "IVa" and
# treat APACHE score == -1 as missing.
full2 <- full %>%
  filter(apacheversion=="IVa" & apachescore!=-1)


######### Select a subset of variables ############
reduced <- full2 %>% dplyr::select(hospitalid,
                          hospitaldischargeyear,
                            region.y,
                           gender.x,
                           age.x,
                           ethnicity,
                           hospitaladmitsource,
                           actualiculos,
                           noninvasivesystolic,
                           noninvasivediastolic,
                           hospitaldischargestatus,
                           apachescore,
                           teachingstatus
                           ) %>%
  mutate(age.x = as.numeric(age.x))

# Keep complete data:
reduced2 <- reduced[complete.cases(reduced),]
# Note: because we used as.numeric on age, induced missingness for ages >89. 
# These individuals are removed in line 75.

### Save the reduced data
saveRDS(reduced2, "./reduced_eICU.rds")
