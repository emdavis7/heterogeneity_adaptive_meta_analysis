#####
## Import vital table
#####


### Reduce vital ####
# reduce vital to first observationoffset for each patientunitstayid
# keep patientunitstayid, observationoffset, noninvasivesystolic, and noninvasivediastolic.
# only keep rows where both bp are not missing and the observation was recorded after arrival in ICU.
# (observationoffset>0)
vital <- read_csv("./vitalAperiodic.csv")
vital2 <- vital %>%
  filter(is.na(noninvasivesystolic)==F &
           is.na(noninvasivediastolic)==F &
           observationoffset >= 0) %>%
  dplyr::select(patientunitstayid, observationoffset, 
                noninvasivesystolic, noninvasivediastolic) %>%
  group_by(patientunitstayid) %>%
  arrange(observationoffset) %>%
  filter(row_number()==1)

saveRDS(vital2, "./vital_reduced.rds")
rm(vital, vital2)
