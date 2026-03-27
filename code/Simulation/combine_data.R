####
# Combine result files from simulation
####

library(purrr)
library(dplyr)

### Study 0 ######################################################################
file_list <- list.files(path = 'data/sim_setting0',
                        pattern = "*.rds", full.names = T)
result_list <- map(file_list, readRDS)

full_results <- bind_rows(result_list)

# Output the full_results
saveRDS(full_results, file = "./data/results_setting0.rds", compress = T)
rm(result_list)




### Study 1 MSE ####################################################################
file_list <- list.files(path = 'data/sim_setting1/mse',
                        pattern = "*.rds", full.names = T)
result_list <- map(file_list, readRDS)

full_results <- bind_rows(result_list)

# Output the full_results
saveRDS(full_results, file = "./data/results_setting1_mse.rds", compress = T)

#### Study 1 CI ####################################################################3

file_list <- list.files(path = 'data/sim_setting1/ci',
                        pattern = "*.rds", full.names = T)
result_list <- map(file_list, readRDS)

full_results <- bind_rows(result_list)

# Output the full_results
saveRDS(full_results, file = "./data/results_setting1_ci.rds", compress = T)





#### Study 2 Rep Level ####################################################################

file_list <- list.files(path = 'data/sim_setting2/rep',
                        pattern = "*.rds", full.names = T)
result_list <- map(file_list, readRDS)

full_results <- bind_rows(result_list)

# Output the full_results
saveRDS(full_results, file = "./data/results_setting2_rep.rds", compress = T)

#### Study 2 Study Level ####################################################################

file_list <- list.files(path = 'data/sim_setting2/study',
                        pattern = "*.rds", full.names = T)
result_list <- map(file_list, readRDS)

full_results <- bind_rows(result_list)

# Output the full_results
saveRDS(full_results, file = "./data/results_setting2_study.rds", compress = T)

#### Study 2 Coeff Level ####################################################################

file_list <- list.files(path = 'data/sim_setting2/coeff',
                        pattern = "*.rds", full.names = T)
result_list <- map(file_list, readRDS)

full_results <- bind_rows(result_list)

# Output the full_results
saveRDS(full_results, file = "./data/results_setting2_coeff.rds", compress = T)



### Study 3 MSE ####################################################################
file_list <- list.files(path = 'data/sim_setting3/mse',
                        pattern = "*.rds", full.names = T)
result_list <- map(file_list, readRDS)

full_results <- bind_rows(result_list)

# Output the full_results
saveRDS(full_results, file = "./data/results_setting3_mse.rds", compress = T)

#### Study3 CI ####################################################################3

file_list <- list.files(path = 'data/sim_setting3/ci',
                        pattern = "*.rds", full.names = T)
result_list <- map(file_list, readRDS)

full_results <- bind_rows(result_list)

# Output the full_results
saveRDS(full_results, file = "./data/results_setting3_ci.rds", compress = T)

