#########
# Sim Setting 1: Analysis
# March 2025
########
library(dplyr)
library(tidyr)
library(stringr)
library(xtable)
library(ggplot2)

okabe_ito <- c("#009E73","#CC79A7","#0072B2", "#000000", "#D55E00", "#56B4E9","#E69F00", "#F0E442")

results <- readRDS("./data/results_setting1_mse.rds")
results_ci <- readRDS("./data/results_setting1_ci.rds")



# Print parameter values at each p.
beta_list <- results_ci %>%
  dplyr::distinct(p, beta, study, parameter)


beta_list_wide <- beta_list %>%
  pivot_wider(id_cols = c(p,parameter), names_from = study, values_from = beta)

xtable(beta_list_wide, digits = c(0, 0, 0, 2, 2, 2)) %>% 
  print(include.rownames = F)



## Summarize performance.
performance <- results %>% group_by(n, p) %>%
  dplyr::summarize(mle_emse = 100*mean(sqerr_mle),
                   fe_emse = 100*mean(sqerr_combined_MLE),
                   ham_emse = 100*mean(sqerr_ham),
                   rl_emse  = 100*mean(sqerr_rl)
                   ) %>%
  mutate(ham_ratio = ham_emse/mle_emse,
         rl_ratio = rl_emse/mle_emse) %>%
  pivot_wider(id_cols = n, names_from = p,
              values_from = c(mle_emse,
                              ham_emse,
                              rl_emse, ham_ratio))


performance <- performance%>%
  dplyr::select(n, mle_emse_2,
                   ham_emse_2,
                   rl_emse_2,
                mle_emse_4,
                ham_emse_4,
                rl_emse_4,
                mle_emse_10,
                ham_emse_10,
                rl_emse_10,
                mle_emse_20,
                ham_emse_20,
                rl_emse_20
                )

# print Table 1.
xtable(performance, digits = 1) %>% print(include.rownames = F)

## Coverage by Dataset
ci_dataset_cov <- results_ci %>% group_by(n, p, study) %>%
  dplyr::summarize(coverage_kld = 100*mean(kld_covers),
            coverage_ham = 100*mean(ham_covers),
            coverage_rl = 100*mean(rl_covers))

ci_perf <- ci_dataset_cov %>%
  pivot_wider(id_cols = n, names_from = c(p, study), values_from = c(coverage_ham),
              names_sort = T)  

# Print table 2.
xtable(ci_perf, digits = 1) %>% print(include.rownames = F)
