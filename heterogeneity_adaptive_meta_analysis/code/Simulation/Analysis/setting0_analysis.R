#########
# Setting 0 Analysis
########
library(dplyr)
library(tidyr)
library(stringr)
library(xtable)
library(ggplot2)
okabe_ito <- c("#009E73","#CC79A7","#0072B2", "#000000", "#D55E00",
               "#56B4E9","#E69F00",
               "#F0E442")

results0 <- readRDS("./data/results_setting0.rds")

## At each sample size, calc empirical MSE for six settings;
emse <- results0 %>% group_by(n) %>%
  dplyr::summarize(emse_mle = mean(sqerr_mle),
            emse_fe = mean(sqerr_combined_MLE),
            emse_kld = mean(sqerr_kld),
            emse_umse = mean(sqerr_umse),
            emse_biased = mean(sqerr_biased),
            emse_ham = mean(sqerr_ham),
            se_emse_mle = sd(sqerr_mle)/sqrt(1000),
            se_emse_fe = sd(sqerr_combined_MLE)/sqrt(1000),
            se_emse_kld = sd(sqerr_kld)/sqrt(1000),
            se_emse_umse = sd(sqerr_umse)/sqrt(1000),
            se_emse_biased = sd(sqerr_biased)/sqrt(1000),
            se_emse_ham = sd(sqerr_ham)/sqrt(1000)
  )

## Analytic MSE
amse <- results0 %>%  group_by(n) %>%
  dplyr::summarize(mle_median = median(mse_mle)*100,
                  kld_median = median(mse_kld_true)*100,
                  kld_worse = mean(mse_ratio_opt_kld>1)*100,
                  kld_ratio_median = median(mse_ratio_opt_kld),
                  umse_median = median(mse_at_min_umse)*100,
                  umse_worse= mean(mse_ratio_umse>1)*100,
                  umse_ratio_median = median(mse_ratio_umse),
                  biased = median(mse_at_biased)*100,
                  biased_worse = mean(mse_ratio_biased>1)*100,
                  biased_ratio_median = median(mse_ratio_biased),
                  ham_median = median(mse_ham)*100,
                  ham_worse= mean(mse_ratio_ham >1)*100,
                  ham_ratio_median = median(mse_ratio_ham))

amse_all <- results0 %>%
  dplyr::summarize(n = "All",

                   kld_worse = mean(mse_ratio_opt_kld>1)*100,
                   kld_ratio_median = median(mse_ratio_opt_kld),

                   umse_worse= mean(mse_ratio_umse>1)*100,
                   umse_ratio_median = median(mse_ratio_umse),

                   biased_worse = mean(mse_ratio_biased>1)*100,
                   biased_ratio_median = median(mse_ratio_biased),

                   ham_worse= mean(mse_ratio_ham >1)*100,
                   ham_ratio_median = median(mse_ratio_ham))

to_print <- bind_rows(amse, amse_all)


xtable(to_print, digits = c(1, 1,
                            2, 2, 2, 2, 2,
                            2, 2, 2, 2, 2,
                            2, 2, 2))


