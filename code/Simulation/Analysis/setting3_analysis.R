#########
# Sim Setting 4: Analysis
########
library(dplyr)
library(tidyr)
library(stringr)
library(xtable)
library(ggplot2)

okabe_ito <- c("#009E73","#CC79A7","#0072B2", "#000000", "#D55E00", 
               "#56B4E9","#E69F00", "#F0E442")

results <- readRDS("./data/results_setting3_mse.rds")
results_ci <- readRDS("./data/results_setting3_ci.rds")


performance_long <- results %>%
  mutate( kld_improvement = sqerr_kld < sqerr_mle,
    ham_improvement = sqerr_ham < sqerr_mle,
    rl_improvement = sqerr_rl < sqerr_mle) %>%
  group_by(design, data_mu, sd, n) %>%
  dplyr::summarize(mle_emse = 100*mean(sqerr_mle),
                   mle_emse_rs = 100*mean(sqerr_mle_rs), 
                   fe_emse = 100*mean(sqerr_combined_MLE), 
                   kld_emse = 100*mean(sqerr_kld),
                   ham_emse = 100*mean(sqerr_ham), 
                   rl_emse  = 100*mean(sqerr_rl),
                   kld_emse_rs = 100*mean(sqerr_kld_rs),
                   ham_emse_rs = 100*mean(sqerr_ham_rs),
                   kld_improvement = mean(kld_improvement),
                   ham_improvement  = mean(ham_improvement),
                   rl_improvement = mean(rl_improvement) 
  ) %>%
  mutate(kld_ratio = kld_emse/mle_emse,
         ham_ratio = ham_emse/mle_emse,
         rl_ratio = rl_emse/mle_emse,
         kld_ratio_rs = kld_emse_rs/mle_emse_rs,
         ham_ratio_rs = ham_emse_rs/mle_emse_rs)


perf_wide <- performance_long %>% arrange(design, n) %>% ungroup() %>%
  dplyr::select(design, mle_emse, ham_emse, rl_emse, mle_emse_rs, ham_emse_rs, n) %>%
  filter(n <= 500) %>%
  pivot_wider(names_from = design, values_from = c(mle_emse, ham_emse, rl_emse, mle_emse_rs, ham_emse_rs))

## Coverage overall
ci_cov <- results_ci %>%  group_by(design, n)  %>%
  filter(n <= 500) %>%
  dplyr::summarize(coverage_ham = 100*mean(ham_covers),
                   coverage_ham_rs = 100*mean(ham_covers_rs)) %>%
  pivot_wider(names_from = design, values_from = c(coverage_ham, coverage_ham_rs) )

## Join to performance and order.
reporting <- perf_wide %>% full_join(ci_cov) %>%
  dplyr::select(n ,
                mle_emse_1, ham_emse_1, rl_emse_1, coverage_ham_1, mle_emse_rs_1, ham_emse_rs_1, coverage_ham_rs_1,
                mle_emse_3, ham_emse_3, rl_emse_3, coverage_ham_3, mle_emse_rs_3, ham_emse_rs_3, coverage_ham_rs_3,
                mle_emse_4, ham_emse_4, rl_emse_4, coverage_ham_4, mle_emse_rs_4, ham_emse_rs_4, coverage_ham_rs_4,
                mle_emse_5, ham_emse_5, rl_emse_5, coverage_ham_5, mle_emse_rs_5, ham_emse_rs_5, coverage_ham_rs_5)


report1 <- reporting %>% dplyr::select(n,
                                       mle_emse_1, ham_emse_1, rl_emse_1, coverage_ham_1,
                                       mle_emse_3, ham_emse_3, rl_emse_3, coverage_ham_3,
                                       mle_emse_4, ham_emse_4, rl_emse_4, coverage_ham_4,
                                       mle_emse_5, ham_emse_5, rl_emse_5, coverage_ham_5)
report2 <- reporting %>% dplyr::select(n,
                                       mle_emse_rs_1, ham_emse_rs_1, coverage_ham_rs_1,
                                       mle_emse_rs_3, ham_emse_rs_3, coverage_ham_rs_3,
                                       mle_emse_rs_4, ham_emse_rs_4, coverage_ham_rs_4,
                                       mle_emse_rs_5, ham_emse_rs_5, coverage_ham_rs_5)

# Print Table 4
xtable(report1, caption = "eMSE and Coverage for Simulation Study 4",
       digits = c(0, 0,
                  1, 1, 1, 1,
                  1, 1, 1, 1,
                  1, 1, 1, 1,
                  1, 1, 1, 1),
       label = "tab:sim4_results1") %>% print(include.rownames = F)

# Print Table 5
xtable(report2, caption = "eMSE, Ratio, and Coverage for Simulation Study 4",
       digits = c(0, 0,
                  1, 1, 1,
                  1, 1, 1,
                  1, 1, 1,
                  1, 1, 1),
       label = "tab:sim4_results2") %>% print(include.rownames = F)

### Visualize MLE CIs across Studies ####
mle_ci <- results_ci %>% group_by(n, design, sd, data_mu, study, parameter) %>%
  summarize(mean_lb = mean(mle_aCI_lb),
            mean_ub = mean(mle_aCI_ub),
            median_lb = median(mle_aCI_lb),
            median_ub = median(mle_aCI_ub),
            var_lb = var(mle_aCI_lb),
            var_ub = var(mle_aCI_ub),
            mean_sd = mean((mle_aCI_ub - mle_aCI_lb)/1.96),
            mean_coverage = mean(mle_covers))

# plot lineranges, group by study, panel by design and parm, restrict to n = 500.
mle_ci$design_f <- factor(mle_ci$design, levels =c(1, 3, 4, 5),
                          labels = c(
                            expression(paste(mu,"=(0, 0, 0)")),
                            expression(paste(mu,"=(10, 10, 10)")),
                            expression(paste(mu,"=(0, 1, 2)")),
                            expression(paste(mu,"=(0, 5, 10)"))
                          ))
mle_ci$parm_f <- factor(mle_ci$parameter,
                        levels = c(1, 2, 3, 4),
                        labels = c(
                          expression(paste(beta[1])),
                          expression(paste(beta[2])),
                          expression(paste(beta[3])),
                          expression(paste(beta[4]))
                        ))
# Figure S1
ggplot(data = filter(mle_ci, n == 500 & design != 2),
       aes(ymin = mean_lb, ymax = mean_ub, x = as.factor(study)))+
  geom_linerange()+
  facet_grid(cols = vars(design_f), rows = vars(parm_f), scales = "free",
             labeller = "label_parsed") +
  labs(x = "Study", y = "Average 95% Confidence Interval from MLE") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"))


ggsave("images/paper_figures/sim4_mleci.png", width = 5.3, height = 3.8)


