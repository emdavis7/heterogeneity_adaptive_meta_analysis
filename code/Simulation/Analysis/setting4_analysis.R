#########
# Sim Setting 4: Analysis
########
library(dplyr)
library(tidyr)
library(stringr)
library(xtable)
library(ggplot2)


okabe_ito <- c("#009E73","#CC79A7",  "#D55E00","#000000","#0072B2",
               "#56B4E9","#E69F00", "#F0E442")

rep <- readRDS("./data/results_setting4_rep_level.rds")
rep2 <- readRDS("./data/results_setting4_rep_level2.rds")
coeff <- readRDS("./data/results_setting4_coeff_level.rds")
coeff2 <- readRDS("./data/results_setting4_coeff_level2.rds")


# Add condition info, then combine.
rep$condition <- 1
rep2$condition <- 2
coeff$condition <- 1
coeff2$condition <- 2
rep_all <- bind_rows(rep, rep2)
coeff_all <- bind_rows(coeff, coeff2)

# eMSE
performance <- rep_all %>%
  group_by(condition) %>%
  dplyr::summarise(
    mle_emse = 100*mean(sqerr_mle),
    rema_emse = 100*mean(sqerr_rema),
    ham_emse = 100*mean(sqerr_ham),
    rl_emse  = 100*mean(sqerr_rl),
    ham_ratio2 = mean(mse_ratio_ham),
    rl_ratio2 = mean(mse_ratio_rl),
    var_bias_ratio = mean(tr_var/bias_sq)
  ) %>%
  mutate(rema_ratio = rema_emse/mle_emse,
         ham_ratio = ham_emse/mle_emse,
         rl_ratio = rl_emse/mle_emse )

performance2 <- performance %>% dplyr::select(mle_emse, ham_emse,
                                              rl_emse, var_bias_ratio)
xtable(performance2)
# Coverage

coverage <- data.frame(
  'HAM Coverage' = c(mean(coeff$ham_covers), mean(coeff2$ham_covers)) )

xtable((coverage))

# Make a table:
to_print <- bind_cols(performance2, coverage)
xtable(to_print, digits = c(0, 1, 1, 1, 1, 3))



# Coverage by study:
coverage_by_study <- coeff_all %>% group_by(study, condition) %>%
  summarize(ham_coverage = mean(ham_covers)*100,
            beta_1 = mean(beta_1),
            ham_pi = mean(pi_ham),
            ham_width_ratio = mean(ham_width/mle_width)) %>%
  ungroup() %>% group_by(condition) %>%
  mutate(q1 = quantile(ham_pi, .25),
    q2 = quantile(ham_pi, .5),
    q3 = quantile(ham_pi, .75),
    q4 = quantile(ham_pi, 1),
    ham_pi_quartile = case_when(ham_pi <= q1 ~ "Q1",
                                ham_pi <=q2~ "Q2",
                                ham_pi <= q3~ "Q3",
                                ham_pi <= q4 ~ "Q4",
                                T~"ERR"))
coverage_by_study$ham_pi_quartile <- factor(coverage_by_study$ham_pi_quartile)
coverage_by_study$beta_alt1 <- ifelse(coverage_by_study$beta_1 < .1, coverage_by_study$beta_1, NA)
coverage_by_study$beta_alt2 <- ifelse(coverage_by_study$beta_1 < .1, NA, coverage_by_study$beta_1)

coverage_by_study$condition <- factor(coverage_by_study$condition, levels = c(1,2),
                                      labels = c(expression(paste(beta[j]," from single distribution")),
                                                 expression(paste(beta[j], " from two distributions"))))


ggplot(data = coverage_by_study, aes(x = beta_1, y = ham_coverage,
                                     color = ham_pi_quartile))+
  geom_point(size = 2)+
  facet_wrap(~condition, scales = "free", labeller = "label_parsed")+
  labs(x = expression(beta[j]), y = "Coverage (95% CI)",
       color = expression(paste(pi["j,HAM"], " Quartile"))) +
  scale_color_viridis_d()+
  theme_bw() +
  theme(# legend.position = "bottom",
        strip.background = element_rect(fill = "white"))
ggsave("images/paper_figures/sim3_beta1.png", width = 5.2, height =2)





