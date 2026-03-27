#########
# Sim Setting 2: Analysis
########
library(dplyr)
library(tidyr)
library(stringr)
library(xtable)
library(ggplot2)
library(scales)

rep <- readRDS("./data/results_setting2_rep.rds")
study <- readRDS("./data/results_setting2_study.rds")
coeff <- readRDS("./data/results_setting2_coeff.rds")

okabe_ito <- c("#009E73","#CC79A7","#0072B2", "#000000", 
               "#D55E00", "#56B4E9","#E69F00", "#F0E442")

## Grab the total number of datasets to add back to the coeff and study
key <- rep %>% dplyr::select(D, seed) %>% rename(K = D)
study <- study %>% left_join(key)
coeff <- coeff %>% left_join(key)

## Summarize replicate level. #####
rep_summary <- rep %>% group_by(condition, D) %>%
  dplyr::summarize(mle_emse = 100*mean(sqerr_mle),
                   mle_se = 100*sd(sqerr_mle)/sqrt(n()),
                   re_emse = 100*mean(sqerr_re_ma),
                   re_se = 100*sd(sqerr_re_ma)/sqrt(n()),
                   fe_emse = 100*mean(sqerr_combined_MLE),
                   fe_se = 100*sd(sqerr_combined_MLE)/sqrt(n()),
                   ham_emse = 100*mean(sqerr_ham),
                   ham_se = 100*sd(sqerr_ham)/sqrt(n()),
                   ham_ratio = mean(mse_ham/mse_mle),
                   rl_emse  = 100*mean(sqerr_rl),
                   rl_se = 100*sd(sqerr_rl)/sqrt(n()),
                   rl_smaller = mean(sqerr_rl < sqerr_ham)
  )

rep %>% group_by(condition, D) %>%
  dplyr::summarize(mean_ham_mse = 100*mean(mse_ham),
                   mean_rl_mse = 100*mean(mse_rl))




 ## Overall CI ####

ci <- coeff %>% group_by(beta_cond, K) %>%
  summarize(ham_covers = mean(ham_aCI_covers)*100) %>%
  rename(condition = beta_cond, D = K)

rep_summary2 <- rep_summary %>% left_join(ci) %>%
  select(condition, D, mle_emse, ham_emse, rl_emse, ham_covers)

# Print table 3
xtable(rep_summary2, digits = c(0, 0, 0, 1, 1, 1, 1),
       caption = "Mean MSEs across conditions",
       label = "tab:sim2_mse" )%>% print(include.rownames = F)


## Pi Analysis #######################
opt_pi <- unlist(rep$opt_pi)
pi_ham <- unlist(rep$pi_ham)
condition <- c()
D <- c()
replicate <- c()
for (i in 1:nrow(rep)){
  condition_tmp <- rep(rep$condition[i], rep$D[i])
  condition <- c(condition, condition_tmp)

  D_tmp <- rep(rep$D[i], rep$D[i])
  D <- c(D, D_tmp)

  replicate_tmp <- rep(rep$seed[i], rep$D[i])
  replicate <- c(replicate, replicate_tmp)
}

pi <- data.frame(replicate = replicate,
                 condition = condition,
                 D = D,
                 euc = study$studywise_eucdist,
                 study = study$D,
                 opt_pi = opt_pi,
                 pi_ham = pi_ham )
pi$cluster <- ifelse(pi$study < 4 & pi$condition == 4, 1,
                     ifelse(pi$condition==4, 0, NA))
pi$cluster <- factor(pi$cluster,
                        labels = c("Heterogeneous",
                                   "Homogeneous"))

pi$condition <- factor(pi$condition,
                       labels = c(expression(paste("All ", beta[j], "\n Equal")),
                                  "'Mild\n Heterogeneity'",
                                  "'Moderate\n Heterogeneity'",
                                  "'Mixture'"))
# within rep, calculate correlation between opt_pi and pi_ham, variance of pi, and IQR
pi_within_rep <- pi %>% group_by(replicate, condition, D) %>% summarize(corr = cor(opt_pi, pi_ham),
                                            variance = var(pi_ham),
                                            iqr = IQR(pi_ham),
                                            range80 = quantile(pi_ham, .9) - quantile(pi_ham, .1),
                                            range = max(pi_ham) - min(pi_ham))


# Mean over reps.
pi_summary <- pi_within_rep %>% group_by(condition, D) %>%
  summarize(mean_corr = mean(corr, na.rm = T),
            mean_variance = mean(variance),
            mean_iqr = mean(iqr),
            mean_range80 = mean(range80),
            mean_range = mean(range))

# Create figure 4
ggplot(data = filter(pi, D==15), aes(x = pi_ham, fill = cluster, group = cluster)) +
  facet_grid(cols = vars(condition),  labeller = label_parsed) +
  geom_histogram(alpha = .7, position = "identity") +
  scale_fill_manual(values = okabe_ito,
                    limits = c("Heterogeneous",
                               "Homogeneous"))+
  theme_bw() +
  labs(x = expression(pi["HAM"]), y = "Count", fill = "Study is:") +
  theme(#legend.position="bottom",
        strip.background = element_rect(fill = "white"),
        axis.text = element_text(size = 8) )
ggsave("images/paper_figures/sim2_pi.png", width = 7.5, height = 2)










