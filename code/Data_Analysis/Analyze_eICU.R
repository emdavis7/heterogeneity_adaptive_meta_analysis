##############
## Generate Estimates with eICU data
##############
gc()

library(dplyr)
library(janitor)
library(tidyr)
library(ggplot2)
library(readr)
library(Matrix)
library(nloptr)
library(gridExtra)
library(cowplot)
library(mvmeta)
library(xtable)
library(stringr)
library(ggokabeito)


## Source functions
source("./code/functions.R")

## Get reduced data
all_data <- readRDS("./reduced_eICU.rds")
hospital <- read.csv("./hospital.csv")
p <- 7

# Add log transform.
all_data$log_los <- log(all_data$actualiculos)


##### Filtering steps: West, discharge Alive, and hospital > 100
data <- all_data %>% filter(region.y=="West")
data <- data %>% filter(hospitaldischargestatus=="Alive")

check_size <- data %>% group_by(hospitalid) %>% 
  dplyr::summarize(count = n())
to_remove <- check_size$hospitalid[check_size$count < 100]
data2 <- data %>% filter(!(hospitalid %in% to_remove))

# count remaining hospitals:
length(unique(data2$hospitalid))
# 29 remaining hospitals in West region with over 100 patients discharged alive.






##### Create necessary matrices for analysis. #####
# List by hospital:
hosp_list <- sort(unique(data2$hospitalid))

# Apply rescaling to data here.
X_list <- lapply(1:length(hosp_list),
                 function(l){
                   hosp <- hosp_list[l]
                   data_sel <- data2 %>% filter(hospitalid==hosp)
                   # Standardize apache score, age, dbp, and sbp.
                   as.matrix(data.frame(int = 1,
                    gender= (data_sel$gender.x=="Male"),
                    admit_from_ED = ifelse(data_sel$hospitaladmitsource == "Emergency Department", 1, 0),
                    apache = scale(data_sel$apachescore, center =T, scale = T),
                    age = scale(as.numeric(data_sel$age.x), center=T, scale = T),
                    dbp = scale(data_sel$noninvasivediastolic, center = T, scale = T),
                    sbp = scale(data_sel$noninvasivesystolic, center = T, scale = T)))
                 })

X <- do.call(rbind, X_list) 

Y_list <-  lapply(1:length(hosp_list),
                  function(l){
                    hosp <- hosp_list[l]
                    data_sel <- data2 %>% filter(hospitalid==hosp)
                    as.numeric(data_sel$log_los)
                  })

# MLE by hospital (on scaled X's):
hosp_mle <- lapply(1:length(hosp_list), function(l){
  hosp <- hosp_list[l]
  X <- X_list[[l]]
  Y <- Y_list[[l]]
  t( solve(t(X) %*% X) %*% t(X) %*% Y )
})


### Estimate sd to generate estimated versions of scaled xtx and xty.
est_sd <-sapply(1:length(hosp_list), function(l){
  y_j <- Y_list[[l]]
  X_j <- X_list[[l]]
  n <- length(y_j)
  mle_j <- hosp_mle[[l]]
  residual_j <- y_j - X_j %*% t(mle_j)
  sse_j <- sum(residual_j^2)
  sqrt(sse_j/(n-p))
})


## Get scaled XTX as list:
scaled_XTX_est_list <- lapply(1:length(hosp_list), function(l){
  1/est_sd[l]^2*
    t(X_list[[l]])%*%
    X_list[[l]]
})
scaled_XTX_est <- as.matrix(bdiag(scaled_XTX_est_list))


### Scaled XTX list is in same order as hosp list.
scaled_XTY_est_list <- lapply(1:length(hosp_list), function(l){
  1/est_sd[l]^2*
    t(X_list[[l]])%*%
    Y_list[[l]]
})
scaled_XTY_est <- unlist(scaled_XTY_est_list)

# Vector version of hosp_mle.
indv_mle <- solve(scaled_XTX_est) %*% scaled_XTY_est

# Make this a dataframe so it's plottable.
mle_df <- do.call(rbind.data.frame, hosp_mle)
mle_df$hospitalid <- hosp_list
mle_df2  <- mle_df %>% left_join(hospital) %>% left_join(check_size)


#### Make 95% CIs for MLE: ####
mle_long <- mle_df %>%
  pivot_longer(cols = c(1:p), 
               names_to = "Covariate", values_to = "Estimate")
mle_long$est_type <- "OLS"
mle_long$std_err <- sqrt(diag(solve(scaled_XTX_est)))
mle_long <- mle_long %>% mutate(lb = Estimate - 1.96*std_err,
                                ub = Estimate + 1.96*std_err)
mle_long <- mle_long %>% mutate(covers_zero = ifelse(lb <= 0 & ub>= 0, 1, 0))


######## Get shrinkage parameter, pi #######################
k <- length(hosp_list)
kron <-  kronecker(rep(1, k), diag(nrow = p))

# combined MLE:
fixed_effect_mle <- solve(t(kron) %*% scaled_XTX_est %*% kron) %*% 
  t(kron) %*% scaled_XTY_est


# Get starting value based on est pi star;
var_part <-  tr( t(kron) %*% solve(scaled_XTX_est) %*% kron -
                   k * solve(t(kron) %*% scaled_XTX_est %*% kron) )
bias_part_est <- t( kron %*% fixed_effect_mle - indv_mle) %*%
  (kron %*% fixed_effect_mle - indv_mle)
pi_star_est <- as.numeric( var_part / (var_part + bias_part_est) )
# pi_star_est = 0.34, which is relatively high.

## This will take a minute or two to run:
pi <- bobyqa(rep(pi_star_est, k),
               fn = mse_hat,
               scaled_XTX = scaled_XTX_est,
               indv_mle = indv_mle, control = list(xtol_rel = 1e-8, maxeval = 2000),
               p = p, D = k, kron = kron,
               lower = rep(0, k),
               upper = rep(1, k))

bPi <- blockLambda_fcn(pi$par, p = 7, D = k)


#### With the selected shrinkage parameters, get beta_ham ######
ham_ests <- kld_linear_alt(pi = pi$par,
                           scaled_XTY = scaled_XTY_est,
                           scaled_XTX = scaled_XTX_est,
                           indv_mle = indv_mle,
                           p=p,
                           D=k,
                           kron=kron)



###### Combined Estimators ############
# Theta + CI:
hbtheta <- ham_ests$theta_hat
hbtheta_var <- solve(t(kron) %*% bPi %*% scaled_XTX_est %*% kron) %*%
  t(kron) %*% bPi %*% scaled_XTX_est %*% bPi %*% kron %*%
  solve(t(kron) %*% bPi %*% scaled_XTX_est %*% kron)
hbtheta_lb <- sapply(1:length(hbtheta), function(l){hbtheta[l] - 1.96*sqrt(hbtheta_var[l,l])})
hbtheta_ub <- sapply(1:length(hbtheta), function(l){hbtheta[l] + 1.96*sqrt(hbtheta_var[l,l])})

# FE + CI:
fe_var <- solve(t(kron) %*% scaled_XTX_est %*% kron)
fe_lb <- sapply(1:length(fixed_effect_mle), function(l){fixed_effect_mle[l] - 1.96*sqrt(fe_var[l,l])})
fe_ub <- sapply(1:length(fixed_effect_mle), function(l){fixed_effect_mle[l] + 1.96*sqrt(fe_var[l,l])})

# REMA + CI
re_ma <- mvmeta(as.matrix(mle_df[1:p]),
                S = lapply(1:length(scaled_XTX_est_list), function(l){
                  solve(scaled_XTX_est_list[[l]])
                }))
re_ma_lb <- sapply(1:p, function(l){re_ma$coefficients[l] - 1.96*sqrt(re_ma$vcov[l,l])})
re_ma_ub <- sapply(1:p, function(l){re_ma$coefficients[l] + 1.96*sqrt(re_ma$vcov[l,l])})

# Coefficient plot for combined estimate;
combined_df <- data.frame(covariate = colnames(mle_df)[1:p],
                         fe_lb = fe_lb, fe_ub = fe_ub,
                         re_ma_lb = re_ma_lb, re_ma_ub = re_ma_ub,
                         hbtheta_lb = hbtheta_lb, hbtheta_ub = hbtheta_ub)



# print Table 6
xtable(combined_df, digits = c(0,1,  3, 3, 3, 3, 3, 3),
       caption = "95% Confidence Intervals for Combined Estimators",
       label = "tab:combined_ci" )%>% print(include.rownames = F)


## Get I^2 from Random-Effects summary
summary(re_ma)

########################################################
#### Create asymptotic confidence intervals from beta_ham #####
ham_aCI <- kld_linear_aCI(pi = pi$par,
                           scaled_XTX =  scaled_XTX_est,
                           indv_mle = indv_mle,
                           p = p, D = k, kron = kron,
                           alpha = .05)

ham_lb <- ham_aCI$lb
ham_ub <- ham_aCI$ub

# Create ham dataframe
long_hosp <- c(sapply(1:length(hosp_list), function(l){
  rep(hosp_list[l], p)
}))
ham_long <- data.frame(hospitalid = long_hosp,
                       Covariate = c("int", "gender", "admit_from_ED",
                                     "apache", "age", "dbp", "sbp"),
                       Estimate = ham_ests$beta_hat,
                       est_type = "HAM",
                       lb = ham_lb,
                       ub = ham_ub)
ham_long$covers_zero <- ifelse(ham_long$lb <=0 & ham_long$ub>=0, 1, 0)

# Merge on hospital and covariate.
test_comparison <- ham_long %>%
  dplyr::select(hospitalid, Covariate, covers_zero, lb, ub) %>%
  rename(covers_zero_ham = covers_zero,
         lb_ham = lb,
         ub_ham = ub) %>%
  left_join(dplyr::select(mle_long, c(hospitalid, Covariate, covers_zero, lb, ub)))


# Create table S4:
to_print <- test_comparison %>% group_by(Covariate) %>%
  dplyr::summarize(both_sig = sum(covers_zero==1 & covers_zero_ham==1),
            ham_sig = sum(covers_zero==0 & covers_zero_ham==1),
            mle_sig = sum(covers_zero==1 & covers_zero_ham==0),
            neither_sig = sum(covers_zero==0 & covers_zero_ham==0))
xtable(to_print) %>% print(include.rownames = F)





#### Make forest plots #####
ham_ests_mat <- as.data.frame(t(matrix(ham_ests$beta_hat, ncol = k)))
colnames(ham_ests_mat) <- c("int", "gender", "admit_from_ED", "apache", 
                            "age", "dbp", "sbp")

ham_ub_mat <- as.data.frame(t(matrix(ham_ub, ncol = k)))
colnames(ham_ub_mat) <- c("int_ub", "gender_ub", "admit_from_ED_ub", 
                          "apache_ub", "age_ub", "dbp_ub", "sbp_ub")

ham_lb_mat <- as.data.frame(t(matrix(ham_lb, ncol = k)))
colnames(ham_lb_mat) <- c("int_lb", "gender_lb", "admit_from_ED_lb", 
                          "apache_lb", "age_lb", "dbp_lb", "sbp_lb")

ham <- cbind(ham_ests_mat, ham_ub_mat, ham_lb_mat)
ham$hospital <- as.factor(hosp_list)

# Format hospital sizes for graph.
smp <- data2 %>% group_by(hospitalid) %>%
  dplyr::summarize(count = n()) %>%
  mutate(hospital = as.factor(hospitalid)) %>%
  dplyr::select(-hospitalid)


###### Add confidence intervals to MLE data frame. ###
mat_ub <- matrix(nrow = k, ncol = p)
mat_lb <- matrix(nrow = k, ncol = p)
for (i in 1:p){
  est_sel <- colnames(mle_df2)[i]

  mle_est <- mle_df2[[eval(est_sel)]]
  mle_ub <- sapply(1:k, function(l){
    est_var <- solve(scaled_XTX_est_list[[l]])
    est_se_admit <-  sqrt(as.numeric(est_var[i,i]))
    mle_est[l] + 1.96*est_se_admit
  })
  mle_lb <- sapply(1:k, function(l){
    est_var <- solve(scaled_XTX_est_list[[l]])
    est_se_admit <- sqrt(as.numeric(est_var[i,i]))
    mle_est[l] - 1.96*est_se_admit
  })
  mat_ub[,i] <- mle_ub
  mat_lb[,i] <- mle_lb

}
# assign colnames
colnames(mat_ub) <- paste0(colnames(mle_df[,1:7]), "_ub")
colnames(mat_lb) <- paste0(colnames(mle_df[,1:7]), "_lb")
# make df in same shape as ham.
mle_to_append <- bind_cols(mle_df2[, 1:7], mat_ub, mat_lb,
                           hospital = as.factor(mle_df2$hospitalid))

all_ests <- bind_rows(ham, mle_to_append)
all_ests$type <- c(rep("HAM", k), rep("MLE", k))
all_ests$pi <- rep(pi$par, 2)

# Join smp;
all_ests <- all_ests %>% left_join(smp)


# Add the CIs from FE, RE, and HBtheta.
combined_point <- rbind(as.numeric(fixed_effect_mle), re_ma$coefficients,
                        as.numeric(hbtheta))
combined_lbs <-rbind(fe_lb, re_ma_lb, hbtheta_lb)
combined_ubs <- rbind(fe_ub, re_ma_ub, hbtheta_ub)
colnames(combined_point) <- colnames(mle_df[,1:7])
colnames(combined_lbs) <- paste0(colnames(mle_df[,1:7]), "_lb")
colnames(combined_ubs) <- paste0(colnames(mle_df[,1:7]), "_ub")
add_meta <- data.frame(type = c("Fixed-Effect", "Random-Effect", "HAM"),
                       hospital = as.factor("Meta-Estimator")) %>%
  bind_cols(combined_point, combined_lbs, combined_ubs)
all_ests <- all_ests %>% bind_rows(add_meta)

## Get hospital list in order of pi.
order <- all_ests %>% dplyr::select(pi, hospital, type) %>% filter(type=="MLE") %>%
  arrange(pi) %>%
  mutate(order = as.ordered(row_number())) %>%
  dplyr::select(hospital, order)

all_ests <- all_ests %>% left_join(order)
all_ests <- all_ests %>% mutate(order = ifelse(is.na(order), k+1, order))


## Panel estimates:
# make long.
ests_long <- all_ests %>% dplyr::select(type, hospital, order, pi,
                                            colnames(mle_df2)[1:7]) %>%
  pivot_longer(cols = c(colnames(mle_df2)[1:7]), names_to = "covariate",
               values_to = "point_est")

ub_long <- all_ests %>% dplyr::select(type, hospital, order,
                                      paste0(colnames(mle_df2)[1:7], "_ub")) %>%
  pivot_longer(cols = c(paste0(colnames(mle_df2)[1:7], "_ub")), names_to = "covariate",
               values_to = "ub") %>%
  mutate(covariate = str_replace(covariate, "_ub", ""))

lb_long <- all_ests %>% dplyr::select(type, hospital, order,
                                      paste0(colnames(mle_df2)[1:7], "_lb")) %>%
  pivot_longer(cols = c(paste0(colnames(mle_df2)[1:7], "_lb")), names_to = "covariate",
               values_to = "lb")%>%
  mutate(covariate = str_replace(covariate, "_lb", ""))

all_ests_long <- ests_long %>% left_join(ub_long) %>% left_join(lb_long)

all_ests_long$type_f <- factor(all_ests_long$type, ordered = T,
                               levels = c("HAM", "MLE", "Fixed-Effect",
                                          "Random-Effect"),
                               labels = c("HAM", "MLE", "Fixed-Effect",
                                          "Random-Effect Meta-Estimate"))
all_ests_long$order_f <- factor(all_ests_long$order, ordered = T,
                                labels = c(1:29, "*"))
all_ests_long$covariate_f <- factor(all_ests_long$covariate,
                                    levels = c("int",
                                               "gender",
                                               "admit_from_ED",
                                               "apache",
                                               "age",
                                               "dbp", "sbp"),
                                    labels = c("Intercept",
                                               "Gender",
                                               "Admit from ED",
                                               "APACHE IV Score",
                                               "Age",
                                               "Diastolic BP",
                                               "Systolic BP"))

## Drop fixed-effect and HAM meta-estimator for graphs.
all_ests_long2 <- all_ests_long %>%
  filter(type!="Fixed-Effect" &
        !(type=="HAM" & hospital=="Meta-Estimator"))


### Make a simpler graph: give HAM estimates and color code for agreement.
rema_cover <- all_ests_long2 %>% filter(type =="Random-Effect")%>%
  mutate(rema_cover = (ub >= 0 & lb <=0)) %>%
  dplyr::select(covariate, rema_cover)
mle_cover <- all_ests_long2 %>% filter(type=="MLE") %>%
  mutate(mle_cover = (ub>=0 & lb <=0)) %>%
  dplyr::select(hospital, covariate, mle_cover)
ham_cover <- all_ests_long2 %>% filter(type=="HAM") %>%
  mutate(ham_cover = (ub >= 0 & lb<=0)) %>%
  left_join(mle_cover) %>% left_join(rema_cover)

ham_cover <- ham_cover %>%
  mutate(summary = case_when(ham_cover == mle_cover & ham_cover == rema_cover ~ "Agreement in HAM, MLE, and RE",
                             ham_cover == mle_cover & ham_cover != rema_cover ~ "Agreement in HAM and MLE; RE Differs",
                             ham_cover != mle_cover & ham_cover == rema_cover ~ "Agreement in HAM and RE; MLE Differs",
                             ham_cover != mle_cover & ham_cover != rema_cover ~ "HAM Differs from MLE and RE",
                             T ~ "ERR"),
         summary2 = case_when(ham_cover == F & mle_cover == F  ~ "Both HAM and MLE Significant",
                             ham_cover == T & mle_cover == F ~ "HAM Not Significant; MLE Significant",
                             ham_cover == F & mle_cover == T ~ "HAM Significant; MLE Not Significant",
                             ham_cover == T & mle_cover == T ~ "Both HAM and MLE Not Significant",
                             T ~ "ERR"))
ham_cover$summary <- factor(ham_cover$summary,
                            levels = c("Agreement in HAM, MLE, and RE",
                                       "Agreement in HAM and MLE; RE Differs",
                                       "Agreement in HAM and RE; MLE Differs",
                                       "HAM Differs from MLE and RE")
)
ham_cover$summary2 <- factor(ham_cover$summary2,
                            levels = c("Both HAM and MLE Significant",
                                       "HAM Not Significant; MLE Significant",
                                       "HAM Significant; MLE Not Significant",
                                       "Both HAM and MLE Not Significant")
)

### Create figure S3
ggplot(data = ham_cover) +
  geom_point(aes(y = point_est, x = order_f, color = summary2),
             size = 1.8)+
  geom_linerange(aes(x = order_f, ymin =  lb, ymax =  ub, color = summary2),
                 linewidth = 1.5) +
  geom_hline(yintercept = 0)+
  scale_x_discrete(breaks = function(x) x[c(T, F)])+
  facet_wrap(~covariate_f, ncol = 1, scales = "free_y")+
  scale_color_okabe_ito(order = c(5, 6, 1, 9)) +
  labs(y = "95% Confidence Intervals", x = "Hospital",
       color = expression(paste("Significant at ",alpha, "=", 0.05))) +
  theme_bw() +
  theme(#axis.text.x = element_blank(), axis.title.x=element_blank(), axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()) +
  guides(color = guide_legend(nrow=2))
ggsave("images/paper_figures/da_coeffs_agreement.png", width = 6.45, height = 8)

ggsave("images/for slides/da_coeffs_agreement.png", width = 5.45, height = 6.5)




########################################
# Create plot for APACHE IV, arranged by ests.
order2 <- all_ests %>% dplyr::select(pi, hospital, type, apache) %>% filter(type=="MLE") %>%
  arrange(apache) %>%
  mutate(order2 = as.ordered(row_number())) %>%
  dplyr::select(hospital, order2)

all_ests <- all_ests %>% left_join(order2)


### Figure 5
gg1 <- ggplot(data = filter(all_ests, is.na(order2)==F)) +
  #coord_flip()+
  geom_point(aes(y = apache, x = order2, color = type), size = 1.3,
             position = position_dodge2(width = .6))+
  geom_linerange(aes(x = order2, ymin =  apache_lb, ymax =  apache_ub, color = type),
                 size = 1.2, position = position_dodge2(width = .6)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c("#CC0000", "#000000")) +
  labs(y = "95% CIs for APACHE IV", x = "Hospital", color = "Estimator") +
  theme_bw() +
  theme(legend.position = "bottom")
gg1
gg2 <- ggplot(data = filter(all_ests, is.na(order2)==F)) +
  geom_point(aes(y = pi, x = order2)) +
  labs(y = expression(paste("Selected Values for ", pi["j,HAM"])), x = "Hospital") +
  theme_bw()
plot3 <- plot_grid(gg1, gg2, ncol = 1, align = "hv")

plot3
ggsave("./images/paper_figures/APACHEIV_95ci.png", plot3,  width = 7, height = 5.5)


