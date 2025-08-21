
library("eRm")
library("ltm")
library(lme4)
library(nlme)
library(lmerTest)
library(mirt)
library(ggmirt)
library(ggpubr)
library(IRTest)
library(robustlmm)
library(oaxaca)


## Fitting 2-parameter IRT model to quality measures
states <- c("ak", "de", "nv", "nm", "nc", "ok")

for (s in 1:6){
  
  dat <- get(load(file = paste0(".../rawdata_", states[s],".RData")))
  
  # Davidian Curve 
  Mod_DC <- IRTest_Dich(data = dat[,c(1:3,5:7)],
                        model = 2,
                        latent_dist = "DC",
                        h=4) 
  
  fscore <- as.data.frame(factor_score(Mod_DC))
  scores<-fscore$theta
  
  dat$z1 <- scores
  mean_z1<- mean(scores)
  sd_z1 <- sd(scores)
  dat$z1_std <- (dat$z1 - mean_z1)/sd_z1
  
  #loglinear smoothing
  Mod_LLS <- IRTest_Dich(data = dat[,c(1:3,5:7)],
                         model = 2,
                         latent_dist = "LLS",
                         h=4) 
  
  fscore <- as.data.frame(factor_score(Mod_DC))
  scores<-fscore$theta
  
  dat$z1 <- scores
  mean_z1<- mean(scores)
  sd_z1 <- sd(scores)
  dat$z1_std <- (dat$z1 - mean_z1)/sd_z1
  
  #Kernel density estimation  
  Mod_KDE <- IRTest_Dich(data = dat[,c(1:3,5:7)],
                         model = 2,
                         latent_dist = "KDE",
                         h=4) 
  fscore <- as.data.frame(factor_score(Mod_KDE))
  scores<-fscore$theta
  
  dat$z1 <- scores
  mean_z1<- mean(scores)
  sd_z1 <- sd(scores)
  dat$z1_std <- (dat$z1 - mean_z1)/sd_z1
  
  #Two-component Gaussian mixture distribution  
  Mod_2NM <- IRTest_Dich(data = dat[,c(1:3,5:7)],
                         model = 2,
                         latent_dist = "2NM",
                         h=4) 
  fscore <- as.data.frame(factor_score(Mod_2NM))
  scores<-fscore$theta
  
  dat$z1 <- scores
  mean_z1<- mean(scores)
  sd_z1 <- sd(scores)
  dat$z1_std <- (dat$z1 - mean_z1)/sd_z1
  
  # check BIC and choose the best fitting IRT model for each state. Merge z1_std with the raw data. 
}


## Fitting robust mixed effect model to measure gaps in quality of dental care with the stanadardized score as the outcome

states <- c("ak", "de", "nv", "nm", "nc", "ok")

for (s in 1:6){
  
  # load data with standardized quality score merged to the data
  model.data <- get(load(file = paste0("../merged_", states[s], ".RData"))) 
  
  robust <- rlmer(z1_std~ race + year + race*year + age_group + sex + high_risk+disability_ + resp_prob + cancer + neuro_prob + 
                    musculo_prob + dvlpmt_prob + infectious_dis + eating_dis + obesity + mental_dis +(1 + year|bene_id), 
                  data = model.data)
  
  save(robust, file = paste0("../reg_output/mixed_robust_", states[s],".RData"))
  
}

## Decomposition analysis 

states <- c("ak", "de", "nv", "nm", "nc", "ok")
state_names <- c("Alaska", "Delaware", "Nevada", "New Mexico", "North Carolina", "Oklahoma")
races <- c("Black", "Hispanic", "Asian", "Other")

for (s in 1:6){
  
  for (r in 1:4){
  # load 2019 data with potential mediators merged
  model.data <- get(load(file = paste0("../merged_data/decomp_", states[s], ".RData")))
  
  ICE_keep <- ifelse(r == 1, "ICE_black", ifelse(r == 2, "ICE_hisp", ifelse(r == 3, "ICE_asian", "ICE_other")))
  dat <- model.data %>% 
    #only look at one race pair at a time
    filter(race %in% c("White", races[r])) %>%
    #only include the variables that you will be using as covariates and mediators
    dplyr::select(colnames(model.data)[colnames(model.data) %in% c("z1_std", "race", "age_", "sex_", "high_risk",
                                                                   "disability_", "resp_prob", "cancer", "neuro_prob", "musculo_prob",
                                                                   "dvlpmt_prob", "infectious_dis", "eating_dis", "obesity", "mental_dis",
                                                                   "pov_rate", "prop_accept_medicaid", "dent_10k", "overall_ranking_SVI",
                                                                   "FQHC", ICE_keep)]) %>%
    dplyr::select(z1_std, race, everything())
  
  dat$race_b<- ifelse(dat$race==races[r],0,1)
  results <- oaxaca(formula = z1_std ~ age_+sex_+high_risk+disability_ + resp_prob + cancer + neuro_prob + musculo_prob + dvlpmt_prob + 
                      infectious_dis +eating_dis+obesity+ mental_dis+overall_ranking_SVI+pov_rate+prop_accept_medicaid+
                      dent_10k+FQHC+ICE_black | race_b , data = dat, R = 100)
  }
}
