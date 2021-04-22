# Code author: AK
# Community spillover effects
# on the spread of covid-19
# updated data 

library(WeightIt)
library(readstata13)

set.seed(1234)

data <- read.dta13('covid_2_24_21_cross_section_cumulative_no_AK_for_causal.dta')
data <- na.omit(data)

# standardization: only applied to a few columns
data$p65_over <- (data$p65_over - mean(data$p65_over))/(2*sd(data$p65_over))
data$ppopurban <- (data$ppopurban - mean(data$ppopurban))/(2*sd(data$ppopurban))
data$disadvantage <- (data$disadvantage - mean(data$disadvantage))/(2*sd(data$disadvantage))
data$n_covid <- (data$n_covid - mean(data$n_covid))/(2*sd(data$n_covid))
data$w_covid <- (data$w_covid - mean(data$w_covid))/(2*sd(data$w_covid))

# binary indicators
# data$abavg_nhb <- data$abavg_nhb - mean(data$abavg_nhb)
# data$abavg_hisp <- data$abavg_hisp - mean(data$abavg_hisp)
# data$abavg_nhw <- data$abavg_nhw - mean(data$abavg_nhw)

# response
#CHANGE BETWEEN DEATHS AND CASES TO SWITCH OUTCOMES
data$deaths <- (data$deaths/data$total_pop)*1000

response <- 'deaths'
all_variables <- colnames(data)
features <- all_variables[all_variables != 'dis_th']
features <- features[features != 'STATEFP']
features <- features[features != 'deaths']
features <- features[features != 'cases']
features <- features[features != 'geoid']
features <- features[features != 'total_pop']
features <- features[features != 'n_delta']
features <- features[features != 'w_delta']
features <- features[features != 'tl_covid']
features <- features[features != 'case_rate']
features <- features[features != 'death_rate']

# result matrix
super_result_standardized <- data.frame()
cbps_result_standardized <- data.frame()
iptw_result_standardized <- data.frame()

# each time consider a feature as treatment
for (treat in features){
  trt_model_covariates <- features[features != treat]
  # Super Learner
  weight_model_super <- WeightIt::weightit(as.formula(paste0(paste0(treat, '~'),paste0(trt_model_covariates, collapse= "+"))),
                                           method = 'super',stop.method = 's.mean',
                                           SL.library = c('SL.gam'), data = data)

  weight_model_super <- trim(weight_model_super, 0.97)

  outcome_model_super_standardized <- glm(as.formula(paste0(paste0(response,'~',treat,'+'),paste0(trt_model_covariates, collapse= "+"))),
                                weights = weight_model_super$weights, data = data)


  super_result_standardized <- rbind (super_result_standardized, summary(outcome_model_super_standardized)$coefficients[2,])

  # CBPS
  weight_model_cbps <- WeightIt::weightit(as.formula(paste0(paste0(treat, '~'),paste0(trt_model_covariates, collapse= "+"))),
                                          method = 'cbps', data = data, over = FALSE)

  weight_model_cbps <- trim(weight_model_cbps, 0.97)

  outcome_model_cbps_standardized <- glm(as.formula(paste0(paste0(response,'~',treat,'+'),paste0(trt_model_covariates, collapse= "+"))),
                            weights = weight_model_cbps$weights, data = data)

  cbps_result_standardized <- rbind(cbps_result_standardized, summary(outcome_model_cbps_standardized)$coefficients[2,])

  # # IPTW
  weight_model_iptw <- WeightIt::weightit(as.formula(paste0(paste0(treat, '~'),paste0(trt_model_covariates, collapse= "+"))),
                                          method = 'ps', data = data, over = FALSE)

  weight_model_iptw <- trim(weight_model_iptw, 0.97)

  outcome_model_iptw_standardized <- glm(as.formula(paste0(paste0(response,'~',treat,'+'),paste0(trt_model_covariates, collapse= "+"))),
                               weights = weight_model_iptw$weights, data = data, maxit = 1000)

  iptw_result_standardized <- rbind(iptw_result_standardized, summary(outcome_model_iptw_standardized)$coefficients[2,])

} # end of for loop

# write the output
#UPDATE FILE NAMES TO SWITCH BETWEEN OUTCOMES
write.csv(super_result_standardized, file = 'deaths_super.csv',
          row.names = c("p65_over","ppopurban","disadvantage","abavg_nhw","abavg_nhb","abavg_hisp","n_covid","w_covid"))
write.csv(cbps_result_standardized, file = 'deaths_cbps.csv',
          row.names = c("p65_over","ppopurban","disadvantage","abavg_nhw","abavg_nhb","abavg_hisp","n_covid","w_covid"))
write.csv(iptw_result_standardized, file = 'deaths_iptw.csv',
          row.names = c("p65_over","ppopurban","disadvantage","abavg_nhw","abavg_nhb","abavg_hisp","n_covid","w_covid"))



