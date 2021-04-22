# Code author: AK
# Community spillover effects
# on the spread of covid-19
# updated data 

library(WeightIt)
library(readstata13)

set.seed(1234)

data <- read.dta13('covid_2_24_21_cross_section_TWO_TIME_cumulative_no_AK_for_causal.dta')
data <- na.omit(data)

# standardization: only applied to a few columns
data$p65_over <- (data$p65_over - mean(data$p65_over))/(2*sd(data$p65_over))
data$ppopurban <- (data$ppopurban - mean(data$ppopurban))/(2*sd(data$ppopurban))
data$disadvantage <- (data$disadvantage - mean(data$disadvantage))/(2*sd(data$disadvantage))
data$n_covid <- (data$n_covid - mean(data$n_covid))/(2*sd(data$n_covid))
data$w_covid <- (data$w_covid - mean(data$w_covid))/(2*sd(data$w_covid))

# response
#SWITCH BETWEEN CASES2 AND DEATHS2 TO CHANGE OUTCOMES
data$deaths2 <- (data$deaths2/data$total_pop)*1000

response <- 'deaths2'
all_variables <- colnames(data)
features <- all_variables[all_variables != 'dis_th']
features <- features[features != 'STATEFP']
features <- features[features != 'deaths2']
features <- features[features != 'cases2']
features <- features[features != 'geoid']
features <- features[features != 'total_pop']
features <- features[features != 'case_rate']
features <- features[features != 'death_rate']
features <- features[features != 'n_delta']
features <- features[features != 'w_delta']
#features <- features[features != 'tl_covid']

# result matrix
super_result<- data.frame()
cbps_result <- data.frame()
iptw_result <- data.frame()

subset_predictors <- c("p65_over","ppopurban","disadvantage","abavg_nhw","abavg_nhb","abavg_hisp","n_covid","w_covid", "tl_covid")

# for loop to go over all variables
for(treat in subset_predictors){
  trt_model_covariates <- features[features != treat]
  
  # super learner
  weight_model_super <- WeightIt::weightit(as.formula(paste0(paste0(treat, '~'),paste0(trt_model_covariates, collapse= "+"))),
                                                   method = 'super',stop.method = 's.mean',
                                                   SL.library = c('SL.gam'), data = data)

  weight_model_super <- trim(weight_model_super, 0.97)

  outcome_model_super <- glm(as.formula(paste0(paste0(response,'~',treat,'+'),paste0(trt_model_covariates, collapse= "+"))),
                                                  weights = weight_model_super$weights, data = data)

  super_result <- rbind(super_result, summary(outcome_model_super)$coefficients[2,])

  # # cbps
  weight_model_cbps <- WeightIt::weightit(as.formula(paste0(paste0(treat, '~'),paste0(trt_model_covariates, collapse= "+"))),
                                          method = 'cbps', data = data, over = FALSE)

  weight_model_cbps <- trim(weight_model_cbps, 0.97)

  outcome_model_cbps <- glm(as.formula(paste0(paste0(response,'~',treat,'+'),paste0(trt_model_covariates, collapse= "+"))),
                                         weights = weight_model_cbps$weights, data = data)

  cbps_result <- rbind(cbps_result, summary(outcome_model_cbps)$coefficients[2,])

  # IPTW
  weight_model_iptw <- WeightIt::weightit(as.formula(paste0(paste0(treat, '~'),paste0(trt_model_covariates, collapse= "+"))),
                                          method = 'ps', data = data, over = FALSE)
  
  weight_model_iptw <- trim(weight_model_iptw, 0.97)
  
  outcome_model_iptw <- glm(as.formula(paste0(paste0(response,'~',treat,'+'),paste0(trt_model_covariates, collapse= "+"))),
                            weights = weight_model_iptw$weights, data = data, maxit = 1000)
  
  iptw_result <- rbind(iptw_result, summary(outcome_model_iptw)$coefficients[2,])
  
}

# write the output
write.csv(super_result, file = 'deaths_super.csv',
          row.names = c("p65_over","ppopurban","disadvantage","abavg_nhw","abavg_nhb","abavg_hisp","n_covid","w_covid", "tl_covid"))
write.csv(cbps_result, file = 'deaths_cbps.csv',
          row.names = c("p65_over","ppopurban","disadvantage","abavg_nhw","abavg_nhb","abavg_hisp","n_covid","w_covid", "tl_covid"))
write.csv(iptw_result, file = 'deaths_iptw.csv',
          row.names = c("p65_over","ppopurban","disadvantage","abavg_nhw","abavg_nhb","abavg_hisp","n_covid","w_covid", "tl_covid"))

