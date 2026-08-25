### Author: Selene Banuelos
### Date: 6/4/2026
### Description: Generalized estimating equations will be used to estimate 
### population average DNAm age deviation given covariates.
### DNAm age ~ chrono age + asthma + chrono age*asthma + covariates

# setup
library(dplyr)
library(gee)
# function to get post-estimate associations from GEE (by Alan Hubbard)
#source('code/gee_post_estimate.R')

# import data
################################################################################
# DNAm age estimates and asthma category
data <- read.csv('data-processed/dnam-age-asthma.csv')

# data wrangling
################################################################################
clean <- data
  # could possibly keep all time points and then just look at ages 9 and 18

# factor variables
clean$asth_lca_3 <- factor(clean$asth_lca_3,
                           levels = c('Never/Infrequent', # reference: N/I
                                      'Intermediate', 
                                      'Persistent'))

# regression using generalized estimating equations
################################################################################


# calculate 95% confidence intervals for all coefficients 
################################################################################
# function that creates 95% CI using robust SE
robust_ci <- function(model, # gee object
                      label # character string describing model
                      ){
  
  # extract coefficients (1) and robust SE (4)
  coef_data <- summary(model)$coefficients[, c(1,4)]
  
  # calculate bounds
  lower_bound <- coef_data[,1] - 1.96 * coef_data[,2]
  upper_bound <- coef_data[,1] + 1.96 * coef_data[,2]
  
  # format 95% CI with estimate and bounds
  ci_table <- data.frame(
    Estimate = round(coef_data[,1], digits = 4),
    Lower_95_CI = round(lower_bound, digits = 4),
    Upper_95_CI = round(upper_bound, digits = 4)
  )
  
  print(paste('95% CIs for', label))
  return(ci_table)
  
}