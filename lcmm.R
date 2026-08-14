# Author: Selene Banuelos
# Date: 8/11/2026
# Description: Repeat measures latent class analysis and latent growth curve
# analysis using lcmm R package

# setup
library(readstata13)
library(dplyr)
library(tidyr)
library(stringr)
library(lcmm)
options(scipen = 999)

# import data ------------------------------------------------------------------
# asthma related data
asthma <- read.dta13('data-raw/de_la_Rosa_07.dta',
                     nonint.factors = TRUE,
                     generate.factors = TRUE)

# data wrangling ---------------------------------------------------------------
# do some reformatting to prepare for LCA
curr_asth_data <- asthma %>%
  # keep asthma related variables only
  dplyr::select(newid, cham, contains('asth_')) %>%
  # remove participants that are missing all data between 9Y-18Y
  filter(!if_all(contains(c('9y', '10Y', '12Y', '14Y', '18Y')), is.na)) %>%
  # make data longer for some data manipulation
  pivot_longer(cols = contains('asth_'),
               names_to = c('.value', 'age'),
               names_pattern = '(asth_.*)_([0-9]+Y)$') %>%
  # remove data from ages 5 and 7
  filter(age != '5Y' & age != '7Y') %>%
  # remove data from participants with only 1 visit (can't create trajectory?)
  group_by(newid) %>%
  filter(n_distinct(age) > 1) %>%
  ungroup() %>%
  # create a 'current asthma' variable for each timepoint
  # current asthma defined as 2/3 of the following: current asthma symptoms, 
  # current asthma medication, ever asthma diagnosis
  # asth_symed = current asthma symptoms OR current asthma medication use
  mutate(
    current_asthma = case_when(
      asth_symed == 'Yes' & asth_diag_ever == 'Yes' ~ 1, # yes
      asth_symed == 'No' | asth_diag_ever == 'No' ~ 0, # no
      is.na(asth_symed) | is.na(asth_diag_ever) ~ NA),
    # # convert to factor to use in LCA
    # current_asthma = factor(current_asthma,
    #                         levels = c(1, 0),
    #                         labels = c('Yes', 'No')),
    age_years = as.numeric(str_remove(age, 'Y')),
    # change subject ID to numeric to use with lcmm()
    newid = as.numeric(newid)
    ) 

# fit latent class mixture model -----------------------------------------------
# define model fitting function
fit_lcmm <- function(k, m){ # k = number of assumed classes
  
  # print status message to console
  message("Currently estimating model with ", k, " classes...")
  
  # fit model
  lcmm(
    fixed = current_asthma ~ age_years,
    mixture = ~ age_years,
    #random = ~ age_years,
    subject = 'newid',
    link = 'thresholds', # binary outcome
    ng = k, # assume k classes
    data = curr_asth_data,
    B = m # 1-class model
  )
  
  }
  
# fit a 1 class model to obtain starting values
set.seed(123)
k_1 <- lcmm(fixed = current_asthma ~ age_years,
            #random = ~ -1, # no within-class random effects (latent class growth analysis)
            subject = 'newid',
            link = 'thresholds', # binary outcome
            ng = 1, # 1 class model
            data = curr_asth_data
            )

# fit 2 class model, with 100 random starting values
set.seed(123)
k_2 <- gridsearch(
  lcmm(fixed = current_asthma ~ age_years,
       mixture = ~ age_years,
       #random = -1, # no within-class random effects
       subject = 'newid',
       link = 'thresholds',
       ng = 2,
       data = curr_asth_data),
  rep = 100, # try 100 different sets of random initial values
  maxiter = 100, # default in lcmm()
  minit = k_1 # 1-class model used to generate random initial values
)

# fit 3 class model, with 100 random starting values
set.seed(123)
k_3 <- gridsearch(
  lcmm(fixed = current_asthma ~ age_years,
       mixture = ~ age_years,
       #random = -1, # no within-class random effects
       subject = 'newid',
       link = 'thresholds',
       ng = 3,
       data = curr_asth_data),
  rep = 100, # try 100 different sets of random initial values
  maxiter = 100, # default in lcmm()
  minit = k_1 # 1-class model used to generate random initial values
)

# fit 4 class model, with 100 random starting values
set.seed(123)
k_4 <- gridsearch(
  lcmm(fixed = current_asthma ~ age_years,
       mixture = ~ age_years,
       #random = -1, # no within-class random effects
       subject = 'newid',
       link = 'thresholds',
       ng = 4,
       data = curr_asth_data),
  rep = 100, # try 100 different sets of random initial values
  maxiter = 100, # default in lcmm()
  minit = k_1 # 1-class model used to generate random initial values
)

# fit 5 class model, with 100 random starting values
set.seed(123)
k_5 <- gridsearch(
  lcmm(fixed = current_asthma ~ age_years,
       mixture = ~ age_years,
       #random = -1, # no within-class random effects
       subject = 'newid',
       link = 'thresholds',
       ng = 5,
       data = curr_asth_data),
  rep = 100, # try 100 different sets of random initial values
  maxiter = 100, # default in lcmm()
  minit = k_1 # 1-class model used to generate random initial values
)

# fit model with 2-5 classes, using starting values from 1-class model
#models <- lapply(2:5, fit_lcmm, k_1)

# save models in list
models <- list(k_1, k_2, k_3, k_4, k_5)

# save fitted models
save.image('data-processed/lcmm.RData')

# post-fit summaries -----------------------------------------------------------
summaryplot(k_1,
            k_2,
            k_3,
            k_4,
            k_5,
            which = c('conv', 'AIC', 'BIC', 'entropy')
            )

summarytable(k_1,
            k_2,
            k_3,
            k_4,
            k_5,
            which = c('G','conv', 'AIC', 'BIC', 'entropy', '%class')
)

# plot trajectories ------------------------------------------------------------
plot_trajectory <- function(model){

  ages <- data.frame(age_years = c(9, 10, 12, 14, 16, 18))

  pred_class <- predictY(model, ages, var.time = 'age_years', draws = TRUE)

  plot(pred_class)

}

#lapply(models, plot_trajectory)

# posterior classification and posterior individual class-membership probabilities
class_probs <- k_3$pprob
postprob(k_3)

# plots
plot_trajectory(k_3)
plot(k_3, which = 'postprob')
plot(k_3, which = 'link')
plot(k_3, which = 'linkfunction')

# below not available for thresholds mixed models
plot(k_3, which = 'residuals')
plot(k_3, which = 'fit', var.time = 'ages_years') 