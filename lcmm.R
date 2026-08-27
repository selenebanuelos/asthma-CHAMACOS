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
library(purrr)
options(scipen = 999)

# import data ------------------------------------------------------------------
# asthma related data
asthma <- read.dta13('data-raw/de_la_Rosa_07.dta',
                     nonint.factors = TRUE,
                     generate.factors = TRUE)

# data wrangling ---------------------------------------------------------------
# create analytic sample
current_asthma <- asthma %>%
  # keep asthma related variables only
  dplyr::select(newid, cham, contains('asth_')) %>%
  # remove participants that are missing all data between 9Y-18Y
  filter(!if_all(contains(c('9y', '10Y', '12Y', '14Y', '18Y')), is.na)) %>%
  # make data longer for some data manipulation & use with lcmm()
  pivot_longer(cols = contains('asth_'),
               names_to = c('.value', 'age'),
               names_pattern = '(asth_.*)_([0-9]+Y)$') %>%
  # remove data from ages 5 and 7
  filter(age != '5Y' & age != '7Y') %>%
  # create a 'current asthma' variable for each timepoint
  # current asthma defined as 2/3 of the following: current asthma symptoms, 
  # current asthma medication, ever asthma diagnosis
  # asth_symed = current asthma symptoms OR current asthma medication use
  mutate(current_asthma = case_when(
    asth_symed == 'Yes' & asth_diag_ever == 'Yes' ~ 1, # yes
    asth_symed == 'No' | asth_diag_ever == 'No' ~ 0, # no
    is.na(asth_symed) | is.na(asth_diag_ever) ~ NA)) %>%
  # remove original asthma variables and keep only current_asthma
  dplyr::select(-contains('asth_')) %>%
  # remove any observations missing a current asthma value
  filter(!is.na(current_asthma)) %>%
  # remove data from participants with only 1 visit
  group_by(newid) %>%
  filter(n_distinct(age) > 1) %>%
  ungroup() %>%
  # create numeric age variable to use with lcmm()
  mutate(age_years = as.numeric(str_remove(age, 'Y')), 
         # change subject ID to numeric to use with lcmm()
         newid = as.numeric(newid))

# fit latent class mixture model -----------------------------------------------
# save any console output to text file (warnings, messages, etc.)
sink('data-processed/lcmm-console.txt')

# fit a 1 class linear model to obtain starting values
set.seed(123)
linear_1 <- lcmm(fixed = current_asthma ~ age_years,
                 # mixture not specified for 1 class models
                 subject = 'newid',
                 link = 'thresholds', # binary outcome
                 ng = 1, # 1 class model
                 data = curr_asth_data
                 )

# fit linear 2-5 class models, using gridsearch() to try random sets of initial 
# values and 1-class model for starting values
set.seed(123)
linear_2 <- gridsearch(
  
  # fit latent class growth model
  lcmm(fixed = current_asthma ~ age_years,
       mixture = ~ age_years,
       #random = ~ age_years, # no within-class random effects (latent class growth analysis)
       subject = 'newid',
       link = 'thresholds', # binary outcome
       ng = 2, # assume k classes
       data = curr_asth_data
  ),
  
  rep = 100, # try 100 different sets of random initial values
  maxiter = 1000, # 1000 iterations max (100 iterations not enough for cubic)
  minit = linear_1 # 1-class model used to generate random initial values
)

set.seed(123)
linear_3 <- gridsearch(
  
  # fit latent class growth model
  lcmm(fixed = current_asthma ~ age_years,
       mixture = ~ age_years,
       subject = 'newid',
       link = 'thresholds', # binary outcome
       ng = 3, # assume k classes
       data = curr_asth_data
  ),
  
  rep = 100, # try 100 different sets of random initial values
  maxiter = 1000, # 1000 iterations max (100 iterations not enough for cubic)
  minit = linear_1 # 1-class model used to generate random initial values
)

set.seed(123)
linear_4 <- gridsearch(
  
  # fit latent class growth model
  lcmm(fixed = current_asthma ~ age_years,
       mixture = ~ age_years,
       subject = 'newid',
       link = 'thresholds', # binary outcome
       ng = 4, # assume k classes
       data = curr_asth_data
  ),
  
  rep = 100, # try 100 different sets of random initial values
  maxiter = 1000, # 1000 iterations max (100 iterations not enough for cubic)
  minit = linear_1 # 1-class model used to generate random initial values
)

set.seed(123)
linear_5 <- gridsearch(
  
  # fit latent class growth model
  lcmm(fixed = current_asthma ~ age_years,
       mixture = ~ age_years,
       subject = 'newid',
       link = 'thresholds', # binary outcome
       ng = 5, # assume k classes
       data = curr_asth_data
  ),
  
  rep = 100, # try 100 different sets of random initial values
  maxiter = 1000, # 1000 iterations max (100 iterations not enough for cubic)
  minit = linear_1 # 1-class model used to generate random initial values
)


# fit a 1-class quadratic model for starting values
set.seed(123)
quadratic_1 <- lcmm(fixed = current_asthma ~ poly(age_years, 2, raw = TRUE),
                    # mixture not specified for 1 class models
                    subject = 'newid',
                    link = 'thresholds', # binary outcome
                    ng = 1, # 1 class model
                    data = curr_asth_data
                    )

set.seed(123)
quadratic_2 <- gridsearch(
  
  # fit latent class growth model
  lcmm(fixed = current_asthma ~ poly(age_years, 2, raw = TRUE),
       mixture = ~ poly(age_years, 2, raw = TRUE),
       subject = 'newid',
       link = 'thresholds', # binary outcome
       ng = 2, # assume k classes
       data = curr_asth_data
  ),
  
  rep = 100, # try 100 different sets of random initial values
  maxiter = 1000, # 1000 iterations max (100 iterations not enough for cubic)
  minit = quadratic_1 # 1-class model used to generate random initial values
)

set.seed(123)
quadratic_3 <- gridsearch(
  
  # fit latent class growth model
  lcmm(fixed = current_asthma ~ poly(age_years, 2, raw = TRUE),
       mixture = ~ poly(age_years, 2, raw = TRUE),
       subject = 'newid',
       link = 'thresholds', # binary outcome
       ng = 3, # assume k classes
       data = curr_asth_data
  ),
  
  rep = 100, # try 100 different sets of random initial values
  maxiter = 1000, # 1000 iterations max (100 iterations not enough for cubic)
  minit = quadratic_1 # 1-class model used to generate random initial values
)

set.seed(123)
quadratic_4 <- gridsearch(
  
  # fit latent class growth model
  lcmm(fixed = current_asthma ~ poly(age_years, 2, raw = TRUE),
       mixture = ~ poly(age_years, 2, raw = TRUE),
       subject = 'newid',
       link = 'thresholds', # binary outcome
       ng = 4, # assume k classes
       data = curr_asth_data
  ),
  
  rep = 100, # try 100 different sets of random initial values
  maxiter = 1000, # 1000 iterations max (100 iterations not enough for cubic)
  minit = quadratic_1 # 1-class model used to generate random initial values
)

set.seed(123)
quadratic_5 <- gridsearch(
  
  # fit latent class growth model
  lcmm(fixed = current_asthma ~ poly(age_years, 2, raw = TRUE),
       mixture = ~ poly(age_years, 2, raw = TRUE),
       subject = 'newid',
       link = 'thresholds', # binary outcome
       ng = 5, # assume k classes
       data = curr_asth_data
  ),
  
  rep = 100, # try 100 different sets of random initial values
  maxiter = 1000, # 1000 iterations max (100 iterations not enough for cubic)
  minit = quadratic_1 # 1-class model used to generate random initial values
)

# fit a 1-class quadratic model for starting values
set.seed(123)
cubic_1 <- lcmm(fixed = current_asthma ~ poly(age_years, 3, raw = TRUE),
                # mixture not specified for 1 class models
                subject = 'newid',
                link = 'thresholds', # binary outcome
                ng = 1, # 1 class model
                data = curr_asth_data,
                maxiter = 1000
)

set.seed(123)
cubic_2 <- gridsearch(
  
  # fit latent class growth model
  lcmm(fixed = current_asthma ~ poly(age_years, 3, raw = TRUE),
       mixture = ~ poly(age_years, 3, raw = TRUE),
       subject = 'newid',
       link = 'thresholds', # binary outcome
       ng = 2, # assume k classes
       data = curr_asth_data
  ),
  
  rep = 100, # try 100 different sets of random initial values
  maxiter = 1000, # 1000 iterations max (100 iterations not enough for cubic)
  minit = cubic_1 # 1-class model used to generate random initial values
)

set.seed(123)
cubic_3 <- gridsearch(
  
  # fit latent class growth model
  lcmm(fixed = current_asthma ~ poly(age_years, 3, raw = TRUE),
       mixture = ~ poly(age_years, 3, raw = TRUE),
       subject = 'newid',
       link = 'thresholds', # binary outcome
       ng = 3, # assume k classes
       data = curr_asth_data
  ),
  
  rep = 100, # try 100 different sets of random initial values
  maxiter = 1000, # 1000 iterations max (100 iterations not enough for cubic)
  minit = cubic_1 # 1-class model used to generate random initial values
)

# stop sinking
sink(NULL)

# refit desired models, having trouble plotting from list of models

# save fitted models
save.image('data-processed/lcmm.RData')

# post-fit summaries -----------------------------------------------------------
summary_table <- function(model) {
  
  # lcmm::summarytable() source code to calculate entropy
  entropy <- function(x)
  {
    z <- log(as.matrix(x$pprob[,c(3:(x$ng+2))]))*as.matrix(x$pprob[,c(3:(x$ng+2))])
    if(any(!is.finite(z)))
    {
      z[which(!is.finite(z))] <- 0
    }
    res <- 1+sum(z)/(x$ns*log(x$ng))
    if(x$ng==1) res <- 1
    return(res)
  }
  
  # calculate entropy
  e <- entropy(model)
  
  # class membership proportions
  class_props <- model$pprob %>%
    group_by(class) %>%
    summarize(prop = n() / nrow(.) * 100) %>%
    pivot_wider(names_from = class,
                names_glue = '%class{class}',
                values_from = prop)
  
  # combine all stats into summary 
  data.frame(k = model$ng,
             entropy = e,
             conv = model$conv,
             AIC = model$AIC,
             BIC = model$BIC) %>%
    cbind(class_props)
  
}

# name each list of models with functional form
linear_models <- set_names(linear_models, 'linear')
quadratic_models <- set_names(quadratic_models, 'quadratic')
cubic_models <- set_names(cubic_models, 'cubic')

# create summary data frame for all models
summaries <- rbind(
  map_df(linear_models, summary_table, .id = 'type'),
  map_df(quadratic_models, summary_table, .id = 'type'),
  map_df(cubic_models, summary_table, .id = 'type')
) %>%
  group_by(type) %>%
  # sort by ascending BIC
  arrange(BIC, .by_group = TRUE)

# plot trajectories ------------------------------------------------------------
plot_trajectory <- function(model){
  
  ages <- data.frame(age_years = c(9, 10, 12, 14, 16, 18))
  
  pred_class <- predictY(model, ages, var.time = 'age_years', draws = TRUE)
  
  plot(pred_class)
  
}

plot_trajectory(linear_3)
plot_trajectory(quadratic_3)
plot_trajectory(cubic_2)
