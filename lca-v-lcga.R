# Author: Selene Banuelos
# Date: 8/17/2026
# Description: Compare participant classifications between latent class analysis
# (LCA) and latent class growth analysis (LCGA)

# setup
library(dplyr)
library(lcmm)

# import data ------------------------------------------------------------------
# load LCA objects
load('data-processed/poLCA.RData')

# load LCGA objects
load('data-processed/lcmm.RData')

# asthma-related data
asthma <- read.dta13('data-raw/de_la_Rosa_07.dta',
                     nonint.factors = TRUE,
                     generate.factors = TRUE)

# move this section to lcmm.R script later -------------------------------------
# plot trajectories
plot_trajectory <- function(model){
  
  ages <- data.frame(age_years = c(9, 10, 12, 14, 16, 18))
  
  pred_class <- predictY(model, ages, var.time = 'age_years', draws = TRUE)
  
  plot(pred_class)
  
}

# linear models
plot_trajectory(linear_2)
plot_trajectory(linear_3) # lowest BIC
plot_trajectory(linear_4)

# quadratic models
plot_trajectory(quadratic_1)
plot_trajectory(quadratic_3)

# cubic models
plot_trajectory(cubic_1)
plot_trajectory(cubic_2)
plot_trajectory(cubic_3)

# get vector of predicted class membership for 3-class linear model



# Data wrangling ---------------------------------------------------------------
# do some reformatting to prepare for LCA
current_asthma <- asthma %>%
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
    # convert to factor to use in LCA
    current_asthma = factor(current_asthma,
                            levels = c(1, 0),
                            labels = c('Yes', 'No')
    )) %>%
  # make data wider to use with poLCA::poLCA()
  pivot_wider(id_cols = c('newid', 'cham'),
              names_from = age,
              values_from = c(asth_sym, 
                              asth_symed, 
                              asth_med, 
                              asth_diag_ever, 
                              current_asthma),
              names_glue = '{.value}_{age}')

# fit statistics ---------------------------------------------------------------
# get summaries of latent class growth models
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

# create summary data frame for all latent class growth models
summaries <- rbind(
  map_df(linear_models, summary_table, .id = 'type'),
  map_df(quadratic_models, summary_table, .id = 'type'),
  map_df(cubic_models, summary_table, .id = 'type')
) %>%
  group_by(type) %>%
  # sort by ascending BIC
  arrange(BIC, .by_group = TRUE)

# get summaries of latent class models
compare_fit <- function(model_list) {
  
  comparison <- data.frame(
    Classes = 2:5,
    Log_Likelihood = sapply(model_list, function(m) round(m$llik, 2)),
    BIC = sapply(model_list, function(m) round(m$bic, 2)),
    AIC = sapply(model_list, function(m) round(m$aic, 2)),
    Smallest_Class_Pct = sapply(model_list, function(m) {round(min(m$P)*100,1)}),
    # higher entropy indicates better class separation
    Entropy = sapply(model_list, poLCA.entropy)
  )
  
  print(comparison)
  
  return(comparison)
  
}

# review comparison table for different k
model_fit <- compare_fit(models)

# compare predicted class membership -------------------------------------------
