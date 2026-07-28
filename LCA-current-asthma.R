# Author: Selene Banuelos
# Date: 7/20/2026
# Description: Latent class analysis to estimate asthma trajectories

# setup
library(readstata13)
library(dplyr)
library(tidyr)
library(ggcorrplot)
library(poLCA) # previously used for LCA
library(purrr)
library(stringr)
library(ggplot2)
options(scipen = 999)

# Import data ------------------------------------------------------------------
# asthma-related data
asthma <- read.dta13('data-raw/de_la_Rosa_07.dta',
                     nonint.factors = TRUE,
                     generate.factors = TRUE)

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

# conditional independence assumption: latent class membership explains all of
# the shared variance among the observed indicators

# Correlation of candidate indicators ------------------------------------------
# correlation between current asthma at different ages
asthma %>%
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
      is.na(asth_symed) | is.na(asth_diag_ever) ~ NA)) %>%
  # make data wider to use with poLCA::poLCA()
  pivot_wider(id_cols = c('newid', 'cham'),
              names_from = age,
              values_from = c(asth_sym, 
                              asth_symed, 
                              asth_med, 
                              asth_diag_ever, 
                              current_asthma),
              names_glue = '{.value}_{age}') %>%
  dplyr::select(contains('current_asthma')) %>%
  # calculate tetrachoric correlations
  sirt::tetrachoric2(.) %>%
  # view correlation matrix
  .$rho
# all correlations are >0.5
# do I need to remove people with only 1 measure??? 
# some resources state that at least three indicators are required for
# identifiability. I guess I could do a sensitvity analysis where I include
# all people (even those with just 1 measure) to see how results differ

# ever asthma diagnosis makes me hesitate due to dependence between time points
# check that ever asthma diagnosis follows the trend no->yes for every person
# maybe change this from a repeated measure to age of asthma diagnosis? have to 
# think this through a bit more

# practitioner's guide says to avoid including clinical diagnosis as indicator

# Check data missingness among analytic sample ---------------------------------
# “FIML approaches in LCA handle any ignorable missing data on the indicators of
# the latent variable, but individuals with missing data on a grouping variable
# or any covariate in the model are deleted from the analysis”
# (Collins and Lanza, 2010)

# Latent class analysis --------------------------------------------------------
# formula: response ~ predictors
current_asthma_ind <- as.formula(
  paste(
    "cbind(", 
    paste(grep('current_asthma', names(current_asthma), value = TRUE),
          collapse = ","),
    ") ~ 1")
)
# not using covariates to estimate latent class membership 
# this means: response ~ 1 in formula above

# estimate models, assuming 2 to 5 classes
set.seed(1234)
models <- lapply(2:5, function(k) 
  {
  
  # print status message to console
  message("Currently estimating model with ", k, " classes...")
  
  # estimate model
  poLCA(current_asthma_ind,
        data = current_asthma,
        nclass = k,
        nrep = 100, # estimate model 100 times to search for global maximum
        maxiter = 5000,
        na.rm = FALSE,
        verbose = FALSE)
  
  })

# Evaluate model fit -----------------------------------------------------------
# function that create the comparison table of models
compare_fit <- function(model_list) {
  
  data.frame(
    Classes = 2:5,
    G2 = sapply(model_list, function(m) round(m$Gsq, 2)),
    Log_Likelihood = sapply(model_list, function(m) round(m$llik, 2)),
    BIC = sapply(model_list, function(m) round(m$bic, 2)),
    AIC = sapply(model_list, function(m) round(m$aic, 2)),
    Smallest_Class_Pct = sapply(model_list, function(m) {round(min(m$P)*100,1)}),
    # higher entropy indicates better class separation
    Entropy = sapply(model_list, poLCA.entropy)
  )
  
}

# review comparison table for different k
compare_fit(models)
# "BIC heavily penalizes the addition of parameters to the model in relation to
# the sample size, where the larger the sample size, the greater the penalty...
# Whereas, as n increases, the AIC has a tendency to select more complex models
# (more classes), as the best fitting, because sample size is not a determining
# factor in its estimation"
# (Sinha et al., 2021)

# In LCA, we're hoping to "find a model for which the null hypothesis is not 
# rejected". The larger value of G^2, "the more evidence
# there is against the null hypothesis" (Collins and Lanza, 2010)
# Therefore, we want a lower G^2?

# function that gets respective latent class sizes from each model
compare_sizes <- function(model) {

  # number of classes assumed
  k <- length(model$P)
  
  # respective size of each latent class
  class_prop <- sort(model$P) %>% round(digits = 2)
  
  # labels for each class
  class_index <- seq(1:k)
  
  # named vector of respective latent class sizes
  setNames(class_prop, class_index)
  
}

# check latent class sizes - don't want small classes
lapply(models, compare_sizes)

# all k produce several small classes
# can also do k-folds cross-validation to select best model 

# need to use multiple random starts to demonstrate sufficient replication
# of the maximum likelihood

# Check local dependency -------------------------------------------------------
# 4 class model
model <- models[[3]]

# get number of classes assumed in the model
k <- length(model$P)

# function that generates correlation matrix for indicators within each

# calculate correlation matrix between indicators
pred_classes <- current_asthma %>%
  # change all 'Yes' to 1, 'No to 0 for use in sirt::tetrachoric2()
  mutate(across(contains('asth'), ~ case_when(. == 'Yes' ~ 1,
                                          . == 'No' ~ 0,
                                          is.na(.) ~ NA))) %>%
  # add predicted class memberships to participant data
  mutate(class = model$predclass)

# visualize correlation matrix of current_asthma var within each class
pred_classes %>%
  filter(class == 1) %>% 
  # keep observations from given class
  dplyr::select(contains('current_asthma')) %>%
  # check correlation between current_asthma_*Y variables within each class
  sirt::tetrachoric2(.) %>%
  # view correlation matrix
  .$rho %>%
  # visualize correlation matrix
  ggcorrplot(., 
             method = 'circle', 
             type = 'upper', 
             lab = TRUE)

pred_classes %>%
  filter(class == 3) %>% 
  # keep observations from given class
  dplyr::select(contains('current_asthma')) %>%
  # check correlation between current_asthma_*Y variables within each class
  sirt::tetrachoric2(.) %>%
  # view correlation matrix
  .$rho %>%
  # visualize correlation matrix
  ggcorrplot(., 
             method = 'circle', 
             type = 'upper', 
             lab = TRUE)

pred_classes %>%
  filter(class == 4) %>% 
  # keep observations from given class
  dplyr::select(contains('current_asthma')) %>%
  # check correlation between current_asthma_*Y variables within each class
  sirt::tetrachoric2(.) %>%
  # view correlation matrix
  .$rho %>%
  # visualize correlation matrix
  ggcorrplot(., 
             method = 'circle', 
             type = 'upper', 
             lab = TRUE)
# class 2 cannot generate correlation matrix  

# visualize correlation matrix of separate asthma vars within each class
pred_classes %>%
  filter(class == 1) %>% 
  # keep observations from given class
  dplyr::select(contains(c('asth_'))) %>%
  # check correlation between current_asthma_*Y variables within each class
  sirt::tetrachoric2(.) %>%
  # view correlation matrix
  .$rho %>%
  # visualize correlation matrix
  ggcorrplot(., 
             method = 'circle', 
             type = 'upper', 
             lab = TRUE)

pred_classes %>%
  filter(class == 2) %>% 
  # keep observations from given class
  dplyr::select(contains(c('asth_'))) %>%
  # check correlation between current_asthma_*Y variables within each class
  sirt::tetrachoric2(.) %>%
  # view correlation matrix
  .$rho %>%
  # visualize correlation matrix
  ggcorrplot(., 
             method = 'circle', 
             type = 'upper', 
             lab = TRUE)

pred_classes %>%
  filter(class == 3) %>% 
  # keep observations from given class
  dplyr::select(contains(c('asth_'))) %>%
  # check correlation between current_asthma_*Y variables within each class
  sirt::tetrachoric2(.) %>%
  # view correlation matrix
  .$rho %>%
  # visualize correlation matrix
  ggcorrplot(., 
             method = 'circle', 
             type = 'upper', 
             lab = TRUE)

pred_classes %>%
  filter(class == 4) %>% 
  # keep observations from given class
  dplyr::select(contains(c('asth_'))) %>%
  # check correlation between current_asthma_*Y variables within each class
  sirt::tetrachoric2(.) %>%
  # view correlation matrix
  .$rho %>%
  # visualize correlation matrix
  ggcorrplot(., 
             method = 'circle', 
             type = 'upper', 
             lab = TRUE)

# Visualize trajectories -------------------------------------------------------
# may want to call them 'profiles' since they're not really modeled as
# longitudinal growth outcomes, and therefore, not continuous trajectories
visualize_classes <- function(model) {

  # class-conditional outcome probabilities
  probs <- model$probs
  
  # reformat data for plotting
  reformat <- map_df(probs, 
         function(m) as.data.frame(m) %>% mutate(class = str_extract(row.names(.), '[0-9]+')), 
         .id = 'age') %>%
    mutate(age = str_extract(age, '[0-9]+'),
           age = factor(age,levels = c('9', '10', '12', '14', '16', '18')))
  
  # plot profiles of class-conditional outcome probabilities over time
  reformat %>%
    ggplot(aes(x = age, 
               y = `Pr(1)`,
               color = class)) +
    geom_line(aes(group = class), size = 1) +
    theme_minimal() +
    labs(y = 'Probability of current asthma')
  
}

lapply(models, visualize_classes)