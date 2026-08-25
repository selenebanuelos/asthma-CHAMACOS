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

# import data ------------------------------------------------------------------
# asthma-related data
asthma <- read.dta13('data-raw/de_la_Rosa_07.dta',
                     nonint.factors = TRUE,
                     generate.factors = TRUE)

# data wrangling ---------------------------------------------------------------
# do some reformatting to prepare for LCA
asthma_sym <- asthma %>%
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
  # make data wider to use with poLCA::poLCA()
  pivot_wider(id_cols = c('newid', 'cham'),
              names_from = age,
              values_from = c(asth_sym), 
                              # asth_symed, 
                              # asth_med, 
                              # asth_diag_ever),
              names_glue = '{.value}_{age}')

# conditional independence assumption: latent class membership explains all of
# the shared variance among the observed indicators

# correlation of candidate indicators ------------------------------------------
# # correlation between current asthma at different ages
# asthma %>%
#   # keep asthma related variables only
#   dplyr::select(newid, cham, contains('asth_')) %>%
#   # remove participants that are missing all data between 9Y-18Y
#   filter(!if_all(contains(c('9y', '10Y', '12Y', '14Y', '18Y')), is.na)) %>%
#   # make data longer for some data manipulation
#   pivot_longer(cols = contains('asth_'),
#                names_to = c('.value', 'age'),
#                names_pattern = '(asth_.*)_([0-9]+Y)$') %>%
#   # remove data from ages 5 and 7
#   filter(age != '5Y' & age != '7Y') %>%
#   # remove data from participants with only 1 visit (can't create trajectory?)
#   group_by(newid) %>%
#   filter(n_distinct(age) > 1) %>%
#   ungroup() %>%
#   # create a 'current asthma' variable for each timepoint
#   # current asthma defined as 2/3 of the following: current asthma symptoms, 
#   # current asthma medication, ever asthma diagnosis
#   # asth_symed = current asthma symptoms OR current asthma medication use
#   mutate(
#     current_asthma = case_when(
#       asth_symed == 'Yes' & asth_diag_ever == 'Yes' ~ 1, # yes
#       asth_symed == 'No' | asth_diag_ever == 'No' ~ 0, # no
#       is.na(asth_symed) | is.na(asth_diag_ever) ~ NA)) %>%
#   # make data wider to use with poLCA::poLCA()
#   pivot_wider(id_cols = c('newid', 'cham'),
#               names_from = age,
#               values_from = c(asth_sym, 
#                               asth_symed, 
#                               asth_med, 
#                               asth_diag_ever, 
#                               current_asthma),
#               names_glue = '{.value}_{age}') %>%
#   dplyr::select(contains('current_asthma')) %>%
#   # calculate tetrachoric correlations
#   sirt::tetrachoric2(.) %>%
#   # view correlation matrix
#   .$rho
# all correlations are >0.5
# do I need to remove people with only 1 measure??? 
# some resources state that at least three indicators are required for
# identifiability. I guess I could do a sensitvity analysis where I include
# all people (even those with just 1 measure) to see how results differ

# correlation between current asthma sym, current med use, ever asthma diagnosis
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
  # make data wider to use with poLCA::poLCA()
  pivot_wider(id_cols = c('newid', 'cham'),
            names_from = age,
            values_from = c(asth_sym, 
                            asth_symed, 
                            asth_med, 
                            asth_diag_ever),
            names_glue = '{.value}_{age}') %>%
  dplyr::select(contains(c('asth_sym_', 'asth_med_', 'asth_symed', 'asth_diag_'))) %>%
  # change all 'Yes' to 1, 'No to 0
  mutate(across(everything(), ~ case_when(. == 'Yes' ~ 1,
                                          . == 'No' ~ 0,
                                          is.na(.) ~ NA))) %>%
  # calculate tetrachoric correlations
  sirt::tetrachoric2(.) %>%
  # view correlation matrix
  .$rho %>%
  # visualize correlation matrix
  ggcorrplot(., 
             method = 'circle', 
             type = 'upper', 
             lab = TRUE)

# ever asthma diagnosis makes me hesitate due to dependence between time points
# check that ever asthma diagnosis follows the trend no->yes for every person
# maybe change this from a repeated measure to age of asthma diagnosis? have to 
# think this through a bit more

# practitioner's guide says to avoid including clinical diagnosis as indicator

# check data missingness among analytic sample ---------------------------------
# “FIML approaches in LCA handle any ignorable missing data on the indicators of
# the latent variable, but individuals with missing data on a grouping variable
# or any covariate in the model are deleted from the analysis”
# (Collins and Lanza, 2010)

# Latent class analysis --------------------------------------------------------
# formula: response ~ predictors
indicators <- as.formula(
  paste(
    "cbind(", 
    paste(grep('asth_sym_', names(asthma_sym), value = TRUE),
          collapse = ","),
    ") ~ 1")
)
# not using covariates to estimate latent class membership 
# this means: response ~ 1 in formula above

# create age at first reported asthma diagnosis

# estimate models, assuming 2 to 5 classes
set.seed(1234)
models <- lapply(2:5, function(k) 
{
  
  # print status message to console
  message("Currently estimating model with ", k, " classes...")
  
  # estimate model
  poLCA(indicators,
        data = asthma_sym,
        nclass = k,
        nrep = 1000, # estimate model 100 times to search for global maximum
        maxiter = 5000,
        na.rm = FALSE,
        verbose = FALSE)
  
})

# Visualize log-likelihood distributions ---------------------------------------
# gives a sense of model identification-ideally looking for one clear global peak

# vector containing the maximum log-likelihood values found in each of the nrep 
# attempts to fit the model
get_loglike <- function(model) {
  
  data.frame(
    # maximum log-likelihood value found at each repeated fit
    log_like = model$attempts,
    # number of classes assumed for model fit
    k = as.character(length(model$P))
  )
}

# create dataframe of all max log-likelihood values obtained for each model
map_df(models, get_loglike) %>%
  # remove model with 2 classes since there is only 1 value obtained
  #filter(k != 2) %>%
  # plot max log-likelihood distributions to inspect frequency of values
  ggplot(aes(x = log_like, fill = k)) +
  geom_histogram(bins = 100) +
  facet_wrap(~k) +
  labs(x = 'Maximum log-likelihood',
       y = 'Count',
       title = 'Frequency of maximum log-likelihood values over repeated starts') +
  theme_minimal()

# Evaluate model fit -----------------------------------------------------------
# function that create the comparison table of models
compare_fit <- function(model_list) {
  
  comparison <- data.frame(
    Classes = 2:5,
    G2 = sapply(model_list, function(m) round(m$Gsq, 2)),
    df = sapply(model_list, function(m) m$npar),
    pvalue = sapply(model_list, function(m) round(pchisq(m$Gsq, df = m$resid.df, lower.tail = FALSE), 4)),
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
# "BIC heavily penalizes the addition of parameters to the model in relation to
# the sample size, where the larger the sample size, the greater the penalty...
# Whereas, as n increases, the AIC has a tendency to select more complex models
# (more classes), as the best fitting, because sample size is not a determining
# factor in its estimation"
# (Sinha et al., 2021)

# In LCA, we're hoping to "find a model for which the null hypothesis is not 
# regected". The larger value of G^2 (the log likelihood), "the more evidence
# there is against the null hypothesis" (Collins and Lanza, 2010)
# Therefore, we want a lower log-likelihood?

# function that gets respective latent class sizes from each model
compare_sizes <- function(model) {
  
  #model = models[[1]]
  
  # counts of predicted class membership, by modal assignment
  table(model$predclass) %>% 
    prop.table(.) %>% 
    round(., digits = 3)
  
}

# check latent class sizes - don't want small classes
lapply(models, compare_sizes)

# can also do k-folds cross-validation to select best model 

# need to use multiple random starts to demonstrate sufficient replication
# of the maximum likelihood

# Visualize trajectories -------------------------------------------------------
# may want to call them 'profiles' since they're not really modeled as
# longitudinal growth outcomes, and therefore, not continuous trajectories

# function that gets class-conditional outcome probabilities from LCA models
get_probs <- function(model) {
  
  # class-conditional outcome probabilities
  probs <- model$probs
  
  # reformat data for plotting
  reformat <- map_df(probs, # class-conditional outcome probabilities
                     function(m) as.data.frame(m) %>% 
                       mutate(class = str_extract(row.names(.), '[0-9]+')), 
                     .id = 'age') %>%
    # create column with number of classes assumed (k) for model
    mutate(k = max(class)) %>%
    # reformat age variable values from "current_asthma_*Y" to just digits
    mutate(age = str_extract(age, '[0-9]+'),
           age = factor(age,levels = c('9', '10', '12', '14', '16', '18')))
  
}

# define custom labels for facets in plot below
facet_names <- c('2' = '2 Classes',
                 '3' = '3 Classes',
                 '4' = '4 Classes',
                 '5' = '5 Classes')

# visualize class profiles of all models in one plot
map_df(models, get_probs) %>% # get class-conditional outcome probabilities
  ggplot(aes(x = age, # show trend over ages
             y = `Pr(2)`, # plot probability of having current asthma symptoms
             color = class # show each class in different color
  )) +
  geom_line(aes(group = class), size = 1) +
  facet_wrap(vars(k), labeller = as_labeller(facet_names)) +
  labs(y = 'Probability of current asthma symptoms') +
  theme_minimal() +
  theme(panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.5))
