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
# create current asthma variable
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
  # remove original asthma variables and keep only current_asthma
  dplyr::select(-contains('asth_')) %>%
  # make data wider to use with poLCA::poLCA()
  pivot_wider(id_cols = c('newid', 'cham'),
              names_from = age,
              values_from = current_asthma,
              names_glue = '{.value}_{age}')

# check if any participants have data at only 1 visit
group_by(newid)

# remove data from participants with only 1 visit
# skeptical of using only 1 data point to evaluate a trajectory
group_by(newid) %>%
  filter(n_distinct(age) > 1) %>%
  ungroup() %>%
  
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
        nrep = 1000, # estimate model 1,000 times to search for global maximum
        maxiter = 5000,
        na.rm = FALSE,
        verbose = FALSE)
  
  })

rm(asthma, current_asthma)
save.image('data-processed/poLCA.RData')

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

# create footnote to be inlcuded in plot
footnote <- 'K = 2 model not included since only one value for maximum log-likelihood was obtained over all repeated starts'

# create dataframe of all max log-likelihood values obtained for each model
map_df(models, get_loglike) %>%
  # remove model with 2 classes since there is only 1 value obtained
  filter(k != 2) %>%
  # plot max log-likelihood distributions to inspect frequency of values
  ggplot(aes(x = log_like, fill = k)) +
  geom_histogram(bins = 100) +
  facet_wrap(~k) +
  labs(x = 'Maximum log-likelihood',
       y = 'Count',
       title = 'Frequency of maximum log-likelihood values over repeated starts',
       caption = footnote) +
  theme_minimal()

# Evaluate model fit -----------------------------------------------------------
# function that create the comparison table of models
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

# "BIC heavily penalizes the addition of parameters to the model in relation to
# the sample size, where the larger the sample size, the greater the penalty...
# Whereas, as n increases, the AIC has a tendency to select more complex models
# (more classes), as the best fitting, because sample size is not a determining
# factor in its estimation"
# (Sinha et al., 2021)

# function that gets respective latent class sizes from each model
compare_sizes <- function(model) {
  
  # counts of predicted class membership, by modal assignment
  table(model$predclass) %>% 
    prop.table(.) %>% 
    round(., digits = 3)
  
}

# check latent class sizes - don't want small classes
lapply(models, compare_sizes)

# all k produce several small classes

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
spaghetti_all <- map_df(models, get_probs) %>% # get class-conditional outcome probabilities
  ggplot(aes(x = age, # show trend over ages
             y = `Pr(1)`, # plot probability of having current asthma as outcome
             color = class # show each class in different color
             )) +
  geom_line(aes(group = class), size = 1) +
  facet_wrap(vars(k), labeller = as_labeller(facet_names)) +
  labs(y = 'Probability of current asthma') +
  theme_minimal() +
  theme(panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.5))

# visualize class profiles of k=3 model only
spaghetti_k3 <- map_df(models, get_probs) %>% # get class-conditional outcome probabilities
  # keep only data from k=4 model
  filter(k == 3) %>%
  # make spaghetti plot
  ggplot(aes(x = age, # show trend over ages
             y = `Pr(1)`, # plot probability of having current asthma as outcome
             color = class # show each class in different color
  )) +
  geom_line(aes(group = class), size = 1) +
  labs(y = 'Probability of current asthma',
       title = 'LCA with 3 Classes') +
  theme_minimal() +
  theme(panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.5))

# define classes ---------------------------------------------------------------
# get vector of predicted class membership from 3-class model
pred_class <- models[[2]]$predclass

# define classes for 3-class model
class_lca_k3 <- data.frame(newid = current_asthma$newid,
                         pred_class = pred_class) %>%
  mutate(class_label = case_when(pred_class == 1 ~ 'late onset',
                                 pred_class == 2 ~ 'never/infrequent',
                                 pred_class == 3 ~ 'persistent'))

# output -----------------------------------------------------------------------
write.csv(model_fit, 'data-processed/curr_sym_model_fit.csv', row.names = FALSE)

# trajectories
ggsave('figures/spaghetti-all-traj.png', spaghetti_all)
ggsave('figures/spaghetti-k3-traj.png', spaghetti_k3)

# remove unecessary objects
rm(asthma, 
   current_asthma, 
   model_fit, 
   current_asthma_ind, 
   facet_names, 
   footnote, 
   pred_class,
   spaghetti_all,
   spaghetti_k3,
   compare_fit,
   compare_sizes,
   get_loglike,
   get_probs
)

# save model fits and classification data
save.image('data-processed/poLCA.RData')