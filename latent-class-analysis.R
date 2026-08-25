# Author: Selene Banuelos
# Date: 7/20/2026
# Description: Latent class analysis to estimate asthma trajectories

# setup
library(readstata13)
library(dplyr)
library(tidyr)
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

# correlation of candidate indicators ------------------------------------------
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
  dplyr::select(contains(c('asth_sym_', 'asth_med_', 'asth_diag_'))) %>%
  # change all 'Yes' to 1, 'No to 0
  mutate(across(everything(), ~ case_when(. == 'Yes' ~ 1,
                                          . == 'No' ~ 0,
                                          is.na(.) ~ NA))) %>%
  # calculate tetrachoric correlations
  sirt::tetrachoric2(.) %>%
  # view correlation matrix
  .$rho

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

# latent class analysis --------------------------------------------------------
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

sym_med_ind <- as.formula(
  paste(
    "cbind(", 
    paste(grep('_sym_|_med_', names(current_asthma), value = TRUE),
          collapse = ","),
    ") ~ 1")
)

#Run models for 2 to 5 classes
set.seed(1234)
model_current_asthma <- lapply(2:5, function(k) {
  message("Currently estimating model with ", k, " classes...")
  poLCA(current_asthma_ind,
        data = current_asthma,
        nclass = k,
        nrep = 100,
        maxiter = 5000,
        na.rm = FALSE,
        verbose = FALSE)
})

set.seed(1234)
model_sym_med <- lapply(2:5, function(k) {
  message("Currently estimating model with ", k, " classes...")
  poLCA(sym_med_ind,
        data = current_asthma,
        nclass = k,
        nrep = 100,
        maxiter = 5000,
        na.rm = FALSE,
        verbose = FALSE)
})

# Create the Comparison Table of Models
create_comparison <- function(lca_object_list) {
  
  data.frame(
    Classes = 2:5,
    Log_Likelihood = sapply(lca_object_list, function(m) round(m$llik, 2)),
    BIC = sapply(lca_object_list, function(m) round(m$bic, 2)),
    AIC = sapply(lca_object_list, function(m) round(m$aic, 2)),
    Smallest_Class_Pct = sapply(lca_object_list, function(m) {round(min(m$P)*100,1)}),
    Entropy = sapply(lca_object_list, poLCA.entropy)
  )
  
}

#Output results
create_comparison(model_current_asthma)
create_comparison(model_sym_med)
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

set.seed(1234)
two_class <- poLCA(function_9_18Y,
                   data = current_asthma,
                   nclass = 2,
                   nrep = 100,
                   maxiter = 5000,
                   na.rm = FALSE,
                   verbose = FALSE)

set.seed(1234)
three_class <- poLCA(function_9_18Y,
                    data = current_asthma,
                    nclass = 3,
                    nrep = 100,
                    maxiter = 5000,
                    na.rm = FALSE,
                    verbose = FALSE)

set.seed(1234)
four_class <- poLCA(function_9_18Y,
                     data = current_asthma,
                     nclass = 4,
                     nrep = 100,
                     maxiter = 5000,
                     na.rm = FALSE,
                     verbose = FALSE)

set.seed(1234)
five_class <- poLCA(function_9_18Y,
                    data = current_asthma,
                    nclass = 5,
                    nrep = 100,
                    maxiter = 5000,
                    na.rm = FALSE,
                    verbose = FALSE)

# check class sizes - don't want small classes
# check entropy - higher entropy indicates better class separation
table(model_current_asthma[[1]]$predclass)

table(model_current_asthma[[2]]$predclass)

table(model_current_asthma[[3]]$predclass)

table(model_current_asthma[[4]]$predclass)
# all k produce several small classes

table(model_sym_med[[1]]$predclass)

table(model_sym_med[[2]]$predclass)

table(model_sym_med[[3]]$predclass)

table(model_sym_med[[4]]$predclass)

# can also do k-folds cross-validation to select best model 

# need to use multiple random starts to demonstrate sufficient replication
# of the maximum likelihood

# Visualize trajectories -------------------------------------------------------
get_probs <- function(age_list) {

  as.data.frame(age_list) %>%
    mutate(class = str_extract(row.names(.), '[0-9]+'))
}

map_df(three_class$probs, get_probs, .id = 'age') %>%
  mutate(age = str_extract(age, '[0-9]+'),
         age = factor(age,levels = c('9', '10', '12', '14', '16', '18'))) %>%
  ggplot(aes(x = age, 
             y = `Pr(1)`,
             color = class)) +
  geom_line(aes(group = class))

map_df(four_class$probs, get_probs, .id = 'age')  %>%
  mutate(age = str_extract(age, '[0-9]+'),
         age = factor(age,levels = c('9', '10', '12', '14', '16', '18'))) %>%
  ggplot(aes(x = age, 
             y = `Pr(1)`,
             color = class)) +
  geom_line(aes(group = class))

map_df(five_class$probs, get_probs, .id = 'age')  %>%
  mutate(age = str_extract(age, '[0-9]+'),
         age = factor(age,levels = c('9', '10', '12', '14', '16', '18'))) %>%
  ggplot(aes(x = age, 
             y = `Pr(1)`,
             color = class)) +
  geom_line(aes(group = class))

