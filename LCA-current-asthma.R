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
curr_asth_only <- asthma %>%
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
  dplyr::select(contains('current_asthma'))

# tetrachoric correlation - incorrect!
curr_asth_only %>%
  # calculate tetrachoric correlations
  sirt::tetrachoric2(.) %>%
  # view correlation matrix
  .$rho %>%
  # visualize correlation matrix
  ggcorrplot(., 
             method = 'circle', 
             type = 'upper', 
             lab = TRUE,
             circle.scale = 2)

# Pearson correlation
# in this case (comparing two binary variables) Pearson coefficient equivalent 
# to phi coefficient/matthews correlation coefficient
curr_asth_plot <- curr_asth_only %>%
  # calculate Pearson correlation coefficients for complete cases
  cor(., 
      method = 'pearson', # equivalent to phi coefficient in binary case
      use = 'complete.obs') %>%
  # visualize correlation matrix
  ggcorrplot(., 
             method = 'circle', 
             type = 'upper', 
             lab = TRUE,
             circle.scale = 2)

ggsave('figures/current_asth_corr.png', curr_asth_plot)
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
        nrep = 1000, # estimate model 1,000 times to search for global maximum
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

write.csv(model_fit, 'data-processed/curr_asth_model_fit.csv', row.names = FALSE)
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

  #model = models[[1]]
  
  # counts of predicted class membership, by modal assignment
  table(model$predclass) %>% 
    prop.table(.) %>% 
    round(., digits = 3)
  
  # # number of classes assumed
  # k <- length(model$P)
  # 
  # # respective size of each latent class
  # class_prop <- sort(model$P) %>% round(digits = 2)
  # 
  # # labels for each class
  # class_index <- seq(1:k)
  # 
  # # named vector of respective latent class sizes
  # setNames(class_prop, class_index)
  
}

# check latent class sizes - don't want small classes
lapply(models, compare_sizes)

# all k produce several small classes
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

# save
ggsave('figures/spaghetti-all-traj.png', spaghetti_all)

# visualize class profiles of k=4 model only
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

# save
ggsave('figures/spaghetti-k3-traj.png', spaghetti_k3)

# visualize class profiles of k=4 model only
spaghetti_k4 <- map_df(models, get_probs) %>% # get class-conditional outcome probabilities
  # keep only data from k=4 model
  filter(k == 4) %>%
  # make spaghetti plot
  ggplot(aes(x = age, # show trend over ages
             y = `Pr(1)`, # plot probability of having current asthma as outcome
             color = class # show each class in different color
  )) +
  geom_line(aes(group = class), size = 1) +
  labs(y = 'Probability of current asthma',
       title = 'LCA with 4 Classes') +
  theme_minimal() +
  theme(panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.5))

# save
ggsave('figures/spaghetti-k4-traj.png', spaghetti_k4)

# Check local dependency -------------------------------------------------------
# Detect local dependence with bivariate residuals
remotes::install_github(repo = "quantmeth/poLCAExtra")
library(poLCAExtra)

# poLCA object from k = 3 model
fit3 <- poLCA(current_asthma_ind,
                  data = current_asthma,
                  nclass = 3,
                  # estimate model 1,000 times to search for global maximum
                  nrep = 1000, 
                  maxiter = 5000,
                  na.rm = FALSE,
                  verbose = FALSE)

# analysis of residual (also known as Tech10 in Mplus) to investiaget local 
# independence
poLCA.residual.pattern(fit3)
# source code: https://github.com/quantmeth/poLCAExtra/blob/main/R/tech10.R
# trying to figure out why this isn't working

# covariance matrix
poLCA.cov(model_k3, nclass = 3)




# Pearson/phi/spearman correlation coefficient is undefined (cov(x,y) = 0, which
# can't be divided by standard deviations to obtain correlation) when comparing 
# two vectors in which one of the vectors has all 0's or 1's. Which is the case 
# in some of the two-way comparison between indicators under predicated class 
# assignments.

# Can also calculate mutual information values to evaluate dependency between
# binary variables (not restricted to linear dependence like in Pearson)
# I think mutual information values can be calculated even when a vector has
# no variability since MI value calculation does not rely on variability?



# function that generates correlation matrix for indicators within each

# calculate correlation matrix between indicators
pred_classes_3 <- current_asthma %>%
  # change all 'Yes' to 1, 'No to 0 for use in sirt::tetrachoric2()
  mutate(across(contains('asth'), ~ case_when(. == 'Yes' ~ 1,
                                              . == 'No' ~ 0,
                                              is.na(.) ~ NA))) %>%
  # add predicted class memberships to participant data
  mutate(class = model_k3$predclass)

# calculate mutual information matrix for current_asthma vars within each class
MI_check <- pred_classes_3 %>%
  filter(class == 3) %>% 
  # keep observations from given class
  dplyr::select(contains('current_asthma')) %>%
  infotheo::mutinformation(., method = 'emp')

# get entropy of each variable
entropy <- pred_classes_3 %>%
  filter(class == 3) %>% 
  # keep observations from given class
  dplyr::select(contains('current_asthma')) %>%
  apply(., 2, infotheo::entropy, method = 'emp')
  
# symmetric normalization of MI value matrix
norm_mi_sym <- 2 * MI_check / outer(entropy, entropy, '+')

split_check <- pred_classes_3 %>%
  group_split(class) %>%
  lapply(., dplyr::select, contains('current_asthma')) %>%
  lapply(., infotheo::mutinformation, method = 'emp')

pred_classes_3 %>%
  group_split(class) %>%
  lapply(., dplyr::select, contains('current_asthma')) %>%

# Pearson correlation
pred_classes_3 %>%
  filter(class == 1) %>% 
  # keep observations from given class
  dplyr::select(contains('current_asthma')) %>%
  # calculate Pearson correlation coefficients for complete cases
  cor(., 
      method = 'pearson', # equivalent to phi coefficient in binary case
      use = 'pairwise.complete.obs') %>%
  # visualize correlation matrix
  ggcorrplot(., 
             method = 'circle', 
             type = 'upper', 
             lab = TRUE,
             circle.scale = 2)

pred_classes_3 %>%
  filter(class == 2) %>% 
  # keep observations from given class
  dplyr::select(contains('current_asthma')) %>%
  # calculate Pearson correlation coefficients for complete cases
  cor(., 
      method = 'pearson', # equivalent to phi coefficient in binary case
      use = 'pairwise.complete.obs') %>%
  # visualize correlation matrix
  ggcorrplot(., 
             method = 'circle', 
             type = 'upper', 
             lab = TRUE,
             circle.scale = 2)
# undefined due to some variables having no variance (warning: SD = 0)

pred_classes_3 %>%
  filter(class == 3) %>% 
  # keep observations from given class
  dplyr::select(contains('current_asthma')) %>%
  # calculate Pearson correlation coefficients for complete cases
  cor(., 
      method = 'pearson', # equivalent to phi coefficient in binary case
      use = 'pairwise.complete.obs') %>%
  # visualize correlation matrix
  ggcorrplot(., 
             method = 'circle', 
             type = 'upper', 
             lab = TRUE,
             circle.scale = 2)
# undefined due to some variables having no variance (warning: SD = 0)