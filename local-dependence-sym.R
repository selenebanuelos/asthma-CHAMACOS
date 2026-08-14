# Author: Selene Banuelos with code adapted from Daniel Oberski
# Date: 8/4/2026
# Description: Evaluating local dependence 

# setup
library(readstata13)
library(dplyr)
library(tidyr)
library(poLCA) # previously used for LCA
library(purrr)
library(stringr)
options(scipen = 999)

# Import data ------------------------------------------------------------------
# asthma-related data
asthma <- read.dta13('data-raw/de_la_Rosa_07.dta',
                     nonint.factors = TRUE,
                     generate.factors = TRUE)

# Data wrangling ---------------------------------------------------------------
# do some reformatting to prepare for LCA
asth_sym <- asthma %>%
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
  mutate(
    # convert to factor to use in LCA
    current_asthma = factor(asth_sym,
                            levels = c('Yes', 'No'),
                            labels = c('Yes', 'No')
    )) %>%
  # make data wider to use with poLCA::poLCA()
  pivot_wider(id_cols = c('newid', 'cham'),
              names_from = age,
              values_from = c(asth_sym), 
                              # asth_symed, 
                              # asth_med, 
                              # asth_diag_ever, 
                              # current_asthma),
              names_glue = '{.value}_{age}')

# formula: response ~ predictors
asthma_sym_ind <- as.formula(
  paste(
    "cbind(", 
    paste(grep('asth_sym', names(asth_sym), value = TRUE),
          collapse = ","),
    ") ~ 1")
)

# poLCA object from k = 3 model
set.seed(1234)
fit3 <- poLCA(asthma_sym_ind,
              data = asth_sym,
              nclass = 3,
              # estimate model 1,000 times to search for global maximum
              nrep = 1000, 
              maxiter = 5000,
              na.rm = FALSE,
              verbose = FALSE)

# Author: Daniel Oberski
# Date: 2017-08-01

# Bivariate residual statistic for latent class analysis
# Calculate the BVR for poLCA objects

# Argument: a poLCA object
# Value: a dist object with BVRs
# Example: bvr(fit)

# adapted from Oberski:
bvr_oberski <- function(fit) {
  stopifnot(inherits(fit, "poLCA"))
  
  ov_names <- names(fit$predcell)[1:(ncol(fit$predcell) - 2)]
  ov_combn <- combn(ov_names, 2)
  
  get_bvr <- function(ov_pair) {
    form_obs <- as.formula(paste0("observed ~ ", ov_pair[1], " + ", ov_pair[2]))
    form_exp <- as.formula(paste0("expected ~ ", ov_pair[1], " + ", ov_pair[2]))
    
    counts_obs <- xtabs(form_obs, data = fit$predcell)
    counts_exp <- xtabs(form_exp, data = fit$predcell)
    
    sum((counts_obs - counts_exp)^2 / counts_exp)
  }
  
  bvr_pairs <- apply(ov_combn, 2, get_bvr)
  names(bvr_pairs) <- apply(ov_combn, 2, paste, collapse = "<->")
  
  attr(bvr_pairs, "class") <- "dist"
  attr(bvr_pairs, "Size") <- length(ov_names)
  attr(bvr_pairs, "Labels") <- ov_names
  attr(bvr_pairs, "Diag") <- FALSE
  attr(bvr_pairs, "Upper") <- FALSE
  
  bvr_pairs
}

# adapted from perplexity:
# ---- helper: sample from fitted poLCA model ----
simulate_poLCA <- function(fit, n = nrow(fit$y)) {
  K <- length(fit$P) # number of classes
  N <- n
  J <- ncol(fit$y)
  
  cl <- sample.int(K, size = N, replace = TRUE, prob = fit$P)
  ysim <- matrix(NA_integer_, nrow = N, ncol = J)
  
  for (i in seq_len(N)) {
    k <- cl[i]
    for (j in seq_len(J)) {
      probs <- fit$probs[[j]][k, ]
      ysim[i, j] <- sample.int(length(probs), size = 1, prob = probs)
    }
  }
  
  ysim <- as.data.frame(ysim)
  names(ysim) <- names(fit$predcell)[1:ncol(ysim)]
  ysim
}

# ---- bootstrap p-values ----
bootstrap_bvr_p <- function(fit, # poLCA object
                            B, # number of bootstrap samples
                            nrep, # number of repeated starts for LCA
                            verbose = FALSE) {
  obs_bvr <- bvr_oberski(fit)
  pair_names <- names(obs_bvr)
  boot_mat <- matrix(NA_real_, nrow = B, ncol = length(obs_bvr))
  colnames(boot_mat) <- pair_names
  
  # Initialize progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = B, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  set.seed(1234)
  # refit LCA model B times
  for (b in seq_len(B)) {
    dat_b <- simulate_poLCA(fit)
    fit_b <- poLCA(
      asthma_sym_ind,
      data = dat_b,
      nclass = length(fit$P),
      nrep = nrep,
      verbose = verbose
    )
    boot_mat[b, ] <- as.numeric(bvr_oberski(fit_b))
    
    # suspend execution briefly to update progress bar
    Sys.sleep(0.1)
    
    # set the progress bar to the current state
    setTxtProgressBar(pb, b)
  }
  
  # close the connection
  close(pb) 
  
  pvals <- (colSums(boot_mat >= rep(obs_bvr, each = B), na.rm = TRUE) + 1) / (B + 1)
  
  data.frame(
    pair = pair_names,
    observed_bvr = as.numeric(obs_bvr),
    p_value = as.numeric(pvals),
    stringsAsFactors = FALSE
  )
}

# bootstrap p-values for observed BVRs
bvr_pvalues <- bootstrap_bvr_p(fit = fit3, B = 1000, nrep = 100)