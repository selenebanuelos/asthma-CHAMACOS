# Author: Selene Banuelos
# Date: 4/28/2026
# Description: Merge asthma classes and DNAm age estimates

# setup
library(readstata13) # work with STATA DTA files
library(dplyr)

# import data ------------------------------------------------------------------
# participant information
demo <- read.dta13('data-raw/de_la_Rosa_07.dta',
                   nonint.factors = TRUE,
                   generate.factors = TRUE)

# estimated epigenetic age from various clocks
dnam_age <- read.dta13("data-raw/de_la_Rosa_epigenetic_07.dta", 
                 nonint.factors=TRUE, 
                 generate.factors=TRUE)

# asthma classification from latent class variable analysis
asthma <- read.csv('data-raw/chamacos.asthma.lca3.csv', sep = ';') %>% select(-X)

# data wrangling 
################################################################################
# specify which method of epigenetic age calculation to use
use_method <- 'Methscore CpG'

# merge epigenetic age and asthma classifications
dnam_asthma <- dnam_age %>%
  # keep DNAm age generated with method of interest
  filter(Method == use_method) %>%
  # change variable type for merging
  mutate(newid = as.integer(newid)) %>%
  # merge asthma class to estimated epi ages
  left_join(., asthma, by = 'newid')

# output
################################################################################
write.csv(dnam_asthma, 'data-processed/dnam-age-asthma.csv', row.names = FALSE)