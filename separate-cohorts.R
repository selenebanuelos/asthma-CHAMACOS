# Author: Selene Banuelos
# Date: 7/15/2026
# Description: NO LONGER NEEDED. Was provided confirmed cohort membership 
# information. Keeping this script as an archive. 

# Original description: Try to separate out CHAM1 from CHAM2 participants

# Notes: I know that CHAM2 recruitment started at time point 9Y. Following this, 
# I can assume that there is a group of participant IDs is only present between 
# time points 9Y-18Y (CHAM2) and another group of IDs that are present across
# all times (CHAM1)

# setup
library(readstata13)
library(dplyr)

# import data ------------------------------------------------------------------
# survey data
asthma <- read.dta13('data-raw/de_la_Rosa_06.dta',
                     nonint.factors = TRUE,
                     generate.factors = TRUE)

# data wrangling ---------------------------------------------------------------
# identify participants with data collected for at least one time point
with_data <- asthma %>%
  # focus on asthma variables 
  select(newid, contains('asth_')) %>%
  # remove rows that have NA in all columns other than newid column
  filter(!if_all(-newid, is.na))

# identify participants that have no data over 5Y-18Y
no_data <- asthma %>%
  # focus on asthma variables 
  select(newid, contains('asth_')) %>%
  # remove rows that have NA in all columns other than newid column
  filter(if_all(-newid, is.na)) %>%
  # create vector of participant IDs
  pull(newid)

# participants with data from 9Y-18Y ONLY (missing data from 5Y & 7Y)
nine_18_only <- with_data %>%
  # remove rows that have NA in all columns other than newid column
  filter(if_all(contains(c('_5Y', '7Y')), is.na)) %>%
  # create vector of participant IDs
  pull(newid)
# I believe these are CHAM2, since data collection started at 9Y

# identify participants that were lost to follow up before 9Y
five_7_only <- with_data %>%
  # filter for IDs in 5Y-7Y but missing entirely from 9Y-18Y
  filter(if_all(
    .cols = contains(c('_9Y','_10Y','_12Y','_14Y','_16Y','_18Y')),
    .fns = is.na))%>%
  # create vector of participant IDs
  pull(newid)
# These are definitely part of CHAM1 since CHAM2 didn't start until 9Y

# participants that have data that spans CHAM1 and CHAM2 time points
spans_all <- with_data %>%
  # at least 1 data point at 5Y or 7Y
  filter(!if_all(contains(c('5Y', '7Y')), is.na) & 
           # AND at least 1 data point between 9Y-18Y
           !if_all(contains(c('_9Y','_10Y','_12Y','_14Y','_16Y','_18Y')), is.na)
         ) %>%
  # create vector of participant IDs
  pull(newid)
# I believe these are CHAM1 participants since they span both cohort timelines

# create new data frame with cohorts assigned
cohorts <- data.frame(newid = asthma$newid) %>%
  mutate(cohort = case_when(newid %in% spans_all ~ 'CHAM1',
                            newid %in% five_7_only ~ 'CHAM1',
                            newid %in% nine_18_only ~ 'CHAM2',
                            newid %in% no_data ~ 'no data'))

# save (unconfirmed) cohort info as csv ----------------------------------------
write.csv(cohorts, 'data-processed/cohorts-unconfirmed.csv', row.names = FALSE)