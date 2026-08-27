# Author: Selene Banuelos
# Date: 7/24/2026
# Description:Investigate response frequencies for ever asthma diagnosis over time
# Since this is an ever/never-type question, want to make sure responses all
# follow no -> yes patterns over time. Checking to make sure no one has any
# yes -> no patterns

# setup
library(readstata13)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggalluvial) # for alluvial plots

# import data ------------------------------------------------------------------
# asthma related variables in CHAMACOS participants
asthma <- read.dta13('data-raw/de_la_Rosa_07.dta',
                     nonint.factors = TRUE,
                     generate.factors = TRUE)

# data wrangling ---------------------------------------------------------------
# clean up data on ever asthma diagnosis by doctor
ever_asthma <- asthma %>%
  # keep asthma related variables only
  select(newid, cham, contains('asth_')) %>%
  # make data longer for easier data manipulation
  pivot_longer(cols = contains('asth_'),
               names_to = c('.value', 'age'),
               names_pattern = '(asth_.*)_([0-9]+Y)$') %>%
  # has_data == 1 if there is at least 1 data point collected at a given age
  # has_data == 0 if all variables at a given age are NA
  mutate(has_data = as.numeric(if_all(contains('asth_'), 
                                      ~ !is.na(.x))
                               )
         ) %>%
  # remove asthma symptoms and asthma med use data
  select(-c(asth_sym, asth_symed, asth_med)) %>%
  # make data wider
  pivot_wider(id_cols = c(newid, cham),
              names_from = age,
              values_from = c(asth_diag_ever, has_data),
              names_glue = '{.value}_{age}')

# identify transitions over time for alluvial graph
transitions <- ever_asthma %>%
  # identify and count all possible transitions that take place in data
  count(asth_diag_ever_5Y, 
        asth_diag_ever_7Y,
        asth_diag_ever_9Y,
        asth_diag_ever_10Y,
        asth_diag_ever_14Y,
        asth_diag_ever_16Y,
        asth_diag_ever_18Y) %>%
  # create IDs for each transition pattern
  mutate(id = row_number()) %>%
  # make data longer for plotting
  pivot_longer(cols = contains('asth'),
               names_to = 'timepoint',
               values_to = 'response') %>%
  # clean up time point values and make them numeric ages
  mutate(age = factor(str_extract(timepoint, '[0-9]+'), 
                      levels = c('5', '7', '9', '10', '14', '16', '18')))

# alluvial graphs --------------------------------------------------------------
# for all participants
transitions %>%
  ggplot(aes(x = age, y = n, stratum = response, fill = response, alluvium = id)) +
  geom_stratum(alpha = 0.5) +
  geom_flow() +
  theme_minimal() +
  # control colors
  scale_fill_manual(values = c('No' = 'lightpink',
                                'Yes' = 'blue'),
                     na.value = '#E9ECEF')
# do not have any Yes -> No transitions in ever asthma diagnosis over time