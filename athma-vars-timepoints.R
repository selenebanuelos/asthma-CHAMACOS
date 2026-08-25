# Author: Selene Banuelos
# Date: 7/15/2026
# Description: Compare asthma-related variables at different time points

# setup
library(readstata13)
library(dplyr)
library(tidyr)
library(gtsummary)
library(gt)
library(ggplot2)

# import data ------------------------------------------------------------------
# asthma variables
asthma <- read.dta13('data-raw/de_la_Rosa_07.dta',
                     nonint.factors = TRUE,
                     generate.factors = TRUE)

# data wrangling ---------------------------------------------------------------
# reformat data to use with gtsummary::tbl_summary()
long <- asthma %>%
  # make table with asthma related variables only
  select(newid, cham, contains('asth_')) %>%
  # make data longer
  pivot_longer(cols = contains('asth_'),
               names_to = c('variable', 'age'),
               names_pattern = 'asth_(.+)_(\\d+Y)$') %>%  
  # format CHAMACOS cohort label
  mutate(cham = paste0('CHAMACOS ', cham))

# spread data wider
wide <- long %>%
  pivot_wider(id_cols = c(newid, cham, age),
              names_from = variable,
              values_from = value)

# factor time point variable to control order displayed in table
long$age <- factor(long$age,
                            levels = c('5Y','7Y','9Y','10Y','12Y','14Y','16Y','18Y'))

wide$age <- factor(wide$age,
                   levels = c('5Y','7Y','9Y','10Y','12Y','14Y','16Y','18Y'))

# table ------------------------------------------------------------------------
# overall table
wide %>%
  tbl_summary(include = c(sym, med, diag_ever), 
              by = age,
              label = list(sym = 'asth_sym_*y',
                           med = 'asth_med_*y',
                           diag_ever = 'asth_diag_ever_*y'),
              # show yes and no n(%)
              type = all_dichotomous() ~ 'categorical',
              missing_text = 'Missing',
              # show missing as n (%)
              missing_stat = '{N_miss} ({p_miss}%)') %>%
  # change header to 'Variable'
  modify_header(label = '**Variable**') %>%
  modify_spanning_header(all_stat_cols() ~ "**Overall**")

# table stratified by CHAMACOS cohort
wide %>%
  tbl_strata(
    strata = cham,
    .tbl_fun = 
      ~ .x %>%
      tbl_summary(include = c(sym, med, diag_ever), 
                  by = age,
                  label = list(sym = 'asth_sym_*y',
                               med = 'asth_med_*y',
                               diag_ever = 'asth_diag_ever_*y'),
                  # show yes and no n(%)
                  type = all_dichotomous() ~ 'categorical',
                  missing_text = 'Missing',
                  # show missing as n (%)
                  missing_stat = '{N_miss} ({p_miss}%)'),
    .header = '**{strata}**'
  ) %>%
  # change header to 'Variable'
  modify_header(label = '**Variable**')

# coutns spaghetti plots -------------------------------------------------------
# plot 'yes' responses over time for each of 3 vars, stratified by cohort
long %>%
  group_by(cham, variable, age) %>%
  summarise(yes_count = sum(value == 'Yes', na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = age, y = yes_count, group = variable)) +
  geom_line(aes(color = variable)) +
  facet_wrap(vars(cham)) +
  theme_minimal()

# plot 'no' responses over time for each of 3 vars, stratified by cohort
long %>%
  group_by(cham, variable, age) %>%
  summarise(no_count = sum(value == 'No', na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = age, y = no_count, group = variable)) +
  geom_line(aes(color = variable)) +
  facet_wrap(vars(cham)) +
  theme_minimal()

# plot 'missing' responses over time for each
long %>%
  group_by(cham, variable, age) %>% 
  summarise(missing_count = sum(is.na(value))) %>%
  ungroup() %>%
  ggplot(aes(x = age, y = missing_count, group = variable)) +
  geom_line(aes(color = variable)) +
  facet_wrap(vars(cham)) +
  theme_minimal()

# output -----------------------------------------------------------------------
# save tables
#gtsave(no_missing, 'figures/asthma-vars-overall.png')