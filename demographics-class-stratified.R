# Author: Selene Banuelos
# Date: 8/19/2026
# Description: Create demographics table, stratified by latent class label

# setup
library(readstata13)
library(dplyr)
library(tidyr)
library(gtsummary)
options(scipen = 999)

# import data ------------------------------------------------------------------
# demographics data
demo <- read.dta13('data-raw/de_la_Rosa_07.dta',
                     nonint.factors = TRUE,
                     generate.factors = TRUE)

# participant classifications with 3-class models
classes <- read.csv('data-processed/classification_k3.csv')

# data wrangling ---------------------------------------------------------------
# recreate sample used in latent class analysis
analytic_sample <- demo %>%
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
  filter(n_distinct(age) > 1) %>% # no participants with <2 measures
  ungroup() %>%
  # make data wide again
  pivot_wider(id_cols = c('newid', 'cham'),
              names_from = age,
              values_from = c(asth_sym, 
                              asth_symed, 
                              asth_med, 
                              asth_diag_ever),
              names_glue = '{.value}_{age}') %>%
  # get participant IDs
  pull(newid)

# manipulate demographics data to match previously generated tables
demo_modified <- demo %>%
  mutate(
    # maternal education
    educcat_mom_mod = case_when(educcat_mom == '>=High School Graduate' ~ '>= High School',
                                educcat_mom == "Don't know" ~ "Don't know",
                                is.na(educcat_mom) ~ NA,
                                # all other responses set to < High school
                                .default = '< High School'),
    # language spoken by mother
    langcat_pg_mod = case_when(langcat_pg == 'Span Most' ~ 'Mostly Spanish',
                               langcat_pg == 'Eng Most' | langcat_pg == 'Both equally' ~ 'Mostly English or Both',
                               langcat_pg == 'Other' ~ 'Other',
                               is.na(langcat_pg) ~ NA),
    # years of maternal U.S. residence at delivery
    USyrcat_dl_mod = case_when(USyrcat_dl == 'Entire life' ~ 'U.S. Born',
                               is.na(USyrcat_dl) ~ NA,
                               # all other responses set to foreign born
                               .default = 'Foreign Born'),
    # maternal adverse childhood experiences
    aces_tot_cat_m_mod = case_when(aces_tot_cat_m == '0' ~ '0',
                                   aces_tot_cat_m == '1' | aces_tot_cat_m == '2' ~ '1-2',
                                   is.na(aces_tot_cat_m) ~ NA,
                                   # all other values set to 3+
                                   .default = '3+'),
    # household income when child was 9 years old
    ipovcat9y_mod = case_when(ipovcat9y == 'At or below Poverty' ~ 'At or below poverty level',
                              is.na(ipovcatbl) ~ NA,
                              # all other values set to above poverty level
                              .default = 'Above poverty level')
  )
  

# merge demographics data with assigned class labels
classified <- demo_modified %>%
  # only keep sample used in LCA
  filter(newid %in% analytic_sample) %>%
  # change newid to integer for joining
  mutate(newid = as.integer(newid)) %>%
  # join demographics data to participant class label data
  left_join(., 
             classes,
             by = 'newid') %>%
  # factor class labels to control order in tables
  mutate(class_label_lca = factor(class_label_lca,
                                  levels = c('never/infrequent',
                                             'late onset',
                                             'persistent')),
         class_label_lcga = factor(class_label_lcga,
                                   levels = c('never/infrequent',
                                              'late onset',
                                              'persistent'))
  )
  
# create table -----------------------------------------------------------------
# characteristics to include in table
vars <- c('momdl_age2', 
          'educcat_mom_mod',
          'marcat_pg',
          'smoke_yn',
          'langcat_pg_mod',
          'USyrcat_dl_mod',
          'yrsusa_dl',
          'aces_tot_cat_m_mod',
          'sex',
          'ipovcat9y_mod'
)

# table stratified on asthma trajectories obtained using latent class analysis
tbl_summary(
  classified,
  include = vars,
  by = class_label_lca,
  type = all_dichotomous() ~ 'categorical', # show stats for yes & no levels
  statistic = list(all_continuous() ~ '{mean} ({sd})'),
  missing_text = 'NA'
  ) %>%
  modify_caption('**Using classes from latent class analysis**')

# table stratified on asthma trajectories obtained using latent class growth models
tbl_summary(
  classified,
  include = vars,
  by = class_label_lcga,
  type = all_dichotomous() ~ 'categorical', # show stats for yes & no levels
  statistic = list(all_continuous() ~ '{mean} ({sd})'),
  missing_text = 'NA'
  ) %>%
  modify_caption('**Using classes from latent class growth model (linear)**')