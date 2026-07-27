# Author: Selene Banuelos
# Date: 7/24/2026
# Description:Investigate progression of ever asthma diagnosis variable over time

# setup
library(readstata13)
library(dplyr)
library(tidyr)
library(ggplot2)

# import data ------------------------------------------------------------------
asthma <- read.dta13('data-raw/de_la_Rosa_07.dta',
                     nonint.factors = TRUE,
                     generate.factors = TRUE)

# data wrangling ---------------------------------------------------------------
ever_asthma <- asthma %>%
  # keep asthma related variables only
  select(newid, cham, contains('asth_')) %>%
  # # remove participants that are missing all data between 9Y-18Y
  # filter(!if_all(contains(c('9y', '10Y', '12Y', '14Y', '18Y')), is.na)) %>%
  # make data longer for some data manipulation
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
              names_glue = '{.value}_{age}') %>%
  
# identify participants with no data from 5Y-18Y
no_data <- ever_asthma %>%
  filter(if_all(contains('asth_'), is.na))

# identify participants with some data from 5Y-18Y
some_data <- ever_asthma %>%
  filter(!if_all(contains('asth_'), is.na)) %>%
  pivot_longer(cols = -c(newid, cham),
               names_to = c('.value', 'age'),
               names_pattern = '(.*)_([0-9]+)') %>%
  mutate(age = factor(age, levels = c(5, 7, 9, 10, 12, 14, 16, 18)))

# identify distinct patterns of responses in participants with some data
# across asth_diag_ever_*Y and has_data*Y
response_groups_2 <- ever_asthma %>%
  filter(!if_all(contains('asth_'), is.na)) %>%
  # remove individual participant IDs
  select(-c(newid, cham)) %>%
  # keep unique/distinct response patterns over time
  distinct() %>%
  # create response group IDs
  mutate(group_id = seq(from=1, to=nrow(.))) %>%
  # make data longer for plotting
  pivot_longer(cols = -c(group_id),
               names_to = c('.value', 'age'),
               names_pattern = '(.*)_([0-9]+)') %>%
  # control order of ages displayed in plot
  mutate(age = factor(age, levels = c(5, 7, 9, 10, 12, 14, 16, 18)))

# identify distinct patterns of responses in participants with some data
# across asth_diag_ever_*Y only
response_groups_1 <- ever_asthma %>%
  filter(!if_all(contains('asth_'), is.na)) %>%
  # remove individual participant IDs
  select(-c(newid, cham, contains('has_data'))) %>%
  # keep unique/distinct response patterns over time
  distinct() %>%
  # create response group IDs
  mutate(group_id = seq(from=1, to=nrow(.))) %>%
  # make data longer for plotting
  pivot_longer(cols = -c(group_id),
               names_to = c('.value', 'age'),
               names_pattern = '(.*)_([0-9]+)') %>%
  # control order of ages displayed in plot
  mutate(age = factor(age, levels = c(5, 7, 9, 10, 12, 14, 16, 18)))

# identify any participants that have a yes/no response but have no data recorded
# at given time point or vice versa


# swimmer plots ----------------------------------------------------------------
# for all individual participants 
some_data %>%
  ggplot(aes(x = age, y = newid, group = newid, col = asth_diag_ever)) +
  geom_line() +
  geom_point(shape = 15)
# can't see anything, too many people in one plot

# for distinct response patterns over time 
# from asth_diag_ever_*Y and has_data*Y data
response_groups_2 %>%
  ggplot(aes(x = age, y = group_id, group = group_id, col = asth_diag_ever)) +
  geom_line(size = 1) +
  geom_point(shape = 15) + 
  theme_minimal() +
  scale_color_manual(values = c('No' = '#90db54',
                                'Yes' = 'red'),
                                na.value = '#E9ECEF')

# for distinct response patterns over time 
# from asth_diag_ever_*Y only
response_groups_1 %>%
  ggplot(aes(x = age, y = group_id, group = group_id, col = asth_diag_ever)) +
  geom_line(size = 1) +
  geom_point(shape = 15) + 
  theme_minimal() +
  scale_color_manual(values = c('No' = '#90db54',
                                'Yes' = 'red'),
                     na.value = '#E9ECEF')

# visually, there doesn't appear to be any Yes -> No ever asthma diagnosis
# looking for red to green transitions

# can also make different type of transision plot 
#https://longitudinalanalysis.com/visualizing-transitions-in-time-using-r-and-alluvial-graphs/
# alluvial graphs --------------------------------------------------------------
