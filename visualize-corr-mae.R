### Author: Selene Banuelos
### Date: 6/4/2026
### Description: Create barplots to visualize correlation and median absolute
### error between DNAm age and chronological age

# setup
library(readstata13) # work with STATA DTA files
library(dplyr)
library(tidyr)
library(purrr)
library(rmcorr) # repeated measures correlation
library(Metrics) # calculate median absolute error
library(ggplot2)

# import data
################################################################################
# estimated epigenetic age from various clocks
dnam_age <- read.dta13("data-raw/de_la_Rosa_epigenetic_06.dta", 
                       nonint.factors=TRUE, 
                       generate.factors=TRUE)

# data wrangling
################################################################################
# specify which method of epigenetic age calculation to use
use_method <- 'Methscore CpG'

# keep method of interest and time point of interest
ages_long <- dnam_age %>%
  filter(Method == use_method,
         # only look at ages 9 and 18
         Timepoint == 'Age 9' | Timepoint == 'Age 18') %>%
  # remove predictors that have NA across all rows
  janitor::remove_empty(which = 'cols') %>%
  # remove any estimates for phenotypes other than chronological age
  select(-DNAmTL) %>%
  # make data longer for calculations and plotting
  pivot_longer(cols = !c(newid, Chrono_Age, Timepoint, Array, Method),
               names_to = 'clock',
               values_to = 'dnam_age')

# calculate DNAm age ~ chrono age correlations 
################################################################################
# within time point correlations
cor <- ages_long %>%
  # calculate correlations at each timepoint separately for each clock
  group_by(Timepoint, clock) %>%
  summarize(corr = round(cor(Chrono_Age, dnam_age), digits = 2)) %>%
  ungroup()

# longitudinal correlations accounting for repeated measures among subjects
long_cor <- ages_long %>%
  # factor participant ID
  mutate(newid = as.factor(newid)) %>%
  # split into multiple data frames based on clock, save in list
  split(f = .$clock) %>%
  # calculate repeated measures correlation coefficient w/ rmcorr()
  lapply(function(df) rmcorr(newid, Chrono_Age, dnam_age, dataset = df)) %>%
  # extract first element (corr coeff) in each list and save as data frame
  map_df(., 1) %>%
  # make longer for plotting
  pivot_longer(cols = everything(),
               names_to = 'clock',
               values_to = 'corr') %>%
  mutate(corr = round(corr, digits = 2))

# calculate DNAm age - chrono age median absolute error (MAE)
################################################################################
# calculate within time point MAE
mae <- ages_long %>%
  # calculate MAE within each timepoint separately for each clock
  group_by(Timepoint, clock) %>%
  summarize(mae = round(mdae(actual = Chrono_Age, predicted = dnam_age), 
                        digits = 2)) %>%
  ungroup()

# calculate MAE across all time points
long_mae <- ages_long %>%
  # calculate MAE separately for each clock
  group_by(clock) %>%
  summarize(mae = round(mdae(actual = Chrono_Age, predicted = dnam_age), 
                        digits = 2)) %>%
  ungroup()

# data visualization
################################################################################
# barplots to visualize correlation
cor_plot <- cor %>%
  ggplot(aes(x = clock,
             y = corr, 
             fill = clock)) +
  # create bar plot
  geom_bar(stat = "Identity") +
  # show Pearson correlation coefficient  above each bar
  geom_text(aes(label = corr), position = position_dodge(width = 0.9), vjust = -0.25) +
  # facet plots
  facet_wrap(vars(Timepoint)) +
  # formatting
  labs(title = 'Correlation between epigenetic and chronological age',
       x = 'Clock',
       y = 'Pearson Correlation Coefficient (r)') +
  theme_light() +
  theme(strip.text = element_text(face = 'bold', size = 12),
        strip.background = element_rect(fill = 'darkgrey'),
        legend.position = 'none')

long_cor_plot <- long_cor %>%
  ggplot(aes(x = clock,
             y = corr, 
             fill = clock)) +
  # create bar plot
  geom_bar(stat = "Identity") +
  # show Pearson correlation coefficient  above each bar
  geom_text(aes(label = corr), position = position_dodge(width = 0.9), vjust = -0.25) +
  # formatting
  labs(title = 'Longitudinal correlation between epigenetic and chronological age',
       x = 'Clock',
       y = 'Pearson Correlation Coefficient (r)') +
  theme_light() +
  theme(strip.text = element_text(face = 'bold', size = 12),
        strip.background = element_rect(fill = 'darkgrey'),
        legend.position = 'none')

# barplots to visualize median absolute error
mae_plot <- mae %>%
  ggplot(aes(x = clock,
             y = mae, 
             fill = clock)) +
  # create bar plot
  geom_bar(stat = "Identity") +
  # show MAE above each bar
  geom_text(aes(label = mae), position = position_dodge(width = 0.9), vjust = -0.25) +
  # facet plots
  facet_wrap(vars(Timepoint)) +
  # formatting
  labs(title = 'Median absolute error between epigenetic and chronological age',
       x = 'Clock',
       y = 'Median absolute error') +
  theme_light() +
  theme(strip.text = element_text(face = 'bold', size = 12),
        strip.background = element_rect(fill = 'darkgrey'),
        legend.position = 'none')

long_mae_plot <- long_mae %>%
  ggplot(aes(x = clock,
             y = mae, 
             fill = clock)) +
  # create bar plot
  geom_bar(stat = "Identity") +
  # show MAE above each bar
  geom_text(aes(label = mae), position = position_dodge(width = 0.9), vjust = -0.25) +
  # formatting
  labs(title = 'Longitudinal median absolute error between epigenetic and chronological age',
       x = 'Clock',
       y = 'Median absolute error') +
  theme_light() +
  theme(strip.text = element_text(face = 'bold', size = 12),
        strip.background = element_rect(fill = 'darkgrey'),
        legend.position = 'none')

# output 
################################################################################
