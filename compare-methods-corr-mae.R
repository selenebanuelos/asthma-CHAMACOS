# Author: Selene Banuelos
# Date: 6/30/2026
# Description: Compare correlation and median absolute error between 
# chronological age and DNAm age estimated with several methods for ages 9-18

# setup
library(readstata13) # work with STATA DTA files
library(dplyr)
library(janitor)
library(tidyr)
library(Metrics) # calculate median absolute error
library(rmcorr) # calculate repeat measure correlations
library(gtsummary)

# import data ------------------------------------------------------------------
# estimated epigenetic age using various methods
dnam_age <- read.dta13("data-raw/de_la_Rosa_epigenetic_06.dta", 
                       nonint.factors=TRUE, 
                       generate.factors=TRUE)

# data wrangling ---------------------------------------------------------------
# create vector of methods we want to compare
methods <- c('methylCIPHER', 
             'Morgan Levine PC', 
             'DunedinPACE',
             'Methscore CpG', 
             'Methscore PC')

# create vector of DNAm biomarker names
clocks <- c('Horvath',
            'SkinBlood',
            'Hannum',
            'PhenoAge',
            'GrimAge',
            'DNAmTL',
            'DunedinPACE')

# create vector of timepoints of interest
timepoints <- c('Age 9', 'Age 12', 'Age 14', 'Age 18')

# reformat data frame
reformat <- dnam_age %>%
  # keep estimates from methods of interest
  filter(Method %in% methods,
         # keep estimates from tiempoints of interest
         Timepoint %in% timepoints) %>%
  # remove empty columns
  remove_empty('cols') %>%
  # make data longer to calculate correlations and MAE
  pivot_longer(cols = clocks,
               names_to = 'clock',
               values_to = 'dnam_age'
               ) %>%
  # remove rows that do not contain DNAm age estimates
  filter(!is.na(dnam_age))

# correlation and MAE ----------------------------------------------------------
# wrapper function to use stats::cor.test() within summarise
get_cross_corr <- function(x, # numeric vector
                       y # numeric vector
                       ){
  
  # Pearson's product-moment correlation
  result <- cor.test(x, y, method = c('pearson'))
  
  # get correlation coefficient, r
  r <- result$estimate %>% round(digits = 2)
  
  # get 95% CI for r
  ci <- result$conf.int %>% round(digits = 2)
  
  # return string containing: "r (95% CI)"
  paste0(r, ' (', ci[1], ', ', ci[2], ')')
  
}

# wrapper function to use rmcorr() within dplyr::do()
get_long_corr <- function(df) {
  
  # repeat measures correlation 
  corr <- rmcorr(newid, Chrono_Age, dnam_age, dataset = df)
  
  # correlation coefficient
  r <- corr$r %>% round(digits = 2)
  
  # 95% CI for r
  ci <- corr$CI %>% round(digits = 2)
  
  # return string containing: "r (95% CI)"
  paste0(r, ' (', ci[1], ', ', ci[2], ')')
  
}
  
# calculate r & MAE within each time point, method, and clock combination
crosssectional <- reformat %>%
  # calculate correlation coefficient within each combination
  group_by(Timepoint, Method, clock) %>%
  # these are cross-sectional correlation coefficients and MAE
  summarise(
    # calculate_corr produces: 'r, 95% CI lower bound, 95% CI upper bound'
    pearson = get_cross_corr(Chrono_Age, dnam_age),
    mae = round(mdae(actual = Chrono_Age, predicted = dnam_age), digits = 2)
    ) %>%
  ungroup(.)

# calculate repeat measure correlations
long_corr <- reformat %>%
  # create groups for each method/clock combination
  group_by(Method, clock) %>%
  # calculate repeat measures correlations for each group with rmcorr()
  do(rm_corr = get_long_corr(.)) %>%
  # save all results in a data frame
  data.frame(.) %>%
  ungroup(.)

# calculate longitudinal MAE and combine with rm correlations
longitudinal <- reformat %>%
  group_by(Method, clock) %>%
  # calculate longitudinal mean absolute error
  summarise(mae = round(mdae(actual = Chrono_Age, predicted = dnam_age), 
                        digits = 2)) %>%
  ungroup(.) %>%
  # join with repeat measures correlation
  full_join(long_corr, by = c('Method', 'clock'))

# reformat data frames to look like desired tables -------------------------------
# cross sectional correlation and MAE
cross_table <- crosssectional %>%
  # make data wider so all methods from each timepoint for each clock in one row
  pivot_wider(id_cols = c(Timepoint, clock),
              names_from = Method,
              values_from = c(pearson, mae),
              names_glue = '{Method}_{.value}') %>%
  # change all columns to character type
  mutate(across(everything(), as.character)) %>%
  # replace all NA with '-'
  replace(is.na(.), '-') %>%
  group_by(Timepoint) %>%
  # control order of clocks within each time point
  slice(order(factor(clock, levels = clocks))) %>%
  ungroup(.) %>%
  # control order of time points
  slice(order(factor(Timepoint, levels = c('Age 9',
                                           'Age 12',
                                           'Age 14',
                                           'Age 18')))) %>%
  # control order of columns
  select(Timepoint, clock, contains(c('methylCIPHER', 
                                      'Methscore CpG', 
                                      'Methscore PC', 
                                      'Morgan Levine PC',
                                      'DunedinPACE')))

# longitudinal correlation and MAE
long_table <- longitudinal %>%
  # change all columns to character type
  mutate(across(everything(), as.character)) %>%
  # make data wider so all methods for each clock are in one row
  pivot_wider(id_cols = clock,
              names_from = Method, 
              values_from = c(rm_corr, mae),
              names_glue = '{Method}_{.value}') %>%
  replace(is.na(.), '-') %>%
  # control order of clocks within each time point
  slice(order(factor(clock, levels = clocks))) %>%
  ungroup(.) %>%
  # control order of columns
  select(clock, contains(c('methylCIPHER',
                           'Methscore CpG',
                           'Methscore PC',
                           'Morgan Levine PC',
                           'DunedinPACE')))
  
# output -----------------------------------------------------------------------
# save as .csv
# cross-sectional correlation and MAE
write.csv(cross_table, 
          'data-processed/crosssectional-corr-mae.csv', 
          row.names = FALSE)

# longitudinal correlation and MAE
write.csv(long_table, 
          'data-processed/longitudinal-corr-mae.csv', 
          row.names = FALSE)