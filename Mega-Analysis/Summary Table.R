######################## This script summarises the multilevel structure of the pooled data for the mega-analysis
########################################################################################################################### 
###########################################################################################################################  

#####################################
### PACKAGE LOADING 

## Load the relevant libraries. 
## If these packages aren't installed on your machine, type in the following line "install.packages("INSERT PACKAGE NAME ONE AT A TIME")"
library(dplyr)
library(ggplot2)
library(brms)
library(flextable)
library(bayesplot)
library(patchwork)
library(tidyverse)
library(stringr)
library(marginaleffects)
library(emmeans)
library(bayestestR)

## Set working directory to current study folder
setwd(file.path(renv::project(), "Mega-Analysis")) 

## Set seed to ensure reproducibility
set.seed(60)

## Loading the custom made table theme to homogenise style across the manuscript
source("Table Theme Function.R")



##########################################
### DATASET LOADING AND SOME QUICK DATA CHECKS

## Load the cleaned datasets from all four studies
d1 <- read.csv("meta_analysis_study1data.csv")
d2 <- read.csv("meta_analysis_study2data.csv")
d3 <- read.csv("meta_analysis_study3data.csv")
d4 <- read.csv("meta_analysis_study4data.csv")

## Combine them into one big dataset
cleaned_dataset <- rbind(d1, d2, d3, d4)


## Some data transformations and checks

## Double-checking if no one has participated twice or whether an anonymised PID has been assigned twice by random chance
## Doesn't seem to be the case
cleaned_dataset %>% count(PID) %>% arrange(desc(n))

## Treat certain variables as factor
cleaned_dataset$chain_position <- as.factor(cleaned_dataset$chain_position)
cleaned_dataset$post_AI_introduction <- as.factor(cleaned_dataset$post_AI_introduction)

## Also specify chain position as an ordered factor in a distinct variable
## This will allow us to model monotonic variants of all our models where we can test trends over the course of transmission
cleaned_dataset$chain_position_trend <- factor(cleaned_dataset$chain_position, levels = c("1", "2", "3", "4", "5", "6", "7"), ordered = TRUE)

## Extract the day on which data was collected
cleaned_dataset$StartDate <- as.POSIXct(cleaned_dataset$StartDate, format="%Y-%m-%d %H:%M:%S")
cleaned_dataset$Day <- weekdays(cleaned_dataset$StartDate)
cleaned_dataset %>% count(study_num, Day)

## Labelling clearly that chains with similar numbers actually belong to distinct experiments
## This would allow for effficient multilevel modelling
cleaned_dataset <- cleaned_dataset %>%
  mutate(chain_ID = case_when(
    # Complexity 1 - chain position 7 observations from experiment 4 can be considered as part of chains from experiment 2 as these observations saw material from experiment 2
    # Thus, assigning them to chains from exp2
    study_num == "exp4" & chain_position == "7" ~ paste0("exp2", chain),
    # All chain position 1 observations from experiment 4 are given a unique random chain
    study_num == "exp4" & is.na(chain)          ~ paste0("exp4_solo_", PID),
    # Rest are labelled as per their experiments
    TRUE                                        ~ paste0(study_num, chain)
  ))

## Checking if some of the experiment 2 chains now have 8 observations (as they should owing to our assignment)
## Worked well - 2 observations were assigned to exp2 chains that initially had their data excluded as some folks in those chains had already participated in the experiments
cleaned_dataset %>% count(study_num)
cleaned_dataset %>% count(chain_ID) %>% arrange(desc(n))





#####################################
### DATA STRUCTURE SUMMARISATION


## Which study's chain does each observation sit in after chain reassignment for experiment 4? (differs from study_num only for reassigned rows)
cleaned_dataset <- cleaned_dataset %>%
  mutate(chain_study = str_sub(chain_ID, 1, 4),
         Day = factor(Day, levels = c("Monday", "Tuesday", "Wednesday", "Thursday",
                                      "Friday", "Saturday", "Sunday")))

## Helper function for total number of rows
add_total_row <- function(df, label_col) {
  bind_rows(df, df %>% summarise(across(where(is.numeric), sum)) %>% mutate({{ label_col }} := "Total"))
}


## Output Flextable Table 1: observations per study by chain position
obs_by_position <- cleaned_dataset %>%
  count(study_num, chain_position) %>%
  pivot_wider(names_from = chain_position, values_from = n,
              names_prefix = "pos", names_sort = TRUE, values_fill = 0) %>%
  mutate(Total = rowSums(across(starts_with("pos"))), .after = study_num) %>%
  add_total_row(study_num)

obs_position_ft <- obs_by_position |>
  mutate(study_num = str_remove(study_num, "^exp")) |>
  flextable() |>
  set_header_labels(study_num = "Study", Total = "N",
                    pos1 = "1", pos2 = "2", pos3 = "3", pos4 = "4",
                    pos5 = "5", pos6 = "6", pos7 = "7") |>
  add_header_row(values = c("", "Observations by chain position"), colwidths = c(2, 7)) |>
  align(align = "center", part = "all") |>
  apa_theme(spanner_cols = list(3:9),
            stripe_rows  = seq(2, nrow(obs_by_position), 2),
            caption = "Observations per study (as collected), broken down by chain position")

obs_position_ft


## Output Flextable Table 2: observations per study, by day of data collection

obs_by_day <- cleaned_dataset %>%
  count(study_num, Day) %>%
  pivot_wider(names_from = Day, values_from = n, names_expand = TRUE, values_fill = 0) %>%
  mutate(Total = rowSums(across(where(is.numeric))), .after = study_num) %>%
  add_total_row(study_num)

obs_day_ft <- obs_by_day |>
  mutate(study_num = str_remove(study_num, "^exp")) |>
  flextable() |>
  set_header_labels(study_num = "Study", Total = "N") |>
  add_header_row(values = c("", "Observations by day of collection"), colwidths = c(2, 7)) |>
  align(align = "center", part = "all") |>
  apa_theme(spanner_cols = list(3:9),
            stripe_rows  = seq(2, nrow(obs_by_day), 2),
            caption = "Observations per study by weekday of data collection")

obs_day_ft


## Output Flextable 3: Chain Structure after reassignment
chain_sizes <- cleaned_dataset %>% count(chain_study, chain_ID, name = "n_obs")

## Detect chain sizes across studies
size_cols <- paste0("size", sort(unique(chain_sizes$n_obs)))   # e.g. size1, size4, size7, size8

chain_structure <- chain_sizes %>%
  group_by(chain_study) %>%
  summarise(Chains = n(), Obs = sum(n_obs), .groups = "drop") %>%
  left_join(chain_sizes %>%
              count(chain_study, n_obs) %>%
              pivot_wider(names_from = n_obs, values_from = n,
                          names_prefix = "size", names_sort = TRUE, values_fill = 0),
            by = "chain_study") %>%
  left_join(cleaned_dataset %>% filter(study_num != chain_study) %>% count(chain_study, name = "Received"),
            by = "chain_study") %>%
  left_join(cleaned_dataset %>% filter(study_num != chain_study) %>% count(study_num,   name = "Donated"),
            by = c("chain_study" = "study_num")) %>%
  mutate(across(c(Received, Donated), ~ replace_na(., 0L))) %>%
  add_total_row(chain_study)

## Sanity check: print the size distribution and eyeball it (4s in one study, 8s, 7s and 1s in exp2, 1s in exp4)
chain_sizes %>% count(chain_study, n_obs)

## Now output as table
n_size <- length(size_cols)
chain_structure_ft <- chain_structure |>
  mutate(chain_study = str_remove(chain_study, "^exp")) |>
  flextable() |>
  set_header_labels(chain_study = "Study", Chains = "Chains", Obs = "N",
                    Received = "Received", Donated = "Donated") |>
  set_header_labels(values = setNames(as.list(sub("^size", "", size_cols)), size_cols)) |>
  add_header_row(values = c("", "Chains by size (observations)", "Reassigned observations"),
                 colwidths = c(3, n_size, 2)) |>
  align(align = "center", part = "all") |>
  apa_theme(spanner_cols = list(4:(3 + n_size), (4 + n_size):(5 + n_size)),
            stripe_rows  = seq(2, nrow(chain_structure), 2),
            caption = "Multilevel structure of the pooled dataset after chain reassignment")

chain_structure_ft






