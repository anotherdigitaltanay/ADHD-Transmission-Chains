######################## This script cleans and analyses the transmission chains data from all four studies
######################## we do this via an individual participant data (IPD) meta analysis or mega analysis
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




##########################################
### DATASET LOADING AND SOME QUICK DATA CHECKS

## Load the cleaned datasets from all four studies
d1 <- read.csv("meta_analysis_study1data.csv")
d2 <- read.csv("meta_analysis_study2data.csv")
d3 <- read.csv("meta_analysis_study3data.csv")
d4 <- read.csv("meta_analysis_Study4data.csv")

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
### DATA ANALYSIS - Model 1 to see the effect of transmission on self-diagnosis

## Specifying the model with our pre-registered priors (priors were justified in the power analysis with a prior predictive simulation)
# Chain_ID random intercept is already calculating the effect per study chain as each is given a unique identifier so specifying study num as a fixed effect to avoid model divergence
h1 <- brm(Self_Diagnosis ~ 1 + chain_position + Day + study_num + (1 | chain_ID),
          data = cleaned_dataset,
          family = cumulative("probit"),
          c(prior(normal(0, 1.5), class = "Intercept"),
            prior(normal(0, 1), class = "b"),
            prior(exponential(1), class = "sd")),
          warmup = 1000,
          iter = 3500,
          cores = 4,
          seed = 92,
          control = list(adapt_delta = 0.95)
)

# Save & load the model output. 
#saveRDS(h1, "MA_mod1.rds")
h1 <- readRDS("MA_mod1.rds")

## Checking the model output and doing a posterior predictive check
# Model predicting existing dataset quite well
summary(h1)
pp_check(h1, type = "bars")


####### HYPOTHESIS TESTS FOR EFFECTS OF INTEREST
### 

# Defining the six contrasts of interests = i.e. chain position contrasts with chain position 1 
position_contrasts <- c(
  pos2 = "chain_position2  > 0",
  pos3 = "chain_position3  > 0",
  pos4 = "chain_position4  > 0",
  pos5 = "chain_position5  > 0",
  pos6 = "chain_position6  > 0",
  pos7 = "chain_position7  > 0"
)

# Testing these contrasts now
hypothesis(h1, position_contrasts)

## Plot predictions
conditional_effects(h1)



#####################################
### DATA ANALYSIS - Model 1 Logit Variant (Odds Ratios) to see the effect of transmission on self-diagnosis

## Specifying the model
h1.1 <- brm(Self_Diagnosis ~ 1 + chain_position + Day + study_num + (1 | chain_ID),
            data = cleaned_dataset,
            family = cumulative("logit"),
            c(prior(normal(0, 1.5), class = "Intercept"),
              prior(normal(0, 1), class = "b"),
              prior(exponential(1), class = "sd")),
            warmup = 1000,
            iter = 3500,
            cores = 4,
            seed = 92,
            control = list(adapt_delta = 0.95)
)

# Save & load the model output. 
#saveRDS(h1.1, "MA_OR_mod1.rds")
h1.1 <- readRDS("MA_OR_mod1.rds")

## Checking the model output and doing a posterior predictive check
# Model predicting existing dataset quite well
summary(h1.1)
pp_check(h1.1, type = "bars")


# Converting the effects in odds-ratio
# NB: that est. error should not be read here as it is being exponentiated as well as the function is applied across the matrix
gen_draws <- exp(as_draws_df(h1.1)[, c("b_chain_position2", "b_chain_position3", "b_chain_position4", "b_chain_position5", "b_chain_position6", "b_chain_position7")])
t(apply(gen_draws, 2, \(x) c(Estimate = mean(x), Est.Error = sd(x), Q2.5 = quantile(x, 0.025), Q975 = quantile(x, 0.975))))




#####################################
### DATA ANALYSIS - Model 1 Monotonic Variant to see the effect of transmission on self-diagnosis

## Specifying the model
## NB: not calculating correlations between varying effects to avoid parameter explosion and convergence issues
h1.2 <- brm(Self_Diagnosis ~ 1 + mo(chain_position_trend) + Day + study_num + (1 + mo(chain_position_trend) || chain_ID),
            data = cleaned_dataset,
            family = cumulative("probit"),
            c(prior(normal(0, 1.5), class = "Intercept"),
              prior(normal(0, 1), class = "b"),
              prior(exponential(1), class = "sd"),
              prior(dirichlet(2), class = "simo", coef = "mochain_position_trend1")),
            warmup = 1000,
            iter = 3500,
            cores = 4,
            seed = 22,
            control = list(adapt_delta = 0.95)
)


## Save and load the model
#saveRDS(h1.2, "MA_monotonic_mod1.RDS")
h1.2 <- readRDS("MA_monotonic_mod1.RDS") 

# Summary of model and posterior predictive check
summary(h1.2)
pp_check(h1.2)


######## HYPOTHESIS TESTS FOR EFFECTS OF INTEREST
hypothesis(h1.2, "mochain_position_trend > 0", class = "bsp")

## Plot predictions
conditional_effects(h1.2)







#####################################
### DATA ANALYSIS - Model 2 to see the effect of transmission on other-diagnosis

## Specifying the model
h2 <- brm(Other_Diagnosis ~ 1 + chain_position + Day + study_num + (1 | chain_ID),
          data = cleaned_dataset,
          family = cumulative("probit"),
          c(prior(normal(0, 1.5), class = "Intercept"),
            prior(normal(0, 1), class = "b"),
            prior(exponential(1), class = "sd")),
          warmup = 1000,
          iter = 3500,
          cores = 4,
          seed = 43,
          control = list(adapt_delta = 0.95)
)


# Save & load the model output. 
# saveRDS(h2, "MA_mod2.rds")
h2 <- readRDS("MA_mod2.rds")

## Checking the model output and doing a posterior predictive check
# Model predicting existing dataset quite well
summary(h2)
pp_check(h2, type = "bars")


####### HYPOTHESIS TESTS FOR EFFECTS OF INTEREST
## 
hypothesis(h2, position_contrasts)




#####################################
### DATA ANALYSIS - Model 2 Logit Variant to see the effect of transmission on other-diagnosis

## Specifying the model
h2.1 <- brm(Other_Diagnosis ~ 1 + chain_position + Day + study_num + (1 | chain_ID),
            data = cleaned_dataset,
            family = cumulative("logit"),
            c(prior(normal(0, 1.5), class = "Intercept"),
              prior(normal(0, 1), class = "b"),
              prior(exponential(1), class = "sd")),
            warmup = 1000,
            iter = 3500,
            cores = 4,
            seed = 13,
            control = list(adapt_delta = 0.95)
)


# Save & load the model output. 
#saveRDS(h2.1, "MA_OR_mod2.rds")
h2.1 <- readRDS("MA_OR_mod2.rds")

## Checking the model output and doing a posterior predictive check
# Model predicting existing dataset quite well
summary(h2.1)
pp_check(h2.1, type = "bars")

# Converting the effects in odds-ratio
# NB: that est. error should not be read here as it is being exponentiated as well as the function is applied across the matrix
gen_draws <- exp(as_draws_df(h2.1)[, c("b_chain_position2", "b_chain_position3", "b_chain_position4", "b_chain_position5", "b_chain_position6", "b_chain_position7")])
t(apply(gen_draws, 2, \(x) c(Estimate = mean(x), Est.Error = sd(x), Q2.5 = quantile(x, 0.025), Q975 = quantile(x, 0.975))))



#####################################
### DATA ANALYSIS - Model 2 Monotonic Variant to see the effect of transmission on other-diagnosis

## Specifying the model
h2.2 <- brm(Other_Diagnosis ~ 1 + mo(chain_position_trend) + Day + study_num + (1 + mo(chain_position_trend) || chain_ID),
            data = cleaned_dataset,
            family = cumulative("probit"),
            c(prior(normal(0, 1.5), class = "Intercept"),
              prior(normal(0, 1), class = "b"),
              prior(exponential(1), class = "sd"),
              prior(dirichlet(2), class = "simo", coef = "mochain_position_trend1")),
            warmup = 1000,
            iter = 3500,
            cores = 4,
            seed = 22,
            control = list(adapt_delta = 0.95)
)

## Save and load the model
#saveRDS(h2.2, "MA_monotonic_mod2.RDS")
h2.2 <- readRDS("MA_monotonic_mod2.RDS") 

# Summary of model and posterior predictive check
summary(h2.2)
pp_check(h2.2)


######## HYPOTHESIS TESTS FOR EFFECTS OF INTEREST
hypothesis(h2.2, "mochain_position_trend > 0", class = "bsp")

## Plot predictions
conditional_effects(h2.2)





#####################################
### DATA ANALYSIS - Model 3 to see the effect of transmission on other help-seeking

## Specifying the model
h3 <- brm(Other_HelpSeeking ~ 1 + chain_position + Day + study_num + (1 | chain_ID),
          data = cleaned_dataset,
          family = cumulative("probit"),
          c(prior(normal(0, 1.5), class = "Intercept"),
            prior(normal(0, 1), class = "b"),
            prior(exponential(1), class = "sd")),
          warmup = 1000,
          iter = 3500,
          cores = 4,
          seed = 24,
          control = list(adapt_delta = 0.95)
)

# Save & load the model output. 
#saveRDS(h3, "MA_mod3.rds")
h3 <- readRDS("MA_mod3.rds")

## Checking the model output and doing a posterior predictive check
# Model predicting existing dataset quite well
summary(h3)
pp_check(h3, type = "bars")


####### HYPOTHESIS TESTS FOR EFFECTS OF INTEREST
### 
hypothesis(h3, position_contrasts)

## Plot predictions
conditional_effects(h3)



#####################################
### DATA ANALYSIS - Model 3 Logit Variant to see the effect of transmission on other help-seeking

## Specifying the model
h3.1 <- brm(Other_HelpSeeking ~ 1 + chain_position + Day + study_num + (1 | chain_ID),
            data = cleaned_dataset,
            family = cumulative("logit"),
            c(prior(normal(0, 1.5), class = "Intercept"),
              prior(normal(0, 1), class = "b"),
              prior(exponential(1), class = "sd")),
            warmup = 1000,
            iter = 3500,
            cores = 4,
            seed = 24,
            control = list(adapt_delta = 0.95)
)

# Save & load the model output. 
#saveRDS(h3.1, "MA_OR_mod3.rds")
h3.1 <- readRDS("MA_OR_mod3.rds")

## Checking the model output and doing a posterior predictive check
# Model predicting existing dataset quite well
summary(h3.1)
pp_check(h3.1, type = "bars")

# Converting the effects in odds-ratio
# NB: that est. error should not be read here as it is being exponentiated as well as the function is applied across the matrix
gen_draws <- exp(as_draws_df(h3.1)[, c("b_chain_position2", "b_chain_position3", "b_chain_position4", "b_chain_position5", "b_chain_position6", "b_chain_position7")])
t(apply(gen_draws, 2, \(x) c(Estimate = mean(x), Est.Error = sd(x), Q2.5 = quantile(x, 0.025), Q975 = quantile(x, 0.975))))



#####################################
### DATA ANALYSIS - Model 3 Monotonic Variant to see the effect of transmission on other help-seeking

## Specifying the model
h3.2 <- brm(Other_HelpSeeking ~ 1 + mo(chain_position_trend) + Day + study_num + (1 + mo(chain_position_trend) || chain_ID),
            data = cleaned_dataset,
            family = cumulative("probit"),
            c(prior(normal(0, 1.5), class = "Intercept"),
              prior(normal(0, 1), class = "b"),
              prior(exponential(1), class = "sd"),
              prior(dirichlet(2), class = "simo", coef = "mochain_position_trend1")),
            warmup = 1000,
            iter = 3500,
            cores = 4,
            seed = 92,
            control = list(adapt_delta = 0.95)
)

## Save and load the model
#saveRDS(h3.2, "MA_monotonic_mod3.RDS")
h3.2 <- readRDS("MA_monotonic_mod3.RDS") 


# Summary of model and posterior predictive check
summary(h3.2)
pp_check(h3.2)


######## HYPOTHESIS TESTS FOR EFFECTS OF INTEREST
hypothesis(h3.2, "mochain_position_trend > 0", class = "bsp")

## Plot predictions
conditional_effects(h3.2)







