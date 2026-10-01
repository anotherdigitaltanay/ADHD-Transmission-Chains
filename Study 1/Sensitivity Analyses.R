######### This script cleans and analyses the transmission chains data from Exploratory Study 1 
######### Here, we treat chain position as an ordered factor which allows us to specify a varying slope on it to account for confounders
########################################################################################################################### 
###########################################################################################################################  

## SOME CONTEXT: Owing to the exploratory nature of this study, its small sample size and differing controls from other pre-registered studies,
## we won't run a proper sensitivity analyses here as this study is probably quite underpowered
## However, to account for confounds, we will run one sensitivity analyses 
## where varying slopes will be specified on chain position to see if the fixed effect still exists after accounting for variability between chains



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
library(psych)

## Set working directory to current study folder
setwd(file.path(renv::project(), "Study 1")) 

set.seed(2334344)



#####################################
### DATASET LOADING AND QUICK DATA CHECKS

## Load the dataset
data <- read.csv("Study1_Data.csv")


### Some quick data-checks

## Let's check if any participant did the experiment twice
## Great, no one did!
data %>% count(Anon_PID, sort = TRUE)

## Now we check if anyone didn't finish the study (i.e. were screened out or didn't provide their consent)
## Everyone provided their consent, but one person was screened out!
data %>% count(Consent)
data %>% count(Q_TerminateFlag)
data %>% count(Progress)

## Let's remove this one individual who got screened out
data <- data %>%
  filter(Q_TerminateFlag == "Complete")

## Check if everyone passed the attention check
## Everyone did
data %>% count(Attention.Check.1)

## Two observations in Study 1 didn't see the stimulus (i.e. mental health information) owing to some programming error. 
data %>% filter(stimulus == "")

## Thankfully, as they are both at the end of the chain, we don't have to exclude the whole chain
# Removing these observations at the end of the chain
data <- data %>% filter(stimulus != "")




#####################################
### DATASET FORMATTING FOR ANALYSES

## Cleaning up column names and converting to appropriate data formats so that analytical pipeline runs smoothly
cleaned_dataset <- data.frame(
  PID = as.factor(data$Anon_PID),
  Age = as.numeric(data$Age.Screener),
  SD_Q1 = data$Self.Diagnosis.Qs_1,
  SD_Q2 = data$Self.Diagnosis.Qs_2,
  SD_Q3 = data$Self.Diagnosis.Qs_3,
  SD_Q4 = data$Self.Diagnosis.Qs_4,
  OD_Q1 = data$Other.Diagnosis.Qs_1,
  OD_Q2 = data$Other.Diagnosis.Qs_2,
  Prevalence = as.numeric(data$Q_Prevalence.Rate),
  Treatment = as.numeric(data$Q_Treatment.Rate),
  Sex = as.factor(data$Q5),
  Prior_medication = as.factor(data$Q11),
  Current_Treatment = as.factor(data$Q12),
  NHS_try = as.factor(data$Q13),
  Prior_knowledge_Q1 = data$Q14_1,
  Prior_knowledge_Q2 = data$Q14_2,
  Prior_knowledge_Q3 = data$Q14_3,
  AI_usage = as.factor(data$Q8),
  AI_reason = data$Q9,
  chain = as.factor(data$chain),
  gen = as.factor(data$gen),
  stimulus = as.factor(data$stimulus),
  prompt = as.factor(data$prompt),
  StartDate = data$StartDate
)

## Now we convert the likert-scale scores to numbers for ordinal analyses

## Define the cutpoints
likert_levels <- c("Strongly Disagree", "Somewhat disagree", "Neither agree nor disagree", "Somewhat agree", "Strongly agree")

## Make the transformation for items sharing similar cutpoints
cleaned_dataset <- cleaned_dataset %>%
  mutate(
    SD_I1 = as.numeric(factor(SD_Q1, levels = likert_levels, ordered = TRUE)),
    SD_I3 = as.numeric(factor(SD_Q3, levels = likert_levels, ordered = TRUE)),
    SD_I4 = as.numeric(factor(SD_Q4, levels = likert_levels, ordered = TRUE)),
    OD_I1 = as.numeric(factor(OD_Q1, levels = likert_levels, ordered = TRUE)),
    OD_I2 = as.numeric(factor(OD_Q2, levels = likert_levels, ordered = TRUE))
  )

## Checking whether the transformations were correctly implemented
## Looks like they were
cleaned_dataset %>% count(SD_Q1)
cleaned_dataset %>% count(SD_I1)

cleaned_dataset %>% count(SD_Q3)
cleaned_dataset %>% count(SD_I3)

cleaned_dataset %>% count(SD_Q4)
cleaned_dataset %>% count(SD_I4)

cleaned_dataset %>% count(OD_Q1)
cleaned_dataset %>% count(OD_I1)

cleaned_dataset %>% count(OD_Q2)
cleaned_dataset %>% count(OD_I2)


## Now note that one of the self-diagnosis items (Q2) is negatively worded. That is, if one strongly agrees with this item, they have a lower tendency to self-diagnose 
## To make this consistent with the other items, we reverse-code it such that scoring higher on this item implies higher self-diagnosis

## Define the cutpoints
rev_likert_levels <- c("Strongly agree",  "Somewhat agree", "Neither agree nor disagree",  "Somewhat disagree", "Strongly Disagree")

## Make the transformation now
cleaned_dataset <- cleaned_dataset %>%
  mutate(
    SD_I2 = as.numeric(factor(SD_Q2, levels = rev_likert_levels, ordered = TRUE))
  )

## Checking whether the transformations were correctly implemented
## Looks like they were
cleaned_dataset %>% count(SD_Q2)
cleaned_dataset %>% count(SD_I2)

## Finally, specify chain position as an ordered factor in a distinct variable
## This will allow us to model monotonic variants of all our models where we can test trends over the course of transmission
cleaned_dataset$gen_trend <- factor(cleaned_dataset$gen, levels = c("1", "2", "3", "4"), ordered = TRUE)



#####################################
### DATA ANALYSIS - Model 1 Monotonic Variant to see the effect of transmission on self-diagnosis

## Specifying the model
## NB: not calculating correlations between varying effects to avoid parameter explosion and convergence issues
mod1 <- brm(SD_I1 ~ 1 + mo(gen_trend) + (1 + mo(gen_trend) || chain),
              data = cleaned_dataset,
              family = cumulative("probit"),
              prior = c(prior(normal(0, 1.5), class = "Intercept"),
                        prior(normal(0, 1), class = "b"),
                        prior(exponential(1), class = "sd"),
                        prior(dirichlet(2), class = "simo", coef = "mogen_trend1")
              ),
              warmup = 1000,
              iter = 3500,
              seed = 42,
              control = list(adapt_delta = 0.95))


## Save and load the model
#saveRDS(mod1, "Study1_monotonic_mod1.RDS")
mod1 <- readRDS("Study1_monotonic_mod1.RDS")  

# Summary of model and posterior predictive check
summary(mod1)
pp_check(mod1)


######## HYPOTHESIS TESTS FOR EFFECTS OF INTEREST
hypothesis(mod1, "mogen_trend > 0", class = "bsp")

## Plot predictions
conditional_effects(mod1)




#####################################
### DATA ANALYSIS - Model 2 monotonic variant to see the effect of transmission on other-diagnosis

## Specifying the model
mod2 <- brm(OD_I1 ~ 1 + mo(gen_trend) + (1 + mo(gen_trend) || chain),
              data = cleaned_dataset,
              family = cumulative("probit"),
              prior = c(prior(normal(0, 1.5), class = "Intercept"),
                        prior(normal(0, 1), class = "b"),
                        prior(exponential(1), class = "sd"),
                        prior(dirichlet(2), class = "simo", coef = "mogen_trend1")
              ),
              warmup = 1000,
              iter = 3500,
              seed = 343,
              control = list(adapt_delta = 0.95))

## Save and load the model
#saveRDS(mod2, "Study1_monotonic_mod2.RDS")
mod2 <- readRDS("Study1_monotonic_mod2.RDS")  

# Summary of model 2 and posterior predictive check
summary(mod2)
pp_check(mod2)


######## HYPOTHESIS TESTS FOR EFFECTS OF INTEREST
hypothesis(mod2, "mogen_trend > 0", class = "bsp")

## Plot predictions
conditional_effects(mod2)


#####################################
### DATA ANALYSIS - Model 3 monotonic variant to see the effect of transmission on other help-seeking

## Specifying the model
mod3 <- brm(OD_I2 ~ 1 + mo(gen_trend) + (1 + mo(gen_trend) || chain),
              data = cleaned_dataset,
              family = cumulative("probit"),
              prior = c(prior(normal(0, 1.5), class = "Intercept"),
                        prior(normal(0, 1), class = "b"),
                        prior(exponential(1), class = "sd"),
                        prior(dirichlet(2), class = "simo", coef = "mogen_trend1")
              ),
              warmup = 1000,
              iter = 3500,
              seed = 453,
              control = list(adapt_delta = 0.95))

## Save and load the model
#saveRDS(mod3, "Study1_monotonic_mod3.RDS")
mod3 <- readRDS("Study1_monotonic_mod3.RDS")  

# Summary of model 3 and posterior predictive check
summary(mod3)
pp_check(mod3)


######## HYPOTHESIS TESTS FOR EFFECTS OF INTEREST
hypothesis(mod3, "mogen_trend > 0", class = "bsp")

## Plot predictions
conditional_effects(mod3)




