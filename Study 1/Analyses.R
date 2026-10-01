######### This script cleans and analyses the transmission chains data from Exploratory Study 1 
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
library(psych)

## Set working directory to current study folder
setwd(file.path(renv::project(), "Study 1")) 

set.seed(6534344)



##########################################
### DATASET LOADING AND SOME QUICK CAVEATS

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


#####################################
### CREATING AI DETECTION VARIABLE

## Note that we are creating this variable for eventual use in our mega-analysis
## Note that in this study, there was only one detection method (i.e. participant self-reporting)
## There was a prompt injection, but as mentioned in the supplement, it was quite crap

cleaned_dataset <- cleaned_dataset %>%
  mutate(chain_position_num = as.numeric(as.character(gen))) %>%                # Convert chain_position/gen to numeric for mathematical comparison
  mutate(                                                                       # Checking observations where participants self-reported AI usage
    used_AI_directly =  (AI_usage == "Yes"),
    used_AI_directly = replace_na(used_AI_directly, FALSE)                      # Catch any potential NAs generated by the logical checks
  ) %>%
  
  group_by(chain) %>%                                                           # Grouping by chain to calculate the temporal exposure of AI explanations to people who used AI + other people in the chain
  mutate(
    any_ai_in_chain = any(used_AI_directly),                                    # Did anyone in this specific chain use AI?
    first_ai_pos = ifelse(any_ai_in_chain,                                                        # If yes, what was the lowest (earliest) chain position where it happened?
                          min(chain_position_num[used_AI_directly == TRUE], na.rm = TRUE),         # If no one used it, assign 'Inf' (infinity) so it never triggers the next step
                          Inf),
    
    post_AI_introduction = ifelse(chain_position_num >= first_ai_pos, 1, 0),    # If the current position is >= the first AI position, they get a 1. Otherwise, 0.
    post_AI_introduction = as.factor(post_AI_introduction)
  ) %>%
  ungroup() %>%
  select(-chain_position_num, -any_ai_in_chain, -first_ai_pos)

#### SANITY CHECK to see if this variable was correctly created. 
## Creating a cross-tabulation of direct AI use vs. the new exposure variable: 

## Direct_AI_Use = TRUE & Post_AI_Status = 0: 
# This number MUST be 0. It is logically impossible for someone to have used AI directly but not be considered "post-AI introduction." 
# If there is a number here, the code broke

## Direct_AI_Use = TRUE & Post_AI_Status = 1: 
# This should perfectly match the total number of people who directly used AI (the rows in the ai_detection_comparison table).

## Direct_AI_Use = FALSE & Post_AI_Status = 0: 
# These are the "pure" observations—people in non-AI chains, plus the people sitting upstream in AI chains before the AI was introduced.
table(Direct_AI_Use = cleaned_dataset$used_AI_directly, 
      Post_AI_Status = cleaned_dataset$post_AI_introduction)


#####################################
### TESTING RELIABILITY OF SELF DIAGNOSIS MEASURE

## Self-diagnosis is a relatively new construct 
# In this study, we assessed whether the item used by Sandra et al. (2025) SD-I1 "I believe I have ADHD"
# tracks well with other items from the SELFI-scale

## Reliability check - Sandra et al.'s item seems to do well
cron_alpha <- alpha(cleaned_dataset %>% select(SD_I1, SD_I2, SD_I3, SD_I4))
cron_alpha






#####################################
### DATA ANALYSIS - Model 1 to see the effect of transmission on self-diagnosis

## Specifying the model
mod1 <- brm(SD_I1 ~ 1 + gen + (1 |chain),
            data = cleaned_dataset,
            family = cumulative("probit"),
            prior = c(prior(normal(0, 1.5), class = "Intercept"),
                      prior(normal(0, 1), class = "b"),
                      prior(exponential(1), class = "sd")
            ),
            warmup = 1000,
            iter = 3500,
            seed = 32,
            control = list(adapt_delta = 0.95))

  
## Save and load the model
#saveRDS(mod1, "Study1_mod1.RDS")
mod1 <- readRDS("Study1_mod1.RDS")  

# Summary of model 1 and posterior predictive check
summary(mod1)
pp_check(mod1)


########### HYPOTHESIS TESTS FOR EFFECTS OF INTEREST
#### 

# Defining the four contrasts of interests = i.e. chain position contrasts with chain position 1 
position_contrasts <- c(
  pos2 = "gen2  > 0",
  pos3 = "gen3  > 0",
  pos4 = "gen4  > 0"
)

# Testing these contrasts now
hypothesis(mod1, position_contrasts)

## Plot predictions
conditional_effects(mod1)



#####################################
### DATA ANALYSIS - Model 1 Logit Variant (Odds Ratio) to see the effect of transmission on self-diagnosis

## Specifying the model
mod1.1 <- brm(SD_I1 ~ 1 + gen + (1 |chain),
            data = cleaned_dataset,
            family = cumulative("logit"),
            prior = c(prior(normal(0, 1.5), class = "Intercept"),
                      prior(normal(0, 1), class = "b"),
                      prior(exponential(1), class = "sd")
            ),
            warmup = 1000,
            iter = 3500,
            seed = 43,
            control = list(adapt_delta = 0.95))


## Save and load the model
#saveRDS(mod1.1, "Study1_OR_mod1.RDS")
mod1.1 <- readRDS("Study1_OR_mod1.RDS")  

# Summary of model and posterior predictive check
summary(mod1.1)
pp_check(mod1.1)

# Converting the effects in odds-ratio
# NB: that est. error should not be read here as it is being exponentiated as well as the function is applied across the matrix
gen_draws <- exp(as_draws_df(mod1.1)[, c("b_gen2", "b_gen3", "b_gen4")])
t(apply(gen_draws, 2, \(x) c(Estimate = mean(x), Est.Error = sd(x), Q2.5 = quantile(x, 0.025), Q975 = quantile(x, 0.975))))





#####################################
### DATA ANALYSIS - Model 2 to see the effect of transmission on other-diagnosis

## Specifying the model
mod2 <- brm(OD_I1 ~ 1 + gen + (1 |chain),
            data = cleaned_dataset,
            family = cumulative("probit"),
            prior = c(prior(normal(0, 1.5), class = "Intercept"),
                      prior(normal(0, 1), class = "b"),
                      prior(exponential(1), class = "sd")
            ),
            warmup = 1000,
            iter = 3500,
            seed = 53,
            control = list(adapt_delta = 0.95))

## Save and load the model
#saveRDS(mod2, "Study1_mod2.RDS")
mod2 <- readRDS("Study1_mod2.RDS")  

# Summary of model 2 and posterior predictive check
summary(mod2)
pp_check(mod2)


####### HYPOTHESIS TESTS FOR EFFECTS OF INTEREST
### 
# Testing these contrasts now
hypothesis(mod2, position_contrasts)

## Plot predictions
conditional_effects(mod2)



#####################################
### DATA ANALYSIS - Model 2 Logit Variant (Odds Ratio) to see the effect of transmission on other-diagnosis

## Specifying the model
mod2.1 <- brm(OD_I1 ~ 1 + gen + (1 |chain),
            data = cleaned_dataset,
            family = cumulative("logit"),
            prior = c(prior(normal(0, 1.5), class = "Intercept"),
                      prior(normal(0, 1), class = "b"),
                      prior(exponential(1), class = "sd")
            ),
            warmup = 1000,
            iter = 3500,
            seed = 23,
            control = list(adapt_delta = 0.95))

## Save and load the model
#saveRDS(mod2.1, "Study1_OR_mod2.RDS")
mod2.1 <- readRDS("Study1_OR_mod2.RDS")  

# Summary of model and posterior predictive check
summary(mod2.1)
pp_check(mod2.1)

# Converting the effects in odds-ratio
gen_draws <- exp(as_draws_df(mod2.1)[, c("b_gen2", "b_gen3", "b_gen4")])
t(apply(gen_draws, 2, \(x) c(Estimate = mean(x), Est.Error = sd(x), Q2.5 = quantile(x, 0.025), Q975 = quantile(x, 0.975))))




#####################################
### DATA ANALYSIS - Model 3 to see the effect of transmission on other help-seeking

## Specifying the model
mod3 <- brm(OD_I2 ~ 1 + gen + (1 |chain),
            data = cleaned_dataset,
            family = cumulative("probit"),
            prior = c(prior(normal(0, 1.5), class = "Intercept"),
                      prior(normal(0, 1), class = "b"),
                      prior(exponential(1), class = "sd")
            ),
            warmup = 1000,
            iter = 3500,
            seed = 270,
            control = list(adapt_delta = 0.95))

## Save and load the model
#saveRDS(mod3, "Study1_mod3.RDS")
mod3 <- readRDS("Study1_mod3.RDS")  

# Summary of model 3 and posterior predictive check
summary(mod3)
pp_check(mod3)


####### HYPOTHESIS TESTS FOR EFFECTS OF INTEREST
#### 
# Testing these contrasts now
hypothesis(mod3, position_contrasts)

## Plot predictions
conditional_effects(mod3)




#####################################
### DATA ANALYSIS - Model 3 Logit Variant to see the effect of transmission on other help-seeking

## Specifying the model
mod3.1 <- brm(OD_I2 ~ 1 + gen + (1 |chain),
            data = cleaned_dataset,
            family = cumulative("logit"),
            prior = c(prior(normal(0, 1.5), class = "Intercept"),
                      prior(normal(0, 1), class = "b"),
                      prior(exponential(1), class = "sd")
            ),
            warmup = 1000,
            iter = 3500,
            seed = 30,
            control = list(adapt_delta = 0.95))

## Save and load the model
#saveRDS(mod3.1, "Study1_OR_mod3.RDS")
mod3.1 <- readRDS("Study1_OR_mod3.RDS")  

# Summary of model and posterior predictive check
summary(mod3.1)
pp_check(mod3.1)

# Converting the effects in odds-ratio
gen_draws <- exp(as_draws_df(mod3.1)[, c("b_gen2", "b_gen3", "b_gen4")])
t(apply(gen_draws, 2, \(x) c(Estimate = mean(x), Est.Error = sd(x), Q2.5 = quantile(x, 0.025), Q975 = quantile(x, 0.975))))




#####################################
### EXPORTING CLEANED DATASET FOR POOLING FOR MEGA-ANALYTIC REGRESSION

## Selecting columns of interest that are vital for merging data across four studies
meta_analytic_dataset <- cleaned_dataset %>%
  select(PID, chain, gen, SD_I1, OD_I1, OD_I2, Prevalence, Treatment, StartDate, post_AI_introduction, stimulus)

## Renaming certain columns for smooth binding of datasets
meta_analytic_dataset <- meta_analytic_dataset %>%
  rename(
    "chain_position" = gen,
    "Self_Diagnosis" = SD_I1,
    "Other_Diagnosis" = OD_I1,
    "Other_HelpSeeking" = OD_I2
  )

## Add experiment number to identify dataset
meta_analytic_dataset$study_num <- "exp1"

## Export dataset
write.csv(meta_analytic_dataset, "meta_analysis_study1data.csv", row.names = FALSE)

