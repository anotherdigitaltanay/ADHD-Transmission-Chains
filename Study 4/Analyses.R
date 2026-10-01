######################## This script cleans and analyses the transmission chains data from Pre-registered Study 4 
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
setwd(file.path(renv::project(), "Study 4")) 

## Set seed to ensure reproducibility
set.seed(45306)




##########################################
### DATASET LOADING AND SOME QUICK DATA CHECKS

## Load the cleaned dataset
d1 <- read.csv("Study4_Data.csv")


################# Some quick data-checks 
###

## Let's remove participants were screened out/didn't finish the survey
# Safe to exclude as the survey was designed to record complete responses and participants weren't allowed to redo the survey (via Qualtrics settings)
## No one did!
d1 <- d1 %>% filter(Progress == "100")
d1 <- d1 %>% filter(Q_TerminateFlag != "Screened")

## Let's check if any participant did the experiment twice
## No one did!
d1 %>% count(Anon_PID) %>% arrange(desc(n))

## Now we check if anyone didn't finish the study (i.e. were screened out or didn't provide their consent)
## Everyone provided their consent and finished the study
d1 %>% count(Consent)
d1 %>% count(Q_TerminateFlag)
d1 %>% count(Progress)

## Check if everyone passed the attention check
## Everyone did
d1 %>% count(Attention.Check.1)

## Now we check the basic pre-screeners (no mental health diagnosis and people only above 18)
# Everyone above 18 and has no current diagnosis
d1 %>% count(Mental.Health.Status)
d1 %>% count(Age.Screener)


## See in chain position 7 whether everyone got a unique stimulus from Study 2 Chain Position 7 assigned to them
## Seems like it
d1 %>% filter(gen == 7) %>% count(Assigned_Chain) %>% 
  arrange(desc(n))

d1 %>% filter(gen == 7) %>% count(Assigned_Text) %>% 
  arrange(desc(n))





##########################################
### DATASET FORMATTING FOR ANALYSIS


## Select columns of interest
cleaned_dataset <- data.frame(
  PID = as.factor(d1$Anon_PID),
  chain = d1$Assigned_Chain,
  chain_position = as.factor(d1$gen),
  stimulus = d1$Assigned_Text,
  self_diagnosis = d1$Self.Diagnosis.Q_1,
  other_diagnosis = d1$Other.Diagnosis.Qs_1,
  Prevalence = d1$Q_Prevalence.Rate,
  Treatment = d1$Q_Treatment.Rate,
  other_helpseeking = d1$Other.Diagnosis.Qs_2,
  age = as.numeric(d1$Age.Screener),
  gender = as.factor(d1$Gender),
  medication = as.factor(d1$Medication),
  treatment = as.factor(d1$Treatment),
  Health_access = as.factor(d1$Health_Access),
  Past_diagnoses = as.factor(d1$Past_MH_diagnosis),
  Priork_friend_help = as.factor(d1$Prior_Know1),
  Priork_friend_helpwait = as.factor(d1$Prior_Know2),
  Priork_socialmedia = as.factor(d1$Prior_Know3),
  AI_self_report = as.factor(d1$AI_Self_Report),
  AI_explanation = d1$AI.freetext,
  StartDate = d1$StartDate,
  prompt_inject = d1$prompt_inject,
  total_chars = d1$total_chars,
  Pasting_examples = as.factor(d1$examples_pasted),
  Example_keys = d1$examples_key_count,
  examples_pasted = d1$examples_pasted
)
  
 


#### Now, we make the final transformations necessary for data analysis
####

## Our aim here is to convert the likert-scale responses to numbers for ordinal analyses. 
# Qualtrics can do this internally while data export, but best to do this ourselves transparently to see if there are any glitches

## Let's first define the verbal responses for self-diagnosis, other-diagnosis and other help-seeking
# Note that all of them share the same response options
likert_levels <- c("Strongly Disagree", "Somewhat disagree", "Neither agree nor disagree", "Somewhat agree", "Strongly agree")

## Making the transformation now
cleaned_dataset <- cleaned_dataset %>% 
  mutate(
    Self_Diagnosis = as.numeric(factor(self_diagnosis, levels = likert_levels, ordered = TRUE)),
    Other_Diagnosis = as.numeric(factor(other_diagnosis, levels = likert_levels, ordered = TRUE)),
    Other_HelpSeeking = as.numeric(factor(other_helpseeking, levels = likert_levels, ordered = TRUE))
  )

## Checking whether the transformations were correctly implemented (counts for transformed and untransformed variables should be the same)
## Looks like they were
cleaned_dataset %>% count(self_diagnosis)
cleaned_dataset %>% count(Self_Diagnosis)

cleaned_dataset %>% count(other_diagnosis)
cleaned_dataset %>% count(Other_Diagnosis)

cleaned_dataset %>% count(other_helpseeking)
cleaned_dataset %>% count(Other_HelpSeeking)

## NB: WE WON"T RUN MONOTONIC MODELS HERE AS THIS IS A SIMPLE BETWEEN-GROUP COMPARISION. 






#####################################
### CREATING AI DETECTION VARIABLE

## Note that we are creating this variable for eventual use in our mega-analysis
## Four methods were used to triangulate AI usage: self-report, prompt-injection, keystroke ratio analysis, copy-paste detector

## NB: For keystroke analysis, the total number of characters used by participants in their written explanation was already computed from their raw response
# The raw responses haven't been shared in the public dataset, owing to some containing PID which we aren't allowed to share. 
# However, all the cleaned responses used during transmission have been shared. 

## Formatting Copy-paste results
cleaned_dataset <- cleaned_dataset %>%
  mutate(
    copy_paste = Pasting_examples
  )

## Method 2: Keystroke analysis
# We compare the number of keystrokes recorded to the total character count of the response.
# A low keystroke-to-character ratio suggests the text was pasted rather than typed.
cleaned_dataset <- cleaned_dataset %>%
  mutate(
    # Extracting the correct keystroke count for the prompt condition
    total_keys = as.numeric(Example_keys),
    
    # Calculating the keystroke-to-character ratio
    # We add 0.001 to the denominator to avoid dividing by zero if total_chars is 0
    key_char_ratio = total_keys / (total_chars + 0.001),
    
    # Flagging massive discrepancies 
    # The ratio should be near 1 or higher for honest typers; can be a bit lower if folks used tablets (e.g. predictive typing); 
    # But below 0.5 seems like potential LLM usage 
    suspected_paste_flag = case_when(
      total_chars > 20 & key_char_ratio < 0.5 ~ TRUE,
      TRUE ~ FALSE
    )
  )


## Compare across the four detection methods
# Here we check whether the four methods agree (self report, copy-paste detection, prompt injection, key-stroke analysis) on which participants are flagged,
# and at which chain position and chain they appear.
ai_detection_comparison <- cleaned_dataset %>%
  mutate(
    flag_copy_paste   = copy_paste == "true",
    flag_self_report  = AI_self_report == "Yes",
    flag_injection    = prompt_inject == TRUE,
    flag_keystroke    = suspected_paste_flag == TRUE
  ) %>%
  filter(flag_copy_paste | flag_self_report | flag_injection | flag_keystroke) %>%
  select(
    PID, chain_position,
    flag_copy_paste, flag_self_report, flag_injection, flag_keystroke
  ) %>%
  arrange(chain_position)


## NB: Note that we won't have to create POST_AI_Introduction variable in the same sense as other studies
# This is because there is no transmission which eliminates the possibility of AI exposure to other participants.
# We still name the variable as post_AI_introduction to homogenise naming convention across scripts and studies
cleaned_dataset <- cleaned_dataset %>%
  mutate(
    post_AI_introduction = (Pasting_examples == "true") | 
      (AI_self_report == "Yes") | 
      (prompt_inject == TRUE) | 
      (suspected_paste_flag == TRUE)
  )






#####################################
### DATA ANALYSIS - Model 1 to see the effect of transmission on self-diagnosis

## Specifying the model with our pre-registered priors (priors were justified in the power analysis with a prior predictive simulation)
h1 <- brm(Self_Diagnosis ~ 1 + chain_position,
          data = cleaned_dataset,
          family = cumulative("probit"),
          c(prior(normal(0, 1.5), class = "Intercept"),
            prior(normal(0, 1), class = "b")),
          warmup = 1000,
          iter = 3500,
          cores = 4,
          seed = 9,
          control = list(adapt_delta = 0.9)
)

# Save & load the model output. 
#saveRDS(h1, "Study4_mod1.rds")
h1 <- readRDS("Study4_mod1.rds")

## Checking the model output and doing a posterior predictive check
# Model predicting existing dataset quite well
summary(h1)
pp_check(h1, type = "bars")


####### HYPOTHESIS TESTS FOR EFFECTS OF INTEREST
### 

# Defining the contrast of interests = i.e. chain position contrasts with chain position 1 
position_contrasts <- c(
  pos7 = "chain_position7  > 0"
)

# Testing these contrasts now
hypothesis(h1, position_contrasts)

## Plot predictions
conditional_effects(h1)




#####################################
### DATA ANALYSIS - Model 1 Logit Variant (Odds Ratios) to see the effect of transmission on self-diagnosis

## Specifying the model
h1.1 <- brm(Self_Diagnosis ~ 1 + chain_position,
            data = cleaned_dataset,
            family = cumulative("logit"),
            c(prior(normal(0, 1.5), class = "Intercept"),
              prior(normal(0, 1), class = "b")),
            warmup = 1000,
            iter = 3500,
            cores = 4,
            seed = 92,
            control = list(adapt_delta = 0.9)
)

# Save & load the model output. 
#saveRDS(h1.1, "Study4_OR_mod1.rds")
h1.1 <- readRDS("Study4_OR_mod1.rds")

## Checking the model output and doing a posterior predictive check
# Model predicting existing dataset quite well
summary(h1.1)
pp_check(h1.1, type = "bars")


# Converting the effects in odds-ratio
# NB: that est. error should not be read here as it is being exponentiated as well as the function is applied across the matrix
gen_draws <- exp(as_draws_df(h1.1)[, c("b_chain_position7")])
t(apply(gen_draws, 2, \(x) c(Estimate = mean(x), Est.Error = sd(x), Q2.5 = quantile(x, 0.025), Q975 = quantile(x, 0.975))))





#####################################
### DATA ANALYSIS - Model 2 to see the effect of transmission on other-diagnosis

## Specifying the model
h2 <- brm(Other_Diagnosis ~ 1 + chain_position,
          data = cleaned_dataset,
          family = cumulative("probit"),
          c(prior(normal(0, 1.5), class = "Intercept"),
            prior(normal(0, 1), class = "b")),
          warmup = 1000,
          iter = 3500,
          cores = 4,
          seed = 43,
          control = list(adapt_delta = 0.9)
)


# Save & load the model output. 
# saveRDS(h2, "Study4_mod2.rds")
h2 <- readRDS("Study4_mod2.rds")

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
h2.1 <- brm(Other_Diagnosis ~ 1 + chain_position,
            data = cleaned_dataset,
            family = cumulative("logit"),
            c(prior(normal(0, 1.5), class = "Intercept"),
              prior(normal(0, 1), class = "b")),
            warmup = 1000,
            iter = 3500,
            cores = 4,
            seed = 13,
            control = list(adapt_delta = 0.9)
)


# Save & load the model output. 
# saveRDS(h2.1, "Study4_OR_mod2.rds")
h2.1 <- readRDS("Study4_OR_mod2.rds")

## Checking the model output and doing a posterior predictive check
# Model predicting existing dataset quite well
summary(h2.1)
pp_check(h2.1, type = "bars")

# Converting the effects in odds-ratio
# NB: that est. error should not be read here as it is being exponentiated as well as the function is applied across the matrix
gen_draws <- exp(as_draws_df(h2.1)[, c("b_chain_position7")])
t(apply(gen_draws, 2, \(x) c(Estimate = mean(x), Est.Error = sd(x), Q2.5 = quantile(x, 0.025), Q975 = quantile(x, 0.975))))





#####################################
### DATA ANALYSIS - Model 3 to see the effect of transmission on other help-seeking

## Specifying the model
h3 <- brm(Other_HelpSeeking ~ 1 + chain_position,
          data = cleaned_dataset,
          family = cumulative("probit"),
          c(prior(normal(0, 1.5), class = "Intercept"),
            prior(normal(0, 1), class = "b")),
          warmup = 1000,
          iter = 3500,
          cores = 4,
          seed = 24,
          control = list(adapt_delta = 0.9)
)

# Save & load the model output. 
#saveRDS(h3, "Study4_mod3.rds")
h3 <- readRDS("Study4_mod3.rds")

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
h3.1 <- brm(Other_HelpSeeking ~ 1 + chain_position,
            data = cleaned_dataset,
            family = cumulative("logit"),
            c(prior(normal(0, 1.5), class = "Intercept"),
              prior(normal(0, 1), class = "b")),
            warmup = 1000,
            iter = 3500,
            cores = 4,
            seed = 24,
            control = list(adapt_delta = 0.9)
)

# Save & load the model output. 
#saveRDS(h3.1, "Study4_OR_mod3.rds")
h3.1 <- readRDS("Study4_OR_mod3.rds")

## Checking the model output and doing a posterior predictive check
# Model predicting existing dataset quite well
summary(h3.1)
pp_check(h3.1, type = "bars")

# Converting the effects in odds-ratio
# NB: that est. error should not be read here as it is being exponentiated as well as the function is applied across the matrix
gen_draws <- exp(as_draws_df(h3.1)[, c("b_chain_position7")])
t(apply(gen_draws, 2, \(x) c(Estimate = mean(x), Est.Error = sd(x), Q2.5 = quantile(x, 0.025), Q975 = quantile(x, 0.975))))



#####################################
### EXPORTING CLEANED DATASET FOR POOLING FOR MEGA-ANALYTIC REGRESSION

## Selecting columns of interest that are vital for merging data across four studies
meta_analytic_dataset <- cleaned_dataset %>%
  select(PID, chain, chain_position, Self_Diagnosis, Other_Diagnosis, Other_HelpSeeking, Prevalence, Treatment, StartDate, post_AI_introduction, stimulus)

## Add experiment number to identify dataset
meta_analytic_dataset$study_num <- "exp4"

## Export dataset
write.csv(meta_analytic_dataset, "meta_analysis_Study4data.csv", row.names = FALSE)














