######################## This script cleans and analyses the transmission chains data from Pre-registered Study 2 
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
library(marginaleffects)
library(emmeans)
library(bayestestR)

## Set working directory to current study folder
setwd(file.path(renv::project(), "Study 2")) 

## Set seed to ensure reproducibility
set.seed(94343)




##########################################
### DATASET LOADING AND SOME QUICK DATA CHECKS

## Load the cleaned dataset
d1 <- read.csv("Study2_Data.csv")


## Some quick data-checks


## Let's check if any participant did the experiment twice
## One person did - we will have to remove one of their whole chain as the experiment is otherwise ruined
d1 %>% count(Anon_PID) %>% arrange(desc(n))
d1 %>% filter(Anon_PID == "5YDtoEvoJO")

## After checking on Prolific, the participant did the experiment the second time in chain 119
# We remove this chain
d1 <- d1 %>% filter(chain != 119)

## Note that two participants who had already participated in Study 1 accidentally participated in this study.
## We removed their chains completely. These exclusions were performed privately as we identified them using Prolific IDs which we can't share here publicly
## Hence, the loaded dataset is missing 14 observations intially

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

## Finally, we check whether all chains have 7 participants each and whether both prompt types have observations divided equally
# Seems to be the case
d1 %>% count(prompt)
d1 %>% count(prompt, chain) %>% arrange(desc(n))



##########################################
### DATASET FORMATTING FOR ANALYSIS


## Select columns of interest
cleaned_dataset <- data.frame(
  PID = as.factor(d1$Anon_PID),
  chain = as.factor(d1$chain),
  chain_position = as.factor(d1$gen),
  stimulus = d1$stimulus,
  self_diagnosis = d1$Self.Diagnosis.Q_1,
  other_diagnosis = d1$Other.Diagnosis.Qs_1,
  Prevalence = d1$Q_Prevalence.Rate,
  Treatment = d1$Q_Treatment.Rate,
  other_helpseeking = d1$Other.Diagnosis.Qs_2,
  transmission_type = d1$prompt,
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
  Pasting_examples = as.factor(d1$examples_pasted),
  Pasting_words = as.factor(d1$words_pasted),
  StartDate = d1$StartDate,
  prompt_inject = d1$prompt_inject,
  total_chars = d1$total_chars,
  Example_keys = d1$examples_key_count,
  Words_keys = d1$words_key_count,
  words_pasted = d1$words_pasted,
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

## Finally, specify chain position as an ordered factor in a distinct variable
## This will allow us to model monotonic variants of all our models where we can test trends over the course of transmission
cleaned_dataset$chain_position_trend <- factor(cleaned_dataset$chain_position, levels = c("1", "2", "3", "4", "5", "6", "7"), ordered = TRUE)



#####################################
### CREATING AI DETECTION VARIABLE

## Note that we are creating this variable for eventual use in our mega-analysis
## Four methods were used to triangulate AI usage: self-report, prompt-injection, keystroke ratio analysis, copy-paste detector

## NB: For keystroke analysis, the total number of characters used by participants in their written explanation was already computed from their raw response
# The raw responses haven't been shared in the public dataset, owing to some containing PID which we aren't allowed to share. 
# However, all the cleaned responses used during transmission have been shared. 


## The copy-paste detection results are currently split by prompt type and are in two separate columns
# This poses issues in summarising the results across chain positions. Let's resolve this by combining these results into one column
cleaned_dataset <- cleaned_dataset %>%
  mutate(
    # Necessary step as misses in both columns are currently empty. Need to convert to NAs so that convenience function that combines these values does it properly
    Pasting_examples = na_if(Pasting_examples, ""),
    Pasting_words = na_if(Pasting_words, ""),
    
    # Now combine the values into one column
    copy_paste = coalesce(Pasting_examples, Pasting_words)
  )


## Compute keystroke ratio
cleaned_dataset <- cleaned_dataset %>%
  mutate(
    
    # Extract the correct keystroke count for the prompt condition
    total_keys = case_when(
      transmission_type == "examples" ~ as.numeric(Example_keys),
      transmission_type == "words" ~ as.numeric(Words_keys),
      TRUE ~ 0
    ),
    
    # Calculate the keystroke-to-character ratio
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
    PID, chain, chain_position,
    flag_copy_paste, flag_self_report, flag_injection, flag_keystroke
  ) %>%
  arrange(chain, chain_position)


## Creating the Post_AI Introduction variable based on these metrics
# Anyone who used AI directly gets a variable value of 1; anyone in a chain who is situated after this AI user also gets a 1 (as they have been exposed to AI output)
# Remaining observations get a value of 0; these are the unpolluted non-AI observation
cleaned_dataset <- cleaned_dataset %>%
  mutate(chain_position_num = as.numeric(as.character(chain_position))) %>%
  mutate(                                                                      # Consolidate all four detection methods into one variable
    used_AI_directly = (copy_paste == "true") | 
      (AI_self_report == "Yes") | 
      (prompt_inject == TRUE) | 
      (suspected_paste_flag == TRUE),
    used_AI_directly = replace_na(used_AI_directly, FALSE)                      # Catch any potential NAs generated by the logical checks
  ) %>%
  
  group_by(chain) %>%
  mutate(
    any_ai_in_chain = any(used_AI_directly),                                                 # Did anyone in this specific chain use AI?
    first_ai_pos = ifelse(any_ai_in_chain,                                                   # If yes, what was the lowest (earliest) chain position where it happened?
                          min(chain_position_num[used_AI_directly == TRUE], na.rm = TRUE),   # If no one used it, assign 'Inf' (infinity) so it never triggers the next step
                          Inf),
    
    post_AI_introduction = ifelse(chain_position_num >= first_ai_pos, 1, 0),
    post_AI_introduction = as.factor(post_AI_introduction)
  ) %>%
  ungroup() %>%
  select(-chain_position_num, -any_ai_in_chain, -first_ai_pos)






#####################################
### DATA ANALYSIS - Model 1 to see the effect of transmission on self-diagnosis

## Specifying the model with our pre-registered priors (priors were justified in the power analysis with a prior predictive simulation)
h1 <- brm(Self_Diagnosis ~ 1 + chain_position + (1 |chain),
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
#saveRDS(h1, "Study2_mod1.rds")
h1 <- readRDS("Study2_mod1.rds")

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
h1.1 <- brm(Self_Diagnosis ~ 1 + chain_position + (1 |chain),
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
#saveRDS(h1.1, "Study2_OR_mod1.rds")
h1.1 <- readRDS("Study2_OR_mod1.rds")

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
h1.2 <- brm(Self_Diagnosis ~ 1 + mo(chain_position_trend) + (1 + mo(chain_position_trend) || chain),
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
#saveRDS(h1.2, "Study2_monotonic_mod1.RDS")
h1.2 <- readRDS("Study2_monotonic_mod1.RDS") 

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
h2 <- brm(Other_Diagnosis ~ 1 + chain_position + (1 |chain),
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
# saveRDS(h2, "Study2_mod2.rds")
h2 <- readRDS("Study2_mod2.rds")

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
h2.1 <- brm(Other_Diagnosis ~ 1 + chain_position + (1 |chain),
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
# saveRDS(h2.1, "Study2_OR_mod2.rds")
h2.1 <- readRDS("Study2_OR_mod2.rds")

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
h2.2 <- brm(Other_Diagnosis ~ 1 + mo(chain_position_trend) + (1 + mo(chain_position_trend) || chain),
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
#saveRDS(h2.2, "Study2_monotonic_mod2.RDS")
h2.2 <- readRDS("Study2_monotonic_mod2.RDS") 

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
h3 <- brm(Other_HelpSeeking ~ 1 + chain_position + (1 |chain),
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
#saveRDS(h3, "Study2_mod3.rds")
h3 <- readRDS("Study2_mod3.rds")

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
h3.1 <- brm(Other_HelpSeeking ~ 1 + chain_position + (1 |chain),
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
#saveRDS(h3.1, "Study2_OR_mod3.rds")
h3.1 <- readRDS("Study2_OR_mod3.rds")

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
h3.2 <- brm(Other_HelpSeeking ~ 1 + mo(chain_position_trend) + (1 + mo(chain_position_trend) || chain),
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
#saveRDS(h3.2, "Study2_monotonic_mod3.RDS")
h3.2 <- readRDS("Study2_monotonic_mod3.RDS") 


# Summary of model and posterior predictive check
summary(h3.2)
pp_check(h3.2)


######## HYPOTHESIS TESTS FOR EFFECTS OF INTEREST
hypothesis(h3.2, "mochain_position_trend > 0", class = "bsp")

## Plot predictions
conditional_effects(h3.2)





################ INVESTIGATING TRANSMISSION FORMAT INTERACTION EFFECTS

#####################################
### DATA ANALYSIS - Model 4 to see whether transmission format moderates the effect of transmission on self-diagnosis

## Specifying the model
h4 <- brm(Self_Diagnosis ~ 1 + chain_position*transmission_type + (1 |chain),
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
#saveRDS(h4, "Study2_mod4.rds")
h4 <- readRDS("Study2_mod4.rds")

## Checking the model output and doing a posterior predictive check
summary(h4)
pp_check(h4, type = "bars")


####### HYPOTHESIS TESTS FOR EFFECTS OF INTEREST
### 

## Test of interaction
hyp_interaction <- c(
  pos2_Int = "chain_position2:transmission_typewords   > 0",
  pos3_Int = "chain_position3:transmission_typewords   > 0",
  pos4_Int = "chain_position4:transmission_typewords   > 0",
  pos5_Int = "chain_position5:transmission_typewords   > 0",
  pos6_Int = "chain_position6:transmission_typewords   > 0",
  pos7_Int = "chain_position7:transmission_typewords   > 0"
)

hypothesis(h4, hyp_interaction)


## Investigating if there is a main effect of transmission format

# Using the marginal effects package
avg_comparisons(h4, variables = "transmission_type", type = "link")

# Using the emmeans package - both yield identical estimates and seems like there is a format main effect
e <- emmeans(h4, ~ transmission_type) |> pairs()
p_direction(e)

#####################################
### DATA ANALYSIS - Model 4 monotonic variant to see whether transmission format moderates the effect of transmission on self-diagnosis

## Specifying the model
h4.1 <- brm(Self_Diagnosis ~ 1 + mo(chain_position_trend)*transmission_type + (1 + mo(chain_position_trend) || chain),
                  data = cleaned_dataset,
                  family = cumulative("probit"),
                  c(prior(normal(0, 1.5), class = "Intercept"),
                    prior(normal(0, 1), class = "b"),
                    prior(exponential(1), class = "sd"),
                    prior(dirichlet(2), class = "simo", coef = "mochain_position_trend1")),
                  warmup = 1000,
                  iter = 3500,
                  cores = 4,
                  seed = 23,
                  control = list(adapt_delta = 0.95)
)


# Save & load the model output. 
#saveRDS(h4.1, "Study2_monotonic_mod4.rds")
h4.1 <- readRDS("Study2_monotonic_mod4.rds")

## Checking the model output and doing a posterior predictive check
summary(h4.1)
pp_check(h4.1, type = "bars")

######## HYPOTHESIS TESTS FOR EFFECTS OF INTEREST

## For interaction
hypothesis(h4.1, "mochain_position_trend:transmission_typewords > 0", class = "bsp")

## For main effect of transmission format
e <- emmeans(h4.1, ~ transmission_type) |> pairs()
p_direction(e)

## NB: We aren't running any odds ratios variants here as latent SDs are very interpretable and this scale makes sensefor comparing a simple between-group difference




#####################################
### DATA ANALYSIS - Model 5 to see whether transmission format moderates the effect of transmission on other-diagnosis

## Specifying the model
h5 <- brm(Other_Diagnosis ~ 1 + chain_position*transmission_type + (1 |chain),
            data = cleaned_dataset,
            family = cumulative("probit"),
            c(prior(normal(0, 1.5), class = "Intercept"),
              prior(normal(0, 1), class = "b"),
              prior(exponential(1), class = "sd")),
            warmup = 1000,
            iter = 3500,
            cores = 4,
            seed = 25,
            control = list(adapt_delta = 0.95)
)

# Save & load the model output
# saveRDS(h5, "Study2_mod5.rds")
h5 <- readRDS("Study2_mod5.rds")

## Checking the model output and doing a posterior predictive check
summary(h5)
pp_check(h5, type = "bars")


####### HYPOTHESIS TESTS FOR EFFECTS OF INTEREST
### 

## Test of interaction
hypothesis(h5, hyp_interaction)

## Main Effect of Transmission Format
e <- emmeans(h5, ~ transmission_type) |> pairs()
p_direction(e)


#####################################
### DATA ANALYSIS - Model 5 Monotonic Variant to see whether transmission format moderates the effect of transmission on other-diagnosis

## Specifying the model
h5.1 <- brm(Other_Diagnosis ~ 1 + mo(chain_position_trend)*transmission_type + (1 + mo(chain_position_trend) || chain),
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

# Save & load the model output. 
#saveRDS(h5.1, "Study2_monotonic_mod5.rds")
h5.1 <- readRDS("Study2_monotonic_mod5.rds")

## Checking the model output and doing a posterior predictive check
summary(h5.1)
pp_check(h5.1, type = "bars")


####### HYPOTHESIS TESTS FOR EFFECTS OF INTEREST
####### 

# Test of interaction 
hypothesis(h5.1, "mochain_position_trend:transmission_typewords < 0", class = "bsp")

# Test of transmission main effect 
e <- emmeans(h5.1, ~ transmission_type) |> pairs()
p_direction(e)




#####################################
### DATA ANALYSIS - Model 6 to see whether transmission format moderates the effect of transmission on other help-seeking

## Specifying the model
h6 <- brm(Other_HelpSeeking ~ 1 + chain_position*transmission_type + (1 |chain),
            data = cleaned_dataset,
            family = cumulative("probit"),
            c(prior(normal(0, 1.5), class = "Intercept"),
              prior(normal(0, 1), class = "b"),
              prior(exponential(1), class = "sd")),
            warmup = 1000,
            iter = 3500,
            cores = 4,
            seed = 23,
            control = list(adapt_delta = 0.95)
)

# Save & load the model output
#saveRDS(h6, "Study2_mod6.rds")
h6 <- readRDS("Study2_mod6.rds")

## Checking the model output and doing a posterior predictive check
summary(h6)
pp_check(h6, type = "bars")


####### HYPOTHESIS TESTS FOR EFFECTS OF INTEREST
####### 

## Test of interaction
hypothesis(h6, hyp_interaction)

# Test of transmission main effect 
e <- emmeans(h6, ~ transmission_type) |> pairs()
p_direction(e)




#####################################
### DATA ANALYSIS - Model 6 to see whether transmission format moderates the effect of transmission on other help-seeking


## Specifying the model
h6.1 <- brm(Other_HelpSeeking ~ 1 + mo(chain_position_trend)*transmission_type + (1 + mo(chain_position_trend) || chain),
                  data = cleaned_dataset,
                  family = cumulative("probit"),
                  c(prior(normal(0, 1.5), class = "Intercept"),
                    prior(normal(0, 1), class = "b"),
                    prior(exponential(1), class = "sd"),
                    prior(dirichlet(2), class = "simo", coef = "mochain_position_trend1")),
                  warmup = 1000,
                  iter = 3500,
                  cores = 4,
                  seed = 56,
                  control = list(adapt_delta = 0.95)
)

# Save & load the model output. 
#saveRDS(h6.1, "Study2_monotonic_mod6.rds")
h6.1 <- readRDS("Study2_monotonic_mod6.rds")

## Checking the model output and doing a posterior predictive check
summary(h6.1)
pp_check(h6.1, type = "bars")


####### HYPOTHESIS TESTS FOR EFFECTS OF INTEREST
###

# Test of interaction 
hypothesis(h6.1, "mochain_position_trend:transmission_typewords > 0", class = "bsp")

# Test of transmission main effect 
e <- emmeans(h6.1, ~ transmission_type) |> pairs()
p_direction(e)



#####################################
### EXPORTING CLEANED DATASET FOR POOLING FOR MEGA-ANALYTIC REGRESSION

## Selecting columns of interest that are vital for merging data across four studies
meta_analytic_dataset <- cleaned_dataset %>%
  select(PID, chain, chain_position, Self_Diagnosis, Other_Diagnosis, Other_HelpSeeking, Prevalence, Treatment, StartDate, post_AI_introduction, stimulus)

## Add experiment number to identify dataset
meta_analytic_dataset$study_num <- "exp2"

## Export dataset
write.csv(meta_analytic_dataset, "meta_analysis_study2data.csv", row.names = FALSE)





