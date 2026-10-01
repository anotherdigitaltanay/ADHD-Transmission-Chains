######################## This script generates a summary table describing the demographics of participants in Study 3
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

## Set working directory to current study folder
setwd(file.path(renv::project(), "Study 4")) 

## Set seed to ensure reproducibility
set.seed(43)

## Loading the custom made table theme to homogenise style across the manuscript
source("Table Theme Function.R")


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
### DATA DEMOGRAPHICS SUMMARY

## Now let's try to summarise these distribution of these variables across each chain position/generation
summary_df <- cleaned_dataset %>%
  group_by(chain_position) %>%
  summarise(
    Observations = n(),
    Age = round(mean(age), 2),
    Gender = round(mean((gender == "Female") * 100), 2),
    Medication = round(mean((medication== "Yes")* 100), 2),
    Treatment = round(mean(treatment == "Yes") * 100, 2),
    NHS_Access = round(mean(Health_access == "Yes") * 100, 2),
    Past_diagnoses = round(mean(Past_diagnoses == "Yes") * 100, 2),
    ADHD_Friend = round(mean(Priork_friend_help == "Yes") * 100, 2),
    ADHD_Friend_Waitlist = round(mean(Priork_friend_helpwait == "Yes") * 100, 2),
    MH_SocialMedia = round(mean(Priork_socialmedia == "Yes") * 100, 2)
  )


## Now let's summarise the AI usage in the sample neatly
ai_usage_summary <- cleaned_dataset %>%
  group_by(post_AI_introduction) %>%
  summarise(Count = n(), .groups = "drop") %>%
  complete(post_AI_introduction = c(FALSE, TRUE), fill = list(Count = 0)) %>%
  mutate(
    Percentage = round(Count / sum(Count) * 100, 2),
    Classification = if_else(
      post_AI_introduction,
      "Participants who used AI to generate their explanation",
      "Pure non-AI observations"
    ),
    post_AI_introduction= if_else(post_AI_introduction, "Yes", "No")
  ) %>%
  select(Classification, post_AI_introduction, Count, Percentage)












## Output the first overall demographics table
descriptive_ft <- summary_df |>
  tibble::add_column(gap = "", .after = "Past_diagnoses") |>
  flextable() |>
  set_header_labels(
    chain_position = "Chain position",
    Observations = "Observations",
    Age = "Age (Mean)",
    Gender = "Female (%)",
    Medication = "MH Medication",
    Treatment = "MH Treatment",
    NHS_Access = "NHS MH Support",
    Past_diagnoses = "Past MH Diagnosis",
    gap = "",
    ADHD_Friend = "Friend Support",
    ADHD_Friend_Waitlist = "Waitlist Support",
    MH_SocialMedia = "Social Media"
  ) |>
  add_header_row(
    values = c("", "Clinical History (%)", "", "ADHD Prior Knowledge (%)"),
    colwidths = c(4, 4, 1, 3)
  ) |>
  align(align = "center", part = "all") |>
  apa_theme(
    spanner_cols = list(5:8, 10:12),
    stripe_rows  = seq(2, nrow(summary_df), 2),
    caption = "Sample characteristics of Study 4 by Chain Position"
  )

descriptive_ft


## Output the AI usage table
ai_usage_ft <- ai_usage_summary %>%
  flextable() %>%
  set_header_labels(
    Classification = "AI Exposure Classification",
    post_AI_introduction = "AI Usage",
    Count = "N",
    Percentage = "% of Total Sample"
  ) %>%
  align(align = "center", part = "all") %>%
  align(j = "Classification", align = "left", part = "all") %>%
  apa_theme(
    caption = "Summary table of direct AI usage",
    stripe_rows  = seq(2, nrow(ai_usage_summary), 2)
  )

ai_usage_ft