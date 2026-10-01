######################## This script generates a summary table describing the demographics of participants in Study 2
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
setwd(file.path(renv::project(), "Study 2")) 

## Set seed to ensure reproducibility
set.seed(643)

## Loading the custom made table theme to homogenise style across the manuscript
source("Table Theme Function.R")


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
  examples_pasted = d1$examples_pasted)             


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


## Okay, now let's present this same distribution breakup by each transmission format category
summary_df2 <- cleaned_dataset %>%
  mutate(
    transmission_type = replace_values(                     # Relabelling the transmission format condition properly
      transmission_type,
      "examples" ~ "Examples",
      "words"    ~ "Narrative Summaries"
    ) 
    ) %>%                                        
  group_by(transmission_type, chain_position) %>%
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
  group_by(
    used_AI_directly = ifelse(used_AI_directly, "Yes", "No"),
    post_AI_introduction = as.character(post_AI_introduction)
  ) %>%
  summarise(Count = n(), .groups = "drop") %>%
  complete(
    used_AI_directly = c("No", "Yes"),
    post_AI_introduction = c("0", "1"),
    fill = list(Count = 0)
  ) %>%
  mutate(
    Percentage = round((Count / sum(Count)) * 100, 2),
    Classification = case_when(
      used_AI_directly == "No"  & post_AI_introduction == "0" ~ "Pure non-AI observations",
      used_AI_directly == "No"  & post_AI_introduction == "1" ~ "Indirect Exposure to an AI Explanation generated by previous participant in the chain",
      used_AI_directly == "Yes" & post_AI_introduction == "1" ~ "Participants who used AI to generate their Explanation",
      used_AI_directly == "Yes" & post_AI_introduction == "0" ~ "Participants who used AI to generate their explanation, but aren't considered to be exposed to AI"
    )
  ) %>%
  select(Classification, used_AI_directly, post_AI_introduction, Count, Percentage)




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
    caption = "Sample characteristics of Study 2 by Chain Position"
  )

descriptive_ft


## Output the second overall demographics table by transmission format
descriptive_ft2 <- summary_df2 |>
  tibble::add_column(gap = "", .after = "Past_diagnoses") |>
  flextable() |>
  set_header_labels(
    chain_position = "Chain position",
    transmission_type = "Transmission Format",
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
    colwidths = c(5, 4, 1, 3)
  ) |>
  align(align = "center", part = "all") |>
  apa_theme(
    spanner_cols = list(6:9, 11:13),
    stripe_rows  = seq(2, nrow(summary_df2), 2),
    caption = "Sample characteristics of Study 2 by Chain Position and Transmission Format"
  )

descriptive_ft2


## Output the AI usage table
ai_usage_ft <- ai_usage_summary %>%
  flextable() %>%
  set_header_labels(
    Classification = "AI Exposure Classification",
    used_AI_directly = "Direct AI User",
    post_AI_introduction = "AI Usage in Transmission Chain",
    Count = "N",
    Percentage = "% of Total Sample"
  ) %>%
  align(align = "center", part = "all") %>%
  align(j = "Classification", align = "left", part = "all") %>%
  apa_theme(
    caption = "Summary table of direct AI usage vs. indirect Exposure to AI generated ADHD explanations",
    stripe_rows  = seq(2, nrow(ai_usage_summary), 2)
  )

ai_usage_ft