######################## This script computes the semantic drift (via cosine similarity) and does a test of mediation
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
library(text)

## Set working directory to current study folder
setwd(file.path(renv::project(), "Mega-Analysis")) 

## Set seed to ensure reproducibility
set.seed(20)

# Initialize the text package environment (run once per session)
textrpp_initialize()




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

## For chain position 1 in Experiment 4, the DSM stimulus was shown in the qualtrics interface directly.
## Hence, it hasn't been recorded as metadata in the initial dataset
## Adding this from one of the experiments so that cosine similarity can be computed smoothly
cleaned_dataset %>% filter(stimulus == "")
dsm_text <- cleaned_dataset$stimulus[cleaned_dataset$study_num == "exp2" & cleaned_dataset$chain_position == 1 & !is.na(cleaned_dataset$stimulus) & cleaned_dataset$stimulus != ""][1]
cleaned_dataset$stimulus[cleaned_dataset$chain_position == 1 & cleaned_dataset$study_num == "exp4"] <- dsm_text
cleaned_dataset %>% filter(stimulus == "")




#####################################
### DATA PREP FOR COSINE SIMILARITY COMPUTATION

## Writing a function to remove the HTML tags in the dataset (inserted so that they would display correctly on Qualtrics) and make the text clean for analysis
clean_html_to_text <- function(text_column) {
  text_clean <- str_replace_all(text_column, "(</p>|</li>)", ". ")   ## Replacing end of paragraph and each example with a period and space
  text_clean <- str_remove_all(text_clean, "<[^>]+>")                ## Removing any remaining tags
  text_clean <- str_replace_all(text_clean, "\\s+", " ")             ## Removing any multiple spaces
  text_clean <- str_replace_all(text_clean, "\\.\\s*\\.", ".")       ## Removing any multiple periods   
  text_clean <- str_trim(text_clean)                                 ## Removing any trailing whitespaces
  return(text_clean)
}


## Cleaning and preparing the dataset for analysis now
cleaned_dataset <- cleaned_dataset %>%
  mutate(
    CleanedStim = clean_html_to_text(stimulus)                             # Applying the cleaning function to the stimulus 
  )



#####################################
### COSINE SIMILARITY COMPUTATION


## Setting the baseline text against which cosine similarity will be computed for all other stimuli
# In our case, it is the DSM which people saw at the beginning of the chain. 
cleaned_dataset$DSM_Baseline <- cleaned_dataset$CleanedStim[cleaned_dataset$study_num == "exp2" & cleaned_dataset$chain_position == 1 & !is.na(cleaned_dataset$CleanedStim) & cleaned_dataset$CleanedStim != ""][1]

## Using the sentence transformer now to extract the embeddings for all of the texts 
print("Extracting embeddings for participant responses...")
embeddings_participant <- textEmbed(
  texts = cleaned_dataset$CleanedStim,
  model = "BAAI/bge-small-en-v1.5"
)

## Ignore the non-ASCII warnings. These refer to usage of curly double quotes, em dashes or non-breaking spaces whose removal is harmless
# Check the texts where these characters are present. 
screen <- textFindNonASCII(tibble(Text = cleaned_dataset$CleanedStim))

nonascii_chars <- cleaned_dataset %>%
  mutate(row = row_number(),
         chars = str_extract_all(CleanedStim, "[^\\x01-\\x7F]")) %>%   # Pulling out every non-ASCII character in each stimulus
  select(row, study_num, chain_position, chars) %>%
  unnest(chars)

# Saving these embeddings locally
#saveRDS(embeddings_participant, "participant_embeddings.rds")
embeddings_participant <- readRDS("participant_embeddings.rds")

## The DSM baseline embeddings have already been extracted 
## So instead of re-running it on the DSM baseline again (which will take a lot of time), we can just extract it
emb_participant <- embeddings_participant$texts[[1]]

dsm_idx <- which(cleaned_dataset$CleanedStim == cleaned_dataset$DSM_Baseline)[1]   # Row index of the first DSM stimulus
emb_dsm <- emb_participant[rep(dsm_idx, nrow(emb_participant)), ]                  # Replicating that single row to match every participant row

## Computing cosine similarity between every stimulus and the DSM baseline
cleaned_dataset$Cosine_Similarity <- textSimilarity(
  x = emb_participant,
  y = emb_dsm,
  method = "cosine"
)

## NB: Note that for experiment 1, chain position 1, cosine similarity rounds to 1 - this is because there was a small typo in the specifiers introductory line
# Should be fine for our analyses



#####################################
### FINAL DATA PREP FOR MEDIATION 

## Standardising the variables to allow for the gaussian models to run smoothly
## Also allows for a smooth comparision of estimates with the probit models
cleaned_dataset$Self_Diagnosis_z <- as.numeric(scale(cleaned_dataset$Self_Diagnosis))
cleaned_dataset$Other_Diagnosis_z <- as.numeric(scale(cleaned_dataset$Other_Diagnosis))
cleaned_dataset$Other_HelpSeeking_z <- as.numeric(scale(cleaned_dataset$Other_HelpSeeking))
cleaned_dataset$Cosine_Similarity_z <- as.numeric(scale(cleaned_dataset$Cosine_Similarity))





#####################################
### DATA ANALYSIS: Model 1 Simple Gaussian Model to see if treating the dependent variables as continuous yields similar estimates to ordinal regressions

## Specifying the model
h1 <- brm(Self_Diagnosis_z ~ 1 + mo(chain_position_trend) + Day + study_num + (1 + mo(chain_position_trend) || chain_ID),
          data = cleaned_dataset,
          family = gaussian(),
          c(prior(normal(0, 1.5), class = "Intercept"),
            prior(normal(0, 1), class = "b"),
            prior(exponential(1), class = "sd"),
            prior(dirichlet(2), class = "simo", coef = "mochain_position_trend1")),
          warmup = 1000,
          iter = 3500,
          cores = 4,
          seed = 20,
          control = list(adapt_delta = 0.95)
)

# Save & load the model output. 
#saveRDS(h1, "MA_gaussian_mod1.rds")
h1 <- readRDS("MA_gaussian_mod1.rds")

## Checking the model output
summary(h1)

## Compare main estimand of interest ("chain position trend") with ordinal results
#  they look similar, so can go ahead and do proper mediation analyses
h1_ordinal <- read_rds("MA_monotonic_mod1.rds")
summary(h1_ordinal)



#####################################
### DATA ANALYSIS: Model 1 Does Semantic Drift mediate transmission effect on self-diagnosis

## Specify mediation model 
m1 <- bf (Cosine_Similarity_z ~ 1 + mo(chain_position_trend) + Day + study_num + (1 + mo(chain_position_trend) || chain_ID))
m2 <- bf (Self_Diagnosis_z ~ 1 + mo(chain_position_trend) + Cosine_Similarity_z + Day + study_num + (1 + mo(chain_position_trend) || chain_ID))

## Fit as a multivariate model
h1.1 <- brm(
  m1 + m2 + set_rescor(FALSE), 
  data = cleaned_dataset, 
  family = gaussian(),
  prior = c(
    # Priors for Cosine_Similarity (m1)
    prior(normal(0, 1), class = "Intercept", resp = "CosineSimilarityz"),
    prior(normal(0, 1), class = "b", resp = "CosineSimilarityz"),
    prior(exponential(1), class = "sd", resp = "CosineSimilarityz"),
    prior(dirichlet(2), class = "simo", coef = "mochain_position_trend1", resp = "CosineSimilarityz"),
    
    # Priors for Self_Diagnosis (m2)
    prior(normal(0, 1), class = "Intercept", resp = "SelfDiagnosisz"),
    prior(normal(0, 1), class = "b", resp = "SelfDiagnosisz"),
    prior(exponential(1), class = "sd", resp = "SelfDiagnosisz"),
    prior(dirichlet(2), class = "simo", coef = "mochain_position_trend1", resp = "SelfDiagnosisz")
  ),
  warmup = 1000,
  iter = 3500,
  cores = 4,
  seed = 8,
  control = list(adapt_delta = 0.95)
)

# Save & load the model output. 
#saveRDS(h1.1, "MA_mediation_mod1.rds")
h1.1 <- readRDS("MA_mediation_mod1.rds")

### Checking evidence for mediation

# Draw from the posteriors
h1_draws <- as_draws_df(h1.1)

# Extract Path coefficients
a        <- h1_draws$bsp_CosineSimilarityz_mochain_position_trend   # chain_position -> Cosine Similarity
b        <- h1_draws$b_SelfDiagnosisz_Cosine_Similarity_z  # Cosine Similarity -> Self Diagnosis
direct   <- h1_draws$bsp_SelfDiagnosisz_mochain_position_trend    # chain_position -> Self Diagnosis (direct)

# Calculate pertinent mediation effects
indirect <- a * b
total    <- direct + indirect
prop_med <- indirect / total

# Summary function to convey the posteriors of these effects
summarise_effect <- function(samples, label) {
  ci  <- quantile(samples, c(0.025, 0.975))
  m <- mean(samples)
  pp  <- mean(samples > 0)
  cat(sprintf("\n%s:\n  Mean = %.4f, 95%% CI [%.4f, %.4f], P(effect > 0) = %.3f\n",
              label, m, ci[1], ci[2], pp))
}

# Effect Summary
summarise_effect(indirect, "Indirect effect (a*b)")
summarise_effect(direct,   "Direct effect")
summarise_effect(total,    "Total effect")
summarise_effect(prop_med, "Proportion mediated")

  


#####################################
### DATA ANALYSIS: Model 2 Simple Gaussian Model to see if treating the dependent variables as continuous yields similar estimates to ordinal regressions

## Specifying the model
h2 <- brm(Other_Diagnosis_z ~ 1 + mo(chain_position_trend) + Day + study_num + (1 + mo(chain_position_trend) || chain_ID),
          data = cleaned_dataset,
          family = gaussian(),
          c(prior(normal(0, 1.5), class = "Intercept"),
            prior(normal(0, 1), class = "b"),
            prior(exponential(1), class = "sd"),
            prior(dirichlet(2), class = "simo", coef = "mochain_position_trend1")),
          warmup = 1000,
          iter = 3500,
          seed = 2,
          cores = 4,
          control = list(adapt_delta = 0.95)
)

# Save & load the model output. 
#saveRDS(h2, "MA_gaussian_mod2.rds")
h2 <- readRDS("MA_gaussian_mod2.rds")

## Checking the model output
summary(h2)

## Compare main estimand of interest ("chain position trend") with ordinal results
#  they look similar, so can go ahead and do proper mediation analyses
h2_ordinal <- read_rds("MA_monotonic_mod2.rds")
summary(h2_ordinal)




#####################################
### DATA ANALYSIS: Model 2 Does Semantic Drift mediate transmission effect on other-diagnosis

## Specify mediation model 
m1 <- bf (Cosine_Similarity_z ~ 1 + mo(chain_position_trend) + Day + study_num + (1 + mo(chain_position_trend) || chain_ID))
m2 <- bf (Other_Diagnosis_z ~ 1 + mo(chain_position_trend) + Cosine_Similarity_z + Day + study_num + (1 + mo(chain_position_trend) || chain_ID))

## Fit as a multivariate model
h2.1 <- brm(
  m1 + m2 + set_rescor(FALSE), 
  data = cleaned_dataset, 
  family = gaussian(),
  prior = c(
    # Priors for Cosine_Similarity (m1)
    prior(normal(0, 1), class = "Intercept", resp = "CosineSimilarityz"),
    prior(normal(0, 1), class = "b", resp = "CosineSimilarityz"),
    prior(exponential(1), class = "sd", resp = "CosineSimilarityz"),
    prior(dirichlet(2), class = "simo", coef = "mochain_position_trend1", resp = "CosineSimilarityz"),
    
    # Priors for Other_Diagnosis (m2)
    prior(normal(0, 1), class = "Intercept", resp = "OtherDiagnosisz"),
    prior(normal(0, 1), class = "b", resp = "OtherDiagnosisz"),
    prior(exponential(1), class = "sd", resp = "OtherDiagnosisz"),
    prior(dirichlet(2), class = "simo", coef = "mochain_position_trend1", resp = "OtherDiagnosisz")
  ),
  warmup = 1000,
  iter = 3500,
  cores = 4,
  seed = 5,
  control = list(adapt_delta = 0.95)
)

# Save & load the model output. 
# saveRDS(h2.1, "MA_mediation_mod2.rds")
h2.1 <- readRDS("MA_mediation_mod2.rds")


### Checking evidence for mediation

# Draw from the posteriors
h2_draws <- as_draws_df(h2.1)

# Extract Path coefficients
a        <- h2_draws$bsp_CosineSimilarityz_mochain_position_trend   # chain_position -> Cosine Similarity
b        <- h2_draws$b_OtherDiagnosisz_Cosine_Similarity_z  # Cosine Similarity -> Other Diagnosis
direct   <- h2_draws$bsp_OtherDiagnosisz_mochain_position_trend    # chain_position -> Other Diagnosis (direct)

# Calculate pertinent mediation effects
indirect <- a * b
total    <- direct + indirect
prop_med <- indirect / total


# Effect Summary
summarise_effect(indirect, "Indirect effect (a*b)")
summarise_effect(direct,   "Direct effect")
summarise_effect(total,    "Total effect")
summarise_effect(prop_med, "Proportion mediated")
  


#####################################
### DATA ANALYSIS: Model 3 Simple Gaussian Model to see if treating the dependent variables as continuous yields similar estimates to ordinal regressions

## Specifying the model
h3 <- brm(Other_HelpSeeking_z ~ 1 + mo(chain_position_trend) + Day + study_num + (1 + mo(chain_position_trend) || chain_ID),
          data = cleaned_dataset,
          family = gaussian(),
          c(prior(normal(0, 1.5), class = "Intercept"),
            prior(normal(0, 1), class = "b"),
            prior(exponential(1), class = "sd"),
            prior(dirichlet(2), class = "simo", coef = "mochain_position_trend1")),
          warmup = 1000,
          iter = 3500,
          cores = 4,
          seed = 3,
          control = list(adapt_delta = 0.95)
)

# Save & load the model output. 
# saveRDS(h3, "MA_gaussian_mod3.rds")
h3 <- readRDS("MA_gaussian_mod3.rds")

## Checking the model output
summary(h3)

## Compare main estimand of interest ("chain position trend") with ordinal results
#  they look similar, so can go ahead and do proper mediation analyses
#h3_ordinal <- read_rds("MA_monotonic_mod3.rds")
summary(h3_ordinal)




#####################################
### DATA ANALYSIS: Model 3 Does Semantic Drift mediate transmission effect on other help-seeking

## Specify mediation model 
m1 <- bf (Cosine_Similarity_z ~ 1 + mo(chain_position_trend) + Day + study_num + (1 + mo(chain_position_trend) || chain_ID))
m2 <- bf (Other_HelpSeeking_z ~ 1 + mo(chain_position_trend) + Cosine_Similarity_z + Day + study_num + (1 + mo(chain_position_trend) || chain_ID))

## Fit as a multivariate model
h3.1 <- brm(
  m1 + m2 + set_rescor(FALSE), 
  data = cleaned_dataset, 
  family = gaussian(),
  prior = c(
    # Priors for Cosine_Similarity (m1)
    prior(normal(0, 1), class = "Intercept", resp = "CosineSimilarityz"),
    prior(normal(0, 1), class = "b", resp = "CosineSimilarityz"),
    prior(exponential(1), class = "sd", resp = "CosineSimilarityz"),
    prior(dirichlet(2), class = "simo", coef = "mochain_position_trend1", resp = "CosineSimilarityz"),
    
    # Priors for Other_HelpSeeking (m2)
    prior(normal(0, 1), class = "Intercept", resp = "OtherHelpSeekingz"),
    prior(normal(0, 1), class = "b", resp = "OtherHelpSeekingz"),
    prior(exponential(1), class = "sd", resp = "OtherHelpSeekingz"),
    prior(dirichlet(2), class = "simo", coef = "mochain_position_trend1", resp = "OtherHelpSeekingz")
  ),
  warmup = 1000,
  iter = 3500,
  cores = 4,
  seed = 20,
  control = list(adapt_delta = 0.95)
)

# Save & load the model output. 
#saveRDS(h3.1, "MA_mediation_mod3.rds")
h3.1 <- readRDS("MA_mediation_mod3.rds")


### Checking evidence for mediation

# Draw from the posteriors
h3_draws <- as_draws_df(h3.1)

# Extract Path coefficients
a        <- h3_draws$bsp_CosineSimilarityz_mochain_position_trend   # chain_position -> Cosine Similarity
b        <- h3_draws$b_OtherHelpSeekingz_Cosine_Similarity_z  # Cosine Similarity -> Other HelpSeeking
direct   <- h3_draws$bsp_OtherHelpSeekingz_mochain_position_trend    # chain_position -> Other HelpSeeking (direct)

# Calculate pertinent mediation effects
indirect <- a * b
total    <- direct + indirect
prop_med <- indirect / total


# Effect Summary
summarise_effect(indirect, "Indirect effect (a*b)")
summarise_effect(direct,   "Direct effect")
summarise_effect(total,    "Total effect")
summarise_effect(prop_med, "Proportion mediated")


