############################################## The following script contains the power analysis for the chains experiment #######################

#### Load the pertinent packages first
library(tidyverse)
library(brms)
library(latent2likert)

### Set a seed for this script
set.seed (250)


###### Now let's give our first pass at data simulation and parameter recovery
###############################################################################

## We will be using the rlikert function from the latent2likert package that allows us to specify the latent effect size for any latent ordinal variable

## Quick explainer on the package's function arguments - 
# n_levels indicates number of likert scale options. In this case, the dependent variable will have 5 options (Strongly Disagree to Strongly Agree)
# n_items indicates number of items. For now, we have one primary dependent variable.
# skew indicates if a group is more likely to respond with options at the tail end of the likert scale question

## Now, while the final experiment would have multiple chains of 7 to 10 people (i.e. 7 to 10 generations), let's build complexity slowly by considering the simplest possible case
## That is, let's compare between generation 1 (i.e. participants in position 1) v/s generation 2 (i.e. participants in position 2). Multilevel structure is being ignored atm, but will be introduced later
# In other words, we are comparing two groups

## Let's assume that data is skewed to the left (as from our pilot sample, we know most folks who don't have any current mental health condition most often select Option 1 = Strongly disagree with "I believe I have ADHD")
## Additionally, let's aim for powering the analysis to detect a medium effect size of 0.3

# Let's start with sample size of 50 initially to see if the simulation works
# Okay, simulate 
gen1 <- replicate(50, rlikert(size = 1, n_levels = 5, n_items = 1, mean = 0, skew = 0.7))
gen2 <- replicate(50, rlikert(size = 1, n_levels = 5, n_items = 1, mean = 0.3, skew = 0.7))

## Quickly inspect the shape of the simulated data
hist(gen1)
hist(gen2)

## Let's combine this into one single dataframe
d1 <- tibble(
  rating = c(gen1, gen2), 
  gen = rep(c("Gen1", "Gen2"), each = 50)
) 

## Let's run the model now to see whether the effect is recovered
m1 <- brm(rating ~ 1 + gen,
          data = d1,
          family = cumulative("probit"),
          prior = c(prior(normal(0, 1.5), class = "Intercept"),
                    prior(normal(0, 1), class = "b")),
          warmup = 1000,
          iter = 3000)

## Save the model output, so that others don't have to re-run the model and waste additional computational resources
## Run only line 56 if you are reading through the script and want the model output directly
saveRDS(m1, "model1.RDS")
m1 <- readRDS("model1.RDS")

## Inspecting model summary
## Model is doing a poor job at recovering the parameter (look at genGen2 coeff which is the difference b/w gen 1 and gen 2 i.e. our effect of interest)
summary(m1)
fixef(m1)


## Re priors, we have chosen what folks called a weakly regularising prior that slightly biases the model towards a null effect
## To ensure that this prior is indeed doing so, let's run a prior predictive simulation and see what predictions the model makes

# In short, before seeing any data, the model should be biased towards making predictions that show no difference between the two groups

# Running the prior predictive simulation
pp_sim_model <- update(                                      ## Using this so that we don't have to recompile m1 again and use additional computational resources
  m1,                                                        ## Specifying which model's compilation to use - this is the regression we specified in line 45
  sample_prior = "only",                                     ## Don't use any data, just sample from the prior
  seed = 23,
  control = c(adapt_delta = 0.9)
)

## Saving model output
saveRDS(pp_sim_model, "prior_predictive_sim.RDS")
pp_sim_model <- readRDS("prior_predictive_sim.RDS")

## Inspecting model summary and predictions
## As you can see, no differences between the two groups (Ignore the continous variable warning for the time-being)
summary(pp_sim_model)
prediction_pp_sim <- conditional_effects(pp_sim_model)
prediction_pp_sim 
prediction_pp_sim_ordinal <- conditional_effects(pp_sim_model, categorical = TRUE)
prediction_pp_sim_ordinal 


###### Now let's reduce code redundancies by writing up the simulation logic in a simulation function
#####################################################################################################

## First, we make a function to simulate the data
data_sim <- function(sim_seed,                     ## Simulation seed
                        n_obs)                     ## Number of observations in each group/generation
                     {                
  
  mean_effect <- 0.3                               ## Effect on the latent ordinal scale              
  skew_ord <- 0.7                                  ## Skew on the latent ordinal scale
  
  ## Set seed for the simulation
  set.seed(sim_seed)
  
  ## Simulate data for both groups
  gen1 <- replicate(n_obs, rlikert(size = 1, n_levels = 5, n_items = 1, mean = 0, skew = skew_ord))
  gen2 <- replicate(n_obs, rlikert(size = 1, n_levels = 5, n_items = 1, mean = mean_effect, skew = skew_ord))
  
  # Combine into the dataframe for modeling
  tibble(
    rating = c(gen1, gen2),
    gen = rep(c("Gen1", "Gen2"), each = n_obs)
  )
}


###### Now, for our first toy power analysis for the simple model
#####################################################################################################

########### General Strategy ###########
########################################
########################################

### We will start by simulating 100 datasets for three sample sizes (n = 120, 160, 220)
### Once we have find a sample size that has acceptable power and precision, we will then do a large scale simulation where we add complexity to the model by adding more generations

### I choose to start with 120 as prior chain experiments have used 120. So that will be our starting point
### For each sample size, we will calculate: (a) Power (how many times out of 100 simulations does our model detect an effect); and (b) Precision (how many times out of 100 simulations do all the confidence intervals have a certain width)

######### Sample Size = 120 ###############
############################################

## Let's define how many simulations we want first
n_sim <- 100

## Defining this to see how long it takes to simulate and analyse
t1 <- Sys.time()

## Simulate and analyse
sim1 <- 
  tibble(seed = 1:n_sim) %>%                                           ## column to identify each simulation via its seed
  mutate(d = map(seed, data_sim, n_obs = 120)) %>%                     ## nested column where all the simulated datasets will be stored. In map(), we first indicate that we will iterate over values of seed. These values are fed to the data simulation function, and the final argument defines the sample size for each simulation
  mutate(b1 = map2(d, seed, ~update(m1, newdata = .x, seed = .y) %>%   ## This passes first two columns to run the model that we have defined above
                     fixef() %>%                                               ## Extracts only the fixed effects of the model and puts in a dataframe
                     data.frame %>%
                     rownames_to_column("parameter") %>%
                     filter(parameter == "genGen2" )))                         ## Only saves the fixed effect of interest

## How long did it take?
## Took 33 mins to run a 100 dataset simulation
t2 <- Sys.time()
t2 - t1

## Save the simulations - If you are reading through, just load line 151
saveRDS(sim1, "simulation1_size120.rds")
sim1 <- readRDS("simulation1_size120.rds")

## Now let's plot the results of the simulation
plot1 <- sim1 %>%
  unnest(b1) %>%                                                                                    ### Since the b1 column in sim1 contains a dataframe in each row, we need to undo this nesting
  ggplot(aes(x = reorder(seed, Q2.5), y = Estimate, ymin = Q2.5, ymax = Q97.5)) +                   ### Plot our parameters and the 95% confidence intervals; we also order the parameter estimates by the lower levels of their credible intervals
  geom_hline(yintercept = c(0, .3), color = "gray50", linetype =2) +                                ### Drawing reference lines at the null effect = 0 and the true effect = 0.3, so that the reader instantly gets how many times we recover the true effect
  geom_pointrange(fatten = 1/2) +
  labs(x = "seed (i.e., simulation index)",
       y = expression(beta[1]))
plot1
 
## Let's calculate power
## Power = 62%
power_sim1 <- sim1 %>%
  select(b1) %>%                                                  ## Selecting the fixed effects
  unnest(b1) %>%                                                  ## Unnesting it
  mutate(check = ifelse(Q2.5 > 0, 1, 0)) %>%                      ## Checking in how many simulations, the confidence interval was above the null effect i.e. 0
  summarise(power = mean(check))                                  ## Calculating power 
  
## Let's calculate precision  
## Precision = where "ALL" credible intervals are of a certain width. All of them are below 0.56
precision_sim1 <- sim1 %>%
  select(b1) %>%                                                  ## Selecting the fixed effects
  unnest(b1) %>%                                                  ## Unnesting it
  mutate(check = ifelse(Q97.5 - Q2.5  < 0.56, 1, 0)) %>%           ## Checking in how many simulations, the width of the confidence interval was less than 0.56
  summarise(precision = mean(check))                              ## Calculating precision 

  
######### Sample Size = 160 ###############
############################################

## Let's define how many simulations we want first
n_sim <- 100

## Defining this to see how long it takes to simulate and analyse
t3 <- Sys.time()

## Simulate and analyse
sim2 <- 
  tibble(seed = 1:n_sim) %>%                                           ## column to identify each simulation via its seed
  mutate(d = map(seed, data_sim, n_obs = 160)) %>%                     ## nested column where all the simulated datasets will be stored. In map(), we first indicate that we will iterate over values of seed. These values are fed to the data simulation function, and the final argument defines the sample size for each simulation
  mutate(b1 = map2(d, seed, ~update(m1, newdata = .x, seed = .y) %>%   ## This passes first two columns to run the model that we have defined above
                     fixef() %>%                                               ## Extracts only the fixed effects of the model and puts in a dataframe
                     data.frame %>%
                     rownames_to_column("parameter") %>%
                     filter(parameter == "genGen2" )))                         ## Only saves the fixed effect of interest

## How long did it take?
## Took 45 mins
t4 <- Sys.time()
t4 - t3

## Save the simulations 
saveRDS(sim2, "simulation2_size160.rds")
sim2 <- readRDS("simulation2_size160.rds")

## Now let's plot the results of the simulation
plot2 <- sim2 %>%
  unnest(b1) %>%                                                                                    ### Since the b1 column in sim1 contains a dataframe in each row, we need to undo this nesting
  ggplot(aes(x = reorder(seed, Q2.5), y = Estimate, ymin = Q2.5, ymax = Q97.5)) +                   ### Plot our parameters and the 95% confidence intervals; we also order the parameter estimates by the lower levels of their credible intervals
  geom_hline(yintercept = c(0, .3), color = "gray50", linetype =2) +                                ### Drawing reference lines at the null effect = 0 and the true effect = 0.3, so that the reader instantly gets how many times we recover the true effect
  geom_pointrange(fatten = 1/2) +
  labs(x = "seed (i.e., simulation index)",
       y = expression(beta[1]))
plot2


## Let's calculate power
## Power = 71%
power_sim2 <- sim2 %>%
  select(b1) %>%                                                  ## Selecting the fixed effects
  unnest(b1) %>%                                                  ## Unnesting it
  mutate(check = ifelse(Q2.5 > 0, 1, 0)) %>%                      ## Checking in how many simulations, the confidence interval was above the null effect i.e. 0
  summarise(power = mean(check))                                  ## Calculating power 

## Let's calculate precision  
## Precision = 1 (i.e. all of the CIs have a confidence interval that are of width 0.48)
precision_sim2 <- sim2 %>%
  select(b1) %>%                                                  ## Selecting the fixed effects
  unnest(b1) %>%                                                  ## Unnesting it
  mutate(check = ifelse(Q97.5 - Q2.5  < 0.48, 1, 0)) %>%           ## Checking in how many simulations, the width of the confidence interval was less than 0.5
  summarise(precision = mean(check))                              ## Calculating precision 


######### Sample Size = 220 ###############
############################################

## Let's define how many simulations we want first
n_sim <- 100

## Defining this to see how long it takes to simulate and analyse
t5 <- Sys.time()

## Simulate and analyse
sim3 <- 
  tibble(seed = 1:n_sim) %>%                                           ## column to identify each simulation via its seed
  mutate(d = map(seed, data_sim, n_obs = 220)) %>%                     ## nested column where all the simulated datasets will be stored. In map(), we first indicate that we will iterate over values of seed. These values are fed to the data simulation function, and the final argument defines the sample size for each simulation
  mutate(b1 = map2(d, seed, ~update(m1, newdata = .x, seed = .y) %>%   ## This passes first two columns to run the model that we have defined above
                     fixef() %>%                                               ## Extracts only the fixed effects of the model and puts in a dataframe
                     data.frame %>%
                     rownames_to_column("parameter") %>%
                     filter(parameter == "genGen2" )))                         ## Only saves the fixed effect of interest

## How long did it take?
## Took 56 mins
t6 <- Sys.time()
t6 - t5

## Save the simulations 
saveRDS(sim3, "simulation3_size220.rds")
sim3 <- readRDS("simulation3_size220.rds")


## Now let's plot the results of the simulation
plot3 <- sim3 %>%
  unnest(b1) %>%                                                                                    ### Since the b1 column in sim1 contains a dataframe in each row, we need to undo this nesting
  ggplot(aes(x = reorder(seed, Q2.5), y = Estimate, ymin = Q2.5, ymax = Q97.5)) +                   ### Plot our parameters and the 95% confidence intervals; we also order the parameter estimates by the lower levels of their credible intervals
  geom_hline(yintercept = c(0, .3), color = "gray50", linetype =2) +                                ### Drawing reference lines at the null effect = 0 and the true effect = 0.3, so that the reader instantly gets how many times we recover the true effect
  geom_pointrange(fatten = 1/2) +
  labs(x = "seed (i.e., simulation index)",
       y = expression(beta[1]))
plot3


## Let's calculate power
## Power = 82%
power_sim3 <- sim3 %>%
  select(b1) %>%                                                  ## Selecting the fixed effects
  unnest(b1) %>%                                                  ## Unnesting it
  mutate(check = ifelse(Q2.5 > 0, 1, 0)) %>%                      ## Checking in how many simulations, the confidence interval was above the null effect i.e. 0
  summarise(power = mean(check))                                  ## Calculating power 

## Let's calculate precision  
## Precision = 1 (i.e. all of the CIs have a confidence interval that are of width 0.41)
precision_sim3 <- sim3 %>%
  select(b1) %>%                                                  ## Selecting the fixed effects
  unnest(b1) %>%                                                  ## Unnesting it
  mutate(check = ifelse(Q97.5 - Q2.5  < 0.41, 1, 0)) %>%           ## Checking in how many simulations, the width of the confidence interval was less than 0.5
  summarise(precision = mean(check))                              ## Calculating precision 

########################################################
### Looks like 220 has both good power and precision
#########################################################

### Plotting our toy simulations power analysis for the three sample sizes 

## First, I combine information from all our simulations into one dataframe so that plotting is easy
toy_sims <- bind_rows(power_sim1, power_sim2, power_sim3) %>%
  mutate(sample_size = c(120, 160, 220),
         precision = c(0.56, 0.48, 0.41))

## Second, I reshape this to a long format
toy_sims_long <- toy_sims %>%
  pivot_longer(
    cols = c(power, precision), 
    names_to = "metric", 
    values_to = "value"
  )

## Plot (Revisit this part of code)
ggplot(toy_sims_long, aes(x = sample_size, y = value)) +
  geom_line(aes(color = metric), size = 1) +  # Add lines
  geom_point(size = 3) +                      # Add points for clarity
  facet_wrap(~ metric, scales = "free_y") +   # Facet by metric, allowing different y-axis scales
  labs(
    x = "Sample Size",
    y = NULL, # Remove generic y-label since facets explain themselves
    title = "Power and Precision vs. Sample Size"
  ) +
  theme_minimal()



###### Now, running a more realistic and serious power analysis for a slightly complicated model
#####################################################################################################
#####################################################################################################

###### In the chains experiment, we know that each chain will atleast be 7 generations long (i.e. 7 groups)
###### Further, it is also plausible that effects go in both direction (+ 0.3 or -0.3)

## So now we will write up a simulation function that incorporates this complexity
## Subsequently, we simulate 1000 datasets of n = 220 and see whether we recover the effects properly


### Let's adapt the simulation function first to accommodate this additional complexities
data_sim_2 <- function(sim_seed,                     ## Simulation seed
                     n_obs)                     ## Number of observations in each group/generation
{                
  
  mean_effect <- 0.3                                   ## Effect on the latent ordinal scale
  mean_effect_2 <- - 0.3                               ## Negative Effect on the latent ordinal scale. Just inserted this to see if the simulation is sensitive to effect direction 
  skew_ord <- 0.7                                      ## Skew on the latent ordinal scale
  
  ## Set seed for the simulation
  set.seed(sim_seed)
  
  ## Simulate data for 7 generations/groups now
  ## Notice the alternative effects - one generation has a positive effect while the next one has a negative effect (all of these are contrasts with generation 1)
  gen1 <- rlikert(size = n_obs, n_levels = 5, n_items = 1, mean = 0, skew = skew_ord)
  gen2 <- rlikert(size = n_obs, n_levels = 5, n_items = 1, mean = mean_effect, skew = skew_ord)
  gen3 <- rlikert(size = n_obs, n_levels = 5, n_items = 1, mean = mean_effect_2, skew = skew_ord)
  gen4 <- rlikert(size = n_obs, n_levels = 5, n_items = 1, mean = mean_effect, skew = skew_ord)
  gen5 <- rlikert(size = n_obs, n_levels = 5, n_items = 1, mean = mean_effect_2, skew = skew_ord)
  gen6 <- rlikert(size = n_obs, n_levels = 5, n_items = 1, mean = mean_effect, skew = skew_ord)
  gen7 <- rlikert(size = n_obs, n_levels = 5, n_items = 1, mean = mean_effect_2, skew = skew_ord)
  
  # Combine into the dataframe for modeling
  tibble(
    rating = c(gen1, gen2, gen3, gen4, gen5, gen6, gen7),
    gen = rep(c("Gen1", "Gen2", "Gen3", "Gen4", "Gen5", "Gen6", "Gen7"), each = n_obs)
  )
}

## Now let's simulate and recover some parameters

######### Sample Size = 220 ###############
############################################

## Let's define how many simulations we want first
n_sim <- 1000

## Defining this to see how long it takes to simulate and analyse
t7 <- Sys.time()

## Simulate and analyse
sim4 <- 
  tibble(seed = 1:n_sim) %>%                                                      ## column to identify each simulation via its seed
  mutate(d = map(seed, data_sim_2, n_obs = 220)) %>%                              ## nested column where all the simulated datasets will be stored. In map(), we first indicate that we will iterate over values of seed. These values are fed to the data simulation function, and the final argument defines the sample size for each simulation
  mutate(b1 = map2(d, seed, ~update(m1, newdata = .x, seed = .y) %>%              ## This passes first two columns to run the model that we have defined above
                     fixef() %>%                                                  ## Extracts only the fixed effects of the model and puts in a dataframe
                     data.frame %>%
                     rownames_to_column("parameter") %>%
                     filter(str_detect(parameter, "^genGen"))))                    ## Now we are saving 6 fixed effects instead of just one

## How long did it take?
## Took roughly a full day on my personal laptop
t8 <- Sys.time()
t8 - t7


## Save the simulations 
saveRDS(sim4, "simulation4_size220_model2.rds")
sim4 <- readRDS("simulation4_size220_model2.rds")




####################### Now let's plot and find out power and precision

## As we are now recovering 6 parameters (instead of just 1 contrast between 2 groups), we need to make some quick changes to make the plots easier to read
## First, for each parameter, we need to highlight whether the effect was positive or negative
results_sim4 <- sim4 %>%
  unnest(b1) %>%
  mutate(true_effect = case_when(
    parameter %in% c("genGen3", "genGen5", "genGen7") ~ -0.3,                      # Assigns the true effect to these generations as -0.3
    parameter %in% c("genGen2", "genGen4", "genGen6") ~ 0.3,                       # Assigns the true effect to these generations as +0.3
    TRUE ~ 0
  ))

## Now let's plot the simulation results for each parameter
plot4 <- results_sim4 %>%
  ggplot(aes(x = reorder(seed, Q2.5), y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_hline(yintercept = 0, color = "gray50", linetype = 2) +
  geom_hline(aes(yintercept = true_effect), color = "gray50", linetype = 2) +   # Now the code from lines 367 -371 will make sense
  geom_pointrange(fatten = 1/2, alpha = 0.5) +
  facet_wrap(~parameter, ncol = 3) +                                            # Facet by parameter to see all 6 contrasts
  labs(x = "Simulation Seed", 
       y = expression(beta),
       title = "Parameter Recovery across 6 Contrasts (N=220)")

print(plot4)


## Calculating power and precision
## Inspecting the data frame, it seems like we recover the parameters with both sufficient power and precision
pow_pre_sim4 <- results_sim4 %>%
  group_by(parameter) %>%
  mutate(
    sig_check = case_when(                        # Determine if we detected the effect based on direction
      true_effect > 0 ~ ifelse(Q2.5 > 0, 1, 0),   # Positive detection
      true_effect < 0 ~ ifelse(Q97.5 < 0, 1, 0),  # Negative detection
      TRUE ~ 0
    ),

    precision_check = ifelse(Q97.5 - Q2.5 < 0.41, 1, 0)  # Check precision (width < 0.41)
  ) %>%
  summarise(
    power = mean(sig_check),
    precision = mean(precision_check),
    mean_estimate = mean(Estimate) # Check for bias
  )



##### Now we simulate for the actual multilevel model that we will use in the chains experiment
###############################################################################################
###############################################################################################

###### Note that the model slightly changes here (i.e. the addition of a random intercept for chain)
## So here are the steps in which we build up the final simulation
## 1. We first setup a slightly different data simulation function that accounts for the final multilevel structure of the experimental data we will observe. 
## 2. We run a 1 model simulation for sample size (n=220) so that we can specify the multilevel regression first. 
## 3. We finish by running a serious 1000 model simulation for the same sample size where we can directly pass the data from the 1000 iterations to the already compiled model


### Starting with Step 1: Setting up the slightly changed simulation function
###############################################################################################

data_sim_3 <- function(sim_seed, n_obs) {
  
  # Parameters
  mean_effect <- 0.3     # Positive effect (Gen 2, 4, 6)
  mean_effect_2 <- -0.3  # Negative effect (Gen 3, 5, 7)
  skew_ord <- 0.7        # Skew on latent scale
  chain_sd <- 0.3        # Variance between chains - thinking of lowest possible variance here (pilot data suggested between-chain variation was around 0.3)
  
  # Setting simulation seed
  set.seed(sim_seed)
  
  # Generate random intercepts (Latent Scale) for the chains 
  # We generate n_obs intercepts (one per chain)
  # This is because the number of observations "in each generation" denotes the number of chains in the experiment
  chain_intercepts <- rnorm(n_obs, mean = 0, sd = chain_sd)
  
  # Calculate "Chain-level" Latent Means for each Generation
  # That is, we add the average fixed effect to the random chain intercepts generated just above for each generation
  # Note: Gen 1 is the reference (0)
  mu_gen1 <- 0 + chain_intercepts
  mu_gen2 <- mean_effect + chain_intercepts
  mu_gen3 <- mean_effect_2 + chain_intercepts
  mu_gen4 <- mean_effect + chain_intercepts
  mu_gen5 <- mean_effect_2 + chain_intercepts
  mu_gen6 <- mean_effect + chain_intercepts
  mu_gen7 <- mean_effect_2 + chain_intercepts
  
  # Generate Data
  # We use sapply to pass the means one-by-one to rlikert.
  # This prevents the "negative probability" error because rlikert 
  # receives a single scalar mean each time, allowing it to handle negative latent means correctly.
  # This logic is written up in the generate_responses function
  
  generate_responses <- function(mu_vector) {
    sapply(mu_vector, function(m) {
      rlikert(size = 1, n_levels = 5, n_items = 1, mean = m, skew = skew_ord)
    })
  }
  
  # Generate data for all 7 generations
  gen1 <- generate_responses(mu_gen1)
  gen2 <- generate_responses(mu_gen2)
  gen3 <- generate_responses(mu_gen3)
  gen4 <- generate_responses(mu_gen4)
  gen5 <- generate_responses(mu_gen5)
  gen6 <- generate_responses(mu_gen6)
  gen7 <- generate_responses(mu_gen7)
  
  # D. Combine into Dataframe
  # We include chain_id so the random effect (1|chain_id) can be estimated
  tibble(
    chain_id = rep(1:n_obs, 7), 
    rating = c(gen1, gen2, gen3, gen4, gen5, gen6, gen7),
    gen = rep(c("Gen1", "Gen2", "Gen3", "Gen4", "Gen5", "Gen6", "Gen7"), each = n_obs)
  )
}


### Step 2: Running the one-model simulation so that the initial multilevel model can be specified
###############################################################################################

# Getting the toy dataset
test_multilevel_data <- data_sim_3(sim_seed = 1, n_obs = 220)

# Specifying the multilevel model now that we will use in the actual experiment
m2 <- brm(rating ~ 1 + gen + (1 | chain_id),
          data = test_multilevel_data,
          family = cumulative("probit"),
          prior = c(prior(normal(0, 1.5), class = "Intercept"),
                    prior(normal(0, 1), class = "b"),
                    prior(exponential(1), class = "sd")),
          warmup = 1000,
          iter = 3500,
          cores = 4,
          control = list(adapt_delta = 0.9))

## Inspect model output
# Seems like variance recovered roughly properly. Effect sizes in the right direction, though the estimate changes
summary(m2)
pp_check(m2)

## Saving model output
saveRDS(m2, "model2.RDS")
m2 <- readRDS("model2.RDS")

## Before proceeding, let's ensure once again whether our priors are indeed defensible and weakly regularising
# Let's go, prior predictive simulation
pps_m2 <- update(
  m2,
  sample_prior = "only",
  seed = 89
)

## Saving model output
saveRDS(pps_m2, "prior_predictive_sim_multilevel_model.RDS")
pps_m2 <- readRDS("prior_predictive_sim_multilevel_model.RDS")

## Inspecting model summary and predictions
## As you can see, no differences between the seven generations (Ignore the continous variable warning for the time-being)
## Our priors seem defensible
summary(pps_m2)
prediction_pp_sim2 <- conditional_effects(pps_m2)
prediction_pp_sim2 
prediction_pp_sim2_ordinal <- conditional_effects(pps_m2, categorical = TRUE)
prediction_pp_sim2_ordinal 

### Step 3: Running the serious simulation i.e. 1000 iterations of model fitting to a 1000 datasets
###############################################################################################

## Let's define how many simulations we want first
n_sim <- 1000

## Defining this to see how long it takes to simulate and analyse
t9 <- Sys.time()

## Simulate and analyse
sim5 <- 
  tibble(seed = 1:n_sim) %>%                                                      ## column to identify each simulation via its seed
  mutate(d = map(seed, data_sim_3, n_obs = 220)) %>%                              ## nested column where all the simulated datasets will be stored. In map(), we first indicate that we will iterate over values of seed. These values are fed to the data simulation function, and the final argument defines the sample size for each simulation
  mutate(b1 = map2(d, seed, ~update(m2, newdata = .x, seed = .y) %>%              ## This passes first two columns to run the model that we have defined above
                     posterior_summary() %>%                                      ## get summary of EVERYTHING
                     data.frame() %>%
                     rownames_to_column("parameter") %>%                          
                     filter(str_detect(parameter, "^b_gen|sd_chain"))))           ## filter out only the fixed effects and the random intercept

## How long did it take?
## Took nearly 2.1 days on the cluster for 1000 modls
t10 <- Sys.time()
t10 - t9

## Save the simulations
saveRDS(sim5, "simulation5_size220_multilevel.RDS")
sim5 <- readRDS("simulation5_size220_multilevel.RDS")


############### Now let's estimate power and precision

## As we are now recovering 7 parameters (6 contrasts + 1 chain level variation term), we need to make some quick changes to make the plots easier to read
## First, for each parameter, we need to highlight whether the effect was positive or negative
results_sim5 <- sim5 %>%
  unnest(b1) %>%
  mutate(true_effect = case_when(
    parameter %in% c("b_genGen3", "b_genGen5", "b_genGen7") ~ -0.3,                      # Assigns the true effect to these generations as -0.3
    parameter %in% c("b_genGen2", "b_genGen4", "b_genGen6") ~ 0.3,                       # Assigns the true effect to these generations as +0.3
    parameter %in% c("sd_chain_id__Intercept") ~ 0.3,
    TRUE ~ 0
  ))

## Now let's plot the simulation results for each parameter
## Seems to roughly recovering the parameters
plot5 <- results_sim5 %>%
  ggplot(aes(x = reorder(seed, Q2.5), y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_hline(yintercept = 0, color = "gray50", linetype = 2) +
  geom_hline(aes(yintercept = true_effect), color = "gray50", linetype = 2) +   # Now the code from lines 367 -371 will make sense
  geom_pointrange(fatten = 1/2, alpha = 0.5) +
  facet_wrap(~parameter, ncol = 3) +                                            # Facet by parameter to see all 6 contrasts
  labs(x = "Simulation Seed", 
       y = expression(beta),
       title = "Parameter Recovery across 6 Contrasts + Between-chain variation (N=220)")
plot5


## Calculating power and precision
## Inspecting the data frame, it seems like we recover the parameters with both sufficient power (~85 - 88%) and precision
pow_pre_sim5 <- results_sim5 %>%
  group_by(parameter) %>%
  mutate(
    sig_check = case_when(                        # Determine if we detected the effect based on direction
      true_effect > 0 ~ ifelse(Q2.5 > 0, 1, 0),   # Positive detection
      true_effect < 0 ~ ifelse(Q97.5 < 0, 1, 0),  # Negative detection
      TRUE ~ 0
    ),
    
    precision_check = ifelse(Q97.5 - Q2.5 < 0.41, 1, 0)  # Check precision (width < 0.41)
  ) %>%
  summarise(
    power = mean(sig_check),
    precision = mean(precision_check),
    mean_estimate = mean(Estimate) # Check for bias
  )

  
