######################## This script makes the regression coefficients table and plots the posteriors from mega mediation analyses
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
library(tidyr)
library(purrr)
library(flextable)
library(officer)
library(ggh4x)
library(emmeans)
library(marginaleffects)

## Set working directory to current study folder
setwd(file.path(renv::project(), "Mega-Analysis")) 

## Load the custom-made function for table Motivationting
source("Table Theme Function.R")

set.seed(42)



#####################################
### STATISTICAL MODEL LOADING

## Load all the pertinent statistical models (Comparing Gaussian Variants of monotonic models)
models <- list(
  "Self-Diagnosis"     = readRDS("MA_gaussian_mod1.rds"),
  "Other-Diagnosis"    = readRDS("MA_gaussian_mod2.rds"),
  "Other Help-Seeking" = readRDS("MA_gaussian_mod3.rds")
)




#####################################
### EXTRACT SUMMARY STATISTICS FOR EACH MODEL

## Function to extract pertinent summary statistics
extract_stats <- function(model, outcome_name) {
  draws <- as_draws_df(model)
  
  # Check if model is on logit (Odds Ratio) scale
  is_or <- grepl("Odds Ratio", outcome_name, ignore.case = TRUE)
  
  # Capture b_ (probit/logit), bsp_ (monotonic main), simo_ (simplexes), and sd_ (random effects)
  target_cols <- grep("^(b_|bsp_|simo_|sd_|sigma$)", names(draws), value = TRUE)
  
  draws[target_cols] |>
    pivot_longer(everything(), names_to = "raw_parameter", values_to = "draw") |>
    mutate(
      is_sd        = grepl("^sd_", raw_parameter),
      is_sigma     = raw_parameter == "sigma",
      is_intercept = grepl("^b_Intercept", raw_parameter),
      is_simplex   = grepl("^simo_", raw_parameter),
      is_main      = grepl("^(b_|bsp_)", raw_parameter) & !is_intercept,
      
      # Exponentiate main effects for logit models
      draw = if_else(is_or & is_main, exp(draw), draw)
    ) |>
    
    # Calculate summary statistics per parameter
    group_by(raw_parameter, is_sd, is_sigma, is_intercept, is_simplex, is_main) |>
    summarise(
      mean_est = mean(draw),
      ci_lower = quantile(draw, 0.025),
      ci_upper = quantile(draw, 0.975),
      pp_pos   = if (is_or) mean(draw > 1) else mean(draw > 0),
      .groups  = "drop"
    ) |>
    
    mutate(
      outcome = outcome_name,
      
      # Clean display labels across all three model variants
      parameter = case_when(
        # Categorical position dummy variables (Probit & Logit)
        raw_parameter == "b_chain_position2" ~ "Chain Position 2",
        raw_parameter == "b_chain_position3" ~ "Chain Position 3",
        raw_parameter == "b_chain_position4" ~ "Chain Position 4",
        raw_parameter == "b_chain_position5" ~ "Chain Position 5",
        raw_parameter == "b_chain_position6" ~ "Chain Position 6",
        raw_parameter == "b_chain_position7" ~ "Chain Position 7",
        
        # Day-of-week covariates (reference = Friday)
        raw_parameter == "b_DayMonday"    ~ "Day: Monday",
        raw_parameter == "b_DayTuesday"   ~ "Day: Tuesday",
        raw_parameter == "b_DayWednesday" ~ "Day: Wednesday",
        raw_parameter == "b_DayThursday"  ~ "Day: Thursday",
        raw_parameter == "b_DaySaturday"  ~ "Day: Saturday",
        
        # Study fixed effects (reference = Study 1)
        raw_parameter == "b_study_numexp2" ~ "Study 2",
        raw_parameter == "b_study_numexp3" ~ "Study 3",
        raw_parameter == "b_study_numexp4" ~ "Study 4",
        
        # Monotonic main effect
        is_main & grepl("mochain_position_trend", raw_parameter) ~ "Chain Position (Monotonic Trend)",
        
        # Monotonic Simplex steps
        is_simplex & grepl("\\[1\\]$", raw_parameter) ~ "Simplex: Step 1 \u2192 2",
        is_simplex & grepl("\\[2\\]$", raw_parameter) ~ "Simplex: Step 2 \u2192 3",
        is_simplex & grepl("\\[3\\]$", raw_parameter) ~ "Simplex: Step 3 \u2192 4",
        is_simplex & grepl("\\[4\\]$", raw_parameter) ~ "Simplex: Step 4 \u2192 5",
        is_simplex & grepl("\\[5\\]$", raw_parameter) ~ "Simplex: Step 5 \u2192 6",
        is_simplex & grepl("\\[6\\]$", raw_parameter) ~ "Simplex: Step 6 \u2192 7",
        is_simplex ~ sub("^simo_", "", raw_parameter),
        
        # Intercept cutpoints
        is_intercept ~ sub("^b_Intercept\\[(\\d+)\\]", "Intercept \\1", raw_parameter),
        
        # Between-chain variation standard deviations
        is_sd & grepl("Intercept", raw_parameter) ~ "SD: Between-Chain Baseline Variation",
        is_sd & grepl("mochain_position_trend", raw_parameter) ~ "SD: Between-Chain Effect Variation",
        is_sd ~ paste0("SD: ", sub("^sd_", "", raw_parameter)),
        
        # Single intercept (Gaussian models)
        raw_parameter == "b_Intercept" ~ "Intercept",
        
        # Residual SD (Gaussian models)
        is_sigma ~ "SD: Residual",
        
        TRUE ~ sub("^b_", "", raw_parameter)
      ),
      
      # Cell Formatting
      # - Intercepts in OR models left blank
      # - Fixed effects display Mean [95% CI] (Pr > 0 / Pr > 1)
      # - Simplexes, Cutpoints, & Random SDs display Mean [95% CI] without direction probability
      cell = case_when(
        is_or & (is_intercept | is_sd) ~ "",
        is_intercept | is_sd | is_sigma | is_simplex ~ sprintf("%.2f [%.2f, %.2f]", mean_est, ci_lower, ci_upper),
        is_main ~ sprintf("%.2f [%.2f, %.2f] (%.2f)", mean_est, ci_lower, ci_upper, pp_pos),
        TRUE ~ sprintf("%.2f [%.2f, %.2f]", mean_est, ci_lower, ci_upper)
      )
    )
}


## Extract the statistics now
results <- imap_dfr(models, extract_stats)




#####################################
### BUILD FORMATTED TABLE FOR EXPORT

## Pivot to wide format
reg_table <- results |>
  select(parameter, outcome, cell) |>
  pivot_wider(names_from = outcome, values_from = cell) |>
  rename(Predictor = parameter) |>
  mutate(across(-Predictor, ~ tidyr::replace_na(.x, ""))) |>   # monotonic-only rows otherwise show NA
  mutate(
    block = case_when(
      grepl("^Intercept",              Predictor) ~ 1,
      grepl("^Chain Position [2-7]$",  Predictor) ~ 2,
      grepl("Monotonic Trend",         Predictor) ~ 3,
      grepl("^Simplex",                Predictor) ~ 4,
      grepl("^Day:",                   Predictor) ~ 5,
      grepl("^Study",                  Predictor) ~ 6,
      grepl("^SD:",                    Predictor) ~ 7,
      TRUE                                        ~ 8
    ),
    day_order = match(
      sub("^Day: (\\w+).*", "\\1", Predictor),
      c("Monday", "Tuesday", "Wednesday", "Thursday", "Saturday")
    )
  ) |>
  arrange(block, day_order, Predictor) |>
  select(-block, -day_order)

print(reg_table)
outcome_cols <- names(models)  


## Build the table
ft <- flextable(reg_table) |>
  
  # Table Header label
  add_header_row(
    values    = c("", "Mean [95% Credible Interval] (Pr(\u03b2 > 0))"),
    colwidths = c(1, length(outcome_cols))
  ) |>
  
  # Column labels
  set_header_labels(
    Predictor = "Predictor",
    .list = setNames(as.list(outcome_cols), outcome_cols)
  ) |>
  
  # Merge the top-left blank cell across the two header rows
  merge_v(part = "header", j = "Predictor") |>
  
  # Alignment
  align(align = "center", part = "header") |>
  align(j = "Predictor", align = "left",   part = "body") |>
  align(j = outcome_cols, align = "center", part = "body") |>
  bold(part = "header") |>
  
  # Apply custom table Motivationting function
  apa_theme(
    caption      = "Mega Analysis: Bayesian regression results across Outcomes.",
    spanner_cols = 2:(length(outcome_cols) + 1),
    stripe_rows  = seq(2, nrow(reg_table), 2)
  )

ft  





#####################################
### MAKING THE MEDIATION PATH COEFFICIENT TABLE 

## Load the mediation models
mediation_models <- list(
  "Self-Diagnosis"     = readRDS("MA_mediation_mod1.rds"),
  "Other-Diagnosis"    = readRDS("MA_mediation_mod2.rds"),
  "Other Help-Seeking" = readRDS("MA_mediation_mod3.rds")
)


## Function to compute the path coefficients and derived mediation quantities from the posterior draws
# a = transmission -> semantic drift; b = semantic drift -> outcome; c' = direct effect of transmission on outcome
extract_paths <- function(model, outcome_name) {
  draws <- as_draws_df(model)
  nm    <- names(draws)
  
  # Locate the three estimated paths without hard-coding the outcome's response label
  a      <- draws[[grep("^bsp_CosineSimilarityz_mochain_position_trend$", nm, value = TRUE)]]
  b      <- draws[[grep("^b_(?!CosineSimilarityz).*_Cosine_Similarity_z$", nm, value = TRUE, perl = TRUE)]]
  direct <- draws[[grep("^bsp_(?!CosineSimilarityz).*_mochain_position_trend$", nm, value = TRUE, perl = TRUE)]]
  
  tibble(
    parameter = c("a: Transmission \u2192 Semantic Drift",
                  "b: Semantic Drift \u2192 Outcome",
                  "c\u2032: Direct Effect (Transmission \u2192 Outcome)",
                  "a \u00d7 b: Indirect Effect",
                  "c: Total Effect",
                  "Proportion Mediated (a \u00d7 b / c)"),
    draw      = list(a, b, direct, a * b, direct + a * b, (a * b) / (direct + a * b))
  ) |>
    mutate(
      outcome  = outcome_name,
      mean_est = map_dbl(draw, mean),
      ci_lower = map_dbl(draw, ~ quantile(.x, 0.025)),
      ci_upper = map_dbl(draw, ~ quantile(.x, 0.975)),
      pp_pos   = map_dbl(draw, ~ mean(.x > 0)),
      cell     = sprintf("%.3f [%.3f, %.3f] (%.2f)", mean_est, ci_lower, ci_upper, pp_pos)
    ) |>
    select(-draw)
}

## Extract and pivot to wide format (rows keep their path order)
path_table <- imap_dfr(mediation_models, extract_paths) |>
  mutate(parameter = factor(parameter, levels = unique(parameter))) |>
  select(parameter, outcome, cell) |>
  pivot_wider(names_from = outcome, values_from = cell) |>
  arrange(parameter) |>
  rename(Path = parameter)

print(path_table)
path_cols <- names(mediation_models)


## Build the table
ft_mediation <- flextable(path_table) |>
  
  # Table Header label
  add_header_row(
    values    = c("", "Mean [95% Credible Interval] (Pr(\u03b2 > 0))"),
    colwidths = c(1, length(path_cols))
  ) |>
  
  # Column labels
  set_header_labels(
    Path = "Path",
    .list = setNames(as.list(path_cols), path_cols)
  ) |>
  
  # Merge the top-left blank cell across the two header rows
  merge_v(part = "header", j = "Path") |>
  
  # Alignment
  align(align = "center", part = "header") |>
  align(j = "Path", align = "left",   part = "body") |>
  align(j = path_cols, align = "center", part = "body") |>
  bold(part = "header") |>
  
  # Apply custom table formatting function
  apa_theme(
    caption      = "Mega Analysis: Mediation of the transmission effect by semantic drift (standardized Gaussian models)",
    spanner_cols = 2:(length(path_cols) + 1),
    stripe_rows  = seq(2, nrow(path_table), 2)
  )

ft_mediation





#####################################
### PLOTTING POSTERIOR DISTRIBUTIONS NOW FOR PLOTTING POSTERIORS OF MEDIATION ANALYSES


#####################################
### DEFINE CUSTOM THEME

## Plot Style Theme
custom_theme <- theme(
  plot.title       = element_text(family = "Helvetica", size = 20, face = "bold", hjust = 0.5, margin = margin(b = 15)),
  axis.title.x     = element_text(family = "Helvetica", size = 15, face = "bold", margin = margin(t = 12)),
  axis.title.y     = element_text(family = "Helvetica", size = 15, face = "bold", margin = margin(r = 12)),
  axis.text.x      = element_text(family = "Helvetica", size = 11, color = "black"),
  axis.text.y      = element_text(family = "Helvetica", size = 11, color = "black"),
  
  # FACET STYLING
  strip.text       = element_text(family = "Helvetica", size = 13, face = "bold", color = "black"),
  strip.background = element_rect(fill = "grey92", color = "black", linewidth = 1),
  strip.placement  = "outside",                                                       
  
  # LEGEND FORMATTING
  legend.position  = "bottom",
  legend.title     = element_text(family = "Helvetica", size = 14, face = "bold"),
  legend.text      = element_text(family = "Helvetica", size = 13, face = "bold"),
  legend.key.size  = unit(1.1, "cm"),
  legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.5),
  
  panel.background = element_rect(fill = "white", color = "black", linewidth = 1),
  panel.grid.major = element_line(color = "grey92"),
  panel.grid.minor = element_blank()
)


#####################################
### EXTRACT POSTERIOR DRAWS & COMPUTE SUMMARY INTERVALS

## Function to pull the indirect (a*b) and direct (c') effect draws from each mediation model
extract_mediation_draws <- function(model, outcome_name) {
  draws  <- as_draws_df(model)
  nm     <- names(draws)
  a      <- draws[[grep("^bsp_CosineSimilarityz_mochain_position_trend$", nm, value = TRUE)]]
  b      <- draws[[grep("^b_(?!CosineSimilarityz).*_Cosine_Similarity_z$", nm, value = TRUE, perl = TRUE)]]
  direct <- draws[[grep("^bsp_(?!CosineSimilarityz).*_mochain_position_trend$", nm, value = TRUE, perl = TRUE)]]
  tibble(
    outcome = outcome_name,
    effect  = rep(c("Indirect Effect (a \u00d7 b)", "Direct Effect (c\u2032)"), each = length(a)),
    draw    = c(a * b, direct)
  )
}

mediation_draws <- imap_dfr(mediation_models, extract_mediation_draws) |>
  mutate(
    outcome = factor(outcome, levels = names(mediation_models)),
    effect  = factor(effect, levels = c("Indirect Effect (a \u00d7 b)", "Direct Effect (c\u2032)"))
  )

## Check the posterior range before fixing the x-axis limits below
range(mediation_draws$draw)

## Summary intervals at stacked vertical offsets (one row per outcome within each panel)
mediation_summary <- mediation_draws |>
  group_by(effect, outcome) |>
  summarise(
    mean_est    = mean(draw),
    ci_lower_95 = quantile(draw, 0.025),
    ci_upper_95 = quantile(draw, 0.975),
    .groups     = "drop"
  ) |>
  mutate(y_pos = -0.06 * as.numeric(outcome))

## Colour legend for outcomes
outcome_colors <- c(
  "Self-Diagnosis"     = "#0072B2",  # Blue
  "Other-Diagnosis"    = "#D55E00",  # Vermilion
  "Other Help-Seeking" = "#009E73"   # Green
)

## Density overlay plot
p_mediation <- ggplot() +
  # Fixed y-range for stacked CIs below the densities
  geom_blank(data = data.frame(y = c(-0.22, 1.05)), aes(y = y)) +
  
  # Null effect reference line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.8) +
  
  # Normalised density curves
  geom_density(
    data = mediation_draws,
    aes(x = draw, y = after_stat(ndensity), fill = outcome),
    color = NA, alpha = 0.40
  ) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.5) +
  
  # Horizontal 95% Credible Intervals and mean points
  geom_segment(
    data = mediation_summary,
    aes(x = ci_lower_95, xend = ci_upper_95, y = y_pos, yend = y_pos, color = outcome),
    linewidth = 1.1, show.legend = FALSE
  ) +
  geom_point(
    data = mediation_summary,
    aes(x = mean_est, y = y_pos, color = outcome),
    size = 2.4, show.legend = FALSE
  ) +
  
  # One panel per effect type, shared x-axis so the decomposition is visually comparable
  facet_wrap(~ effect, ncol = 2) +
  
  scale_x_continuous(limits = c(-0.1, 0.1), breaks = seq(-0.1, 0.1, by = 0.05), oob = scales::oob_keep) +
  scale_color_manual(values = outcome_colors, name = "Outcome") +
  scale_fill_manual(values = outcome_colors, name = "Outcome") +
  guides(fill = guide_legend(override.aes = list(alpha = 0.6))) +
  
  labs(
    x = "Posterior Estimate",
    y = "Density",
    title = "Posterior Distributions of Indirect and Direct Effects with 95% CIs and Means"
  ) +
  custom_theme

print(p_mediation)

# Export publication-ready plot
ggsave(
  filename = "Mediation Effect Posteriors.png",
  plot     = p_mediation,
  width    = 13,
  height   = 7,
  dpi      = 300
)


