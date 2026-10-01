######################## This script makes the regression coefficients table and plots the posteriors for the sensitivity analyses from Study 1 data 
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

## Set working directory to current study folder
setwd(file.path(renv::project(), "Study 1")) 

## Load the custom-made function for table formatting
source("Table Theme Function.R")

set.seed(239)



#####################################
### STATISTICAL MODEL LOADING

## Load all the pertinent statistical models
models <- list(
  "Self-Diagnosis" = readRDS("Study1_monotonic_mod1.rds"),
  "Other-Diagnosis" = readRDS("Study1_monotonic_mod2.rds"),
  "Other Help-Seeking" = readRDS("Study1_monotonic_mod3.rds")
)

 


#####################################
### EXTRACT SUMMARY STATISTICS FOR EACH MODEL

## Function to extract these summary statistics
extract_stats <- function(model, outcome_name) {
  draws <- as_draws_df(model)
  
  # Check if the model is a logit model and thus has co-effs on the log-odds scale
  is_or <- grepl("Odds Ratio", outcome_name, ignore.case = TRUE)
  
  # Select columns of interest to summarise
  target_cols <- grep("^(b_|bsp_|simo_|sd_)", names(draws), value = TRUE)
  
  draws[target_cols] |>
    pivot_longer(everything(), names_to = "raw_parameter", values_to = "draw") |>
    mutate(
      is_sd        = grepl("^sd_", raw_parameter),
      is_intercept = grepl("^b_Intercept", raw_parameter),
      is_simplex   = grepl("^simo_", raw_parameter),
      is_main      = grepl("^(b_|bsp_)", raw_parameter) & !is_intercept,
      
      # Exponentiate main effect if model is on Odds Ratio scale
      draw = if_else(is_or & is_main, exp(draw), draw)
    ) |>
    
    # Calculate summary statistics per parameter
    group_by(raw_parameter, is_sd, is_intercept, is_simplex, is_main) |>
    summarise(
      mean_est = mean(draw),
      ci_lower = quantile(draw, 0.025),
      ci_upper = quantile(draw, 0.975),
      pp_pos   = if (is_or) mean(draw > 1) else mean(draw > 0),
      .groups  = "drop"
    ) |>
    
    mutate(
      outcome = outcome_name,
      
      # Clean Parameter Display Labels
      parameter = case_when(
        is_main ~ "Chain Position (Monotonic Trend)",
        
        is_simplex & grepl("1\\[1\\]$", raw_parameter) ~ "Simplex: Step 1 \u2192 2",
        is_simplex & grepl("1\\[2\\]$", raw_parameter) ~ "Simplex: Step 2 \u2192 3",
        is_simplex & grepl("1\\[3\\]$", raw_parameter) ~ "Simplex: Step 3 \u2192 4",
        is_simplex ~ sub("^simo_", "", raw_parameter),
        
        is_intercept ~ sub("^b_Intercept\\[(\\d+)\\]", "Intercept \\1", raw_parameter),
        
        is_sd & grepl("Intercept", raw_parameter) ~ "SD: Between-Chain Baseline Variation",
        is_sd & grepl("mogen_trend|gen_trend", raw_parameter) ~ "SD: Between-Chain Effect Variation",
        is_sd ~ paste0("SD: ", raw_parameter),
        
        TRUE ~ raw_parameter
      ),
      
      # Table Formatting Logic:
      # 1. Main Predictors: Mean [95% CI] (Pr > 0 or Pr > 1)
      # 2. Simplexes, Cutpoints, & Random SDs: Mean [95% CI] (No Pr calculation)
      cell = case_when(
        is_main ~ sprintf("%.2f [%.2f, %.2f] (%.2f)", mean_est, ci_lower, ci_upper, pp_pos),
        TRUE    ~ sprintf("%.2f [%.2f, %.2f]", mean_est, ci_lower, ci_upper)
      ),
    )
}

## Extract the statistics now
results <- imap_dfr(models, extract_stats)




#####################################
### BUILD FORMATTED TABLE FOR EXPORT

## Building the table
reg_table <- results |>
  select(parameter, outcome, cell) |>
  pivot_wider(names_from = outcome, values_from = cell) |>
  rename(Predictor = parameter) |>
  arrange(
    case_when(
      grepl("^Intercept", Predictor) ~ 1,
      grepl("^Simplex", Predictor)   ~ 3,
      grepl("^SD:", Predictor)       ~ 4,
      TRUE                           ~ 2
    )
  )
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
  
  # Apply custom table formatting function
  apa_theme(
    caption      = "Study 1: Bayesian Monotonic regression results across Outcomes.",
    spanner_cols = 2:(length(outcome_cols) + 1),
    stripe_rows  = seq(2, nrow(reg_table), 2)
  )

ft  



#####################################
### PLOTTING POSTERIOR DISTRIBUTIONS NOW



#####################################
### DEFINE COLOR PALETTE & CUSTOM THEME

## Color legend for different outcomes
outcome_colors <- c(
  "Self-Diagnosis"     = "#0072B2",  # Blue
  "Other-Diagnosis"    = "#D55E00",  # Vermilion
  "Other Help-Seeking" = "#009E73"   # Green
)

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
  legend.position  = "none",
  
  panel.background = element_rect(fill = "white", color = "black", linewidth = 1),
  panel.grid.major = element_line(color = "grey92"),
  panel.grid.minor = element_blank(),

  
  # Defining aspect ratio of proportions
  aspect.ratio     = 0.55)


#####################################
### EXTRACT POSTERIOR DRAWS & COMPUTE SUMMARY INTERVALS

## Function to extract the pertinent posterior draws
extract_full_draws <- function(model, model_name) {
  draws <- as_draws_df(model)

  gen_cols <- grep("^bsp_", names(draws), value = TRUE)
  
  draws[gen_cols] |>
    pivot_longer(everything(), names_to = "raw_parameter", values_to = "draw") |>
    mutate(
      parameter = "Chain Position (Monotonic Trend)",
      outcome   = model_name
    )
}

## Extracting the posteriors from the model now
posterior_draws <- imap_dfr(models, extract_full_draws) |>
  mutate(
    outcome = factor(outcome, levels = c("Self-Diagnosis", "Other-Diagnosis", "Other Help-Seeking"))
  )

# Compute Mean and 95% Credible Intervals at fixed vertical offsets
posterior_summary <- posterior_draws |>
  group_by(outcome, parameter) |>
  summarise(
    mean_est    = mean(draw),
    ci_lower_95 = quantile(draw, 0.025),
    ci_upper_95 = quantile(draw, 0.975),
    .groups     = "drop"
  ) |>
  mutate(
    # Exact fixed vertical positions across all panels for perfectly uniform spacing between the stacked credible intervals
    y_pos = -0.05
  )


#####################################
### CONFIGURE DUMMY BOUNDS FOR UNIFORM BREATHING SPACE

# Anchor x-axis endpoints and y-axis baseline/peak space
dummy_x_limits <- data.frame(x = c(-0.4, 0.4))



#####################################
### DENSITY OVERLAY PLOT GENERATION

p_density <- ggplot() +
  # Null reference line at zero
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.8) +
  geom_blank(data = dummy_x_limits, aes(x = x)) +

  
  # Density curves per outcome
  geom_density(
    data = posterior_draws,
    aes(x = draw, y = after_stat(ndensity), fill = outcome),
    color = NA,
    alpha = 0.45
  ) +
  
  # Solid baseline reference line at y = 0
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.5) +
  
  # 95% Credible Interval bar
  geom_segment(
    data = posterior_summary,
    aes(x = ci_lower_95, xend = ci_upper_95, y = y_pos, yend = y_pos, color = outcome),
    linewidth = 1.1
  ) +
  
  # Posterior Mean point
  geom_point(
    data = posterior_summary,
    aes(x = mean_est, y = y_pos, color = outcome),
    size = 2.6
  ) +
  
  # Facet panel by Outcome
  facet_wrap(~ outcome, ncol = 3) +
  
  # Set x-axis limits
  scale_x_continuous(
    limits = c(-0.4, 0.4),
    breaks = round(seq(-0.4, 0.4, by = 0.1), 1)
  ) +
  
  scale_color_manual(values = outcome_colors) +
  scale_fill_manual(values = outcome_colors) +
  
  labs(
    x = "Posterior Estimate (\u03b2)",
    y = "Normalized Density",
    title = "Posterior Distributions of Chain Position (Monotonic Variant)"
  ) +
  custom_theme

print(p_density)

# Export publication-ready plot
ggsave(
  filename = "Posterior Distributions_SA.png",
  plot     = p_density,
  width    = 13,
  height   = 5.5,
  dpi      = 300
)





