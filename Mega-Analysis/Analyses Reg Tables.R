######################## This script makes the regression coefficients table and plots the posteriors from mega-analyses
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

set.seed(22)



#####################################
### STATISTICAL MODEL LOADING

## Load all the pertinent statistical models (Simple Outcome Analyses)
models <- list(
  "Self-Diagnosis"                   = readRDS("MA_mod1.rds"),     
  "Self-Diagnosis (Odds Ratio)"      = readRDS("MA_OR_mod1.rds"),
  "Self-Diagnosis (Monotonic)"      = readRDS("MA_monotonic_mod1.rds"),

  "Other-Diagnosis"                  = readRDS("MA_mod2.rds"),
  "Other-Diagnosis (Odds Ratio)"     = readRDS("MA_OR_mod2.rds"),
  "Other-Diagnosis (Monotonic)"      = readRDS("MA_monotonic_mod2.rds"),
  
  "Other Help-Seeking"               = readRDS("MA_mod3.rds"),
  "Other Help-Seeking (Odds Ratio)"  = readRDS("MA_OR_mod3.rds"),
  "Other Help-Seeking (Monotonic)"      = readRDS("MA_monotonic_mod3.rds")
)





#####################################
### EXTRACT SUMMARY STATISTICS FOR EACH MODEL

## Function to extract pertinent summary statistics
extract_stats <- function(model, outcome_name) {
  draws <- as_draws_df(model)
  
  # Check if model is on logit (Odds Ratio) scale
  is_or <- grepl("Odds Ratio", outcome_name, ignore.case = TRUE)
  
  # Capture b_ (probit/logit), bsp_ (monotonic main), simo_ (simplexes), and sd_ (random effects)
  target_cols <- grep("^(b_|bsp_|simo_|sd_)", names(draws), value = TRUE)
  
  draws[target_cols] |>
    pivot_longer(everything(), names_to = "raw_parameter", values_to = "draw") |>
    mutate(
      is_sd        = grepl("^sd_", raw_parameter),
      is_intercept = grepl("^b_Intercept", raw_parameter),
      is_simplex   = grepl("^simo_", raw_parameter),
      is_main      = grepl("^(b_|bsp_)", raw_parameter) & !is_intercept,
      
      # Exponentiate main effects for logit models
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
        
        TRUE ~ sub("^b_", "", raw_parameter)
      ),
      
      # Cell Motivationting
      # - Intercepts in OR models left blank
      # - Fixed effects display Mean [95% CI] (Pr > 0 / Pr > 1)
      # - Simplexes, Cutpoints, & Random SDs display Mean [95% CI] without direction probability
      cell = case_when(
        is_or & (is_intercept | is_sd) ~ "",
        is_intercept | is_sd | is_simplex ~ sprintf("%.2f [%.2f, %.2f]", mean_est, ci_lower, ci_upper),
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
### PLOTTING POSTERIOR DISTRIBUTIONS NOW (FIRST FOR MAIN ANALYSES)


#####################################
### DEFINE COLOR PALETTE & CUSTOM THEME

## Color legend for different predictors
chain_pos_colors <- c(
  "Chain Position 2" = "#FF1744",  # Vivid Pink/Coral
  "Chain Position 3" = "#00B0FF",  # Bright Cyan
  "Chain Position 4" = "#00E676",   # Neon Green
  "Chain Position 5" = "#FF9100",  # Amber / Orange
  "Chain Position 6" = "#AA00FF",  # Vivid Purple
  "Chain Position 7" = "#FF007F",  # Deep Magenta
  "Chain Position (Monotonic Trend)" = "#2CD"   
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

## Function to extract the pertinent posterior draws across distinct model types
extract_full_draws <- function(model, model_name) {
  draws    <- as_draws_df(model)
  is_or    <- grepl("Odds Ratio", model_name, ignore.case = TRUE)
  is_mono  <- grepl("Monotonic", model_name, ignore.case = TRUE)
  all_cols <- names(draws)
  
  if (is_mono) {
    mono_cols <- grep("^bsp_.*mochain_position_trend", all_cols, value = TRUE)
    
    draws[mono_cols] |>
      pivot_longer(everything(), names_to = "raw_param", values_to = "draw") |>
      mutate(
        parameter = "Chain Position (Monotonic Trend)",
        draw      = draw,
        outcome   = sub(" \\(Monotonic\\)", "", model_name),
        scale     = "Monotonic Model (\u03b2)"
      )
  } else {
    cat_cols <- grep("^(b_chain_position|b_gen)[2-7]$", all_cols, value = TRUE)
    
    draws[cat_cols] |>
      pivot_longer(everything(), names_to = "raw_param", values_to = "draw") |>
      mutate(
        pos_num   = sub("^(b_chain_position|b_gen)", "", raw_param),
        parameter = paste("Chain Position", pos_num),
        draw      = if (is_or) exp(draw) else draw,
        outcome   = sub(" \\(Odds Ratio\\)", "", model_name),
        scale     = if (is_or) "Logit Model (Odds Ratio)" else "Probit Model (Standardized Latent \u03b2)"
      )
  }
}

## Extracting the posteriors from the model now
posterior_draws <- imap_dfr(models, extract_full_draws) |>
  mutate(
    outcome   = factor(outcome, levels = c("Self-Diagnosis", "Other-Diagnosis", "Other Help-Seeking")),
    parameter = factor(parameter,levels = c(
      "Chain Position 2", "Chain Position 3", "Chain Position 4",
      "Chain Position 5", "Chain Position 6", "Chain Position 7",
      "Chain Position (Monotonic Trend)"
    )
    ),
    scale     = factor(
      scale,
      levels = c(
        "Probit Model (Standardized Latent \u03b2)",
        "Logit Model (Odds Ratio)",
        "Monotonic Model (\u03b2)"
      )
    )
  )

# Compute Mean and 95% Credible Intervals at fixed vertical offsets
posterior_summary <- posterior_draws |>
  group_by(scale, outcome, parameter) |>
  summarise(
    mean_est    = mean(draw),
    ci_lower_95 = quantile(draw, 0.025),
    ci_upper_95 = quantile(draw, 0.975),
    .groups     = "drop"
  ) |>
  mutate(
    y_pos = case_when(
      parameter == "Chain Position 2"                 ~ -0.05,
      parameter == "Chain Position 3"                 ~ -0.10,
      parameter == "Chain Position 4"                 ~ -0.15,
      parameter == "Chain Position 5"                 ~ -0.20,
      parameter == "Chain Position 6"                 ~ -0.25,
      parameter == "Chain Position 7"                 ~ -0.30,
      parameter == "Chain Position (Monotonic Trend)" ~ -0.10
    )
  )

# Null reference lines across all three model scales
ref_lines <- data.frame(
  scale = factor(
    c("Probit Model (Standardized Latent \u03b2)", 
      "Logit Model (Odds Ratio)", 
      "Monotonic Model (\u03b2)"),
    levels = c("Probit Model (Standardized Latent \u03b2)", 
               "Logit Model (Odds Ratio)", 
               "Monotonic Model (\u03b2)")
  ),
  xintercept = c(0, 1, 0)
)


#####################################
### CONFIGURE X-AXIS LIMITS 

# Define fixed x-axis limits per scale column
probit_x_limits <- c(-1.0, 1.0)  # Fixed scale for Standardized Latent Betas
or_x_limits     <- c(0.0, 7.0)   # Fixed scale for Odds Ratios
mono_x_limits   <- c(-0.4, 0.4)

dummy_x_limits <- data.frame(
  scale = factor(
    c("Probit Model (Standardized Latent \u03b2)", "Probit Model (Standardized Latent \u03b2)",
      "Logit Model (Odds Ratio)", "Logit Model (Odds Ratio)",
      "Monotonic Model (\u03b2)", "Monotonic Model (\u03b2)"),
    levels = c("Probit Model (Standardized Latent \u03b2)", 
               "Logit Model (Odds Ratio)", 
               "Monotonic Model (\u03b2)")
  ),
  x = c(probit_x_limits, or_x_limits, mono_x_limits)
)

# Dummy y-limits to standardize panel vertical bounds (-0.35 baseline space up (for more stacked credible intervals) to 1.05 peak height)
dummy_y_limits <- data.frame(y = c(-0.35, 1.05))



#####################################
### DENSITY OVERLAY PLOT GENERATION

p_density <- ggplot() +
  # Force uniform x-limits per column scale & fixed y-range
  geom_blank(data = dummy_x_limits, aes(x = x)) +
  geom_blank(data = dummy_y_limits, aes(y = y)) +
  
  # Null effect reference lines (0 for Probit & Monotonic, 1 for Odds Ratio)
  geom_vline(
    data = ref_lines,
    aes(xintercept = xintercept),
    linetype = "dashed", color = "grey35", linewidth = 0.8
  ) +
  
  # Normalized density curves
  geom_density(
    data = posterior_draws,
    aes(x = draw, y = after_stat(ndensity), fill = parameter),
    color = NA,
    alpha = 0.40
  ) +
  
  # Baseline reference line at y = 0
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.5) +
  
  # Horizontal 95% Credible Intervals
  geom_segment(
    data = posterior_summary,
    aes(
      x = ci_lower_95, xend = ci_upper_95,
      y = y_pos, yend = y_pos,
      color = parameter
    ),
    linewidth = 1.1,
    show.legend = FALSE
  ) +
  
  # Mean posterior estimate points
  geom_point(
    data = posterior_summary,
    aes(x = mean_est, y = y_pos, color = parameter),
    size = 2.4,
    show.legend = FALSE
  ) +
  
  # 3x3 Facet Grid
  facet_grid(outcome ~ scale, scales = "free", switch = "y") +
  
  # Custom X-Axis Breaks per Facet Scale using ggh4x
  facetted_pos_scales(
    x = list(
      scale == "Probit Model (Standardized Latent \u03b2)" ~ scale_x_continuous(
        limits = c(-1.0, 1.0),
        breaks = seq(-1.0, 1.0, by = 0.2),
        oob    = scales::oob_keep
      ),
      scale == "Logit Model (Odds Ratio)" ~ scale_x_continuous(
        limits = c(0.0, 7.0),
        breaks = seq(0.0, 7.0, by = 1.0),
        oob    = scales::oob_keep
      ),
      scale == "Monotonic Model (\u03b2)" ~ scale_x_continuous(
        limits = c(-0.4, 0.4),
        breaks = seq(-0.4, 0.4, by = 0.1),
        oob    = scales::oob_keep
      )
    )
  ) +
  
  # Scale Mappings
  scale_color_manual(values = chain_pos_colors, name = "Predictor") +
  scale_fill_manual(values = chain_pos_colors, name = "Predictor") +
  
  guides(
    fill = guide_legend(override.aes = list(alpha = 0.6))
  ) +
  
  labs(
    x = "Posterior Estimate",
    y = "Density",
    title = "Full Posterior Distributions with Stacked 95% CIs and Means"
  ) +
  custom_theme

print(p_density)


# Export publication-ready plot
ggsave(
  filename = "Posterior Distributions.png",
  plot     = p_density,
  width    = 13,
  height   = 10,
  dpi      = 300
)





