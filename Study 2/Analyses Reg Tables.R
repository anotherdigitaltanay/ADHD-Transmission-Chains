######################## This script makes the regression coefficients table and plots the posteriors from Study 2 data
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
setwd(file.path(renv::project(), "Study 2")) 

## Load the custom-made function for table formatting
source("Table Theme Function.R")

set.seed(39)



#####################################
### STATISTICAL MODEL LOADING

## Load all the pertinent statistical models (Simple Outcome Analyses)
models <- list(
  "Self-Diagnosis" = readRDS("Study2_mod1.rds"),
  "Self-Diagnosis (Odds Ratio)" = readRDS("Study2_OR_mod1.rds"),
  "Self-Diagnosis (Monotonic)" = readRDS("Study2_monotonic_mod1.rds"),
  
  "Other-Diagnosis" = readRDS("Study2_mod2.rds"),
  "Other-Diagnosis (Odds Ratio)" = readRDS("Study2_OR_mod2.rds"),
  "Other-Diagnosis (Monotonic)" = readRDS("Study2_monotonic_mod2.rds"),
  
  "Other Help-Seeking" = readRDS("Study2_mod3.rds"),
  "Other Help-Seeking (Odds Ratio)" = readRDS("Study2_OR_mod3.rds"),
  "Other Help-Seeking (Monotonic)" = readRDS("Study2_monotonic_mod3.rds")
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
      
      # Cell Formatting
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
  arrange(
  case_when(
    grepl("^Intercept", Predictor) ~ 1,
    grepl("^SD:", Predictor)       ~ 3,
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
    caption      = "Study 2: Bayesian regression results across Outcomes.",
    spanner_cols = 2:(length(outcome_cols) + 1),
    stripe_rows  = seq(2, nrow(reg_table), 2)
  )

ft  




#####################################
### STATISTICAL MODEL LOADING (REPEATING THE TABLE GENERATING PROCESS FOR THE INTERACTION MODELS)

## Load all the pertinent statistical models (Simple Outcome Analyses)
int_models <- list(
  "Self-Diagnosis" = readRDS("Study2_mod4.rds"),
  "Self-Diagnosis (Monotonic)" = readRDS("Study2_monotonic_mod4.rds"),
  
  "Other-Diagnosis" = readRDS("Study2_mod5.rds"),
  "Other-Diagnosis (Monotonic)" = readRDS("Study2_monotonic_mod5.rds"),
  
  "Other Help-Seeking" = readRDS("Study2_mod6.rds"),
  "Other Help-Seeking (Monotonic)" = readRDS("Study2_monotonic_mod6.rds")
)



#####################################
### EXTRACT SUMMARY STATISTICS FOR EACH INTERACTION MODEL 

## Extracting pertinent summary statistics
extract_stats_interaction <- function(model, outcome_name) {
  
  # Calculate transmission type main effect using the marginal effects package
  cmp <- avg_comparisons(model, variables = "transmission_type", type = "link")
  
  # Extract the numeric vector of posterior draws from the 'draw' column
  main_tx_draws <- get_draws(cmp)$draw
  tx_main_df <- tibble(
    raw_parameter = "main_transmission_type",
    draw          = main_tx_draws
  )
  
  # Extract Other Posterior Draws now
  draws <- as_draws_df(model)
  all_cols <- names(draws)
  is_monotonic <- any(grepl("mochain_position_trend", all_cols))
  
  if (!is_monotonic) {
    # Categorical model interaction terms (Positions 1-7 Contrasts)
    int_cols <- c(
      "b_transmission_typewords",
      all_cols[startsWith(all_cols, "b_chain_position") & endsWith(all_cols, ":transmission_typewords")]
    )
  } else {
    # Monotonic model interaction term
    int_cols <- all_cols[startsWith(all_cols, "bsp_mochain_position_trend") & endsWith(all_cols, ":transmission_typewords")]
  }
  
  sd_cols <- all_cols[startsWith(all_cols, "sd_")]
  target_cols <- unique(c(int_cols, sd_cols))
  
  # Calculate summary statistics
  param_summary <- draws |>
    select(all_of(target_cols)) |>
    pivot_longer(everything(), names_to = "raw_parameter", values_to = "draw") |>
    bind_rows(tx_main_df) |>
    mutate(
      is_sd      = startsWith(raw_parameter, "sd_"),
      is_main_tx = raw_parameter == "main_transmission_type"
    ) |>
    group_by(raw_parameter, is_sd, is_main_tx) |>
    summarise(
      mean_est = mean(draw),
      ci_lower = quantile(draw, 0.025),
      ci_upper = quantile(draw, 0.975),
      pp_pos   = mean(draw > 0),
      .groups  = "drop"
    )
  
  # Clean up variable names
  param_summary |>
    mutate(
      outcome = outcome_name,
      parameter = case_when(
        is_main_tx ~ "Transmission Type (Main Effect: Narrative Summaries - Examples)",
        
        # Random Effects Standard Deviations
        is_sd & grepl("Intercept", raw_parameter)              ~ "SD: Between-Chain Baseline Variation",
        is_sd & grepl("mochain_position_trend", raw_parameter) ~ "SD: Between-Chain Effect Variation",
        is_sd                                                  ~ paste0("SD: ", sub("^sd_", "", raw_parameter)),
        
        # Categorical Interaction Terms
        raw_parameter == "b_transmission_typewords" ~ "Chain Position 1 \u00d7 Transmission Format",
        grepl("chain_position2", raw_parameter)     ~ "Chain Position 2 \u00d7 Transmission Format",
        grepl("chain_position3", raw_parameter)     ~ "Chain Position 3 \u00d7 Transmission Format",
        grepl("chain_position4", raw_parameter)     ~ "Chain Position 4 \u00d7 Transmission Format",
        grepl("chain_position5", raw_parameter)     ~ "Chain Position 5 \u00d7 Transmission Format",
        grepl("chain_position6", raw_parameter)     ~ "Chain Position 6 \u00d7 Transmission Format",
        grepl("chain_position7", raw_parameter)     ~ "Chain Position 7 \u00d7 Transmission Format",
        
        # Monotonic Interaction Term
        startsWith(raw_parameter, "bsp_mochain_position_trend") ~ "Chain Position (Monotonic Trend) \u00d7 Transmission Format",
        
      
        TRUE ~ raw_parameter
      ),
      cell = case_when(
        is_sd ~ sprintf("%.2f [%.2f, %.2f]", mean_est, ci_lower, ci_upper),
        TRUE  ~ sprintf("%.2f [%.2f, %.2f] (%.2f)", mean_est, ci_lower, ci_upper, pp_pos)
      )
    )
}

int_results <- imap_dfr(int_models, extract_stats_interaction)



#####################################
### BUILD FORMATTED TABLE FOR EXPORT


## Pivot to a wide table formal
reg_table_int <- int_results |>
  select(parameter, outcome, cell) |>
  distinct(parameter, outcome, .keep_all = TRUE) |>
  pivot_wider(names_from = outcome, values_from = cell) |>
  rename(Predictor = parameter) |>
  mutate(across(-Predictor, ~ tidyr::replace_na(.x, ""))) |>
  arrange(
    case_when(
      grepl("Main Effect", Predictor)          ~ 1,
      grepl("Chain Position [1-7]", Predictor) ~ 2,
      grepl("Monotonic Trend", Predictor)      ~ 3,
      grepl("^SD:", Predictor)                  ~ 4,
      TRUE                                     ~ 5
    )
  )

outcome_cols <- names(int_models)


# Build the flextable
ft_int <- flextable(reg_table_int) |>
  add_header_row(
    values    = c("Predictor", "Mean [95% Credible Interval] (Pr(\u03b2 > 0))"),
    colwidths = c(1, length(outcome_cols))
  ) |>

  # Column Labels
  set_header_labels(
    Predictor = "Predictor",
    .list     = setNames(as.list(outcome_cols), outcome_cols)
  ) |>
  
  # Merge vertical header for Predictor column
  merge_v(part = "header", j = "Predictor") |>
  
  # Alignment & Formatting
  align(align = "center", part = "header") |>
  align(j = "Predictor", align = "left",   part = "body") |>
  align(j = outcome_cols, align = "center", part = "body") |>
  bold(part = "header") |>
  
  apa_theme(
    caption      = "Study 2: Bayesian interaction regression models (Chain Position \u00d7 Transmission Format).",
    spanner_cols = 2:(length(outcome_cols) + 1),
    stripe_rows  = seq(2, nrow(reg_table_int), 2)
  )

ft_int




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




#####################################
### PLOTTING POSTERIOR DISTRIBUTIONS NOW (FOR MODERATION EFFECTS)

#####################################
### DEFINE COLOR PALETTE MATCHING MAIN EFFECTS PLOT

int_colors <- c(
  "Transmission Format (Main Effect)"             = "#24f",  
  "Chain Position 1"    = "#D50000",  # Crimson 
  "Chain Position 2"    = "#FF1744",  # Vivid Pink/Coral 
  "Chain Position 3"    = "#00B0FF",  # Bright Cyan 
  "Chain Position 4"    = "#00E676",  # Neon Green 
  "Chain Position 5"    = "#FF9100",  # Amber / Orange 
  "Chain Position 6"    = "#AA00FF",  # Vivid Purple 
  "Chain Position 7"    = "#FF007F",  # Deep Magenta 
  "Monotonic Trend"     = "#2CD"      # Cyan/Teal 
)



#####################################
### COLUMN HEADER LABELS (CONTRAST DIRECTION)

col_label_main  <- "Format Main Effect\n(Narrative Summaries \u2212 Examples)"
col_label_cat   <- "Categorical Model Contrasts\n(Narrative Summaries \u2212 Examples)"
col_label_mono  <- "Monotonic Model Contrast\n(Narrative Summaries \u2212 Examples)"


#####################################
### EXTRACT POSTERIOR DRAWS FOR INTERACTION MODELS

base_outcomes <- c("Self-Diagnosis", "Other-Diagnosis", "Other Help-Seeking")

extract_int_full_draws <- function(outcomes_list, int_models) {
  map_dfr(outcomes_list, function(out_name) {
    cat_model  <- int_models[[out_name]]
    mono_model <- int_models[[paste0(out_name, " (Monotonic)")]]
    
    # Extract Main Effect
    cmp <- avg_comparisons(cat_model, variables = "transmission_type", type = "link")
    main_draws <- tibble(
      draw      = get_draws(cmp)$draw,
      parameter = "Transmission Format (Main Effect)",
      scale     = col_label_main,
      outcome   = out_name
    )
    
    # Extract Categorical Interaction terms
    cat_df   <- as_draws_df(cat_model)
    cat_cols <- c(
      grep("^b_transmission_typewords$", names(cat_df), value = TRUE),
      grep("^b_chain_position.*:transmission_typewords$", names(cat_df), value = TRUE)
    )
    
    cat_draws <- cat_df |>
      select(all_of(cat_cols)) |>
      pivot_longer(everything(), names_to = "raw_param", values_to = "draw") |>
      mutate(
        pos_num = case_when(
          grepl("chain_position2", raw_param) ~ "2",
          grepl("chain_position3", raw_param) ~ "3",
          grepl("chain_position4", raw_param) ~ "4",
          grepl("chain_position5", raw_param) ~ "5",
          grepl("chain_position6", raw_param) ~ "6",
          grepl("chain_position7", raw_param) ~ "7",
          TRUE ~ "1"
        ),
        parameter = paste0("Chain Position ", pos_num),
        scale     = col_label_cat,
        outcome   = out_name
      ) |>
      select(draw, parameter, scale, outcome)
    
    # Extract Monotonic Contrast
    mono_df   <- as_draws_df(mono_model)
    mono_cols <- grep("^bsp_mochain_position_trend.*:transmission_typewords$", names(mono_df), value = TRUE)
    
    mono_draws <- mono_df |>
      select(all_of(mono_cols)) |>
      pivot_longer(everything(), names_to = "raw_param", values_to = "draw") |>
      mutate(
        parameter = "Monotonic Trend",
        scale     = col_label_mono,
        outcome   = out_name
      ) |>
      select(draw, parameter, scale, outcome)
    
    bind_rows(main_draws, cat_draws, mono_draws)
  })
}

# Extract draws & format factor levels
int_posterior_draws <- extract_int_full_draws(base_outcomes, int_models) |>
  mutate(
    outcome   = factor(outcome, levels = base_outcomes),
    scale     = factor(
      scale,
      levels = c(col_label_main, col_label_cat, col_label_mono)
    ),
    parameter = factor(
      parameter,
      levels = c(
        "Transmission Format (Main Effect)",
        "Chain Position 1",
        "Chain Position 2",
        "Chain Position 3",
        "Chain Position 4",
        "Chain Position 5",
        "Chain Position 6",
        "Chain Position 7",
        "Monotonic Trend"
      )
    )
  )


#####################################
### COMPUTE SUMMARY INTERVALS & STRICT VERTICAL OFFSETS

int_posterior_summary <- int_posterior_draws |>
  group_by(scale, outcome, parameter) |>
  summarise(
    mean_est    = mean(draw),
    ci_lower_95 = quantile(draw, 0.025),
    ci_upper_95 = quantile(draw, 0.975),
    .groups     = "drop"
  ) |>
  mutate(
    y_pos = case_when(
      grepl("Position 1", parameter) ~ -0.05,
      grepl("Position 2", parameter) ~ -0.10,
      grepl("Position 3", parameter) ~ -0.15,
      grepl("Position 4", parameter) ~ -0.20,
      grepl("Position 5", parameter) ~ -0.25,
      grepl("Position 6", parameter) ~ -0.30,
      grepl("Position 7", parameter) ~ -0.35,
      TRUE                           ~ -0.10
    )
  )

# Null reference lines at 0 across all interaction columns
ref_lines_int <- data.frame(
  scale = factor(
    c(col_label_main, col_label_cat, col_label_mono),
    levels = c(col_label_main, col_label_cat, col_label_mono)
  ),
  xintercept = c(0, 0, 0)
)


#####################################
### CONFIGURE AXIS LIMITS

# Unified -1.0 to 1.0 range for Column 1 & Column 2
dummy_x_limits_int <- data.frame(
  scale = factor(
    c(col_label_main, col_label_main,
      col_label_cat,  col_label_cat,
      col_label_mono, col_label_mono),
    levels = c(col_label_main, col_label_cat, col_label_mono)
  ),
  x = c(-1.0, 1.0, -1.0, 1.0, -0.4, 0.4)
)

# Y-axis bounds to accommodate stacked 95% CIs 
dummy_y_limits_int <- data.frame(y = c(-0.42, 1.05))


#####################################
### DENSITY OVERLAY PLOT GENERATION

p_density_int <- ggplot() +
  # Force uniform x-limits per column scale & fixed vertical bounds
  geom_blank(data = dummy_x_limits_int, aes(x = x)) +
  geom_blank(data = dummy_y_limits_int, aes(y = y)) +
  
  # Null effect reference line at 0
  geom_vline(
    data = ref_lines_int,
    aes(xintercept = xintercept),
    linetype = "dashed", color = "grey35", linewidth = 0.8
  ) +
  
  # Normalized density curves
  geom_density(
    data = int_posterior_draws,
    aes(x = draw, y = after_stat(ndensity), fill = parameter),
    color = NA,
    alpha = 0.40
  ) +
  
  # Baseline reference line at y = 0
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.5) +
  
  # Horizontal 95% Credible Intervals
  geom_segment(
    data = int_posterior_summary,
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
    data = int_posterior_summary,
    aes(x = mean_est, y = y_pos, color = parameter),
    size = 2.4,
    show.legend = FALSE
  ) +
  
  # 3x3 Facet Grid
  facet_grid(outcome ~ scale, scales = "free", switch = "y") +
  
  # Custom X-Axis Scales per Column
  facetted_pos_scales(
    x = list(
      scale == col_label_main ~ scale_x_continuous(
        limits = c(-1.0, 1.0),
        breaks = seq(-1.0, 1.0, by = 0.2),
        oob    = scales::oob_keep
      ),
      scale == col_label_cat ~ scale_x_continuous(
        limits = c(-1.0, 1.0),
        breaks = seq(-1.0, 1.0, by = 0.2),
        oob    = scales::oob_keep
      ),
      scale == col_label_mono ~ scale_x_continuous(
        limits = c(-0.4, 0.4),
        breaks = seq(-0.4, 0.4, by = 0.1),
        oob    = scales::oob_keep
      )
    )
  ) +
  
  # Color & Fill Mappings
  scale_color_manual(values = int_colors, name = "Predictor") +
  scale_fill_manual(values = int_colors, name = "Predictor") +
  
  guides(
    fill = guide_legend(override.aes = list(alpha = 0.6), nrow = 3, byrow = TRUE)
  ) +
  
  labs(
    x = "Posterior Estimate (\u03b2)",
    y = "Density",
    title = "Interaction Models Posterior Distributions with Stacked 95% CIs and Means"
  ) +
  custom_theme

print(p_density_int)

# Export plot
ggsave(
  filename = "Posterior Distributions_Moderation.png",
  plot     = p_density_int,
  width    = 13,
  height   = 10,
  dpi      = 300
)
