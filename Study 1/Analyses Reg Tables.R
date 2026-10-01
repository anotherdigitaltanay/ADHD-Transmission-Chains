######################## This script makes the regression coefficients table and plots the posteriors from Study 1 data
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

## Set working directory to current study folder
setwd(file.path(renv::project(), "Study 1")) 

## Load the custom-made function for table formatting
source("Table Theme Function.R")

set.seed(43539)




#####################################
### STATISTICAL MODEL LOADING

## Load all the pertinent statistical models
models <- list(
  "Self-Diagnosis" = readRDS("Study1_mod1.RDS"),
  "Self-Diagnosis (Odds Ratio)" = readRDS("Study1_OR_mod1.RDS"),
  
  "Other-Diagnosis" = readRDS("Study1_mod2.RDS"),
  "Other-Diagnosis (Odds Ratio)" = readRDS("Study1_OR_mod2.RDS"),
  
  "Other Help-Seeking" = readRDS("Study1_mod3.RDS"),
  "Other Help-Seeking (Odds Ratio)" = readRDS("Study1_OR_mod3.RDS")
)

 


#####################################
### EXTRACT SUMMARY STATISTICS FOR EACH MODEL

## Function to extract these summary statistics
extract_stats <- function(model, outcome_name) {
  draws <- as_draws_df(model)
  
  # Check if the model is a logit model and thus has co-effs on the log-odds scale
  is_or <- grepl("Odds Ratio", outcome_name, ignore.case = TRUE)
  
  # Capture fixed effects (b_) AND random effect standard deviations (sd_)
  target_cols <- grep("^(b_|sd_)", names(draws), value = TRUE)
  
  draws[target_cols] |>
    pivot_longer(everything(), names_to = "parameter", values_to = "draw") |>
    mutate(
      is_sd        = grepl("^sd_", parameter),
      is_intercept = grepl("Intercept", parameter),
      
      # Clean parameter names (e.g., "sd_group__Intercept" -> "SD: group (Intercept)")
      parameter = case_when(
        is_sd ~ paste0(gsub("__", " (", sub("^sd_", "SD: ", parameter)), ")"),
        TRUE  ~ sub("^b_", "", parameter)
      ),
      
      # Convert fixed effects to odds ratios for logit models (leave Intercepts and SDs on link scale)
      draw = if_else(is_or & !is_intercept & !is_sd, exp(draw), draw)
    ) |>
    
    # Calculate posterior summaries
    group_by(parameter, is_sd, is_intercept) |>
    summarise(
      mean_est  = mean(draw),
      ci_lower  = quantile(draw, 0.025),
      ci_upper  = quantile(draw, 0.975),
      pp_pos    = if (is_or) mean(draw > 1) else mean(draw > 0),
      .groups   = "drop"
    ) |>
    mutate(
      parameter = case_when(
        is_sd & grepl("Intercept", parameter) ~ "SD: Between-Chain Baseline Variation",
        is_sd                                 ~ paste0("SD: ", parameter),
        is_intercept                          ~ sub("^Intercept\\[(\\d+)\\]", "Intercept \\1", parameter),
        parameter == "gen2"                   ~ "Chain Position 2",
        parameter == "gen3"                   ~ "Chain Position 3",
        parameter == "gen4"                   ~ "Chain Position 4",
        TRUE                                  ~ parameter
      ),   
      outcome = outcome_name,
      
      # Format table cells:
      # 1. Blank for Intercepts in Odds Ratio models
      # 2. Mean [95% CI] without Pr(>0) for Intercepts and SDs
      # 3. Mean [95% CI] (Pr > 0) for fixed effect predictors
      cell = case_when(
        is_or & is_intercept ~ "",
        is_intercept | is_sd ~ sprintf("%.2f [%.2f, %.2f]", mean_est, ci_lower, ci_upper),
        TRUE                 ~ sprintf("%.2f [%.2f, %.2f] (%.2f)", mean_est, ci_lower, ci_upper, pp_pos)
      )
    ) |>
    select(-is_intercept, -is_sd)
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
    caption      = "Study 1: Bayesian regression results across Outcomes.",
    spanner_cols = 2:(length(outcome_cols) + 1),
    stripe_rows  = seq(2, nrow(reg_table), 2)
  )

ft  




#####################################
### PLOTTING POSTERIOR DISTRIBUTIONS NOW



#####################################
### DEFINE COLOR PALETTE & CUSTOM THEME

## Color legend for different predictors
chain_pos_colors <- c(
  "Chain Position 2" = "#FF1744",  # Vivid Pink/Coral
  "Chain Position 3" = "#00B0FF",  # Bright Cyan
  "Chain Position 4" = "#00E676"   # Neon Green
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

## Function to extract the pertinent posterior draws
extract_full_draws <- function(model, model_name) {
  draws <- as_draws_df(model)
  is_or <- grepl("Odds Ratio", model_name, ignore.case = TRUE)
  
  gen_cols <- grep("^b_gen[234]$", names(draws), value = TRUE)
  
  draws[gen_cols] |>
    pivot_longer(everything(), names_to = "parameter", values_to = "draw") |>
    mutate(
      parameter = case_when(
        parameter == "b_gen2" ~ "Chain Position 2",
        parameter == "b_gen3" ~ "Chain Position 3",
        parameter == "b_gen4" ~ "Chain Position 4"
      ),
      draw    = if (is_or) exp(draw) else draw,
      outcome = sub(" \\(Odds Ratio\\)", "", model_name),
      scale   = if (is_or) "Logit Model (Odds Ratio)" else "Probit Model (Standardized Latent \u03b2)"
    )
}

## Extracting the posteriors from the model now
posterior_draws <- imap_dfr(models, extract_full_draws) |>
  mutate(
    outcome   = factor(outcome, levels = c("Self-Diagnosis", "Other-Diagnosis", "Other Help-Seeking")),
    parameter = factor(parameter, levels = c("Chain Position 2", "Chain Position 3", "Chain Position 4")),
    scale     = factor(scale, levels = c("Probit Model (Standardized Latent \u03b2)", "Logit Model (Odds Ratio)"))
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
    # Exact fixed vertical positions across all panels for perfectly uniform spacing between the stacked credible intervals
    y_pos = case_when(
      parameter == "Chain Position 2" ~ -0.05,
      parameter == "Chain Position 3" ~ -0.10,
      parameter == "Chain Position 4" ~ -0.15
    )
  )

# Null reference lines
ref_lines <- data.frame(
  scale      = factor(c("Probit Model (Standardized Latent \u03b2)", "Logit Model (Odds Ratio)"),
                      levels = c("Probit Model (Standardized Latent \u03b2)", "Logit Model (Odds Ratio)")),
  xintercept = c(0, 1)
)


#####################################
### CONFIGURE X-AXIS LIMITS 

# Define fixed x-axis limits per scale column
probit_x_limits <- c(-1.0, 1.0)  # Fixed scale for Standardized Latent Betas
or_x_limits     <- c(0.0, 7.0)   # Fixed scale for Odds Ratios

dummy_x_limits <- data.frame(
  scale = factor(
    c("Probit Model (Standardized Latent \u03b2)", "Probit Model (Standardized Latent \u03b2)",
      "Logit Model (Odds Ratio)", "Logit Model (Odds Ratio)"),
    levels = c("Probit Model (Standardized Latent \u03b2)", "Logit Model (Odds Ratio)")
  ),
  x = c(probit_x_limits, or_x_limits)
)

# Dummy y-limits to standardize panel vertical bounds (-0.20 baseline space up to 1.05 peak height)
dummy_y_limits <- data.frame(y = c(-0.20, 1.05))



#####################################
### DENSITY OVERLAY PLOT GENERATION

p_density <- ggplot() +
  # Force uniform x-limits per column scale & fixed y-range
  geom_blank(data = dummy_x_limits, aes(x = x)) +
  geom_blank(data = dummy_y_limits, aes(y = y)) +
  
  # Null effect reference line
  geom_vline(
    data = ref_lines,
    aes(xintercept = xintercept),
    linetype = "dashed", color = "grey35", linewidth = 0.8
  ) +
  
  # Normalized density curves (peaks scaled to exactly 1.0 in every panel)
  geom_density(
    data = posterior_draws,
    aes(x = draw, y = after_stat(ndensity), fill = parameter),
    color = NA,
    alpha = 0.40
  ) +
  
  # Solid baseline reference line at y = 0
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.5) +
  
  # Horizontal 95% Credible Intervals at fixed vertical offsets
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
  
  # Facet Grid with left outcome headers
  facet_grid(outcome ~ scale, scales = "free", switch = "y") +
  
  # Custom X-Axis Breaks Per Facet Scale
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
      )
    )
  ) +
  
  # Scale Mappings
  scale_color_manual(values = chain_pos_colors, name = "Chain Position") +
  scale_fill_manual(values = chain_pos_colors, name = "Chain Position") +
  
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