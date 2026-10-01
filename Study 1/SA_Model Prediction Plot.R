######################## This script plots the model predicted effects from Study 1 (both as continuous effects and across likert scale options)
########################################################################################################################### 
###########################################################################################################################  

#####################################
### PACKAGE LOADING 

## Load the relevant libraries. 
## If these packages aren't installed on your machine, type in the following line "install.packages("INSERT PACKAGE NAME ONE AT A TIME")"
library(dplyr)
library(ggplot2)
library(brms)
library(patchwork)
library(tidyverse)
library(scales)

## Set working directory to current study folder
setwd(file.path(renv::project(), "Study 1")) 

set.seed(344)


#####################################
### STATISTICAL MODEL LOADING

mod1 <- readRDS("Study1_monotonic_mod1.rds")  
mod2 <- readRDS("Study1_monotonic_mod2.rds")  
mod3 <- readRDS("Study1_monotonic_mod3.rds") 

models <- list(
  "Self-Diagnosis"     = mod1,
  "Other-Diagnosis"    = mod2,
  "Other Help-Seeking" = mod3
)


#####################################
### PLOT THEMES

## Define colors for each outcome
outcome_colors <- c(
  "Self-Diagnosis"     = "#0072B2",  # Blue
  "Other-Diagnosis"    = "#D55E00",  # Vermilion
  "Other Help-Seeking" = "#009E73"   # Green
)

## Define color scheme for each likert scale option (1 = Strongly Disagree that one has ADHD; 5 = Strongly Agree that one has ADHD)
likert_colors <- c(
  "1" = "#FFCDD2",  # Soft Rose 
  "2" = "#EF9A9A",  # Light Coral
  "3" = "#E57373",  # Medium Red
  "4" = "#E53935",  # Bright Crimson
  "5" = "#B71C1C"   # Deep Burgundy / Maroon (Highest Intensity)
)

## Defining custom plot theme 
custom_theme <- theme(
  plot.title         = element_text(family = "Helvetica", size = 15, face = "bold", hjust = 0.5, margin = margin(b = 10)),
  axis.title.x       = element_text(family = "Helvetica", size = 13, face = "bold", margin = margin(t = 8)),
  axis.title.y       = element_text(family = "Helvetica", size = 13, face = "bold", margin = margin(r = 8)),
  axis.text.x        = element_text(family = "Helvetica", size = 11, color = "black"),
  axis.text.y        = element_text(family = "Helvetica", size = 11, color = "black"),
  
  legend.position    = "bottom",
  legend.title       = element_text(family = "Helvetica", size = 12, face = "bold"),
  legend.text        = element_text(family = "Helvetica", size = 11),
  
  panel.background   = element_rect(fill = "white", color = "black", linewidth = 0.8),
  panel.grid.major.y = element_line(color = "grey92"),
  panel.grid.major.x = element_blank(),
  panel.grid.minor   = element_blank()
)


#####################################
### FUNCTIONS FOR EXTRACTING PLOT DATA FROM MODEL PREDICTIONS

## Function 1 for continuous predictions
make_continuous_plot <- function(model, panel_title, y_label, point_color, show_x_lab = FALSE) {
  ce_data <- conditional_effects(model, effects = "gen_trend", categorical = FALSE, re_formula = NULL)$gen_trend
  ref_val <- ce_data$estimate__[1]
  
  p <- ggplot(ce_data, aes(x = gen_trend, y = estimate__, group = 1)) +
    geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = point_color, alpha = 0.20) +
    geom_line(color = point_color, linewidth = 1.1) +
    geom_point(color = point_color, size = 3.2) +
    
    coord_cartesian(ylim = c(1, 3.5)) +
    labs(
      title = panel_title,
      y     = y_label,
      x     = if (show_x_lab) "Chain Position" else NULL
    ) +
    custom_theme +
    theme(plot.margin = margin(t = 8, r = 28, b = 8, l = 12))
  
  return(p)
}


## Function 2 for decomposing effects as probabilities of choosing each likert scale option
make_categorical_plot <- function(model, panel_title, show_x_lab = FALSE) {
  ce_cat <- conditional_effects(model, effects = "gen_trend", categorical = TRUE, re_formula = NULL)$gen_trend
  cat_col <- if ("cats__" %in% names(ce_cat)) "cats__" else "response__"
  
  plot_data <- ce_cat |>
    mutate(Likert = factor(.data[[cat_col]], levels = c("1", "2", "3", "4", "5")))
  
  p <- ggplot(plot_data, aes(x = gen_trend, y = estimate__, group = Likert, color = Likert, fill = Likert)) +
   
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.5) +
    
    scale_y_continuous(labels = scales::percent, limits = c(0, 1.0), expand = expansion(mult = c(0, 0.02))) +
    scale_color_manual(values = likert_colors, name = "Likert Rating Option") +
    scale_fill_manual(values = likert_colors, name = "Likert Rating Option") +
    labs(
      title = panel_title,
      y     = "Predicted Probability",
      x     = if (show_x_lab) "Chain Position" else NULL
    ) +
    custom_theme +
    theme(plot.margin = margin(t = 8, r = 28, b = 8, l = 12))
  
  return(p)
}


#####################################
### GENERATE PLOTS FOR EACH MODEL

# Outcome 1: Self-Diagnosis
p1_cont <- make_continuous_plot(mod1, "(a) Self-Diagnosis", "Self-Diagnosis Rating", outcome_colors["Self-Diagnosis"], show_x_lab = FALSE)
p1_cat  <- make_categorical_plot(mod1, "(b) Likert Scale Choice Breakdown (Self-Diagnosis)", show_x_lab = FALSE)

# Outcome 2: Other-Diagnosis
p2_cont <- make_continuous_plot(mod2, "(c) Other-Diagnosis", "Other-Diagnosis Rating", outcome_colors["Other-Diagnosis"], show_x_lab = FALSE)
p2_cat  <- make_categorical_plot(mod2, "(d) Likert Scale Choice Breakdown (Other-Diagnosis)", show_x_lab = FALSE)

# Outcome 3: Other Help-Seeking
p3_cont <- make_continuous_plot(mod3, "(e) Other Help-Seeking", "Other Help-Seeking Rating", outcome_colors["Other Help-Seeking"], show_x_lab = TRUE)
p3_cat  <- make_categorical_plot(mod3, "(f) Likert Scale Choice Breakdown (Other Help-Seeking)", show_x_lab = TRUE)


#####################################
### PREPARE FINAL PLOT FOR EXPORT

## Combine all the plots
combined_plot_3x2 <- (p1_cont | p1_cat) /
  (p2_cont | p2_cat) /
  (p3_cont | p3_cat) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(combined_plot_3x2)

## Export the plot
ggsave(
  filename = "Continuous_vs_Categorical_SA.png",
  plot     = combined_plot_3x2,
  width    = 12,
  height   = 13,
  dpi      = 300
)



