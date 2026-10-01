######################## This script plots the model predicted effects of chain position (both as continuous effects and across likert scale options)
######################## It does this for the PROBIT models Study 4, producing one figure for each.
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
library(ggh4x)      # for facet_nested() in the supplementary breakdown figure

## Set working directory to current study folder
setwd(file.path(renv::project(), "Study 4")) 

set.seed(4)


#####################################
### STATISTICAL MODEL LOADING

## The name of the main predictor variable in each type of model.
## If either name ever changes in future study datasets, these are the only lines that need updating.
predictor_names <- c(
  "probit"    = "chain_position",
  "monotonic" = "chain_position_trend"
)

## Probit models (chain position treated as a categorical predictor)
probit_models <- list(
  "Self-Diagnosis"     = readRDS("Study4_mod1.rds"),
  "Other-Diagnosis"    = readRDS("Study4_mod2.rds"),
  "Other Help-Seeking" = readRDS("Study4_mod3.rds")
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

## Every function below takes a "model_type" argument, which can be either "probit" or "monotonic".
## This single argument controls:
##   (1) how the predictions are pulled out of the model, and
##   (2) how they are drawn (probit = separate points/bars per chain position; monotonic = connected lines with ribbons).


## Function 0 (helper): pull the model predictions for the predictor out of a brms model
## - categorical = FALSE gives the predicted mean rating at each chain position
## - categorical = TRUE  gives the predicted probability of choosing each likert option at each chain position
## - interaction_with (optional) is the name of a second variable; if given, predictions are returned separately for each level of it
## - re_formula controls whether group-level (random) effects are included; by default this follows the model type (see below)
get_predictions <- function(model, model_type, categorical = FALSE, interaction_with = NULL,
                            re_formula = if (model_type == "monotonic") NULL else NA) {
  
  ## Look up which variable name to ask brms for, based on the model type (see predictor_names at the top of the script)
  predictor_name <- predictor_names[[model_type]]
  
  ## Note: by default the monotonic models include group-level (random) effects in the predictions (re_formula = NULL),
  ## whereas the probit models use brms' default of population-level effects only (re_formula = NA).
  
  if (is.null(interaction_with)) {
    ## No second variable: ask brms for the main effect of chain position
    effect_name <- predictor_name
    ce_list <- conditional_effects(model, effects = effect_name, categorical = categorical, re_formula = re_formula)
    
  } else if (!categorical) {
    ## Second variable, mean ratings: ask brms for the interaction directly (e.g. "chain_position:motivation")
    effect_name <- paste0(predictor_name, ":", interaction_with)
    ce_list <- conditional_effects(model, effects = effect_name, categorical = categorical, re_formula = re_formula)
    
  } else {
    ## Second variable, likert probabilities: brms doesn't allow interactions when categorical = TRUE,
    ## so instead we ask for the chain position effect separately at each level of the second variable ("conditions")
    effect_name  <- predictor_name
    group_levels <- levels(factor(model$data[[interaction_with]]))
    conds        <- data.frame(group_levels); names(conds) <- interaction_with; rownames(conds) <- group_levels
    ce_list <- conditional_effects(model, effects = effect_name, categorical = categorical, re_formula = re_formula, conditions = conds)
  }
  
  ## Stop with a helpful message if brms returned nothing for this predictor (e.g. if the variable name is misspelt)
  if (length(ce_list) == 0) {
    stop("conditional_effects() returned no effects for '", effect_name, "' in the ", model_type,
         " model. Check that predictor_names[['", model_type, "']] matches the variable in the model.")
  }
  
  ## We only asked for one effect, so take the first (and only) element of the list.
  ## This avoids problems where brms labels the effect slightly differently to the variable name (e.g. for mo() terms).
  ce_data <- ce_list[[1]]
  
  ## brms always stores the x-axis variable in a column called "effect1__".
  ## When a second variable is involved, brms also carries that variable through as its own column (e.g. "motivation"),
  ## whichever route we took above, so we read the grouping from there.
  ## Rename these to fixed names ("x" and "group") so the plotting code below never needs to know what the variables are called.
  ce_data <- ce_data |> rename(x = effect1__)
  if (!is.null(interaction_with)) {
    ce_data$group <- factor(as.character(ce_data[[interaction_with]]), levels = levels(factor(model$data[[interaction_with]])))
  }
  
  return(ce_data)
}


## Function 1 for continuous predictions
make_continuous_plot <- function(model, model_type, panel_title, y_label, point_color, show_x_lab = FALSE) {
  ce_data <- get_predictions(model, model_type, categorical = FALSE)
  ref_val <- ce_data$estimate__[1]   # predicted rating at the first chain position (used as the dashed reference line)
  
  ## Start the plot with the parts that are the same for both model types
  p <- ggplot(ce_data, aes(x = x, y = estimate__, group = 1))
  
  ## Add the geoms that differ by model type
  if (model_type == "probit") {
    ## Probit: dashed reference line, error bars, and a point per chain position
    p <- p +
      geom_hline(yintercept = ref_val, color = point_color, linetype = "dashed", linewidth = 0.9, alpha = 0.5) +
      geom_errorbar(aes(ymin = lower__, ymax = upper__), width = 0.18, color = point_color, linewidth = 0.9) +
      geom_point(color = point_color, size = 3.8)
  } else {
    ## Monotonic: shaded uncertainty ribbon, a connecting line, and a point per chain position
    p <- p +
      geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = point_color, alpha = 0.20) +
      geom_line(color = point_color, linewidth = 1.1) +
      geom_point(color = point_color, size = 3.2)
  }
  
  ## Finish the plot with the parts that are the same for both model types
  p <- p +
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
## - split_by / keep_group (optional): if given, predictions are computed separately for each level of "split_by" (e.g. Motivation)
##   and only the level named in "keep_group" is drawn. This is used for the supplementary breakdown figure at the end of the script.
make_categorical_plot <- function(model, model_type, panel_title, show_x_lab = FALSE, split_by = NULL, keep_group = NULL) {
  ce_cat  <- get_predictions(model, model_type, categorical = TRUE, interaction_with = split_by)
  if (!is.null(keep_group)) ce_cat <- ce_cat |> filter(group == keep_group)
  
  ## Stop with a helpful message if there is nothing to plot (e.g. keep_group doesn't match a level in the data)
  if (nrow(ce_cat) == 0) {
    stop("No predictions left to plot for '", panel_title, "'. Check that keep_group ('", keep_group,
         "') matches a level of '", split_by, "' in the model data: ", paste(levels(factor(model$data[[split_by]])), collapse = ", "))
  }
  
  cat_col <- if ("cats__" %in% names(ce_cat)) "cats__" else "response__"   # brms versions name this column differently
  
  ## The likert options as brms labels them (normally "1" to "5"), in order. Colours are matched to these by position,
  ## so the plot still works even if brms labels the categories slightly differently.
  likert_levels <- levels(factor(ce_cat[[cat_col]]))
  panel_colors  <- setNames(unname(likert_colors)[seq_along(likert_levels)], likert_levels)
  
  plot_data <- ce_cat |>
    mutate(Likert = factor(.data[[cat_col]], levels = likert_levels))
  
  ## Add the geoms that differ by model type
  if (model_type == "probit") {
    ## Probit: side-by-side bars, one per likert option, at each chain position
    p <- ggplot(plot_data, aes(x = x, y = estimate__, fill = Likert)) +
      geom_col(position = position_dodge(width = 0.78), width = 0.70, color = "grey30", linewidth = 0.3) +
      scale_fill_manual(values = panel_colors, name = "Likert Rating Option")
  } else {
    ## Monotonic: one line (with ribbon and points) per likert option across chain positions
    p <- ggplot(plot_data, aes(x = x, y = estimate__, group = Likert, color = Likert, fill = Likert)) +
      geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.15, color = NA) +
      geom_line(linewidth = 1.0) +
      geom_point(size = 2.5) +
      scale_color_manual(values = panel_colors, name = "Likert Rating Option") +
      scale_fill_manual(values = panel_colors, name = "Likert Rating Option")
  }
  
  ## Finish the plot with the parts that are the same for both model types
  p <- p +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1.0), expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = panel_title,
      y     = "Predicted Probability",
      x     = if (show_x_lab) "Chain Position" else NULL
    ) +
    custom_theme +
    theme(plot.margin = margin(t = 8, r = 28, b = 8, l = 12))
  
  return(p)
}


## Function 3: build the full 3x2 figure (one row per outcome; continuous plot on the left, likert breakdown on the right) and export it
make_combined_figure <- function(model_list, model_type, filename, width = 12, height = 13, dpi = 300) {
  
  ## Panel letters run (a), (b) down the rows: (a)/(b) for outcome 1, (c)/(d) for outcome 2, (e)/(f) for outcome 3
  panel_letters <- letters[seq_len(2 * length(model_list))]
  outcome_names <- names(model_list)
  n_outcomes    <- length(model_list)
  
  ## Build the pair of plots for each outcome
  rows <- lapply(seq_len(n_outcomes), function(i) {
    outcome    <- outcome_names[i]
    is_last    <- (i == n_outcomes)                      # only the bottom row shows the x-axis label
    letter_con <- panel_letters[2 * i - 1]               # letter for the continuous panel
    letter_cat <- panel_letters[2 * i]                   # letter for the likert breakdown panel
    
    p_cont <- make_continuous_plot(
      model       = model_list[[outcome]],
      model_type  = model_type,
      panel_title = paste0("(", letter_con, ") ", outcome),
      y_label     = paste(outcome, "Rating"),
      point_color = outcome_colors[[outcome]],
      show_x_lab  = is_last
    )
    
    p_cat <- make_categorical_plot(
      model       = model_list[[outcome]],
      model_type  = model_type,
      panel_title = paste0("(", letter_cat, ") Likert Scale Choice Breakdown (", outcome, ")"),
      show_x_lab  = is_last
    )
    
    p_cont | p_cat   # place the two panels side by side
  })
  
  ## Stack the rows on top of each other and share one legend at the bottom
  combined_plot <- wrap_plots(rows, ncol = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  print(combined_plot)
  
  ## Export the plot
  ggsave(filename = filename, plot = combined_plot, width = width, height = height, dpi = dpi)
  
  return(invisible(combined_plot))
}


#####################################
### GENERATE AND EXPORT THE FIGURES

## Figure 1: model predictions from the PROBIT models
probit_figure <- make_combined_figure(
  model_list = probit_models,
  model_type = "probit",
  filename   = "Continuous_vs_Categorical_Predictions_Probit.png"
)






