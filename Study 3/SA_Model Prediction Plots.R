######################## This script plots the model predicted effects of chain position (both as continuous effects and across likert scale options)
######################## It does this first for the PROBIT models and then for the MONOTONIC models in Study 3 which adjust for various covariates (SENSITIVITY ANALYSES), producing one figure for each. 
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
setwd(file.path(renv::project(), "Study 3")) 

set.seed(4)


#####################################
### STATISTICAL MODEL LOADING

## The name of the main predictor variable in each type of model.
## If either name ever changes in future study datasets, these are the only lines that need updating.
predictor_names <- c(
  "probit"    = "chain_position",
  "monotonic" = "chain_position_trend"
)

## Defining the non-AI usage sample levels, given that the model is a three-way-interaction model
# All figures will be drawn at this level
ai_variable <- "post_AI_introduction"
ai_hold_level <- "0"

sa_hold_at <- setNames(list(ai_hold_level), ai_variable)

## Some variables have been scaled using the scale() which can output matrices in the dataset. To avoid this format to cause
## further troubles down the line, this function changes it to numeric format
fix_matrix_columns <- function(model) {
  for (v in names(model$data)) {
    if (is.matrix(model$data[[v]]) && ncol(model$data[[v]]) == 1) model$data[[v]] <- as.numeric(model$data[[v]])
  }
  model
}



## Probit models (chain position treated as a categorical predictor)
probit_models <- list(
  "Self-Diagnosis"     = fix_matrix_columns(readRDS("SA_Study3_mod1.rds")),
  "Other-Diagnosis"    = fix_matrix_columns(readRDS("SA_Study3_mod2.rds")),
  "Other Help-Seeking" = fix_matrix_columns(readRDS("SA_Study3_mod3.rds"))
)

## Monotonic models (chain position treated as an ordered predictor via mo())
monotonic_models <- list(
  "Self-Diagnosis"     = fix_matrix_columns(readRDS("SA_Study3_monotonic_mod1.RDS")),
  "Other-Diagnosis"    = fix_matrix_columns(readRDS("SA_Study3_monotonic_mod2.RDS")),
  "Other Help-Seeking" = fix_matrix_columns(readRDS("SA_Study3_monotonic_mod3.RDS"))
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
## Sensitivity Analyses CHANGE: new "hold_at" argument: a named list of values for predictors that should be FIXED while chain position
## >>> (and, optionally, interaction_with) varies. Defaults to sa_hold_at (post_AI_introduction = "0"), i.e. every figure
## >>> in this script is drawn for the unpolluted stratum only.
get_predictions <- function(model, model_type, categorical = FALSE, interaction_with = NULL,
                            re_formula = if (model_type == "monotonic") NULL else NA,
                            hold_at = sa_hold_at) {
  
  ## Look up which variable name to ask brms for, based on the model type (see predictor_names at the top of the script)
  predictor_name <- predictor_names[[model_type]]
  
  ## Note: by default the monotonic models include group-level (random) effects in the predictions (re_formula = NULL),
  ## whereas the probit models use brms' default of population-level effects only (re_formula = NA).
  
  ## Specifying the conditions at which the other variable values should be fixed at
  hold_at <- hold_at[setdiff(names(hold_at), c(predictor_name, interaction_with))]
  hold_df <- NULL
  if (length(hold_at) > 0) {
    hold_df <- as.data.frame(hold_at, stringsAsFactors = FALSE)
    for (v in names(hold_df)) {
      if (is.factor(model$data[[v]])) hold_df[[v]] <- factor(hold_df[[v]], levels = levels(model$data[[v]]))
    }
  }
  
  
  
  if (is.null(interaction_with)) {
    ## No second variable: ask brms for the main effect of chain position
    effect_name <- predictor_name
    ce_list <- conditional_effects(model, effects = effect_name, categorical = categorical, re_formula = re_formula,
                                   conditions = hold_df)
    
  } else if (!categorical) {
    ## Second variable, mean ratings: ask brms for the interaction directly (e.g. "chain_position:motivation")
    effect_name <- paste0(predictor_name, ":", interaction_with)
    ce_list <- conditional_effects(model, effects = effect_name, categorical = categorical, re_formula = re_formula,
                                   conditions = hold_df)
    
  } else {
    ## Second variable, likert probabilities: brms doesn't allow interactions when categorical = TRUE,
    ## so instead we ask for the chain position effect separately at each level of the second variable ("conditions")
    effect_name  <- predictor_name
    group_levels <- levels(factor(model$data[[interaction_with]]))
    conds        <- data.frame(group_levels); names(conds) <- interaction_with
    if (!is.null(hold_df)) conds <- cbind(conds, hold_df[rep(1, nrow(conds)), , drop = FALSE])
    rownames(conds) <- group_levels
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
## - split_by / keep_group (optional): if given, predictions are computed separately for each level of "split_by" (e.g. transmission type)
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
  filename   = "SA_Continuous_vs_Categorical_Predictions_Probit.png"
)

## Figure 2: model predictions from the MONOTONIC models
monotonic_figure <- make_combined_figure(
  model_list = monotonic_models,
  model_type = "monotonic",
  filename   = "SA_Continuous_vs_Categorical_Predictions_Monotonic.png"
)


#####################################
### MANUAL CHECKS OF PLOT OUTPUT WITH CONDITIONAL EFFECTS

## Select the model
mod        <- probit_models[["Self-Diagnosis"]]      # swap outcome as needed
gender_ref <- levels(mod$data$gender)[1]  

## See plot
ce <- conditional_effects(
  mod, effects = "chain_position", re_formula = NULL,
  categorical = TRUE,
  conditions = data.frame(post_AI_introduction = factor("0", levels = c("0", "1")))
)

plot(ce)



###########################################################################################################################
###########################################################################################################################
######################## SUPPLEMENT FIGURE: TRANSMISSION TYPE x CHAIN POSITION INTERACTIONS
######################## Top row    = probit interaction models (points + error bars, one colour per transmission type)
######################## Bottom row = monotonic "trend" interaction models (line + shaded ribbon per transmission type)
###########################################################################################################################
###########################################################################################################################

#####################################
### STATISTICAL MODEL LOADING (INTERACTION MODELS)

## Probit interaction models (chain position x transmission type)
interaction_probit_models <- list(
  "Self-Diagnosis"     = fix_matrix_columns(readRDS("SA_Study3_mod4.rds")),
  "Other-Diagnosis"    = fix_matrix_columns(readRDS("SA_Study3_mod5.rds")),
  "Other Help-Seeking" = fix_matrix_columns(readRDS("SA_Study3_mod6.rds"))
)

## Monotonic "trend" interaction models (rate of change in the transmission effect across chain positions)
interaction_trend_models <- list(
  "Self-Diagnosis"     = fix_matrix_columns(readRDS("SA_Study3_monotonic_mod4.rds")),
  "Other-Diagnosis"    = fix_matrix_columns(readRDS("SA_Study3_monotonic_mod5.rds")),
  "Other Help-Seeking" = fix_matrix_columns(readRDS("SA_Study3_monotonic_mod6.rds"))
)


#####################################
### PLOT THEMES (AS PER MAIN MANUSCRIPT)

## The name of the variable that chain position interacts with
interaction_variable <- "motivation"

## How the transmission type levels are stored in the data (left) vs. how they should appear in the figure titles and legends (right)
transmission_labels <- c(
  "virality" = "Virality",
  "accuracy"    = "Accuracy"
)

## Line/point colours for each outcome: a darker and lighter shade of that outcome's colour, one per transmission type
transmission_colors <- list(
  "Self-Diagnosis"     = c("#0072B2", "#56B4E9"),
  "Other-Diagnosis"    = c("#D55E00", "#E69F00"),
  "Other Help-Seeking" = c("#009E73", "#52C4A3")
)

## Ribbon (shaded area) colours used in the trend plots; the same for every outcome
trend_ribbon_colors <- c("#1B3A6B", "#F0A500")

## Theme for the manuscript figure (larger text than the plots above)
manuscript_theme <- theme(
  plot.title         = element_text(family = "Helvetica", size = 24, face = "bold", hjust = 0.5, margin = margin(b = 20)), # Margin argument pushes the plot title away from the graph
  axis.title.x       = element_text(family = "Helvetica", size = 18, margin = margin(t = 15)),               
  axis.title.y       = element_text(family = "Helvetica", size = 18, margin = margin(r = 15)),
  axis.text.x        = element_text(family = "Helvetica", size = 16),
  axis.text.y        = element_text(family = "Helvetica", size = 18),
  
  panel.background   = element_rect(fill = "white", color = "black"),
  panel.grid.major.y = element_line(color = "grey90"),
  panel.grid.major.x = element_blank(),
  panel.grid.minor   = element_blank()
)


#####################################
### FUNCTIONS FOR THE MANUSCRIPT FIGURE

## Function 4: plot the chain position x transmission type interaction for one outcome
## - model_type = "probit"    draws points and error bars per chain position, plus a dotted reference line at chain position 1 for each transmission type
## - model_type = "monotonic" draws a line with a shaded uncertainty ribbon per transmission type
make_interaction_plot <- function(model, model_type, outcome, panel_letter) {
  
  ## Pull the predictions, split by transmission type. re_formula = NA (population-level effects only) for BOTH model types here.
  ## The interaction is a population-level question (does the trend differ between transmission types?), and the random slope
  ## in the monotonic models doesn't vary by transmission type, so adding chain-level variability here would widen the ribbons
  ## without telling us anything extra about the interaction itself.
  plot_data <- get_predictions(model, model_type, categorical = FALSE, interaction_with = interaction_variable, re_formula = NA)
  
  line_colors <- transmission_colors[[outcome]]
  
  ## Start the plot with the parts that are the same for both model types
  p <- ggplot(plot_data, aes(x = x, y = estimate__, group = group, color = group))
  
  ## Add the geoms that differ by model type
  if (model_type == "probit") {
    ## Chain position 1 estimate for each transmission type, used as a visual reference line
    ref_lines <- plot_data |>
      filter(x == levels(factor(x))[1]) |>
      select(group, ref_value = estimate__)
    
    ## Defining a slight dodge so the error bars don't overlap on top of each other
    pd <- position_dodge(width = 0.2)
    
    p <- p +
      geom_hline(data = ref_lines, aes(yintercept = ref_value, color = group),                    ## Adding reference lines
                 linetype = "dotted", linewidth = 1, alpha = 0.5, show.legend = FALSE) +
      geom_errorbar(aes(ymin = lower__, ymax = upper__), width = 0.2, linewidth = 1, position = pd) +   ## Drawing the error bars
      geom_point(size = 4, position = pd) +                                                        ## Drawing the means
      scale_color_manual(values = line_colors, labels = transmission_labels, name = "Transmission Type")
    
  } else {
    p <- p +
      geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = group), color = NA, alpha = 0.2) +    ## Drawing the shaded uncertainty ribbon
      geom_line(linewidth = 1.2) +                                                                 ## Drawing the trend line
      scale_color_manual(values = line_colors,        labels = transmission_labels, name = "Transmission Type") +
      scale_fill_manual( values = trend_ribbon_colors, labels = transmission_labels, name = "Transmission Type")
  }
  
  ## Finish the plot with the parts that are the same for both model types
  p <- p +
    coord_cartesian(ylim = c(1, 3.5)) +                                                            ## Setting the axis limits
    labs(
      y     = paste(outcome, "Rating"),                                                            ## Labelling the axes
      x     = "Chain Position",
      title = paste0("(", panel_letter, ")")
    ) +
    manuscript_theme
  
  return(p)
}


## Function 5: build the 3x2 manuscript figure (probit interaction plots on the top row, trend interaction plots on the bottom row) and export it
make_manuscript_figure <- function(probit_list, trend_list, filename, width = 18, height = 12, dpi = 300) {
  
  outcome_names <- names(probit_list)
  
  ## Top row: one probit interaction plot per outcome, lettered (a), (b), (c)
  top_row <- lapply(seq_along(outcome_names), function(i) {
    make_interaction_plot(probit_list[[i]], "probit", outcome_names[i], panel_letter = letters[i])
  })
  
  ## Bottom row: one trend interaction plot per outcome, lettered (d), (e), (f)
  bottom_row <- lapply(seq_along(outcome_names), function(i) {
    make_interaction_plot(trend_list[[i]], "monotonic", outcome_names[i], panel_letter = letters[i + length(outcome_names)])
  })
  
  ## Arrange the two rows, share the axis titles, and place a small legend inside each panel
  combined_plot <- (wrap_plots(top_row, nrow = 1) / wrap_plots(bottom_row, nrow = 1)) +
    plot_layout(axis_titles = "collect") &
    theme(
      legend.position   = c(0.15, 0.9),                                                                          # Legend position
      legend.text       = element_text(size = 7.5),                                                              # Legend size
      legend.title      = element_text(size = 8, face = "bold"),
      legend.background = element_rect(fill = ggplot2::alpha("white", 0.8), color = "black", linewidth = 0.3)    # Making legend transparent
    )
  
  print(combined_plot)
  
  ## Export the plot
  ggsave(filename = filename, plot = combined_plot, width = width, height = height, dpi = dpi)
  
  return(invisible(combined_plot))
}


#####################################
### GENERATE AND EXPORT THE SUPPLEMENT MANUSCRIPT FIGURE

manuscript_figure <- make_manuscript_figure(
  probit_list = interaction_probit_models,
  trend_list  = interaction_trend_models,
  filename    = "SA_Study_3_Combined_Plots_Main Text.png"
)



###########################################################################################################################
###########################################################################################################################
######################## ADDITIONAL SUPPLEMENTARY FIGURES: LIKERT SCALE CHOICE BREAKDOWN WITHIN EACH TRANSMISSION TYPE
######################## One faceted 6x2 figure. Rows = outcome and model type (probit, then monotonic); columns = transmission type (Example / Narrative Summary).
######################## Each panel shows the predicted probability of choosing each likert option across chain positions.
###########################################################################################################################
###########################################################################################################################


## Function 6: build the full likert breakdown figure as ONE faceted plot and export it
## Rows = outcome (outer strip) and model type (inner strip); columns = transmission type.
## Probit panels are drawn as bars, monotonic panels as lines with ribbons, all inside the same ggplot.
make_interaction_breakdown_figure <- function(probit_list, trend_list, filename, width = 12, height = 22, dpi = 300) {
  
  outcome_names <- names(probit_list)
  
  ## Facet strip styling for this figure (borrowed from the posterior-distribution plots so the supplement looks uniform)
  facet_theme <- theme(
    strip.text       = element_text(family = "Helvetica", size = 13, face = "bold", color = "black"),
    strip.background = element_rect(fill = "grey92", color = "black", linewidth = 1),
    strip.placement  = "outside",
    panel.spacing    = unit(0.8, "lines")
  )
  
  ## Step 1: gather the predictions for every outcome x model type into one long data frame.
  ## Chain position comes back as a plain factor from the probit models and an ordered factor from the monotonic models,
  ## which can't be stacked directly, so it (and the likert column) are converted to plain text first.
  to_text <- function(d) d |> mutate(across(any_of(c("x", "group", "cats__", "response__")), as.character))
  
  breakdown_data <- bind_rows(lapply(outcome_names, function(outcome) {
    bind_rows(
      get_predictions(probit_list[[outcome]], "probit",    categorical = TRUE, interaction_with = interaction_variable) |>
        to_text() |> mutate(model_type = "Probit"),
      get_predictions(trend_list[[outcome]],  "monotonic", categorical = TRUE, interaction_with = interaction_variable) |>
        to_text() |> mutate(model_type = "Monotonic")
    ) |>
      mutate(outcome = outcome)
  }))
  
  cat_col <- if ("cats__" %in% names(breakdown_data)) "cats__" else "response__"   # brms versions name this column differently
  
  ## Step 2: tidy up the labelling so the facets and legend read nicely
  likert_levels <- sort(unique(breakdown_data[[cat_col]]))
  panel_colors  <- setNames(unname(likert_colors)[seq_along(likert_levels)], likert_levels)
  
  breakdown_data <- breakdown_data |>
    mutate(
      x            = factor(x, levels = sort(unique(x))),                                        # chain position, same axis type for both model types
      Likert       = factor(.data[[cat_col]], levels = likert_levels),
      outcome      = factor(outcome, levels = outcome_names),
      model_type   = factor(model_type, levels = c("Probit", "Monotonic")),
      transmission = factor(transmission_labels[group], levels = unname(transmission_labels))
    )
  
  probit_data    <- breakdown_data |> filter(model_type == "Probit")
  monotonic_data <- breakdown_data |> filter(model_type == "Monotonic")
  
  ## Step 3: draw everything in one plot
  combined_plot <- ggplot() +
    
    ## Probit panels: side-by-side bars, one per likert option, at each chain position
    geom_col(data = probit_data, aes(x = x, y = estimate__, fill = Likert),
             position = position_dodge(width = 0.78), width = 0.70, color = "grey30", linewidth = 0.3) +
    
    ## Monotonic panels: one line (with ribbon and points) per likert option across chain positions
    geom_ribbon(data = monotonic_data, aes(x = x, ymin = lower__, ymax = upper__, fill = Likert, group = Likert),
                alpha = 0.15, color = NA) +
    geom_line(data = monotonic_data, aes(x = x, y = estimate__, color = Likert, group = Likert), linewidth = 1.0) +
    geom_point(data = monotonic_data, aes(x = x, y = estimate__, color = Likert), size = 2.5) +
    
    ## Facet grid: outcome (outer) and model type (inner) down the rows, transmission type across the columns.
    ## facet_nested() (from ggh4x) merges the outcome strip across its two rows; if ggh4x isn't available we fall back to
    ## ggplot2's facet_grid(), which draws the same layout but repeats the outcome strip on each row.
    (if (requireNamespace("ggh4x", quietly = TRUE)) ggh4x::facet_nested(outcome + model_type ~ transmission, switch = "y")
     else facet_grid(outcome + model_type ~ transmission, switch = "y")) +
    
    scale_y_continuous(labels = scales::percent, limits = c(0, 1.0), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual( values = panel_colors, name = "Likert Rating Option") +
    scale_color_manual(values = panel_colors, name = "Likert Rating Option") +
    labs(x = "Chain Position", y = "Predicted Probability") +
    custom_theme +
    facet_theme
  
  print(combined_plot)
  
  ## Export the plot
  ggsave(filename = filename, plot = combined_plot, width = width, height = height, dpi = dpi)
  
  return(invisible(combined_plot))
}


#####################################
### GENERATE AND EXPORT THE SUPPLEMENTARY FIGURE

breakdown_figure <- make_interaction_breakdown_figure(
  probit_list = interaction_probit_models,
  trend_list  = interaction_trend_models,
  filename    = "SA_Study_3_Likert_Breakdown_by_Transmission.png"
)



####################################
### QUICK MANUAL CHECK OF THE FIGURES ABOVE (optional)

####################################
### QUICK MANUAL CHECK OF THE INTERACTION LIKERT SCALE FIGURES ABOVE (optional)

## This asks brms to draw its own default plot of the same predictions, with none of the custom plotting code involved,
## so the supplementary figure can be checked against it by eye.
## Swap the model in the first line for whichever one you want to inspect (e.g. interaction_trend_models[["Other-Diagnosis"]]),
## and swap "chain_position" for "chain_position_trend" when checking a monotonic model.
## brms draws one panel per transmission type, with one line per likert option across chain positions.

check_model <- interaction_trend_models[["Self-Diagnosis"]]

check_ce <- conditional_effects(check_model, effects = "chain_position_trend", categorical = TRUE,
                                conditions = data.frame(motivation = c("virality", "accuracy"),
                                                        post_AI_introduction = factor("0", levels = c("0", "1"))))

plot(check_ce, ask = FALSE)

