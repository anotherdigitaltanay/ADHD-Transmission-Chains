######################## This script plots the model predicted effects of chain position (both as continuous effects and across likert scale options)
######################## It does this first for the PROBIT models and then for the MONOTONIC models in the mega-analyses, producing one figure for each.
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
library(ggh4x)      
library(marginaleffects)   ## Using this to marginalise predictions over other controls, as conditional effects doesn't do this
library(collapse)

## Set working directory to current study folder
setwd(file.path(renv::project(), "Mega-Analysis")) 

set.seed(23)


#####################################
### STATISTICAL MODEL LOADING

## The name of the main predictor variable in each type of model.
## If either name ever changes in future study datasets, these are the only lines that need updating.
predictor_names <- c(
  "probit"    = "chain_position",
  "monotonic" = "chain_position_trend"
)

## The variables that the plotted predictions are AVERAGED over.
## Rather than showing predictions for one study on one day (the first/reference level of each, which is what
## conditional_effects() does by default), predictions are made for every study/day combination seen in the data and
## then averaged, with each combination counting in proportion to the number of participants in it.
average_over <- c("study_num", "Day")

## The name of the transmission chain identifier (the grouping variable for the random effects)
chain_variable <- "chain_ID"

## Probit models (chain position treated as a categorical predictor)
probit_models <- list(
  "Self-Diagnosis"     = readRDS("MA_mod1.rds"),
  "Other-Diagnosis"    = readRDS("MA_mod2.rds"),
  "Other Help-Seeking" = readRDS("MA_mod3.rds")
)

## Monotonic models (chain position treated as an ordered predictor via mo())
monotonic_models <- list(
  "Self-Diagnosis"     = readRDS("MA_monotonic_mod1.RDS"),
  "Other-Diagnosis"    = readRDS("MA_monotonic_mod2.RDS"),
  "Other Help-Seeking" = readRDS("MA_monotonic_mod3.RDS")
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
get_predictions <- function(model, model_type, categorical = FALSE, interaction_with = NULL,
                            re_formula = if (model_type == "monotonic") NULL else NA, resp = NULL) {
  
  ## Look up which variable name to use, based on the model type (see predictor_names at the top of the script)
  predictor_name <- predictor_names[[model_type]]
  model_data     <- model$data
  
  ## Stop with a helpful message if the predictor isn't in the model (e.g. if the variable name is misspelt)
  if (!predictor_name %in% names(model_data)) {
    stop("'", predictor_name, "' was not found in the ", model_type,
         " model. Check that predictor_names[['", model_type, "']] matches the variable in the model.")
  }
  
  ## Note: by default the monotonic models include group-level (random) effects in the predictions (re_formula = NULL),
  ## whereas the probit models use population-level effects only (re_formula = NA).
  
  ## MARGINALISATION STEPS: STEP 1: List every study/day combination seen in the data, and count how many participants are in each (the weights)
  avg_vars <- intersect(average_over, names(model_data))
  if (length(avg_vars) < length(average_over)) {
    warning("Not averaging over '", paste(setdiff(average_over, avg_vars), collapse = "', '"), "' as it is not a variable in this model.")
  }
  if (length(avg_vars) > 0) {
    combos <- model_data |> count(across(all_of(avg_vars)), name = "n_obs")
  } else {
    combos <- data.frame(n_obs = nrow(model_data))
  }
  
  ## STEP 2: build the prediction grid = every chain position (x every level of the second variable, if given) x every study/day combination
  ## The chain position variable is kept as the same type of factor as in the model (needed for the mo() terms)
  x_var    <- model_data[[predictor_name]]
  x_values <- if (is.factor(x_var)) factor(levels(x_var), levels = levels(x_var), ordered = is.ordered(x_var)) else sort(unique(x_var))
  grid      <- setNames(data.frame(x_values), predictor_name)
  plot_vars <- predictor_name   # the variable(s) that define one plotted point
  
  if (!is.null(interaction_with)) {
    group_values <- setNames(data.frame(sort(unique(model_data[[interaction_with]]))), interaction_with)
    grid         <- tidyr::expand_grid(grid, group_values)
    plot_vars    <- c(predictor_name, interaction_with)
  }
  
  grid <- tidyr::expand_grid(grid, combos) |> as.data.frame()
  
  ## Predictions are made for a new (unobserved) chain. This only matters when group-level effects are included (re_formula = NULL):
  ## for each posterior draw, brms then borrows the chain-level effects of a randomly chosen existing chain, so the
  ## intervals include chain-to-chain variation. This is the same thing conditional_effects() does.
  grid[[chain_variable]] <- "new_chain"
  
  ## The outcome variable(s) are not used to make predictions, but the columns are expected to exist, so give them a placeholder value
  model_formulas <- if (!is.null(model$formula$forms)) model$formula$forms else list(model$formula)
  outcome_vars   <- unique(unlist(lapply(model_formulas, function(f) all.vars(f$formula[[2]]))))
  for (v in setdiff(intersect(outcome_vars, names(model_data)), names(grid))) grid[[v]] <- model_data[[v]][1]
  
  ## Safety net: stop if the model contains a predictor that is neither plotted nor averaged over (so nothing is silently held at one level)
  unhandled_vars <- setdiff(names(model_data), c(names(grid), outcome_vars))
  if (length(unhandled_vars) > 0) {
    stop("The model also contains '", paste(unhandled_vars, collapse = "', '"),
         "'. Add it to average_over at the top of the script so that the predictions are averaged over it too.")
  }
  
  ## STEP 3: let marginaleffects predict every row of the grid and average over the study/day combinations
  ## - by   = the variable(s) that define one plotted point; rows sharing these are averaged together.
  ##          For the likert (ordinal) models marginaleffects returns one prediction per likert option in a column called "group",
  ##          so "group" is added to keep the five options separate.
  ## - wts  = the column holding the weights (number of participants in each study/day combination)
  ## - the remaining arguments are passed straight on to brms 
  is_ordinal <- is.null(resp) && isTRUE(model$family$family %in% c("cumulative", "sratio", "cratio", "acat"))
  by_vars    <- if (is_ordinal) c("group", plot_vars) else plot_vars
  
  if (is.null(resp)) {
    avg <- avg_predictions(model, newdata = grid, by = by_vars, wts = "n_obs", type = "response",
                           re_formula = re_formula, allow_new_levels = TRUE, sample_new_levels = "uncertainty")
  } else {
    avg <- avg_predictions(model, newdata = grid, by = by_vars, wts = "n_obs", type = "response",
                           re_formula = re_formula, allow_new_levels = TRUE, sample_new_levels = "uncertainty", resp = resp)
  }
  
  ## STEP 4: put the result into the column names the plotting functions expect (estimate__, lower__, upper__, cats__)
  if (categorical) {
    if (!is_ordinal) stop("categorical = TRUE can only be used with the likert (ordinal) models.")
    
    ## Probability of choosing each likert option: marginaleffects already reports the posterior median and 95% interval
    avg     <- as.data.frame(avg)
    ce_data <- data.frame(avg[plot_vars], cats__ = avg$group,
                          estimate__ = avg$estimate, lower__ = avg$conf.low, upper__ = avg$conf.high)
    
  } else if (is_ordinal) {
    ## Predicted mean rating = each likert option (1 to 5) multiplied by its predicted probability, added up
    ## (this is how conditional_effects() treats ordinal models when categorical = FALSE).
    ## This is done within each posterior draw, and then summarised as a posterior median and 95% interval.
    draws <- as.data.frame(posterior_draws(avg))   # one row per posterior draw x chain position x likert option ("draw" = the probability)
    draws$likert_value <- as.numeric(as.character(draws$group))
    if (anyNA(draws$likert_value)) stop("The likert options are not labelled as numbers, so a mean rating cannot be computed.")
    
    ce_data <- draws |>
      group_by(drawid, across(all_of(plot_vars))) |>
      summarise(rating = sum(likert_value * draw), .groups = "drop") |>
      group_by(across(all_of(plot_vars))) |>
      summarise(estimate__ = median(rating),
                lower__    = quantile(rating, probs = 0.025),
                upper__    = quantile(rating, probs = 0.975),
                .groups = "drop") |>
      as.data.frame()
    
  } else {
    ## Other models (e.g. cosine similarity): marginaleffects already reports the posterior median and 95% interval
    avg     <- as.data.frame(avg)
    ce_data <- data.frame(avg[plot_vars], estimate__ = avg$estimate, lower__ = avg$conf.low, upper__ = avg$conf.high)
  }
  
  ## Use fixed names ("x" and "group") so the plotting code below never needs to know what the variables are called,
  ## and sort the rows so that the first chain position comes first (the plotting code relies on this for the reference line).
  ce_data$x <- factor(as.character(ce_data[[predictor_name]]), levels = as.character(x_values), ordered = is.ordered(x_var))
  if (!is.null(interaction_with)) {
    ce_data$group <- factor(as.character(ce_data[[interaction_with]]), levels = levels(factor(model$data[[interaction_with]])))
  }
  ce_data <- ce_data[order(ce_data$x), ]
  rownames(ce_data) <- NULL
  
  return(ce_data)
}



## Function 1 for continuous predictions
make_continuous_plot <- function(model, model_type, panel_title, y_label, point_color, show_x_lab = FALSE, show_points = TRUE) {
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
      geom_line(color = point_color, linewidth = 1.1)
    if (show_points) p <- p + geom_point(color = point_color, size = 3.2)
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

## Figure 2: model predictions from the MONOTONIC models
monotonic_figure <- make_combined_figure(
  model_list = monotonic_models,
  model_type = "monotonic",
  filename   = "Continuous_vs_Categorical_Predictions_Monotonic.png"   
)



#####################################
### GENERATE MAIN MANUSCRIPT FIGURE

# Using one of the mediation models as the the mediator equation (chain position -> cosine similarity) is identical in all of them
mediation_model <- readRDS("MA_mediation_mod1.rds")

## The mediation models were fitted on standardised cosine similarity, so predictions come back in SD units.
## To plot on the raw cosine scale, we need the mean and SD used for standardisation (from "Mediation Analyses.R":
## mean(cleaned_dataset$Cosine_Similarity) and sd(cleaned_dataset$Cosine_Similarity)). Update these if the data change.
cosine_mean <- 0.74
cosine_sd   <- 0.133


## Function 4: plot the model-predicted cosine similarity to the DSM at each chain position (monotonic style: line + ribbon + points)
make_drift_plot <- function(model, panel_title, line_color = "#7B3FA0") {
  
  
  ## Pull the mediator-equation predictions (averaged over studies and days, and including chain random effects,
  ## as for the other monotonic panels), then convert from the standardised scale back to raw cosine similarity
  ce_data <- get_predictions(model, model_type = "monotonic", resp = "CosineSimilarityz") |>
    mutate(across(c(estimate__, lower__, upper__), ~ .x * cosine_sd + cosine_mean))

  
  ggplot(ce_data, aes(x = x, y = estimate__, group = 1)) +
    geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = line_color, alpha = 0.20) +
    geom_line(color = line_color, linewidth = 1.1) +
    coord_cartesian(ylim = c(0, 1.02)) +
    labs(
      title = panel_title,
      y     = "Cosine Similarity",
      x     = "Chain Position"
    ) +
    custom_theme +
    theme(plot.margin = margin(t = 8, r = 28, b = 8, l = 12))
}


## Function 5: build the 3-row figure and export it
## Row 1: probit models (point estimates per chain position); Row 2: monotonic models (trend + ribbon); Row 3: semantic drift, centred
make_main_figure <- function(probit_list, monotonic_list, mediation_model, filename, width = 18, height = 18, dpi = 300) {
  
  outcome_names <- names(probit_list)
  n_outcomes    <- length(probit_list)
  
  ## Row 1 (a-c): probit models
  row1 <- lapply(seq_len(n_outcomes), function(i) {
    make_continuous_plot(
      model       = probit_list[[i]],
      model_type  = "probit",
      panel_title = paste0("(", letters[i], ") ", outcome_names[i]),
      y_label     = paste(outcome_names[i], "Rating"),
      point_color = outcome_colors[[outcome_names[i]]],
      show_x_lab  = TRUE
    )
  })
  
  ## Row 2 (d-f): monotonic models
  row2 <- lapply(seq_len(n_outcomes), function(i) {
    make_continuous_plot(
      model       = monotonic_list[[i]],
      model_type  = "monotonic",
      panel_title = paste0("(", letters[n_outcomes + i], ") ", outcome_names[i]),
      y_label     = paste(outcome_names[i], "Rating"),
      point_color = outcome_colors[[outcome_names[i]]],
      show_x_lab  = TRUE,
      show_points = FALSE
    )
  })
  
  ## Row 3 (g): semantic drift, centred with spacers so it aligns with the 3-column grid above
  drift_panel <- make_drift_plot(mediation_model, panel_title = paste0("(", letters[2 * n_outcomes + 1], ") Semantic Drift"))
  row3        <- list(plot_spacer(), drift_panel, plot_spacer())
  
  ## Assemble: each row is wrapped so the three rows stack cleanly with equal column widths
  combined_plot <- wrap_elements(wrap_plots(row1, ncol = 3)) /
    wrap_elements(wrap_plots(row2, ncol = 3)) /
    wrap_elements(wrap_plots(row3, ncol = 3, widths = c(1, 1, 1)))
  
  print(combined_plot)
  ggsave(filename = filename, plot = combined_plot, width = width, height = height, dpi = dpi)
  return(invisible(combined_plot))
}


## Generate and export the figure
main_figure <- make_main_figure(
  probit_list     = probit_models,
  monotonic_list  = monotonic_models,
  mediation_model = mediation_model,
  filename        = "MegaAnalysis_Combined_Plots.png"  
)






