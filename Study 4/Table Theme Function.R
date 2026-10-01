######################## This script builds a custom function for formatting tables for output in the supplemental materials
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

set.seed(43539)



#####################################
### WRITING FUNCTION

apa_theme <- function(ft, caption, stripe_rows = NULL,
                      spanner_cols = NULL, gap_cols = NULL, gap_width = 0.15) {
  ft <- ft |>
    border_remove() |>
    ## Thick rule above the header
    hline_top(part = "header", border = officer::fp_border(width = 1.5))
  
  ## Thin rule under each spanner; pass one range per spanner, e.g. list(4:8, 10:12)
  if (!is.null(spanner_cols)) {
    if (!is.list(spanner_cols)) spanner_cols <- list(spanner_cols)
    for (cols in spanner_cols) {
      ft <- hline(ft, i = 1, j = cols, part = "header",
                  border = officer::fp_border(width = 0.5))
    }
  }
  
  ft <- ft |>
    ## Rule separating header from body
    hline_bottom(part = "header", border = officer::fp_border(width = 0.75)) |>
    ## Thick rule closing the table
    hline_bottom(part = "body", border = officer::fp_border(width = 1.5))
  
  if (!is.null(stripe_rows)) {
    ft <- bg(ft, i = stripe_rows, bg = "#F5F5F5", part = "body")
  }
  
  ft <- ft |>
    set_caption(caption) |>
    font(fontname = "Times New Roman", part = "all") |>
    fontsize(size = 10, part = "all") |>
    padding(padding.top = 3, padding.bottom = 3, part = "all") |>
    autofit()
  
  ## Narrow the gap columns last, so autofit() can't undo it
  if (!is.null(gap_cols)) {
    ft <- ft |>
      padding(j = gap_cols, padding.left = 0, padding.right = 0, part = "all") |>
      width(j = gap_cols, width = gap_width)
  }
  
  ft
}
