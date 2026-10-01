###### This script generates the agreement metrics between the two coders who screened participant ADHD explanations for various language errors in Study 2
##########################################################################################################################################
##########################################################################################################################################


#####################################
### PACKAGE LOADING

## Load the necessary packages
library(readxl)
library(dplyr)
library(flextable) 

## Set working directory to current study folder
setwd(file.path(renv::project(), "Explanation Screening Agreement", "Study 2"))

## Setting seed
set.seed(435)


#####################################
### HOW AGREEMENT BETWEEN CODERS IS MEASURED HERE

### Each participant ADHD explanation was screened independently by two coders, who marked a 1
### in a column whenever a given type of language error was present, and left the cell
### empty otherwise. Agreement between coders can be calculated at two levels:

##  Metric 1 -- did the two coders agree that the participant's ADHD explanation contained ANY error at all?
##              This collapses the nine error columns into a single yes/no judgement per
##              coder, so two coders who both saw a problem still count as agreeing here
##              even if they disagreed about which error it was.

##  Metric 2 -- did the two coders agree on THE SAME TYPES of error? This requires all
##              nine columns to match, so it is necessarily the stricter of the two and
##              will always be the lower number.

### The gap between metric 1 and metric 2 is the quantity of interest: it is the share of
### responses where both coders saw a problem but described it differently. The usage
### table produced further down breaks that gap apart by error category, showing which
### categories the coders applied at similar rates and which they did not.


#####################################
### DATA CLEANING FUNCTION

### readxl assigns each column a type based on its contents, so an error column reads as
### numeric only if every cell in it is numeric or empty. A single stray character in one
### cell types the entire column as text, and the same column may then be numeric in one
### workbook and text in another. This helper function normalises an error column to 0/1 regardless
### of how it arrived:

##   - a blank cell becomes 0, since the coders left cells empty to mean "error absent"
##   - "1" and 1 both become 1
##   - anything else halts the script, naming the file, column and offending value

### That last behaviour is deliberate. The alternative -- coercing unrecognised text to NA
### and then to 0 -- would record a cell the coders did fill in as an absent error, which
### is a silent miscount rather than a visible failure.

to_binary <- function(x, colname, file) {
  v <- trimws(as.character(x))                        
  v[v == ""] <- NA_character_
  
  unexpected <- !is.na(v) & !(v %in% c("0", "0.0", "1", "1.0"))                    ## Flagging if any other value apart from 1 or 0 is present           
  if (any(unexpected)) {
    stop(sprintf("Unexpected value(s) in column '%s' of %s: %s",
                 colname, basename(file),
                 paste(sprintf("'%s'", unique(v[unexpected])), collapse = ", ")),
         call. = FALSE)
  }
  
  out <- rep(0L, length(v))
  out[v %in% c("1", "1.0")] <- 1L
  out
}


#####################################
### THE MAIN FUNCTION

### Takes one chain position's workbook and returns three data frames:

##   $agreement  -- one row: metrics 1 and 2, with the denominators they rest on
##   $usage      -- one row per error category: how often each coder applied it, and how
##                  often they agreed about it
##   $exclusions -- one row per response excluded from the analysis, with the reason

### Arguments:
##   file          path to the workbook for this chain position
##   position      the chain position number, carried through into the output tables
##   codes         positions of the nine error columns (columns C to K)
##   response_col  accepted spellings of the response column header;
##   sheets        names of the two coders' sheets within the workbook
##   on_one_sided  what to do when a response is blank on one sheet but not the other

screening_agreement <- function(file,
                                position,
                                codes        = 3:11,
                                response_col = c("Raw Response"),
                                sheets       = c("Coder1", "Coder2"),
                                on_one_sided = c("stop", "exclude")) {
  
  on_one_sided <- match.arg(on_one_sided)
  
  ## --- Load the sheets from the two coders -----------------------------------------------------------------
  ## One sheet per coder, both from the same workbook
  d1 <- read_excel(file, sheet = sheets[1])
  d2 <- read_excel(file, sheet = sheets[2])
  
  ## --- Locate the response column -------------------------------------------
  ## Match the header against the column name that contains the participant explanations
  ## This is checked explicitly because of how R behaves if the column is missing: d1[["x"]]
  ## returns NULL, NULL == NULL evaluates to logical(0), and all(logical(0)) is TRUE. The
  ## ordering check further down would therefore report success having compared nothing.
  resp1 <- intersect(response_col, names(d1))
  resp2 <- intersect(response_col, names(d2))
  if (length(resp1) != 1 || length(resp2) != 1) {
    stop("Could not identify exactly one response column in ", basename(file),
         ". Columns found: ", paste(names(d1), collapse = ", "), call. = FALSE)
  }
  
  ## --- Check the two sheets are comparable ----------------------------------
  ## Every metric below compares row i of one sheet with row i of the other, so the sheets
  ## must have the same number of rows before anything else is attempted.
  if (nrow(d1) != nrow(d2)) {
    stop("Sheets have different row counts in ", basename(file),
         " (", nrow(d1), " vs ", nrow(d2), ")", call. = FALSE)
  }
  
  n_in_sheet <- nrow(d1)
  
  ## --- Identify rows with no response to code -------------------------------
  ## Some rows might be empty because the response was collected late during data collection. These
  ## might pop up, but they cannot simply be left in place: a blank compares as NA rather
  ## than TRUE, and two blanks would otherwise count as a matching row and inflate every
  ## percentage below. Each excluded row is logged so the exclusions can be reported.
  blank1 <- is.na(d1[[resp1]]) | trimws(as.character(d1[[resp1]])) == ""
  blank2 <- is.na(d2[[resp2]]) | trimws(as.character(d2[[resp2]])) == ""
  
  ## Blank on both sheets: no response exists, so there was nothing for either coder to do
  drop_row <- blank1 & blank2
  reason   <- ifelse(drop_row, "blank on both sheets", NA_character_)
  
  ## Blank on one sheet only: the two sheets disagree about whether a response exists at
  ## all, which is a different situation and is treated as an error unless told otherwise
  one_sided <- xor(blank1, blank2)
  if (any(one_sided)) {
    if (on_one_sided == "stop") {
      stop("Response blank on one sheet but not the other in ", basename(file),
           " at Excel row(s): ", paste(which(one_sided) + 1, collapse = ", "),
           "\n  If these are replaced responses, re-run with on_one_sided = \"exclude\".",
           call. = FALSE)
    }
    drop_row <- drop_row | one_sided
    reason[one_sided] <- ifelse(blank1[one_sided],
                                "blank on Coder1 sheet only",
                                "blank on Coder2 sheet only")
  }
  
  ## --- Log the exclusions ----------------------------------------------------
  ## Excel row numbers (index + 1 for the header) so rows can be found in the workbook.
  ## Where one sheet still holds the response text, keep it for identification.
  exclusions <- if (any(drop_row)) {
    idx  <- which(drop_row)
    text <- ifelse(!blank1[idx], as.character(d1[[resp1]])[idx],
                   ifelse(!blank2[idx], as.character(d2[[resp2]])[idx], NA_character_))
    data.frame(chain_position = position,
               file           = basename(file),
               excel_row      = idx + 1,
               reason         = reason[idx],
               response_text  = substr(text, 1, 80),
               stringsAsFactors = FALSE)
  } else {
    
    ## An empty frame with the same columns, so binding across files still works
    data.frame(chain_position = integer(), file = character(), excel_row = integer(),
               reason = character(), response_text = character(),
               stringsAsFactors = FALSE)
  }
  
  if (any(drop_row)) {
    message(basename(file), ": excluding ", sum(drop_row), " of ", n_in_sheet,
            " rows with no response to code (Excel row(s): ",
            paste(which(drop_row) + 1, collapse = ", "), ")")
    d1 <- d1[!drop_row, ]
    d2 <- d2[!drop_row, ]
  }
  
  if (nrow(d1) == 0) {
    stop("No analysable responses left in ", basename(file), call. = FALSE)
  }
  
  ## --- Normalise the error columns -------------------------------------------
  ## Done after the exclusions so dropped rows cannot contribute codes to the totals
  d1[codes] <- Map(to_binary, d1[codes], names(d1)[codes], file)
  d2[codes] <- Map(to_binary, d2[codes], names(d2)[codes], file)
  
  ## --- Confirm the rows correspond -------------------------------------------
  ## Row i of each sheet must be the same participant response. If either sheet were ever
  ## sorted or filtered independently, the comparisons below would pair the wrong rows and
  ## produce plausible-looking but meaningless numbers.
  if (!all(d1[[resp1]] == d2[[resp2]])) {
    stop("Responses are not in the same order in ", basename(file), call. = FALSE)
  }
  
  ## --- Compute the metrics ----------------------------------------------------
  ## Matrices of 0s and 1s, one row per response and one column per error category
  c1 <- as.matrix(d1[, codes])
  c2 <- as.matrix(d2[, codes])
  
  ## Metric 1: collapses the nine columns to a single "did this need editing" judgement,
  ## then ask how often the two coders reached the same judgement
  any1 <- rowSums(c1) > 0
  any2 <- rowSums(c2) > 0
  
  ## Metric 2: a row counts as agreement only if all nine columns match
  same_errors <- rowSums(c1 != c2) == 0
  
  ## Per-category agreement, used in the usage table to explain the gap between the two
  by_code <- colMeans(c1 == c2) * 100
  
  list(
    ## Headline metrics, with the denominators visible alongside them
    agreement = data.frame(
      chain_position      = position,
      n_rows_in_sheet     = n_in_sheet,
      n_excluded          = sum(drop_row),
      n_responses         = nrow(c1),
      metric1_any_error   = mean(any1 == any2) * 100,
      metric2_same_errors = mean(same_errors) * 100
    ),
    
    ## Per-category detail. pct_coder1 and pct_coder2 are the share of responses each
    ## coder marked with that error; 
    usage = data.frame(
      chain_position = position,
      error_type     = colnames(c1),
      pct_agreement  = by_code,
      pct_coder1     = colMeans(c1) * 100,
      pct_coder2     = colMeans(c2) * 100,
      row.names      = NULL
    ),
    ## Every row excluded from the metrics above, and why
    exclusions = exclusions
  )
}


#####################################
### RUNNING THIS FUNCTION ACROSS ALL CHAIN POSITIONS

## The seven files share a naming convention and an identical structure
files <- sprintf("ChainPosition%d_Responses.xlsx", 1:7)

## Fail immediately on a missing or misnamed file, rather than part-way through the loop
stopifnot(all(file.exists(files)))

## Apply the function to each file, pairing every file with its chain position
results <- Map(screening_agreement, files, seq_along(files))

## Metrics 1 and 2, one row per chain position
agreement_table <- bind_rows(lapply(results, `[[`, "agreement")) |>
  mutate(across(where(is.numeric), \(x) round(x, 1)))

## Per-category detail, one row per chain position per error type
usage_table <- bind_rows(lapply(results, `[[`, "usage")) |>
  mutate(across(where(is.numeric), \(x) round(x, 1)))

## Every response excluded from the analysis, and the reason for it
exclusions_table <- bind_rows(lapply(results, `[[`, "exclusions"))



#####################################
### GENERATING TABLES FOR EXPORTING TO MANUSCRIPT

## Cleaning Column Names for Agreement Table
agreement_table <- agreement_table %>% 
  select(-c(n_excluded, n_rows_in_sheet)) %>%             # Removing this here as we didn't find any exclusions and this will only confuse a lay-reader
  rename(
    "Chain Position" = chain_position,
    "# of Participant Explanations" = n_responses,
    "Any Error Present" = metric1_any_error,
    "Exact Same Types of Error" = metric2_same_errors
  )

## Defining a table theme that can be used across both tables
apa_theme <- function(ft, caption, stripe_rows, spanner_cols) {
  ft |>
    border_remove() |>
    
    ## Thick rule above the header
    hline_top(part = "header", border = officer::fp_border(width = 1.5)) |>
    
    ## Thin rule under the spanner, spanning only the columns the spanner covers
    hline(i = 1, j = spanner_cols, part = "header",
          border = officer::fp_border(width = 0.5)) |>
   
    ## Rule separating header from body
    hline_bottom(part = "header", border = officer::fp_border(width = 0.75)) |>
   
     ## Thick rule closing the table
    hline_bottom(part = "body", border = officer::fp_border(width = 1.5)) |>
    bg(i = stripe_rows, bg = "#F5F5F5", part = "body") |>
    set_caption(caption) |>
    font(fontname = "Times New Roman", part = "all") |>
    fontsize(size = 10, part = "all") |>
    padding(padding.top = 3, padding.bottom = 3, part = "all") |>
    autofit()
}


## Export the agreement table as a Flextable
agreement_ft <- agreement_table |>
  flextable() |>
  add_header_row(values = c("", "", "Coder agreement (%)"),
                 colwidths = c(1, 1, 2)) |>
  align(align = "center", part = "all") |>
  align(j = 1, align = "left", part = "all") |>
  colformat_double(j = c("Any Error Present", "Exact Same Types of Error"), digits = 1) |>
  apa_theme(
    caption = paste("Table 1. Agreement between two independent coders at each chain position",
                    "for whether a participant explanation contained any error and for whether the",
                    "coders identified exactly the same errors."),
    stripe_rows  = seq(2, nrow(agreement_table), 2),
    spanner_cols = 3:4
  )


## Now doing the same for Usage table
usage_ft <- usage_table |>
  flextable() |>
  set_header_labels(
    chain_position = "Chain position",
    error_type     = "Type of error",
    pct_agreement  = "Coder Agreement (%)",
    pct_coder1     = "Coder 1",
    pct_coder2     = "Coder 2"
  ) |>
  add_header_row(values = c("", "", "", "Explanations marked with this error (%)"),
                 colwidths = c(1, 1, 1, 2)) |>
  merge_v(j = "chain_position") |>
  valign(j = 1, valign = "top", part = "body") |>
  align(align = "center", part = "all") |>
  align(j = 1:2, align = "left", part = "all") |>
  colformat_double(j = c("pct_agreement", "pct_coder1", "pct_coder2"), digits = 1) |>
  apa_theme(
    caption = paste("Table 2. Agreement between the two coders within each error",
                    "category, alongside the percentage of explanations each coder",
                    "marked with that error."),
    stripe_rows  = which(usage_table$chain_position %% 2 == 0),
    spanner_cols = 4:5
  ) |>
  ## Thin rules between chain position blocks, added after the theme so they are not
  ## removed by border_remove()
  hline(i = head(cumsum(rle(usage_table$chain_position)$lengths), -1),
        part = "body", border = officer::fp_border(width = 0.5))

usage_ft

