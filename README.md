# ADHD-Transmission-Chains

# Transmission of Mental Health Information Increases ADHD Self-Diagnosis

Data and analysis code for four transmission-chain studies and a pooled mega-analysis, together with the supplementary information for the paper.

- **Authors:** Tanay Katiyar, George Gillett, Adam Hunt, Amanda Ferguson, Alessia Pascale, Alberto Acerbi, Nikhil Chaudhary*, Amy Orben*
- **Preprint:** [link or DOI]
- **Complete archive (code, data and fitted models):** [Zenodo DOI]

## What is where

| Location | Contains |
|---|---|
| This GitHub repository | Analysis scripts, the supplementary information for the paper, and the `renv` files that record package versions |
| Zenodo archive ([DOI]) | Everything in this repository **plus the fitted model files (`.rds`) and the anonymised data** |

The fitted model files are too large for GitHub, and the anonymised data are provided only in the Zenodo archive. The scripts load both, so most of them will not run from the GitHub copy alone. You have two options:

- **Recommended:** download the complete archive from Zenodo and unzip it. Everything is already in place.
- **Alternative:** clone this repository, download the Zenodo archive, and copy the `.rds` and `.csv` files from each folder in the archive into the folder with the same name here.

## Getting started

1. **Open the project.** Double-click `Chains Repo.Rproj` in the top-level folder. Do not open the script files directly: each script finds its own folder through the project, and the package setup only activates when the project is open.
2. **Restore the packages (first time only).** In the R console, run:

   ```r
   renv::restore()
   ```

   Answer `y` when asked. This installs the package versions used for the paper. `renv` installs itself the first time the project opens, and the first restore can take a while.
3. **Run a script.** Open it from the Files pane inside RStudio and run it.

## Folder structure

```
Chains Repo/
├── Chains Repo.Rproj
├── renv.lock, .Rprofile, renv/     package versions (do not edit)
├── Supplementary Materials/
├── Power Analysis/
├── Explanation Screening Agreement/
├── Study 1/
├── Study 2/
├── Study 3/
├── Study 4/
└── Mega-Analysis/
```

### Inside `Supplementary Materials`

The supplementary information for the associated paper.

### Inside `Power Analysis`

| File | What it does |
|---|---|
| `Power Analysis_Chains.R` | Runs a Bayesian simulation-based power and precision analysis, assuming ordinal data |
| `*.rds` | Saved simulation models and results (Zenodo archive only) |

**This analysis takes a very long time to run.** It fits a large number of simulated models and was run on a computing cluster. Running it on a personal computer is not recommended. To inspect or reproduce the reported power and precision results, load the saved models and results from the Zenodo archive instead of re-running the simulation.

### Inside `Explanation Screening Agreement`

| File | What it does |
|---|---|
| `Explanation Screening Coder Agreement.R` | Computes raw percentage agreement between two independent coders on minor linguistic error recognition coding in participant ADHD explanations |

### Inside each study folder

| File | What it does |
|---|---|
| `Analyses.R` | Cleans the data, fits the main models, runs the hypothesis tests, and exports the dataset used in the mega-analysis |
| `Sensitivity Analyses.R` (Study 1), `SA_Analyses.R` (Studies 2–4) | Fits the sensitivity-analysis models |
| `Demographics Table.R` | Builds the sample characteristics tables |
| `Analyses Reg Tables.R` | Builds the regression tables and posterior plots for the main models |
| `SA_Reg Tables.R` (Study 1), `SA_Analyses Reg Tables.R` (Studies 2–4) | The same for the sensitivity-analysis models |
| `Model Prediction Plots.R` | Plots the model-predicted effects for the main models |
| `SA_Model Prediction Plot.R` (Study 1), `SA_Model Prediction Plots.R` (Studies 2–3) | The same for the sensitivity-analysis models |
| `Table Theme Function.R` | Table formatting helper, loaded by the other scripts; not run on its own |
| `Study[N]_Data.csv` | Anonymised data (Zenodo archive only) |
| `*.rds` | Fitted models (Zenodo archive only) |

### Inside `Mega-Analysis`

The mega-analysis is equivalent to a one-stage individual participant data (IPD) meta-analysis: the participant-level data from all four studies are pooled and analysed together in a single set of models, instead of combining summary estimates from each study.

| File | What it does |
|---|---|
| `Analyses.R` | Pools the four studies and fits the mega-analysis models |
| `Mediation Analyses.R` | Computes semantic drift and fits the mediation models |
| `Summary Table.R` | Summarises the structure of the pooled dataset |
| `Analyses Reg Tables.R`, `Mediation Analyses Reg Tables.R` | Build the regression and mediation tables and posterior plots |
| `Model Prediction Plots.R` | Plots the model-predicted effects, including the main manuscript figure |
| `Table Theme Function.R` | Table formatting helper |
| `*.rds` | Fitted models and saved text embeddings (Zenodo archive only) |

## Order to run things

1. Within a study folder, run `Analyses.R` first. The table and plot scripts depend on its models.
2. Run the sensitivity-analysis script before the `SA_` table and plot scripts.
3. Run the `Mega-Analysis` scripts last. They use the datasets exported by each study's `Analyses.R`. Within that folder, run `Analyses.R`, then `Mediation Analyses.R`, then the table and plot scripts.

The `Power Analysis` and `Explanation Screening Agreement` scripts are independent of the study folders and can be run at any point.

## Reproducing results without refitting

In the analysis scripts, each model has a `brm()` call that fits it, followed by a `readRDS()` line that loads the saved fit. Fitting takes a long time, and refitted models can differ very slightly from the saved ones. To reproduce the reported numbers, skip the `brm()` call and run the `readRDS()` line beneath it. This applies most of all to the power analysis, which should be reproduced from the saved files.

## Software requirements

- R 4.5.2 and RStudio. Package versions are recorded in `renv.lock`.
- Refitting models requires a C++ toolchain (Rtools on Windows, Xcode Command Line Tools on Mac). Loading the saved models does not.
- Recomputing the text embeddings in `Mediation Analyses.R` requires a Python environment, set up with `text::textrpp_install()`. This is not needed if you use the saved embeddings from the Zenodo archive.

## Citation

[How to cite the paper and the Zenodo archive]

## Contact

Tanay Katiyar: tanay.katiyar20@gmail.com
