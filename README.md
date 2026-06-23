# Longitudinal Semi-Continuous Data Analysis

This repository contains simulation studies and a real-data application comparing statistical models for longitudinal semi-continuous outcomes. Semi-continuous data contain many zero values along with positive continuous measurements, creating an important modeling question: should zeros be treated as true absence, nondetection, or left-censoring below a detection threshold?

The analyses focus on two main modeling frameworks:

- **Tobit mixed-effects models**, which treat zeros as left-censored values from an underlying latent continuous process.
- **Two-part hurdle mixed-effects models**, which separately model whether an outcome is positive and the magnitude of the positive outcome.

The repository combines methodological simulations with an applied clonal hematopoiesis analysis using longitudinal mosaic chromosomal alteration data from the Prostate, Lung, Colorectal, and Ovarian Cancer Screening Trial.

## Repository Structure

```text
Longitudinal-Semi-Continuous-Data-Analysis/
├── Clonal Hematopoiesis Application/
│   ├── README.md
│   ├── mLOY Slope Analysis.Rmd
│   ├── mLOX Slope Analysis.Rmd
│   └── Autosomal Slope Analysis.Rmd
│
└── Simulations/
    ├── README.md
    ├── Correctly Specified Model Simulations/
    ├── Cross Model Fitting Simulations/
    ├── Excess Zero Simulations/
    └── Partial Equivalence Simulations/
```

## Project Overview

This project has two connected goals:

1. **Simulation studies**  
   Evaluate when Tobit and hurdle models recover parameters accurately, when they are robust to model misspecification, and when excess zeros create bias.

2. **Clonal hematopoiesis application**  
   Reanalyze longitudinal clonal fraction data from PLCO to assess whether previously reported mosaic chromosomal alteration trends remain consistent when semi-continuous outcomes are modeled explicitly.

Together, these analyses help determine when a Tobit model is reasonable, when a hurdle model is preferred, and how model choice affects interpretation in longitudinal biomedical data.

## Clonal Hematopoiesis Application

The `Clonal Hematopoiesis Application` folder reanalyzes longitudinal mosaic chromosomal alteration data from the PLCO Cancer Screening Trial.

The analyses are performed separately for:

- **mLOY**: mosaic loss of chromosome Y
- **mLOX**: mosaic loss of chromosome X
- **Autosomal mCAs**: autosomal gains, losses, and copy-neutral loss of heterozygosity

The outcome is clonal fraction, which is semi-continuous because many observations are zero while detected clones have positive continuous values.

### Models Used

The application compares several modeling approaches:

- Linear mixed-effects models for positive clonal fractions
- Logistic/probit mixed models for detection probability
- Bayesian two-part hurdle models
- Bayesian Tobit models for left-censored outcomes

For mLOY and mLOX, models include age, smoking status, and age-by-smoking interactions. Autosomal analyses also account for multiple mCA events within the same participant using hierarchical random effects.

## Simulation Studies

The `Simulations` folder evaluates Tobit and hurdle models under controlled data-generating mechanisms.

Simulation scenarios include:

- Correctly specified hurdle and Tobit models
- Cross-model fitting, where data are generated from one model and fit with another
- Excess-zero scenarios, where zeros arise from an additional selection process
- Partial-equivalence scenarios where hurdle and Tobit formulations become more similar

Across most simulation settings, datasets include:

- 1,000 subjects
- 5 repeated measurements per subject
- Equally spaced time points
- Subject-level random effects
- Bayesian model fitting using `brms`

Simulation outputs are saved as CSV files and summarized using:

- Mean parameter estimates
- Monte Carlo standard deviations
- Average posterior standard errors
- Bias
- Relative bias

## Data Requirements

The clonal hematopoiesis application expects the following input file:

```text
PLCO_longitudinal_mCAs_0s.xlsx
```

Required variables include:

| Variable | Description |
|---|---|
| `plco_id` | Participant identifier |
| `mCA_code3` | mCA event identifier |
| `type_FINAL` | mCA type |
| `cf_FINAL` | Clonal fraction |
| `STUDY_YR` | Years since baseline |
| `age_at_collection` | Age at blood collection |
| `smoking_status` | Smoking category |

The PLCO data file is not included in this repository and must be supplied separately according to data-use requirements.

## Software Requirements

The analyses are written in R. Main packages include:

```r
install.packages(c(
  "readxl",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "lme4",
  "glmmTMB",
  "brms",
  "MASS",
  "readr",
  "truncnorm"
))
```

## Running the Application Analyses

From the `Clonal Hematopoiesis Application` folder, render the R Markdown workflows:

```r
rmarkdown::render("mLOY Slope Analysis.Rmd")
rmarkdown::render("mLOX Slope Analysis.Rmd")
rmarkdown::render("Autosomal Slope Analysis.Rmd")
```

## Running the Simulations

Simulation scripts are designed for high-performance computing using SLURM job arrays. Example:

```bash
sbatch run_sim_hurd_hurd.txt
```

Each SLURM array task runs one Monte Carlo replication and saves results to a scenario-specific output directory.

## Outputs

The repository produces:

- Model estimates for longitudinal clonal fraction growth
- Detection probability estimates
- Random-effect variance components
- Random-effect correlations
- Bayesian posterior summaries
- Simulation CSV files
- Simulation performance summaries

## Purpose

This project supports methodological comparison and applied interpretation for longitudinal semi-continuous biomedical outcomes. The simulations clarify when Tobit and hurdle models perform well or fail, while the PLCO application shows how these models affect inference for clonal hematopoiesis and mosaic chromosomal alteration growth.
