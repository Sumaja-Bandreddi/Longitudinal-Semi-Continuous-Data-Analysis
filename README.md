# Longitudinal Semi-Continuous Data Analysis

This repository contains simulation studies and a real-data application comparing statistical models for longitudinal semi-continuous outcomes. Semi-continuous data contain many zero values along with positive continuous measurements, creating an important modeling question: should zeros be treated as true absence, nondetection, or left-censoring below a detection threshold?

The analyses focus on two main modeling frameworks:

- **Tobit mixed-effects models**, which treat zeros as left-censored values from an underlying latent continuous process.
- **Two-part hurdle mixed-effects models**, which separately model whether an outcome is positive and the magnitude of the positive outcome.

The repository combines methodological simulations with an applied clonal hematopoiesis analysis using longitudinal mosaic chromosomal alteration data from the Prostate, Lung, Colorectal, and Ovarian Cancer Screening Trial.

## Usage Examples

The example below generates one longitudinal semi-continuous dataset from a Tobit model and then fits two models to the same simulated data:

1. a correctly specified lognormal Tobit model;
2. a two-part hurdle model.

### Generate data from a Tobit model

The data-generating model uses:

- `n_id = 1000` subjects;
- `m = 5` repeated measurements per subject;
- time scaled from 0 to 1;
- fixed effects `beta_0 = 0.5` and `beta_t = 0.4`;
- random-intercept standard deviation `sigma_b = 0.8`;
- residual standard deviation `sigma_eps = 0.6`;
  
```r

# Required packages 
library(brms)
library(dplyr)

# The following function simulates longitudinal left-censored data from a Tobit mixed-effects model for illustrating model fitting in `brms`. It generates repeated measurements for 1,000 subjects, applies left censoring at a specified threshold, and returns a data frame containing the subject identifier, observation time, observed outcome, and censoring indicator.

generate_tobit_data <- function(seed = 1) {

 # Set random seed for reproducibility
  set.seed(seed)

  #----------------------------------------------------------
  # 1) Study design
  #----------------------------------------------------------
  
  # Number of subjects
  n_id <- 1000
  
  # Number of repeated measurements per subject
  m <- 5
  
  # Observation times, equally spaced between 0 and 1
  t_grid <- seq(0, 1, length.out = m)

  # Create subject IDs (each subject appears m times)
  id <- rep(seq_len(n_id), each = m)

  # Assign observation times to each subject
  tij <- rep(t_grid, times = n_id)

  # Total number of observations
  n <- length(id)

  #----------------------------------------------------------
  # 2) True model parameters
  #----------------------------------------------------------

  # Fixed intercept
  B0 <- 0.5

  # Fixed effect of time
  B1 <- 0.4

  # Standard deviation of the subject-specific random intercept
  sigma_b <- 0.8

  # Residual (within-subject) standard deviation
  sigma_eps <- 0.6

  # Left-censoring threshold (around 15% censored values) 
  c <- 0.7

  #----------------------------------------------------------
  # 3) Generate subject-specific random intercepts
  #----------------------------------------------------------

  # Draw one random intercept for each subject
  b_i0 <- rnorm(n_id, mean = 0, sd = sigma_b)

  # Expand the subject-level random intercepts to every observation
  b0 <- b_i0[id]

  #----------------------------------------------------------
  # 4) Generate the latent (uncensored) outcome
  #----------------------------------------------------------

  # Compute the mean of the latent log outcome
  mu_log <- B0 + B1 * tij + b0

  # Generate latent log-scale observations
  logY_star <- rnorm(
    n,
    mean = mu_log,
    sd = sigma_eps
  )

  # Transform back to the original outcome scale (this guarantees the latent outcome is positive)
  Y_star <- exp(logY_star)

  #----------------------------------------------------------
  # 5) Apply left censoring
  #----------------------------------------------------------

  # Observations below the detection limit are recorded as zero
  Y_obs <- ifelse(Y_star <= c, 0, Y_star)

  # Create the censoring indicator required by brms
  cens <- ifelse(Y_star <= c, "left", "none")

  #----------------------------------------------------------
  # 6) Create the final analysis dataset
  #----------------------------------------------------------
  
  dat <- data.frame(
    id   = id,      # Subject identifier (1000 subjects, each measured at 5 time points)
    t    = tij,     # Observation time (scaled from 0 to 1)
    Y    = Y_obs,   # Observed semi-continuous outcome (0 if left-censored; otherwise the observed positive value)
    cens = cens     # Censoring indicator ("left" = left-censored, "none" = uncensored)
  )

return(dat)
}
```

In this example, `Y` is the observed semi-continuous outcome, where censored values are recorded as zero. 

A lognormal Tobit model is used because the latent response is generated on the log scale and exponentiated, so the uncensored observations are lognormally distributed. The observed zeros result exclusively from left censoring below the detection limit rather than from a distinct zero-generating process, making the lognormal Tobit model a special case of a hurdle model.

### Example 1: Fit a lognormal Tobit model

```r
#The following function fits a Bayesian Tobit mixed-effects model to the simulated data using `brms`. The model includes a fixed effect for time, a subject-specific random intercept, and a lognormal response distribution with left censoring.

fit_tobit_model <- function(dat) {
  brm(
    bf(Y | cens(cens) ~ t + (1 | id)),
    data = dat,
    family = lognormal(),
    chains = 4,
    cores = 4,
    iter = 4000,
    warmup = 2000,
    control = list(adapt_delta = 0.95),
    backend = "rstan"
  )
}

#The following function extracts the posterior means and posterior standard deviations for the fixed effects, random-intercept standard deviation, and residual standard deviation from the fitted Tobit model.

summarize_tobit_fit <- function(fit) {
  res <- as.data.frame(
    posterior_summary(fit)[
      c(
        "b_Intercept",
        "b_t",
        "sd_id__Intercept",
        "sigma"
      ),
      c("Estimate", "Est.Error")
    ]
  )

  res$parameter <- rownames(res)
  rownames(res) <- NULL
  res
}

tobit_fit <- fit_tobit_model(dat)
tobit_results <- summarize_tobit_fit(tobit_fit)

```

The Tobit model estimates the fixed intercept and time effect, the random-intercept standard deviation, and the residual log-scale standard deviation. 

### Example 2: Fit a hurdle model to the same Tobit-generated data

```r
fit_hurdle_model <- function(dat) {
  brm(
    bf(Ipos ~ t + (1 | p | id), family = bernoulli(link = "probit")) +
      bf(Y | subset(Ipos == 1) + trunc(lb = 0.7) ~ t + (1 | p | id),
         family = lognormal()),
    data = dat,
    chains = 4,
    cores = 4,
    iter = 4000,
    warmup = 2000,
    control = list(adapt_delta = 0.95),
    backend = "rstan"
  )
}

summarize_hurdle_fit <- function(fit) {
  res <- as.data.frame(
    posterior_summary(fit)[
      c(
        "b_Y_Intercept",
        "b_Y_t",
        "b_Ipos_Intercept",
        "b_Ipos_t",
        "sd_id__Y_Intercept",
        "sd_id__Ipos_Intercept",
        "cor_id__Ipos_Intercept__Y_Intercept",
        "sigma_Y"
      ),
      c("Estimate", "Est.Error")
    ]
  )

  res$parameter <- rownames(res)
  rownames(res) <- NULL
  res
}

hurdle_fit <- fit_hurdle_model(dat)
hurdle_results <- summarize_hurdle_fit(hurdle_fit)

print(hurdle_results)
```

The hurdle model separates the probability of being above the censoring limit from the positive outcome distribution. A probit link is used for the hurdle component because the Tobit model can be viewed as a special case of the hurdle model under a probit link, where both the binary and continuous components arise from the same latent normal process

### Running the Provided Simulation Scripts

Each simulation script is written to run one simulation replicate at a time. The replicate seed is read from the SLURM array task ID:

```r
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
res <- simulate_and_fit_tobit(seed = task_id)
```

For example, to run one replicate interactively:

```r
Sys.setenv(SLURM_ARRAY_TASK_ID = 1)
source("Correctly Specified Model Simulations/Simulation_tob_tob.R")
```

or:

```r
Sys.setenv(SLURM_ARRAY_TASK_ID = 1)
source("Cross Model Fitting Simulations/Simulation_tob_hurd.R")
```

The scripts save one CSV file per replicate:

- `results_tob_tob/sim_<task_id>.csv` for the Tobit fit;
- `results_tob_hurd/sim_<task_id>.csv` for the hurdle fit.

### Output Parameters

The Tobit script returns posterior summaries for:

- `b_Intercept`: fixed intercept on the log scale;
- `b_t`: fixed time effect on the log scale;
- `sd_id__Intercept`: subject-level random-intercept standard deviation;
- `sigma`: residual standard deviation on the log scale.

The hurdle script returns posterior summaries for:

- `b_Y_Intercept`: fixed intercept for the positive lognormal outcome model;
- `b_Y_t`: fixed time effect for the positive lognormal outcome model;
- `b_Ipos_Intercept`: fixed intercept for the binary positive-outcome model;
- `b_Ipos_t`: fixed time effect for the binary positive-outcome model;
- `sd_id__Y_Intercept`: subject-level random-intercept standard deviation for the positive outcome model;
- `sd_id__Ipos_Intercept`: subject-level random-intercept standard deviation for the binary model;
- `cor_id__Ipos_Intercept__Y_Intercept`: correlation between subject-level random intercepts across the binary and positive outcome parts;
- `sigma_Y`: residual standard deviation for the positive lognormal outcome model.


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

## References 

Su, L. and Tom, B. D. and Farewell, V. T. (2009) Bias in 2-part mixed models for longitudinal semicontinuous data. Biostatistics, 10(2): 374-389.

Albert, P. S. and Shen, J. (2005) Modelling longitudinal semicontinuous emesis volume data with serial correlation in an acupuncture clinical trial. Journal of the Royal Statistical Society Series C: Applied Statistics, 54(4): 707-720.

Buerkner, Paul-Christian (2017) brms: An R package for Bayesian multilevel models using Stan. Journal of Statistical Software, 80(1): 1-28.

Amemiya, Takeshi (1984) Tobit models: A survey. Journal of Econometrics, 24(1-2): 3-61.

Kelly, R. L. and Brown, D. W. and Zhou, W. and Hubbard, A. K. and Young, C. D. and Barnao, K. M. and Klein, A. and Dutta, D. and Vogt, A. and Liu, J. and Wang, J. and Huang, W.-Y. and Freedman, N. D. and Chanock, S. J. and Albert, P. S. and Machiela, M. J. (2026) Longitudinal characterization of mosaic chromosomal alterations identifies factors influencing clonal dynamics of leukocytes. Nature Communications. Under revision.
