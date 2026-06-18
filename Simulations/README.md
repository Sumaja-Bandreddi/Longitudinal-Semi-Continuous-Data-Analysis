# Simulations

This directory contains simulation studies used to evaluate longitudinal semi-continuous outcome models under a variety of data-generating mechanisms and model-fitting scenarios.

## Overview

The simulations examine the performance of Hurdle and Tobit modeling approaches for longitudinal semi-continuous outcomes with repeated measurements. Scenarios include:

1. Correct model specification
2. Cross-model fitting (model misspecification)
3. Excess-zero data generation
4. Partial equivalence between hurdle and Tobit processes

All simulations were implemented in R using the `brms` package and executed on an HPC cluster using SLURM job arrays.

---

## Folder Structure

### Correctly Specified Model Simulations

Evaluates parameter recovery when the fitted model matches the true data-generating mechanism.

| File | Description |
|--------|------------|
| `Simulation_hurd_hurd.R` | Generates data from a longitudinal hurdle model and fits a hurdle model. |
| `Simulation_tob_tob.R` | Generates data from a longitudinal Tobit model and fits a Tobit model. |
| `run_sim_hurd_hurd.txt` | SLURM submission script for hurdle → hurdle simulations. |
| `run_sim_tob_tob.txt` | SLURM submission script for Tobit → Tobit simulations. |

---

### Cross Model Fitting Simulations

Evaluates robustness to model misspecification.

| File | Description |
|--------|------------|
| `Simulation_hurd_tob.R` | Generates hurdle data and fits a Tobit model. |
| `Simulation_tob_hurd.R` | Generates Tobit data and fits a hurdle model. |
| `run_sim_hurd_tob.txt` | SLURM submission script for hurdle → Tobit simulations. |
| `run_sim_tob_hurd.txt` | SLURM submission script for Tobit → hurdle simulations. |

---

### Excess Zero Simulations

Examines situations where there are excess zeros in the data.

| File | Description |
|--------|------------|
| `Simulation_zero_hurd.R` | Excess-zero data generated using a two-stage selection process and analyzed with a hurdle model. |
| `Simulation_zero_tobit.R` | Excess-zero data generated using a two-stage selection process and analyzed with a Tobit model. |
| `run_zero_hurd.txt` | SLURM submission script for excess-zero hurdle simulations. |
| `run_zero_tobit.txt` | SLURM submission script for excess-zero Tobit simulations. |

---

### Partial Equivalence Simulations

Investigates settings where hurdle and Tobit formulations become partially equivalent (alpha1 = 0).

| File | Description |
|--------|------------|
| `hurd_hurd.R` | Generates and fits a hurdle model under partial-equivalence conditions. |
| `hurd_tob.R` | Generates hurdle-type data and fits a Tobit model under partial-equivalence conditions. |
| `run_hurd.txt` | SLURM submission script for hurdle simulations. |
| `run_tob.txt` | SLURM submission script for Tobit simulations. |

---

## Simulation Design

Across most scenarios:

- Number of subjects: **1,000**
- Repeated measurements per subject: **5**
- Time points: equally spaced on [0,1]
- Random intercept structure included
- Bayesian estimation performed using `brms`
- Posterior summaries extracted for:
  - Fixed effects
  - Random-effect variances
  - Random-effect correlations
  - Residual variance parameters

Simulation results are saved as individual CSV files and subsequently combined to compute:

- Mean parameter estimates
- Monte Carlo standard deviations
- Average posterior standard errors
- Bias
- Relative bias

---

## Software Requirements

Required R packages include:

```r
library(brms)
library(MASS)
library(dplyr)
library(readr)
library(truncnorm)
```

---

## HPC Execution

Simulations were designed for execution on NIH's Biowulf cluster using SLURM array jobs.

Example:

```bash
sbatch run_sim_hurd_hurd.txt
```

Each array task corresponds to a single Monte Carlo replication and saves results to a scenario-specific output directory. 1,000 simulations were run for each file.

---

## Purpose

These simulations were developed to compare the statistical performance of longitudinal hurdle and Tobit models under varying assumptions regarding:

- Structural zeros
- Left censoring
- Model misspecification
- Partial model equivalence
- Random-effect structures

The results inform model selection and interpretation for longitudinal semi-continuous outcomes.
