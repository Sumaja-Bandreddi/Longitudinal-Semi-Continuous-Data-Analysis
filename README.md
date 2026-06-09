## Background
---

This repository reproduces and extends the longitudinal mosaic chromosomal alteration (mCA) analyses reported by Kelly et al. (2025) using data from the Prostate, Lung, Colorectal, and Ovarian (PLCO) Cancer Screening Trial.

To evaluate the robustness of previously reported findings, we reanalyze the PLCO longitudinal mCA data using statistical methods that explicitly accommodate semi-continuous outcomes, including:

- Tobit mixed-effects models
- Two-part (hurdle) mixed-effects models
- Linear and logistic models for comparison

The analyses are performed separately for:

- Mosaic loss of chromosome Y (mLOY)
- Mosaic loss of chromosome X (mLOX)
- Autosomal mCAs (gains, losses, and copy-neutral loss of heterozygosity)
Our goal is to compare covariate effects and longitudinal growth estimates across modeling frameworks and assess whether the associations reported by Kelly et al. remain consistent when the semi-continuous nature of clonal fraction measurements is explicitly modeled.
---

## Repository Structure

| File | Description |
|--------|------------|
| `mLOY Slope Analysis.Rmd` | Longitudinal analysis of mosaic Loss of Y (mLOY) |
| `mLOX Slope Analysis.Rmd` | Longitudinal analysis of mosaic Loss of X (mLOX) |
| `Autosomal Slope Analysis.Rmd` | Longitudinal analysis of autosomal mCAs (Gain, Loss, CN-LOH) |

---

## Data Requirements

The scripts expect an input file:

```text
PLCO_longitudinal_mCAs_0s.xlsx
```

Required variables include:

| Variable | Description |
|-----------|-------------|
| `plco_id` | Participant identifier |
| `mCA_code3` | mCA event identifier |
| `type_FINAL` | mCA type |
| `cf_FINAL` | Clonal fraction |
| `STUDY_YR` | Years since baseline |
| `age_at_collection` | Age at blood collection |
| `smoking_status` | Smoking category |

---

## Statistical Models

### 1. Linear Mixed-Effects Models (LME)

Used for positive clonal fractions after log transformation.

Implemented using:

```r
lme4::lmer()
```

---

### 2. Gamma Mixed Models

Used for strictly positive clonal fractions without transformation.

Implemented using:

```r
glmmTMB::glmmTMB(
  family = Gamma(link = "log")
)
```

---

### 3. Logistic Mixed Models (GLMM)

Models the probability that an mCA is detectable (`CF > 0`).

Implemented using:

```r
lme4::glmer()
```

with a probit link.

---

### 4. Bayesian Joint Models

Jointly model:

- Presence/absence of an mCA
- Magnitude of clonal fraction among positive observations

Implemented using:

```r
brms
cmdstanr
```

Two continuous outcome distributions are considered:

#### Lognormal Joint Model

```r
family = lognormal()
```

#### Gamma Joint Model

```r
family = Gamma(link = "log")
```

Shared random effects allow correlation between:

- Detection probability
- Clonal fraction growth

---

### 5. Tobit Models

Used to account for left-censoring caused by observations with clonal fraction equal to zero.

Implemented using:

```r
brms
family = gaussian()
cens()
```

---

## Autosomal mCA Modeling

Autosomal analyses account for multiple mCA events per individual using hierarchical random effects:

```r
(1 + STUDY_YR | plco_id / mca_num)
```

Bayesian joint models include random effects at:

1. Participant level
2. Event-within-participant level

This structure captures heterogeneity across both individuals and distinct autosomal mCA events.

---

## Additional Analyses

### Smoking and Age Models

The mLOY and mLOX workflows include extensions evaluating:

- Age effects
- Smoking status effects
- Age × Smoking interactions

Example model:

```r
cf_FINAL ~ age_65 * smoking_status +
           (1 + age_65 | plco_id)
```

---

## Model Comparisons

For each mCA category, the scripts compare:

- Separate logistic and continuous models
- Joint Bayesian models

Reported summaries include:

- Fixed effects
- Random-effect variance components
- Random-effect correlations
- Standard errors

These comparisons assess whether joint modeling improves parameter estimation and captures dependence between mCA detection and clonal fraction growth.

---

## Software Requirements

### Required R Packages

```r
install.packages(c(
  "readxl",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "lme4",
  "glmmTMB",
  "brms"
))
```

For Bayesian models, install and configure CmdStanR:

```r
install.packages("cmdstanr")
```

See:

<https://mc-stan.org/cmdstanr/>

---

## Running the Analyses

Render any workflow directly from R:

```r
rmarkdown::render("mLOY Slope Analysis.Rmd")
```

```r
rmarkdown::render("mLOX Slope Analysis.Rmd")
```

```r
rmarkdown::render("Autosomal Slope Analysis.Rmd")
```

---

## Outputs

The analyses produce:

- Fixed-effect estimates
- Random-effect variance components
- Random-effect correlations
- Bayesian posterior summaries
- Comparisons between joint and separate modeling approaches

---
