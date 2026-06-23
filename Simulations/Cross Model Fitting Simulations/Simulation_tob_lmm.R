library(lme4)
library(dplyr)
library(readr)

simulate_and_fit_lmer <- function(seed) {
  
  set.seed(seed)
  
  # --------------------------
  # 1) Design
  # --------------------------
  n_id   <- 1000
  m      <- 5
  t_grid <- seq(0, 1, length.out = m)
  
  id  <- rep(seq_len(n_id), each = m)
  tij <- rep(t_grid, times = n_id)
  n   <- length(id)
  
  # --------------------------
  # 2) True parameters
  # --------------------------
  B0 <- 0.5
  B1 <- 0.4
  
  sigma_b   <- 0.8
  sigma_eps <- 0.6
  
  c <- 0.7   # censoring threshold
  
  # --------------------------
  # 3) Random intercepts
  # --------------------------
  b_i0 <- rnorm(n_id, 0, sigma_b)
  b0   <- b_i0[id]
  
  # --------------------------
  # 4) Latent log-scale model
  # --------------------------
  mu_log   <- B0 + B1 * tij + b0
  logY_star <- rnorm(n, mean = mu_log, sd = sigma_eps)
  Y_star    <- exp(logY_star)
  
  # --------------------------
  # 5) Left censoring
  # --------------------------
  Y_obs <- ifelse(Y_star <= c, 0, Y_star)
  cens  <- ifelse(Y_star <= c, "left", "none")
  
  dat <- data.frame(
    id = factor(id),
    t  = tij,
    Y  = Y_obs,
    cens = cens
  )
  
  # Keep only positive observations
  dat_pos <- subset(dat, Y > 0)
  
  # Fit LMM on positive values only
  fit <- lmer(
    log(Y) ~ t + (1 | id),
    data = dat_pos,
    REML = FALSE
  )
  
  # --------------------------
  # 7) Extract summaries
  # --------------------------
  summ <- summary(fit)
  
  fe <- as.data.frame(summ$coefficients)
  re_sd <- as.numeric(attr(VarCorr(fit)$id, "stddev"))
  resid_sd <- sigma(fit)
  
  res <- data.frame(
    parameter = c("b_Intercept", "b_t", "sd_id__Intercept", "sigma"),
    Estimate = c(
      fe["(Intercept)", "Estimate"],
      fe["t", "Estimate"],
      re_sd,
      resid_sd
    ),
    Est.Error = c(
      fe["(Intercept)", "Std. Error"],
      fe["t", "Std. Error"],
      NA,
      NA
    )
  )
  
  return(res)
}

# -------------------------------
# SLURM array task ID
# -------------------------------
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# -------------------------------
# Run simulation
# -------------------------------
res <- simulate_and_fit_lmer(seed = task_id)
print(res)

# -------------------------------
# Save results
# -------------------------------
dir.create("results_lmer", showWarnings = FALSE)

write.csv(
  res,
  file = paste0("results_lmer/sim_", task_id, ".csv"),
  row.names = FALSE
)

# ==========================================================
# COMBINE RESULTS
# ==========================================================

files <- list.files(
  path = "results_lmer",
  pattern = "^sim_.*\\.csv$",
  full.names = TRUE
)

all_sim <- lapply(files, read_csv) %>% bind_rows()

true_values <- data.frame(
  parameter = c("b_Intercept", "b_t", "sd_id__Intercept", "sigma"),
  true_value = c(0.5, 0.4, 0.8, 0.6)
)

summary_table <- all_sim %>%
  group_by(parameter) %>%
  summarise(
    mean_estimate = mean(Estimate, na.rm = TRUE),
    mean_error    = mean(Est.Error, na.rm = TRUE),
    mc_sd         = sd(Estimate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(true_values, by = "parameter") %>%
  mutate(
    bias = mean_estimate - true_value,
    relative_bias = bias / true_value
  ) %>%
  select(
    parameter,
    true_value,
    mean_estimate,
    bias,
    relative_bias,
    mean_error,
    mc_sd
  ) %>%
  arrange(factor(parameter, levels = c(
    "b_Intercept",
    "b_t",
    "sd_id__Intercept",
    "sigma"
  )))

summary_table$parameter <- c("β_Intercept", "β_t", "σ_b0", "σ_ε")

print(summary_table)
