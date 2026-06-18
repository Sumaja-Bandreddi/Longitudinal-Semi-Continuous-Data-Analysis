library(brms)
library(MASS)
library(dplyr)
library(readr)

simulate_and_fit_tobit <- function(seed) {
  
  set.seed(seed)
  
  # --------------------------
  # 1) Design
  # --------------------------
  n_id  <- 1000
  m     <- 5
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
  mu_log <- B0 + B1 * tij + b0
  
  logY_star <- rnorm(n, mean = mu_log, sd = sigma_eps)
  Y_star    <- exp(logY_star)
  
  # --------------------------
  # 5) Left censoring 
  # --------------------------
  Y_obs <- ifelse(Y_star <= c, c, Y_star)
  cens  <- ifelse(Y_star <= c, "left", "none")
  
  dat <- data.frame(
    id = id,
    t  = tij,
    Y  = Y_obs,
    cens = cens
  )
  
  # message("Proportion censored: ",
  #         round(mean(dat$cens == "left"), 3))
  
  # --------------------------
  # 6) Fit lognormal Tobit
  # --------------------------
  fit <- brm(
    bf(Y | cens(cens) ~ t + (1 | id)),
    data   = dat,
    family = lognormal(),
    chains = 4,
    cores  = 4,
    iter   = 4000,
    warmup = 2000,
    control = list(adapt_delta = 0.95),
    backend = "rstan"
  )
  
  # ---------------------------
  # 9) Extract summaries (TOBIT params)
  # ---------------------------
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
  
  return(res)
}

# -------------------------------
# SLURM array task ID
# -------------------------------
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# -------------------------------
# Run simulation
# -------------------------------
res <- simulate_and_fit_tobit(seed = task_id)
print(res)

# -------------------------------
# Save results
# -------------------------------
dir.create("results_tob_tob", showWarnings = FALSE)

write.csv(
  res,
  file = paste0("results_tob_tob/sim_", task_id, ".csv"),
  row.names = FALSE
)


# ==========================================================
# COMBINE RESULTS
# ==========================================================

files <- list.files(
  pattern = "^sim_.*\\.csv$",
  full.names = TRUE
)

all_sim <- lapply(files, read_csv) %>% bind_rows()

# ---------------------------
# True values
# ---------------------------

true_values <- data.frame(
  parameter = c(
    "b_Intercept",
    "b_t",
    "sd_id__Intercept",
    "sigma"
  ),
  true_value = c(
    0.5,
    0.4,
    0.8,
    0.6
  )
)

# ---------------------------
# Summarize across simulations
# ---------------------------
summary_table <- all_sim %>%
  group_by(parameter) %>%
  summarise(
    mean_estimate = mean(Estimate),
    mean_error    = mean(Est.Error),
    mc_sd         = sd(Estimate),
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
  )

# -------------------------------
# Reorder for better output
# -------------------------------
summary_table <- summary_table %>%
  arrange(factor(parameter, levels = c(
    "b_Intercept",
    "b_t",
    "sd_id__Intercept",
    "sigma"
  )))

# -------------------------------
# Clean parameter names
# -------------------------------
summary_table$parameter <- c(
  "β_Intercept",
  "β_t",
  "σ_b0",
  "σ_ε"
)

print(summary_table)

