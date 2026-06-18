library(brms)
library(MASS)
library(dplyr)
library(readr)

simulate_tobit_probit <- function(seed) {
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
  # Selection equation: probit P(D_ij = 1)
  A0 <- -0.2
  A1 <-  1.0
  sigma_u <- 0.7
  
  # Outcome equation for latent positive Y*
  B0 <- 0.5
  B1 <- 0.4
  
  sigma_b   <- 0.8
  sigma_eps <- 0.6
  
  c <- 0.7   # left-censoring threshold for the Tobit outcome
  
  # --------------------------
  # 3) Random intercepts
  # --------------------------
  u_i <- rnorm(n_id, 0, sigma_u)   # selection random intercept
  u0  <- u_i[id]
  
  b_i0 <- rnorm(n_id, 0, sigma_b)   # outcome random intercept
  b0   <- b_i0[id]
  
  # --------------------------
  # 4) Selection process
  #    pi_ij = P(D_ij = 1) via probit
  # --------------------------
  eta_d <- A0 + A1 * tij + u0
  pi_ij <- pnorm(eta_d)
  D_ij  <- rbinom(n, size = 1, prob = pi_ij)
  
  # --------------------------
  # 5) Latent outcome process
  # --------------------------
  mu    <- B0 + B1 * tij + b0
  Y_star <- rnorm(n, mean = mu, sd = sigma_eps)

  # Tobit-observed outcome
  # If latent outcome falls below c, it is 0
  Y_ij    <- ifelse(Y_star <= c, 0, Y_star)
  cens_ij <- ifelse(Y_star <= c, "left", "none")
  
  # --------------------------
  # 6) Create Y_ij_star = all zeros
  # --------------------------
  Y_ij_star <- rep(0, n)
  
  # --------------------------
  # 7) Final observed outcome
  #    If D_ij = 1 use Y_ij, otherwise 0
  # --------------------------
  Y_ij_final <- ifelse(D_ij == 1, Y_ij, Y_ij_star)
  
  # Censoring indicator for final observed outcome
  # - zeros from D_ij = 0 are left-censored at 0
  # - positive Tobit-censored values stay left-censored at c
  cens_final <- ifelse(D_ij == 0, "left", cens_ij)
  
  dat <- data.frame(
    id         = id,
    t          = tij,
    pi_ij      = pi_ij,
    D_ij       = D_ij,
    Y_star     = Y_star,
    Y_ij       = Y_ij,
    Y_ij_final = Y_ij_final,
    cens_final = cens_final
  )
  
  # --------------------------
  # 8) Fit Tobit model in brms
  # --------------------------
  fit <- brm(
    bf(Y_ij_final | cens(cens_final) ~ t + (1 | id)),
    data   = dat,
    family = gaussian(),
    chains = 4,
    cores  = 4,
    iter   = 4000,
    warmup = 2000,
    control = list(adapt_delta = 0.95),
    backend = "rstan"
  )
  
  # --------------------------
  # 9) Extract summaries
  # --------------------------
  res <- as.data.frame(
    posterior_summary(fit)[
      c("b_Intercept", 
        "b_t", 
        "sd_id__Intercept", 
        "sigma"),
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
res <- simulate_tobit_probit(seed = task_id)
print(res)

# -------------------------------
# Save results
# -------------------------------
dir.create("results_zero_tob", showWarnings = FALSE)

write.csv(
  res,
  file = paste0("results_zero_tob/sim_", task_id, ".csv"),
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
  dplyr::select(
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



