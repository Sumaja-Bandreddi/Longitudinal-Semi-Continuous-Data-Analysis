library(brms)
library(MASS)
library(dplyr)
library(readr)
library(truncnorm)

# around 12% below censored value 

simulate_hurdle_fit_tobit <- function(seed) {
  
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
  # 2) Positive-part parameters
  # --------------------------
  B0 <- 0.5
  B1 <- 0.4
  
  sigma_b   <- 0.8
  sigma_eps <- 0.6
  
  c = 0.7
  
  # --------------------------
  # 3) Draw random intercepts
  # --------------------------
  b_i0 <- rnorm(n_id, 0, sigma_b)
  b0   <- b_i0[id]
  
  # --------------------------
  # 4) Construct implied probit parameters
  # --------------------------
  alpha0 <- (B0 - log(c)) / sigma_eps
  alpha1 <- 0
  c0     <- b0 / sigma_eps   
  
  # --------------------------
  # 5) Probit hurdle
  # --------------------------
  eta_hu <- alpha0 + alpha1 * tij + c0
  p_pos  <- pnorm(eta_hu)
  I_pos  <- rbinom(n, 1, p_pos)
  
  # --------------------------
  # 6) Positive part
  # --------------------------
  mu_pos <- B0 + B1 * tij + b0
  
  log_Y <- rep(NA_real_, n) 
  idx <- which(I_pos == 1)
  if (length(idx) > 0) {
    log_Y[idx] <- rtruncnorm(
      n = length(idx),
      a = log(c),
      b = Inf,
      mean = mu_pos[idx],
      sd = sigma_eps
    )
  }
  
  Y <- rep(NA_real_, n)
  Y[idx] <- exp(log_Y[idx])
  
  dat <- data.frame(id = id, t = tij, Y = Y, Ipos = I_pos, log_Y = log_Y)
  
  # # diagonostics
  # # 1) Histogram of log_Y for positives
  # logY_pos <- log_Y[!is.na(log_Y)]
  # 
  # if (length(logY_pos) > 0) {
  #   hist(logY_pos,
  #        breaks = 40,
  #        main = "Histogram of log(Y) for positives",
  #        xlab = "log(Y)")
  #   
  #   # show truncation point
  #   abline(v = log(c), col = "red", lwd = 2)
  # }
  
  
  # --------------------------
  # 7) Tobit prep: censor at c
  # --------------------------
  dat$logY_obs <- ifelse(!is.na(dat$log_Y), dat$log_Y, log(c))
  dat$cens_log  <- ifelse(!is.na(dat$log_Y), 0L, -1L) 
  
  # --------------------------
  # 8) Fit TOBIT model in brms
  # --------------------------
  fit <- brm(
    bf(logY_obs | cens(cens_log) ~ t + (1 | id)),
    data   = dat,
    family = gaussian(),
    chains = 4,
    cores  = 4,
    iter   = 4000,
    warmup = 2000,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    backend = "rstan"
  )
  
  # ---------------------------
  # 9) Extract summaries (TOBIT params)
  # ---------------------------
  res <- as.data.frame(
    posterior_summary(fit)[
      c(
        "b_Intercept",          # fixed intercept (log scale)
        "b_t",                  # fixed slope (log scale)
        "sd_id__Intercept",     # SD random intercept
        "sigma"                 # residual SD (log scale)
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
res <- simulate_hurdle_fit_tobit(seed = task_id)
print(res)

# -------------------------------
# Save results
# -------------------------------
dir.create("results_tob", showWarnings = FALSE)

write.csv(
  res,
  file = paste0("results_tob/sim_", task_id, ".csv"),
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

