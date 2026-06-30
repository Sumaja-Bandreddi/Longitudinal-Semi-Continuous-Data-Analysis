library(brms)
library(MASS)
library(dplyr)
library(readr)

Sys.setenv(TMPDIR = Sys.getenv("LSCRATCH"))

simulate_and_fit <- function(seed) {
  
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
  # 2) Positive-part parameters (ONLY these are fixed)
  # --------------------------
  B0 <- 0.5
  B1 <- 0.4
  
  sigma_b   <- 0.8
  sigma_eps <- 0.6
  
  c <- 0.7
  
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
  
  Y <- numeric(n)
  idx <- which(I_pos == 1)
  Y[idx] <- exp(mu_pos[idx] + rnorm(length(idx), 0, sigma_eps))
  
  dat <- data.frame(id = id, t = tij, Y = Y, Ipos = I_pos)
  
  # # check if c value is good
  # table(dat$Ipos)
  # prop_censored <- mean(dat$Y <= c)    # fraction at or below threshold
  # prop_censored
  
  # --------------------------
  # 6) Fit model
  # --------------------------
  fit <- brm(
    bf(Ipos ~ t + (1 | p | id), family = bernoulli(link = "probit")) +
      bf(Y | subset(Ipos == 1) ~ t + (1 | p | id), family = lognormal()),
    data   = dat,
    chains = 4,
    cores  = 4,
    iter   = 4000,
    warmup = 2000,
    control = list(adapt_delta = 0.95),
    backend = "rstan"
  )
  
  # ---------------------------
  # Extract summaries
  # ---------------------------
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
  
  return(res)
}

# -------------------------------
# SLURM array task ID
# -------------------------------
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# -------------------------------
# Run simulation
# -------------------------------
res <- simulate_and_fit(seed = task_id)
print(res)

# -------------------------------
# Save results
# -------------------------------
dir.create("results_hurd", showWarnings = FALSE)

write.csv(
  res,
  file = paste0("results_hurd/sim_", task_id, ".csv"),
  row.names = FALSE
)


#############################################################################################################################
# -------------------------------
# Combine simulation results
# -------------------------------
files <- list.files(
  pattern = "^sim_.*\\.csv$",
  full.names = TRUE
)

all_sim <- lapply(files, read_csv) %>% bind_rows()

# -------------------------------
# True values
# -------------------------------
true_values <- data.frame(
  parameter = c(
    "b_Y_Intercept",            # fixed effect lognormal intercept
    "b_Y_t",                     # fixed effect lognormal slope
    "b_Ipos_Intercept",          # fixed effect probit intercept
    "b_Ipos_t",                  # fixed effect probit slope
    "sd_id__Y_Intercept",        # SD of lognormal intercept RE
    "sd_id__Ipos_Intercept",     # SD of probit intercept RE
    "cor_id__Ipos_Intercept__Y_Intercept",  # intercept correlation
    "sigma_Y"                    # residual SD lognormal
  ),
  true_value = c(
    0.5,            
    0.4,         
    (0.5 - log(0.7)) / 0.6,
    0.0,      
    0.8,            
    0.8 / 0.6,      
    1.0,            
    0.6           
  )
)


# -------------------------------
# Summarize across simulations
# -------------------------------
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
    "b_Y_Intercept",                  # β_Intercept
    "b_Y_t",                           # β_t
    "b_Ipos_Intercept",                # α_Intercept
    "b_Ipos_t",                        # α_t
    "sd_id__Y_Intercept",              # σ_b0
    "sd_id__Ipos_Intercept",           # σ_c0
    "cor_id__Ipos_Intercept__Y_Intercept", # cor(b0, c0)
    "sigma_Y"                          # σ_ε
  )))

# -------------------------------
# Clean parameter names
# -------------------------------
summary_table$parameter <- c(
  "β_Intercept",
  "β_t",
  "α_Intercept",
  "α_t",
  "σ_b0",
  "σ_c0",
  "cor(b0, c0)",
  "σ_ε"
)

print(summary_table)

