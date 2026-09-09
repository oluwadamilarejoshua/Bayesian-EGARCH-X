# ============================================================
# run_baseline_dispersed.R
# Diagnostic: refit the ACTUAL baseline model (tight informative
# priors -- betaGDP ~ N(-0.5,0.4), betaDebt ~ N(0.5,0.4), the model
# that "cleanly converged" before) but with 4 chains DELIBERATELY
# seeded at very different starting points:
#   Chain 1: the standard "safe" start used in every prior run
#   Chain 2: seeded near the baseline's OWN previously-found answer
#            (phi~0.91, betaDebt~+0.11)
#   Chain 3: seeded near the 2016-exclusion run's answer
#            (phi~0.80, betaDebt~-0.90)
#   Chain 4: a wildcard start (low persistence, wrong-signed betas)
#
# If the tight prior genuinely dominates, all 4 chains should
# reconverge to the same place despite wildly different starts.
# If they don't (Rhat blows up here too), the baseline's earlier
# "clean" convergence was an artifact of always using the same
# shared initialization, not genuine robustness.
# Reduced iterations (diagnostic purposes, not final estimates).
# ============================================================

if (.Platform$OS.type == "windows" && nchar(Sys.which("make")) == 0) {
  Sys.setenv(PATH = paste("C:/rtools45/usr/bin", Sys.getenv("PATH"), sep = ";"))
}

suppressMessages(library(rstan))
rstan_options(auto_write = TRUE)
options(mc.cores = 1)

data_path <- "c:/Users/joshua.ogundairo/Downloads/My Research Works/Impacts_of_GDP_and_Debt_Servicing_on_Exchange_Rate_Volatility_of_the_Nigerian_Naira__A_Bayesian_EGARCH_Approach/Analysis/Data/nigeria_monthly_stan_input_post1999.csv"
stan_file <- "c:/Users/joshua.ogundairo/Downloads/My Research Works/Impacts_of_GDP_and_Debt_Servicing_on_Exchange_Rate_Volatility_of_the_Nigerian_Naira__A_Bayesian_EGARCH_Approach/Analysis/Bayesian-EGARCH/egarch_x_model.stan"
out_csv   <- "C:/Users/JOSHUA~1.OGU/AppData/Local/Temp/claude/c--Users-joshua-ogundairo-Downloads-My-Research-Works-Impacts-of-GDP-and-Debt-Servicing-on-Exchange-Rate-Volatility-of-the-Nigerian-Naira--A-Bayesian-EGARCH-Approach/89ca73a5-75a5-4283-8ef8-9cb30222e7e7/scratchpad/sensitivity/baseline_dispersed_perchain.csv"

df <- read.csv(data_path, stringsAsFactors = FALSE)
df$date <- as.Date(df$date)
T_full <- nrow(df)
y      <- df$y[2:T_full]
y_lag  <- df$y[1:(T_full - 1)]
gdp    <- as.numeric(scale(df$gdp[2:T_full]))
debt   <- as.numeric(scale(df$debt[2:T_full]))
T_obs  <- length(y)
h1_init <- var(df$y)

stan_data <- list(T = T_obs, y = y, y_lag = y_lag, gdp = gdp, debt = debt, h1_init = h1_init)

logit <- function(p) log(p / (1 - p))

init_list <- list(
  # Chain 1: standard safe start (phi = 0.50)
  list(beta0 = 0.006, betaExR = 0.20, omega = -0.50, alpha = 0.15,
       gamma = 0.10, phi_raw = 0, betaGDP = -0.30, betaDebt = 0.30),
  # Chain 2: near the baseline's own previously-found mode (phi ~ 0.91)
  list(beta0 = 0.0004, betaExR = 0.26, omega = -0.52, alpha = 0.74,
       gamma = 0.19, phi_raw = logit(0.91), betaGDP = -0.59, betaDebt = 0.11),
  # Chain 3: near the 2016-exclusion run's mode (phi ~ 0.80)
  list(beta0 = 0.0001, betaExR = 0.02, omega = -1.50, alpha = 1.20,
       gamma = 0.34, phi_raw = logit(0.80), betaGDP = -0.70, betaDebt = -0.85),
  # Chain 4: wildcard -- low persistence, wrong-signed betas
  list(beta0 = 0.001, betaExR = 0.10, omega = -2.00, alpha = 0.30,
       gamma = -0.10, phi_raw = logit(0.27), betaGDP = 0.20, betaDebt = -0.30)
)

cat("Fitting baseline (tight prior) with DISPERSED inits (reduced iter)...\n")
t0 <- Sys.time()
fit <- stan(
  file    = stan_file,
  data    = stan_data,
  iter    = 3000,
  warmup  = 1000,
  chains  = 4,
  seed    = 42,
  init    = init_list,
  control = list(adapt_delta = 0.92, max_treedepth = 15)
)
cat("Fit time:", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")

pars_all <- c("omega", "alpha", "gamma", "phi", "betaGDP", "betaDebt")
csum <- summary(fit, pars = pars_all)$c_summary

rows <- list()
for (p in pars_all) {
  for (ch in 1:4) {
    rows[[length(rows) + 1]] <- data.frame(
      parameter = p, chain = ch,
      mean = csum[p, "mean", ch], sd = csum[p, "sd", ch],
      q2.5 = csum[p, "2.5%", ch], q97.5 = csum[p, "97.5%", ch]
    )
  }
}
per_chain <- do.call(rbind, rows)
write.csv(per_chain, out_csv, row.names = FALSE)
cat("Saved:", out_csv, "\n\n")

cat("=== PER-CHAIN MEANS (this is the whole point of this run) ===\n")
print(per_chain, row.names = FALSE)

pooled <- summary(fit, pars = pars_all)$summary
cat("\n=== Pooled Rhat ===\n")
print(round(pooled[, "Rhat"], 3))
cat("\n=== Pooled summary ===\n")
print(round(pooled[, c("mean", "2.5%", "97.5%", "n_eff", "Rhat")], 4))

cat("\nDone: baseline_dispersed\n")
