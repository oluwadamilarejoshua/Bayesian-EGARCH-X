# ============================================================
# run_flat_perchain.R
# Diagnostic rerun of the flat-prior variant, retaining PER-CHAIN
# summaries (not just the pooled summary) to determine whether the
# non-convergence we saw is genuine multimodality (chains settle
# individually into different stable places) or a sampler failure
# (chains never settle anywhere at all).
# Reduced iterations (diagnostic purposes, not final estimates).
# ============================================================

if (.Platform$OS.type == "windows" && nchar(Sys.which("make")) == 0) {
  Sys.setenv(PATH = paste("C:/rtools45/usr/bin", Sys.getenv("PATH"), sep = ";"))
}

suppressMessages(library(rstan))
rstan_options(auto_write = TRUE)
options(mc.cores = 1)

data_path <- "c:/Users/joshua.ogundairo/Downloads/My Research Works/Impacts_of_GDP_and_Debt_Servicing_on_Exchange_Rate_Volatility_of_the_Nigerian_Naira__A_Bayesian_EGARCH_Approach/Analysis/Data/nigeria_monthly_stan_input_post1999.csv"
stan_file <- "C:/Users/JOSHUA~1.OGU/AppData/Local/Temp/claude/c--Users-joshua-ogundairo-Downloads-My-Research-Works-Impacts-of-GDP-and-Debt-Servicing-on-Exchange-Rate-Volatility-of-the-Nigerian-Naira--A-Bayesian-EGARCH-Approach/89ca73a5-75a5-4283-8ef8-9cb30222e7e7/scratchpad/sensitivity/egarch_x_flat.stan"
out_csv   <- "C:/Users/JOSHUA~1.OGU/AppData/Local/Temp/claude/c--Users-joshua-ogundairo-Downloads-My-Research-Works-Impacts-of-GDP-and-Debt-Servicing-on-Exchange-Rate-Volatility-of-the-Nigerian-Naira--A-Bayesian-EGARCH-Approach/89ca73a5-75a5-4283-8ef8-9cb30222e7e7/scratchpad/sensitivity/flat_perchain.csv"

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

# Same common "safe" init as every prior run -- the point here is only
# to see where each chain DRIFTS TO under the flat prior, starting
# from the same place they always have.
init_list <- lapply(1:4, function(i) list(
  beta0 = 0.006, betaExR = 0.20, omega = -0.50, alpha = 0.15,
  gamma = 0.10, phi_raw = 0.00, betaGDP = -0.30, betaDebt = 0.30
))

cat("Fitting flat-prior variant (per-chain diagnostic, reduced iter)...\n")
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
csum <- summary(fit, pars = pars_all)$c_summary   # [param, stat, chain]

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
cat("\n=== Pooled Rhat (for reference) ===\n")
print(round(pooled[, "Rhat"], 3))

cat("\nDone: flat_perchain\n")
