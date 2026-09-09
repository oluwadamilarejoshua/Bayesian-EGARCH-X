# ============================================================
# run_variant.R
# Prior-sensitivity refit driver. Fits egarch_x_model.stan (or a
# prior-variant copy of it) on the same post-1999 data used by the
# main analysis, and prints/saves posterior summaries for betaGDP
# and betaDebt so they can be compared against the baseline fit.
#
# Usage: Rscript run_variant.R <stan_file> <label> <out_csv>
# ============================================================

if (.Platform$OS.type == "windows" && nchar(Sys.which("make")) == 0) {
  rtools_candidates <- c("C:/rtools45/usr/bin", "C:/rtools44/usr/bin", "C:/Rtools/bin")
  found <- rtools_candidates[dir.exists(rtools_candidates)]
  if (length(found) > 0) {
    Sys.setenv(PATH = paste(found[1], Sys.getenv("PATH"), sep = ";"))
  } else {
    stop("Rtools not found.")
  }
}

suppressMessages({
  library(rstan)
})

rstan_options(auto_write = TRUE)
options(mc.cores = 1)

args      <- commandArgs(trailingOnly = TRUE)
stan_file <- args[1]
label     <- args[2]
out_csv   <- args[3]

data_path <- "c:/Users/joshua.ogundairo/Downloads/My Research Works/Impacts_of_GDP_and_Debt_Servicing_on_Exchange_Rate_Volatility_of_the_Nigerian_Naira__A_Bayesian_EGARCH_Approach/Analysis/Data/nigeria_monthly_stan_input_post1999.csv"

df <- read.csv(data_path, stringsAsFactors = FALSE)
df$date <- as.Date(df$date)
cat("[", label, "] Loaded:", nrow(df), "obs,", format(min(df$date)), "to", format(max(df$date)), "\n")

T_full <- nrow(df)
y      <- df$y[2:T_full]
y_lag  <- df$y[1:(T_full - 1)]
gdp    <- as.numeric(scale(df$gdp[2:T_full]))
debt   <- as.numeric(scale(df$debt[2:T_full]))
T_obs  <- length(y)
h1_init <- var(df$y)

stan_data <- list(T = T_obs, y = y, y_lag = y_lag, gdp = gdp, debt = debt, h1_init = h1_init)

init_list <- lapply(1:4, function(i) list(
  beta0 = 0.006, betaExR = 0.20, omega = -0.50, alpha = 0.15,
  gamma = 0.10, phi_raw = 0.00, betaGDP = -0.30, betaDebt = 0.30
))

cat("[", label, "] Fitting:", stan_file, "\n")
t0 <- Sys.time()
fit <- stan(
  file    = stan_file,
  data    = stan_data,
  iter    = 4000,
  warmup  = 1000,
  chains  = 4,
  seed    = 42,
  init    = init_list,
  control = list(adapt_delta = 0.92, max_treedepth = 15)
)
cat("[", label, "] Fit time:", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")

rhat_vals <- summary(fit)$summary[, "Rhat"]
rhat_bad  <- rhat_vals[!is.na(rhat_vals) & rhat_vals > 1.01]
cat("[", label, "] R-hat > 1.01 count:", length(rhat_bad), "\n")

pars_all <- c("beta0", "betaExR", "omega", "alpha", "gamma", "phi", "betaGDP", "betaDebt")
post_summary <- as.data.frame(summary(fit, pars = pars_all)$summary)
post_summary$parameter <- rownames(post_summary)
post_summary$variant   <- label
post_summary <- post_summary[, c("variant", "parameter", "mean", "sd",
                                  "2.5%", "50%", "97.5%", "n_eff", "Rhat")]
colnames(post_summary) <- c("variant", "parameter", "mean", "sd",
                             "ci_2.5", "ci_50", "ci_97.5", "n_eff", "rhat")

write.csv(post_summary, out_csv, row.names = FALSE)
cat("[", label, "] Saved:", out_csv, "\n")

cat("\n[", label, "] === betaGDP / betaDebt summary ===\n")
print(post_summary[post_summary$parameter %in% c("betaGDP", "betaDebt"), ], row.names = FALSE)
cat("\nDone:", label, "\n")
