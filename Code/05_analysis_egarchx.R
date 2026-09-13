# ============================================================
# analysis_egarchx.R
# Bayesian EGARCH(1,1)-X estimation for:
#   "Impacts of GDP and Debt Servicing on Exchange Rate
#    Volatility of the Nigerian Naira"
#
# Model  : egarch_x_model.stan
# Data   : ../Data/processed/nigeria_monthly_stan_input_post1999.csv
#          (Mar 1999 – Dec 2020, 262 obs, post-1999 managed-float era)
#          produced by 01_data_fetch_monthly.R + 02_make_post1999.R
#          y    = log return of official NGN/USD exchange rate
#          gdp  = first-differenced interpolated monthly GDP growth
#          debt = first-differenced interpolated monthly debt service
#
# This is step 05 of the pipeline — run after 01-04 (data fetch,
# post-1999 restriction, stationarity tests, structural break tests).
# GDP and debt service enter the VARIANCE equation (EGARCH-X),
# not the mean equation.
# ============================================================

# ---- Rtools PATH guard (Windows) ---------------------------
# Stan compiles C++ and requires Rtools. If 'make' is not found
# despite Rtools being installed, add its bin directory to PATH here.
if (.Platform$OS.type == "windows" && nchar(Sys.which("make")) == 0) {
  rtools_candidates <- c(
    "C:/rtools45/usr/bin",
    "C:/rtools44/usr/bin",
    "C:/Rtools/bin"
  )
  found <- rtools_candidates[dir.exists(rtools_candidates)]
  if (length(found) > 0) {
    Sys.setenv(PATH = paste(found[1], Sys.getenv("PATH"), sep = ";"))
    message("Added Rtools to PATH: ", found[1])
  } else {
    stop("Rtools not found. Install from https://cran.r-project.org/bin/windows/Rtools/ and restart RStudio.")
  }
}

library(rstan)
library(bayesplot)
library(loo)
library(dplyr)
library(ggplot2)

rstan_options(auto_write = TRUE)
options(mc.cores = 1)  # run chains sequentially on Windows to avoid socket overhead

# ---- Load data ---------------------------------------------
data_path <- tryCatch({
  d <- dirname(rstudioapi::getSourceEditorContext()$path)
  file.path(d, "../Data/processed/nigeria_monthly_stan_input_post1999.csv")
}, error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  flag <- grep("--file=", args, value = TRUE)
  if (length(flag) > 0)
    file.path(dirname(normalizePath(sub("--file=", "", flag[1]))),
              "../Data/processed/nigeria_monthly_stan_input_post1999.csv")
  else
    file.path(getwd(), "../Data/processed/nigeria_monthly_stan_input_post1999.csv")
})

if (!file.exists(data_path))
  stop("Data file not found: ", data_path,
       "\nRun 01_data_fetch_monthly.R and 02_make_post1999.R first.")

df <- read.csv(data_path, stringsAsFactors = FALSE)
df$date <- as.Date(df$date)
cat("Loaded:", nrow(df), "obs,",
    format(min(df$date)), "to", format(max(df$date)), "\n")

# ---- Build Stan data list ----------------------------------
# y[t] = current log return; y_lag[t] = lag for AR(1) mean
# gdp[t] and debt[t] = contemporaneous variance regressors
T_full <- nrow(df)

y      <- df$y[2:T_full]
y_lag  <- df$y[1:(T_full - 1)]

# Standardise macro regressors to mean 0, SD 1 so that betaGDP and
# betaDebt are interpretable as "change in log conditional variance
# per 1-SD change in the regressor" and informative priors are on
# a common scale regardless of the original units.
gdp_raw  <- df$gdp[2:T_full]
debt_raw <- df$debt[2:T_full]
gdp      <- as.numeric(scale(gdp_raw))
debt     <- as.numeric(scale(debt_raw))

cat("Standardisation check:\n")
cat("  gdp  — mean:", round(mean(gdp), 6), " sd:", round(sd(gdp), 6), "\n")
cat("  debt — mean:", round(mean(debt), 6), " sd:", round(sd(debt), 6), "\n\n")

T_obs  <- length(y)

# Initial conditional variance: sample variance of the full y series
h1_init <- var(df$y)

cat("Estimation sample: T =", T_obs, "observations\n")
cat("  y    : mean =", round(mean(y), 5), "  sd =", round(sd(y), 5), "\n")
cat("  gdp  : mean =", round(mean(gdp), 5),  "  sd =", round(sd(gdp), 5), "\n")
cat("  debt : mean =", round(mean(debt), 5), "  sd =", round(sd(debt), 5), "\n")
cat("  h1_init =", round(h1_init, 6), "\n\n")

stan_data <- list(
  T       = T_obs,
  y       = y,
  y_lag   = y_lag,
  gdp     = gdp,
  debt    = debt,
  h1_init = h1_init
)

# ---- Fit model ---------------------------------------------
stan_file <- tryCatch({
  d <- dirname(rstudioapi::getSourceEditorContext()$path)
  file.path(d, "egarch_x_model.stan")
}, error = function(e) file.path(getwd(), "egarch_x_model.stan"))

cat("Fitting EGARCH(1,1)-X model...\n")
cat("Stan file:", stan_file, "\n\n")

# Explicit initialization at phi = 0.50 (phi_raw = 0) to keep all chains in
# the numerically stable region and avoid the high-phi explosive mode.
init_list <- lapply(1:4, function(i) list(
  beta0    =  0.006,
  betaExR  =  0.20,
  omega    = -0.50,
  alpha    =  0.15,
  gamma    =  0.10,
  phi_raw  =  0.00,   # inv_logit(0) = 0.50
  betaGDP  = -0.30,
  betaDebt =  0.30
))

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

# ---- Posterior summaries -----------------------------------
pars_mean <- c("beta0", "betaExR")
pars_var  <- c("omega", "alpha", "gamma", "phi", "betaGDP", "betaDebt")
# phi is a transformed parameter (inv_logit(phi_raw)); phi_raw is sampled directly
pars_var_raw <- c("omega", "alpha", "gamma", "phi_raw", "betaGDP", "betaDebt")

cat("\n=== MEAN EQUATION (AR(1)) ===\n")
print(fit, pars = pars_mean, digits_summary = 4)

cat("\n=== VARIANCE EQUATION (EGARCH-X) ===\n")
# phi is reported on the original (0,1) scale via transformed parameters
print(fit, pars = pars_var, digits_summary = 4)

cat("\n(phi_raw on unconstrained scale for diagnostics)\n")
print(fit, pars = pars_var_raw, digits_summary = 4)

# ---- Convergence diagnostics -------------------------------
cat("\n=== CONVERGENCE DIAGNOSTICS ===\n")
sp <- get_sampler_params(fit, inc_warmup = FALSE)
accept_rates <- sapply(sp, function(x) mean(x[, "accept_stat__"]))
cat("Mean acceptance rates per chain:",
    paste(round(accept_rates, 3), collapse = "  "), "\n")
cat("Overall mean acceptance rate:", round(mean(accept_rates), 3), "\n")

# R-hat check
rhat_vals <- summary(fit)$summary[, "Rhat"]
rhat_bad  <- rhat_vals[!is.na(rhat_vals) & rhat_vals > 1.01]
if (length(rhat_bad) > 0) {
  cat("\nWARNING: R-hat > 1.01 for", length(rhat_bad), "parameters:\n")
  print(round(rhat_bad, 4))
} else {
  cat("All R-hat values <= 1.01 — convergence looks good.\n")
}

# ---- LOO-CV for model comparison ---------------------------
cat("\n=== LOO-CV (for model comparison) ===\n")
log_lik <- extract_log_lik(fit, parameter_name = "log_lik", merge_chains = FALSE)
loo_result <- loo(log_lik, r_eff = relative_eff(exp(log_lik)))
print(loo_result)

# ---- Output directories -------------------------------------
# Figures (PDFs) and Results (CSVs) are kept separate, both one level
# up from this script's own location (Code/), regardless of caller cwd.
script_dir <- tryCatch(dirname(rstudioapi::getSourceEditorContext()$path),
                       error = function(e) getwd())
figures_dir <- file.path(script_dir, "../Figures")
results_dir <- file.path(script_dir, "../Results")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# ---- Plots -------------------------------------------------
color_scheme_set("gray")

# Trace plots — variance equation parameters (use raw sampled params)
p_trace <- stan_trace(fit, pars = pars_var_raw)
ggsave(file.path(figures_dir, "trace_variance_eq.pdf"), plot = p_trace, width = 12, height = 8)
cat("\nSaved: trace_variance_eq.pdf\n")

# Posterior density — all parameters
p_dens <- stan_dens(fit, pars = c(pars_mean, pars_var), separate_chains = FALSE)
ggsave(file.path(figures_dir, "posterior_densities.pdf"), plot = p_dens, width = 10, height = 7)
cat("Saved: posterior_densities.pdf\n")

# Posterior intervals
p_int <- stan_plot(fit, pars = c(pars_mean, pars_var),
                   point_est = "median", show_density = TRUE)
ggsave(file.path(figures_dir, "posterior_intervals.pdf"), plot = p_int, width = 8, height = 6)
cat("Saved: posterior_intervals.pdf\n")

# ---- Extract and save conditional volatility ---------------
log_h_post <- extract(fit, pars = "log_h")$log_h  # draws x T matrix
h_median   <- apply(exp(log_h_post), 2, median)
h_lo       <- apply(exp(log_h_post), 2, quantile, 0.025)
h_hi       <- apply(exp(log_h_post), 2, quantile, 0.975)

vol_df <- data.frame(
  date     = df$date[2:T_full],
  h_median = h_median,
  h_lo     = h_lo,
  h_hi     = h_hi,
  vol_median = sqrt(h_median),
  vol_lo   = sqrt(h_lo),
  vol_hi   = sqrt(h_hi)
)

vol_path <- file.path(results_dir, "conditional_volatility.csv")
write.csv(vol_df, vol_path, row.names = FALSE)
cat("Saved conditional volatility estimates:", vol_path, "\n")

# ---- Volatility plot ---------------------------------------
pdf(file.path(figures_dir, "conditional_volatility_plot.pdf"), width = 10, height = 5)
plot(vol_df$date, vol_df$vol_median,
     type = "l", col = "black", lwd = 1.2,
     xlab = "Date", ylab = "Conditional std dev (sqrt h_t)",
     main = "Posterior median conditional volatility — NGN/USD\nBayesian EGARCH(1,1)-X, Mar 1999 – Dec 2020")
lines(vol_df$date, vol_df$vol_lo, lty = 2, col = "grey50")
lines(vol_df$date, vol_df$vol_hi, lty = 2, col = "grey50")
abline(v = as.Date("2016-06-01"), col = "red", lty = 3)
text(as.Date("2016-06-01"), max(vol_df$vol_hi) * 0.9, "Jun 2016\nCBN float",
     col = "red", cex = 0.75, pos = 4)
legend("topleft", legend = c("Posterior median", "95% credible interval"),
       lty = c(1, 2), col = c("black", "grey50"), bty = "n", cex = 0.8)
dev.off()
cat("Saved: conditional_volatility_plot.pdf\n")

# ---- Save posterior results as CSV -------------------------
# 1. Posterior summary table (parameters of interest)
pars_all <- c(pars_mean, pars_var)   # beta0, betaExR, omega, alpha, gamma, phi, betaGDP, betaDebt
post_summary <- as.data.frame(summary(fit, pars = pars_all)$summary)
post_summary$parameter <- rownames(post_summary)
post_summary <- post_summary[, c("parameter", "mean", "sd",
                                  "2.5%", "25%", "50%", "75%", "97.5%",
                                  "n_eff", "Rhat")]
colnames(post_summary) <- c("parameter", "mean", "sd",
                             "ci_2.5", "ci_25", "ci_50", "ci_75", "ci_97.5",
                             "n_eff", "rhat")

summ_path <- file.path(results_dir, "posterior_summary.csv")
write.csv(post_summary, summ_path, row.names = FALSE)
cat("\nSaved posterior summary:", summ_path, "\n")

# 2. phi_raw summary (for diagnostics)
phi_raw_summ <- as.data.frame(summary(fit, pars = "phi_raw")$summary)
phi_raw_summ$parameter <- "phi_raw"
phi_raw_summ <- phi_raw_summ[, c("parameter", "mean", "sd", "2.5%", "50%", "97.5%", "n_eff", "Rhat")]
colnames(phi_raw_summ) <- c("parameter", "mean", "sd", "ci_2.5", "ci_50", "ci_97.5", "n_eff", "rhat")
phi_raw_path <- file.path(results_dir, "phi_raw_summary.csv")
write.csv(phi_raw_summ, phi_raw_path, row.names = FALSE)

# 3. LOO-CV metrics
loo_csv <- data.frame(
  metric    = c("elpd_loo", "p_loo", "looic"),
  estimate  = c(loo_result$estimates["elpd_loo", "Estimate"],
                loo_result$estimates["p_loo",    "Estimate"],
                loo_result$estimates["looic",    "Estimate"]),
  se        = c(loo_result$estimates["elpd_loo", "SE"],
                loo_result$estimates["p_loo",    "SE"],
                loo_result$estimates["looic",    "SE"])
)
loo_path <- file.path(results_dir, "loo_cv_metrics.csv")
write.csv(loo_csv, loo_path, row.names = FALSE)
cat("Saved LOO-CV metrics:   ", loo_path, "\n")

# 4. Individual posterior draws for the 8 model parameters
#    (12000 rows × 8 columns; useful for secondary analysis or plots)
draws_mat <- as.data.frame(extract(fit, pars = pars_all))
draws_path <- file.path(results_dir, "posterior_draws.csv")
write.csv(draws_mat, draws_path, row.names = FALSE)
cat("Saved posterior draws:  ", draws_path, "\n")

cat("\nAll CSV files saved in:", results_dir, "\n")
cat("All figure PDFs saved in:", figures_dir, "\n")

cat("\nDone.\n")
