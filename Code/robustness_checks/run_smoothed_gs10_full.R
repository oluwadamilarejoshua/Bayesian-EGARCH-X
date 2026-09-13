# ============================================================
# run_smoothed_gs10_full.R
# Refits the baseline EGARCH-X model (tight, theory-motivated priors
# -- unchanged from the main model) on the SMOOTHED-GS10-indicator
# data, with full diagnostic plots (trace, posterior densities,
# intervals, conditional volatility), saved permanently into
# robustness_checks/05_smoothed_gs10/.
# ============================================================

if (.Platform$OS.type == "windows" && nchar(Sys.which("make")) == 0) {
  Sys.setenv(PATH = paste("C:/rtools45/usr/bin", Sys.getenv("PATH"), sep = ";"))
}

suppressMessages({ library(rstan); library(bayesplot); library(loo); library(ggplot2) })
rstan_options(auto_write = TRUE)
options(mc.cores = 1)

data_path <- "c:/Users/joshua.ogundairo/Downloads/My Research Works/Impacts_of_GDP_and_Debt_Servicing_on_Exchange_Rate_Volatility_of_the_Nigerian_Naira__A_Bayesian_EGARCH_Approach/Analysis/Data/smoothed_gs10_test/nigeria_monthly_stan_input_post1999_smoothed.csv"
stan_file <- "c:/Users/joshua.ogundairo/Downloads/My Research Works/Impacts_of_GDP_and_Debt_Servicing_on_Exchange_Rate_Volatility_of_the_Nigerian_Naira__A_Bayesian_EGARCH_Approach/Analysis/Bayesian-EGARCH/egarch_x_model.stan"
OUT_DIR <- "c:/Users/joshua.ogundairo/Downloads/My Research Works/Impacts_of_GDP_and_Debt_Servicing_on_Exchange_Rate_Volatility_of_the_Nigerian_Naira__A_Bayesian_EGARCH_Approach/Analysis/Bayesian-EGARCH/robustness_checks/05_smoothed_gs10"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(data_path, stringsAsFactors = FALSE)
df$date <- as.Date(df$date)
cat("Loaded:", nrow(df), "obs,", format(min(df$date)), "to", format(max(df$date)), "\n")

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

cat("Fitting EGARCH-X (smoothed GS10 indicator) with baseline tight priors...\n")
t0 <- Sys.time()
fit <- stan(
  file = stan_file, data = stan_data, iter = 4000, warmup = 1000, chains = 4,
  seed = 42, init = init_list, control = list(adapt_delta = 0.92, max_treedepth = 15)
)
cat("Fit time:", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")

rhat_vals <- summary(fit)$summary[, "Rhat"]
rhat_bad  <- rhat_vals[!is.na(rhat_vals) & rhat_vals > 1.01]
cat("R-hat > 1.01 count:", length(rhat_bad), "\n")
if (length(rhat_bad) > 0) print(round(rhat_bad, 3))

pars_mean <- c("beta0", "betaExR")
pars_var  <- c("omega", "alpha", "gamma", "phi", "betaGDP", "betaDebt")
pars_all  <- c(pars_mean, pars_var)

cat("\n=== POSTERIOR SUMMARY ===\n")
print(fit, pars = pars_all, digits_summary = 4)

post_summary <- as.data.frame(summary(fit, pars = pars_all)$summary)
post_summary$parameter <- rownames(post_summary)
post_summary <- post_summary[, c("parameter", "mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%", "n_eff", "Rhat")]
colnames(post_summary) <- c("parameter", "mean", "sd", "ci_2.5", "ci_25", "ci_50", "ci_75", "ci_97.5", "n_eff", "rhat")
write.csv(post_summary, file.path(OUT_DIR, "posterior_summary_smoothed.csv"), row.names = FALSE)
cat("\nSaved:", file.path(OUT_DIR, "posterior_summary_smoothed.csv"), "\n")

# ---- LOO-CV ----
log_lik <- extract_log_lik(fit, parameter_name = "log_lik", merge_chains = FALSE)
loo_result <- loo(log_lik, r_eff = relative_eff(exp(log_lik)))
print(loo_result)
loo_csv <- data.frame(metric = c("elpd_loo", "p_loo", "looic"),
                       estimate = c(loo_result$estimates["elpd_loo","Estimate"], loo_result$estimates["p_loo","Estimate"], loo_result$estimates["looic","Estimate"]),
                       se = c(loo_result$estimates["elpd_loo","SE"], loo_result$estimates["p_loo","SE"], loo_result$estimates["looic","SE"]))
write.csv(loo_csv, file.path(OUT_DIR, "loo_cv_metrics_smoothed.csv"), row.names = FALSE)

# ---- Plots ----
color_scheme_set("gray")
pars_var_raw <- c("omega", "alpha", "gamma", "phi_raw", "betaGDP", "betaDebt")

p_trace <- stan_trace(fit, pars = pars_var_raw)
ggsave(file.path(OUT_DIR, "trace_variance_eq_smoothed.pdf"), p_trace, width = 12, height = 8)

p_dens <- stan_dens(fit, pars = pars_all, separate_chains = FALSE)
ggsave(file.path(OUT_DIR, "posterior_densities_smoothed.pdf"), p_dens, width = 10, height = 7)

p_int <- stan_plot(fit, pars = pars_all, point_est = "median", show_density = TRUE)
ggsave(file.path(OUT_DIR, "posterior_intervals_smoothed.pdf"), p_int, width = 8, height = 6)

# posterior densities SPLIT BY CHAIN too -- important given the bimodality history
p_dens_bychain <- stan_dens(fit, pars = pars_var, separate_chains = TRUE)
ggsave(file.path(OUT_DIR, "posterior_densities_by_chain_smoothed.pdf"), p_dens_bychain, width = 10, height = 7)

log_h_post <- extract(fit, pars = "log_h")$log_h
h_median <- apply(exp(log_h_post), 2, median)
h_lo     <- apply(exp(log_h_post), 2, quantile, 0.025)
h_hi     <- apply(exp(log_h_post), 2, quantile, 0.975)
vol_df <- data.frame(date = df$date[2:T_full], h_median = h_median, h_lo = h_lo, h_hi = h_hi,
                      vol_median = sqrt(h_median), vol_lo = sqrt(h_lo), vol_hi = sqrt(h_hi))
write.csv(vol_df, file.path(OUT_DIR, "conditional_volatility_smoothed.csv"), row.names = FALSE)

pdf(file.path(OUT_DIR, "conditional_volatility_plot_smoothed.pdf"), width = 10, height = 5)
plot(vol_df$date, vol_df$vol_median, type = "l", col = "black", lwd = 1.2,
     xlab = "Date", ylab = "Conditional std dev (sqrt h_t)",
     main = "Posterior median conditional volatility -- smoothed GS10 indicator")
lines(vol_df$date, vol_df$vol_lo, lty = 2, col = "grey50")
lines(vol_df$date, vol_df$vol_hi, lty = 2, col = "grey50")
abline(v = as.Date("2016-06-01"), col = "red", lty = 3)
text(as.Date("2016-06-01"), max(vol_df$vol_hi) * 0.9, "Jun 2016\nCBN float", col = "red", cex = 0.75, pos = 4)
legend("topleft", legend = c("Posterior median", "95% credible interval"), lty = c(1,2), col = c("black","grey50"), bty = "n", cex = 0.8)
dev.off()

# per-chain summary too, for direct comparability with the bimodality diagnostics
csum <- summary(fit, pars = pars_var)$c_summary
rows <- list()
for (p in pars_var) for (ch in 1:4)
  rows[[length(rows)+1]] <- data.frame(parameter = p, chain = ch, mean = csum[p,"mean",ch],
                                        sd = csum[p,"sd",ch], q2.5 = csum[p,"2.5%",ch], q97.5 = csum[p,"97.5%",ch])
write.csv(do.call(rbind, rows), file.path(OUT_DIR, "per_chain_summary_smoothed.csv"), row.names = FALSE)

cat("\nAll outputs saved in:", OUT_DIR, "\n")
cat("Done: smoothed_gs10_full\n")
