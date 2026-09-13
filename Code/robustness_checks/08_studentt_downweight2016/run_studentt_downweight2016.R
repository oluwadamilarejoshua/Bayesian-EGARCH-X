# ============================================================
# run_studentt_downweight2016.R
# Robustness refit of the ADOPTED baseline (AR(1) + Student-t,
# smoothed-GS10 data) with June & July 2016 excluded from the
# likelihood (weight = 0), to test whether betaDebt (and the other
# structural parameters) survive without the Naira-float months.
# Same iter/warmup as the baseline runs, given current system speed.
# ============================================================

if (.Platform$OS.type == "windows" && nchar(Sys.which("make")) == 0) {
  Sys.setenv(PATH = paste("C:/rtools45/usr/bin", Sys.getenv("PATH"), sep = ";"))
}

suppressMessages({ library(rstan); library(loo) })
rstan_options(auto_write = TRUE)
options(mc.cores = 1)

data_path <- "c:/Users/joshua.ogundairo/Downloads/My Research Works/Impacts_of_GDP_and_Debt_Servicing_on_Exchange_Rate_Volatility_of_the_Nigerian_Naira__A_Bayesian_EGARCH_Approach/Analysis/Data/smoothed_gs10_test/nigeria_monthly_stan_input_post1999_smoothed.csv"
stan_file <- "C:/Users/JOSHUA~1.OGU/AppData/Local/Temp/claude/c--Users-joshua-ogundairo-Downloads-My-Research-Works-Impacts-of-GDP-and-Debt-Servicing-on-Exchange-Rate-Volatility-of-the-Nigerian-Naira--A-Bayesian-EGARCH-Approach/89ca73a5-75a5-4283-8ef8-9cb30222e7e7/scratchpad/sensitivity/egarch_x_studentt_downweight2016.stan"
OUT_DIR <- "c:/Users/joshua.ogundairo/Downloads/My Research Works/Impacts_of_GDP_and_Debt_Servicing_on_Exchange_Rate_Volatility_of_the_Nigerian_Naira__A_Bayesian_EGARCH_Approach/Analysis/Bayesian-EGARCH/robustness_checks/08_studentt_downweight2016"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(data_path, stringsAsFactors = FALSE)
df$date <- as.Date(df$date)
T_full <- nrow(df)
y      <- df$y[2:T_full]
y_lag  <- df$y[1:(T_full - 1)]
gdp    <- as.numeric(scale(df$gdp[2:T_full]))
debt   <- as.numeric(scale(df$debt[2:T_full]))
T_obs  <- length(y)
h1_init <- var(df$y)
dates_for_y <- df$date[2:T_full]

weight <- rep(1, T_obs)
excl_idx <- which(dates_for_y %in% as.Date(c("2016-06-01", "2016-07-01")))
cat("Excluding indices:", excl_idx, "-> dates:", format(dates_for_y[excl_idx]), "\n")
weight[excl_idx] <- 0

stan_data <- list(T = T_obs, y = y, y_lag = y_lag, gdp = gdp, debt = debt,
                   h1_init = h1_init, weight = weight)

init_list <- lapply(1:4, function(i) list(
  beta0 = 0.006, betaExR = 0.20, omega = -0.50, alpha = 0.15,
  gamma = 0.10, phi_raw = 0.00, betaGDP = -0.30, betaDebt = 0.30, nu = 10
))

cat("Fitting Student-t downweight-2016 variant (smoothed-GS10 data)...\n")
t0 <- Sys.time()
fit <- stan(
  file = stan_file, data = stan_data, iter = 1500, warmup = 600, chains = 4,
  seed = 42, init = init_list, control = list(adapt_delta = 0.95, max_treedepth = 15)
)
cat("Fit time:", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")

rhat_vals <- summary(fit)$summary[, "Rhat"]
rhat_bad  <- rhat_vals[!is.na(rhat_vals) & rhat_vals > 1.01]
cat("R-hat > 1.01 count:", length(rhat_bad), "\n")
if (length(rhat_bad) > 0) print(round(rhat_bad, 3))

pars_all <- c("beta0", "betaExR", "omega", "alpha", "gamma", "phi", "betaGDP", "betaDebt", "nu")
cat("\n=== POSTERIOR SUMMARY (2016-06 & 2016-07 excluded from likelihood) ===\n")
print(fit, pars = pars_all, digits_summary = 4)

post_summary <- as.data.frame(summary(fit, pars = pars_all)$summary)
post_summary$parameter <- rownames(post_summary)
post_summary$variant   <- "studentt_downweight2016"
post_summary <- post_summary[, c("variant", "parameter", "mean", "sd", "2.5%", "50%", "97.5%", "n_eff", "Rhat")]
colnames(post_summary) <- c("variant", "parameter", "mean", "sd", "ci_2.5", "ci_50", "ci_97.5", "n_eff", "rhat")
write.csv(post_summary, file.path(OUT_DIR, "posterior_summary_studentt_downweight2016.csv"), row.names = FALSE)
cat("\nSaved:", file.path(OUT_DIR, "posterior_summary_studentt_downweight2016.csv"), "\n")

# LOO-CV computed only over the 259 INCLUDED observations (excluded months
# contribute log_lik under the fitted params but were given zero weight in
# training; drop them from the LOO comparison so it's apples-to-apples).
log_lik_all <- extract_log_lik(fit, parameter_name = "log_lik", merge_chains = FALSE)
incl_idx <- which(weight == 1)
log_lik_incl <- log_lik_all[, , incl_idx]
loo_dw <- loo(log_lik_incl, r_eff = relative_eff(exp(log_lik_incl)))
print(loo_dw)
loo_csv <- data.frame(metric = c("elpd_loo", "p_loo", "looic"),
                       estimate = c(loo_dw$estimates["elpd_loo","Estimate"], loo_dw$estimates["p_loo","Estimate"], loo_dw$estimates["looic","Estimate"]),
                       se = c(loo_dw$estimates["elpd_loo","SE"], loo_dw$estimates["p_loo","SE"], loo_dw$estimates["looic","SE"]))
write.csv(loo_csv, file.path(OUT_DIR, "loo_cv_metrics_studentt_downweight2016.csv"), row.names = FALSE)
pk <- loo_dw$diagnostics$pareto_k
cat("Pareto-k (259 included obs): good(<=0.7)=", sum(pk <= 0.7), " bad(0.7-1]=", sum(pk > 0.7 & pk <= 1), " very bad(>1)=", sum(pk > 1), " of", length(pk), "\n")

# Standardized residuals over the full 261 months (weight only zeroed the
# likelihood contribution during fitting; h[t] and epsilon[t] are still
# defined for every t via the transformed-parameters recursion).
post_mean <- setNames(post_summary$mean, post_summary$parameter)
h_post <- extract(fit, pars = "h")$h
h_median <- apply(h_post, 2, median)
epsilon <- y - (post_mean["beta0"] + post_mean["betaExR"] * y_lag)
z <- epsilon / sqrt(h_median)

cat("\n=== Standardized residuals, full sample incl. the 2 excluded months ===\n")
cat(sprintf("mean=%.4f sd=%.4f\n", mean(z), sd(z)))
idx_sorted <- order(-abs(z))
cat("Top 6 |z|:\n")
for (i in idx_sorted[1:6]) cat(sprintf("  %s  z=%.2f\n", format(dates_for_y[i]), z[i]))

write.csv(data.frame(date = dates_for_y, epsilon = epsilon, h_median = h_median, z = z, weight = weight),
          file.path(OUT_DIR, "standardized_residuals_studentt_downweight2016.csv"), row.names = FALSE)

cat("\nAll outputs saved in:", OUT_DIR, "\n")
cat("Done: studentt_downweight2016\n")
