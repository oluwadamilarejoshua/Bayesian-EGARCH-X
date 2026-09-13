# ============================================================
# run_ar2_studentt.R
# Refits EGARCH-X on the smoothed-GS10 data with:
#   - AR(2) mean equation (was AR(1))
#   - Student-t mean-equation innovations (was Gaussian)
# Extending to AR(2) costs one more observation at the start of
# the sample (needs y[t-2]), so T drops from 261 to 260.
# ============================================================

if (.Platform$OS.type == "windows" && nchar(Sys.which("make")) == 0) {
  Sys.setenv(PATH = paste("C:/rtools45/usr/bin", Sys.getenv("PATH"), sep = ";"))
}

suppressMessages({ library(rstan); library(loo) })
rstan_options(auto_write = TRUE)
options(mc.cores = 1)

data_path <- "c:/Users/joshua.ogundairo/Downloads/My Research Works/Impacts_of_GDP_and_Debt_Servicing_on_Exchange_Rate_Volatility_of_the_Nigerian_Naira__A_Bayesian_EGARCH_Approach/Analysis/Data/smoothed_gs10_test/nigeria_monthly_stan_input_post1999_smoothed.csv"
stan_file <- "C:/Users/JOSHUA~1.OGU/AppData/Local/Temp/claude/c--Users-joshua-ogundairo-Downloads-My-Research-Works-Impacts-of-GDP-and-Debt-Servicing-on-Exchange-Rate-Volatility-of-the-Nigerian-Naira--A-Bayesian-EGARCH-Approach/89ca73a5-75a5-4283-8ef8-9cb30222e7e7/scratchpad/sensitivity/egarch_x_ar2_studentt.stan"
OUT_DIR <- "c:/Users/joshua.ogundairo/Downloads/My Research Works/Impacts_of_GDP_and_Debt_Servicing_on_Exchange_Rate_Volatility_of_the_Nigerian_Naira__A_Bayesian_EGARCH_Approach/Analysis/Bayesian-EGARCH/robustness_checks/07_ar2_studentt"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(data_path, stringsAsFactors = FALSE)
df$date <- as.Date(df$date)
T_full <- nrow(df)

y      <- df$y[3:T_full]
y_lag  <- df$y[2:(T_full - 1)]
y_lag2 <- df$y[1:(T_full - 2)]
gdp    <- as.numeric(scale(df$gdp[3:T_full]))
debt   <- as.numeric(scale(df$debt[3:T_full]))
dates  <- df$date[3:T_full]
T_obs  <- length(y)
h1_init <- var(df$y)

cat("T_obs =", T_obs, " (was 261 for the AR(1) model; -1 for the extra lag)\n")

stan_data <- list(T = T_obs, y = y, y_lag = y_lag, y_lag2 = y_lag2, gdp = gdp, debt = debt, h1_init = h1_init)

init_list <- lapply(1:4, function(i) list(
  beta0 = 0.006, betaExR = 0.20, beta2 = 0.0, omega = -0.50, alpha = 0.15,
  gamma = 0.10, phi_raw = 0.00, betaGDP = -0.30, betaDebt = 0.30, nu = 10
))

cat("Fitting AR(2) + Student-t EGARCH-X (smoothed-GS10 data)...\n")
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

pars_all <- c("beta0", "betaExR", "beta2", "omega", "alpha", "gamma", "phi", "betaGDP", "betaDebt", "nu")
cat("\n=== POSTERIOR SUMMARY ===\n")
print(fit, pars = pars_all, digits_summary = 4)

post_summary <- as.data.frame(summary(fit, pars = pars_all)$summary)
post_summary$parameter <- rownames(post_summary)
post_summary <- post_summary[, c("parameter", "mean", "sd", "2.5%", "50%", "97.5%", "n_eff", "Rhat")]
colnames(post_summary) <- c("parameter", "mean", "sd", "ci_2.5", "ci_50", "ci_97.5", "n_eff", "rhat")
write.csv(post_summary, file.path(OUT_DIR, "posterior_summary_ar2_studentt.csv"), row.names = FALSE)

csum <- summary(fit, pars = c("omega","alpha","gamma","phi","betaGDP","betaDebt","nu","beta2"))$c_summary
rows <- list()
for (p in c("omega","alpha","gamma","phi","betaGDP","betaDebt","nu","beta2")) for (ch in 1:4)
  rows[[length(rows)+1]] <- data.frame(parameter = p, chain = ch, mean = csum[p,"mean",ch],
                                        q2.5 = csum[p,"2.5%",ch], q97.5 = csum[p,"97.5%",ch])
write.csv(do.call(rbind, rows), file.path(OUT_DIR, "per_chain_summary_ar2_studentt.csv"), row.names = FALSE)

log_lik_t <- extract_log_lik(fit, parameter_name = "log_lik", merge_chains = FALSE)
loo_t <- loo(log_lik_t, r_eff = relative_eff(exp(log_lik_t)))
print(loo_t)
loo_csv <- data.frame(metric = c("elpd_loo", "p_loo", "looic"),
                       estimate = c(loo_t$estimates["elpd_loo","Estimate"], loo_t$estimates["p_loo","Estimate"], loo_t$estimates["looic","Estimate"]),
                       se = c(loo_t$estimates["elpd_loo","SE"], loo_t$estimates["p_loo","SE"], loo_t$estimates["looic","SE"]))
write.csv(loo_csv, file.path(OUT_DIR, "loo_cv_metrics_ar2_studentt.csv"), row.names = FALSE)
pk <- loo_t$diagnostics$pareto_k
cat("Pareto-k: good(<=0.7)=", sum(pk <= 0.7), " bad(0.7-1]=", sum(pk > 0.7 & pk <= 1), " very bad(>1)=", sum(pk > 1), " of", length(pk), "\n")

post_mean <- setNames(post_summary$mean, post_summary$parameter)
h_post <- extract(fit, pars = "h")$h
h_median <- apply(h_post, 2, median)
epsilon <- y - (post_mean["beta0"] + post_mean["betaExR"] * y_lag + post_mean["beta2"] * y_lag2)
z <- epsilon / sqrt(h_median)

cat("\n=== Re-diagnostics on AR(2)+Student-t model's standardized residuals ===\n")
skew <- mean((z-mean(z))^3)/sd(z)^3
kurt <- mean((z-mean(z))^4)/sd(z)^4
n <- length(z)
jb_stat <- (n/6)*(skew^2 + ((kurt-3)^2)/4)
jb_p <- 1 - pchisq(jb_stat, df=2)
cat(sprintf("skewness=%.3f excess kurtosis=%.3f  JB=%.2f p=%.4f [%s]\n",
            skew, kurt-3, jb_stat, jb_p, ifelse(jb_p<0.05,"REJECT normality","fail to reject normality")))
for (lag in c(6,12,24)) {
  lb <- Box.test(z, lag=lag, type="Ljung-Box", fitdf=3)
  cat(sprintf("Ljung-Box on z, lag=%2d: Q=%.4f p=%.4f [%s]\n", lag, lb$statistic, lb$p.value,
              ifelse(lb$p.value<0.05,"REJECT","clean")))
}
lb2 <- Box.test(z^2, lag=12, type="Ljung-Box", fitdf=7)
cat(sprintf("Ljung-Box on z^2, lag=12: Q=%.4f p=%.4f\n", lb2$statistic, lb2$p.value))

write.csv(data.frame(date=dates, epsilon=epsilon, h_median=h_median, z=z),
          file.path(OUT_DIR, "standardized_residuals_ar2_studentt.csv"), row.names=FALSE)

cat("\nAll outputs saved in:", OUT_DIR, "\n")
cat("Done: ar2_studentt\n")
