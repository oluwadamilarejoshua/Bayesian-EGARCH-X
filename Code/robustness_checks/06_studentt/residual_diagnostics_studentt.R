# ============================================================
# residual_diagnostics_studentt.R
# Standalone, reproducible diagnostics on the ADOPTED model's
# (06_studentt) standardized residuals -- no refitting required,
# reads the already-saved residuals CSV.
#
# Produces:
#   (1) A genuine ARCH-LM test (Engle 1982: regress z^2 on its own
#       lags up to 12, LM = n*R^2 ~ chisq(12)) -- this was cited in
#       the manuscript but never actually computed for this model;
#       this script closes that gap.
#   (2) Ljung-Box on z and z^2, unadjusted (for reference).
#   (3) Ljung-Box on z and z^2 after winsorizing the 5 most extreme
#       |z| observations to the 99th percentile of the remainder --
#       reproduces the p=0.112 figure reported in the manuscript's
#       Model Verification section, which previously existed only
#       as an ad hoc calculation with no saved script.
# ============================================================

data_path <- "c:/Users/joshua.ogundairo/Downloads/My Research Works/Impacts_of_GDP_and_Debt_Servicing_on_Exchange_Rate_Volatility_of_the_Nigerian_Naira__A_Bayesian_EGARCH_Approach/Analysis/Bayesian-EGARCH/robustness_checks/06_studentt/standardized_residuals_studentt.csv"
OUT_DIR <- "c:/Users/joshua.ogundairo/Downloads/My Research Works/Impacts_of_GDP_and_Debt_Servicing_on_Exchange_Rate_Volatility_of_the_Nigerian_Naira__A_Bayesian_EGARCH_Approach/Analysis/Bayesian-EGARCH/robustness_checks/06_studentt"

df <- read.csv(data_path, stringsAsFactors = FALSE)
z <- df$z
n <- length(z)
cat("Loaded", n, "standardized residuals from:", data_path, "\n\n")

# ---------------------------------------------------------------
# (1) ARCH-LM test (Engle 1982), 12 lags, on the RAW (unadjusted) z
# ---------------------------------------------------------------
arch_lm_test <- function(resid, lags = 12) {
  e2 <- resid^2
  n  <- length(e2)
  X  <- sapply(1:lags, function(k) c(rep(NA, k), e2[1:(n - k)]))
  dat <- data.frame(y = e2, X)
  dat <- na.omit(dat)
  fit <- lm(y ~ ., data = dat)
  r2  <- summary(fit)$r.squared
  n_eff <- nrow(dat)
  LM  <- n_eff * r2
  p   <- 1 - pchisq(LM, df = lags)
  list(LM = LM, p = p, df = lags, n_eff = n_eff)
}

arch_raw <- arch_lm_test(z, lags = 12)
cat("=== ARCH-LM test (12 lags), raw z ===\n")
cat(sprintf("LM = %.4f, df = %d, p = %.4f\n\n", arch_raw$LM, arch_raw$df, arch_raw$p))

# ---------------------------------------------------------------
# (2) Ljung-Box on raw z and z^2 (unadjusted, for reference)
# ---------------------------------------------------------------
cat("=== Ljung-Box, raw z (unadjusted) ===\n")
for (lag in c(6, 12, 24)) {
  lb <- Box.test(z, lag = lag, type = "Ljung-Box", fitdf = 2)
  cat(sprintf("lag=%2d: Q=%.4f p=%.4f\n", lag, lb$statistic, lb$p.value))
}
cat("\n=== Ljung-Box, raw z^2 (unadjusted) ===\n")
for (lag in c(12, 24)) {
  lb2 <- Box.test(z^2, lag = lag, type = "Ljung-Box", fitdf = 7)
  cat(sprintf("lag=%2d: Q=%.4f p=%.4f\n", lag, lb2$statistic, lb2$p.value))
}

# ---------------------------------------------------------------
# (3) Winsorized versions: cap the 5 most extreme |z| at the 99th
#     percentile of the remaining 256 observations, sign preserved.
# ---------------------------------------------------------------
k <- 5
idx_sorted <- order(-abs(z))
top_idx <- idx_sorted[1:k]
cap <- quantile(abs(z[-top_idx]), 0.99)
z_wins <- z
z_wins[top_idx] <- sign(z[top_idx]) * cap

cat("\n=== Top", k, "|z| observations, capped at +/-", round(cap, 3), "===\n")
for (i in top_idx) cat(sprintf("  %s  z=%.2f -> %.2f\n", df$date[i], z[i], z_wins[i]))

cat("\n=== ARCH-LM test (12 lags), winsorized z ===\n")
arch_wins <- arch_lm_test(z_wins, lags = 12)
cat(sprintf("LM = %.4f, df = %d, p = %.4f\n\n", arch_wins$LM, arch_wins$df, arch_wins$p))

cat("=== Ljung-Box, winsorized z ===\n")
lb_w  <- Box.test(z_wins,    lag = 12, type = "Ljung-Box", fitdf = 2)
cat(sprintf("lag=12: Q=%.4f p=%.4f\n", lb_w$statistic, lb_w$p.value))
cat("=== Ljung-Box, winsorized z^2 ===\n")
lb2_w <- Box.test(z_wins^2, lag = 12, type = "Ljung-Box", fitdf = 7)
cat(sprintf("lag=12: Q=%.4f p=%.4f\n", lb2_w$statistic, lb2_w$p.value))

# ---------------------------------------------------------------
# Save everything to a CSV for the manuscript / repo
# ---------------------------------------------------------------
out <- data.frame(
  test      = c("ARCH-LM (raw)", "ARCH-LM (winsorized)",
                 "Ljung-Box z (raw, lag12)", "Ljung-Box z^2 (raw, lag12)",
                 "Ljung-Box z (winsorized, lag12)", "Ljung-Box z^2 (winsorized, lag12)"),
  statistic = c(arch_raw$LM, arch_wins$LM,
                Box.test(z, 12, "Ljung-Box", fitdf=2)$statistic,
                Box.test(z^2, 12, "Ljung-Box", fitdf=7)$statistic,
                lb_w$statistic, lb2_w$statistic),
  p_value   = c(arch_raw$p, arch_wins$p,
                Box.test(z, 12, "Ljung-Box", fitdf=2)$p.value,
                Box.test(z^2, 12, "Ljung-Box", fitdf=7)$p.value,
                lb_w$p.value, lb2_w$p.value)
)
write.csv(out, file.path(OUT_DIR, "residual_diagnostics_full.csv"), row.names = FALSE)
cat("\nSaved:", file.path(OUT_DIR, "residual_diagnostics_full.csv"), "\n")
