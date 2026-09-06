# ============================================================
# Script: stationarity_tests.R
# Purpose: Unit root tests on the post-1999 monthly series
#          before Bayesian EGARCH-X estimation
#
# Tests applied (all from urca; tseries used as cross-check):
#   ADF  — Augmented Dickey-Fuller (urca::ur.df)
#   PP   — Phillips-Perron        (urca::ur.pp)
#   KPSS — Kwiatkowski-Phillips-Schmidt-Shin (urca::ur.kpss)
#
# Series tested:
#   y    — log return of official NGN/USD exchange rate
#   gdp  — Denton-Cholette interpolated monthly GDP growth
#   debt — Denton-Cholette interpolated monthly debt service (% exports)
#
# Input:  Analysis/Data/nigeria_monthly_stan_input_post1999.csv
# Output: Analysis/Data/stationarity_results.txt  (formatted table)
# ============================================================

required <- c("urca", "tseries", "dplyr", "lubridate")
new_pkg  <- required[!(required %in% installed.packages()[, "Package"])]
if (length(new_pkg)) install.packages(new_pkg)

library(urca)
library(tseries)
library(dplyr)
library(lubridate)

# ---- Locate input file -------------------------------------
input_csv <- tryCatch({
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
  file.path(script_dir, "../Data/processed/nigeria_monthly_stan_input_post1999.csv")
}, error = function(e) {
  args      <- commandArgs(trailingOnly = FALSE)
  file_flag <- grep("--file=", args, value = TRUE)
  if (length(file_flag) > 0) {
    file.path(dirname(normalizePath(sub("--file=", "", file_flag[1]))),
              "../Data/processed/nigeria_monthly_stan_input_post1999.csv")
  } else {
    file.path(getwd(), "../Data/processed/nigeria_monthly_stan_input_post1999.csv")
  }
})

if (!file.exists(input_csv)) {
  stop("Input file not found: ", input_csv,
       "\nRun data_fetch_monthly.R first to generate it.")
}

df <- read.csv(input_csv, stringsAsFactors = FALSE) %>%
  mutate(date = as.Date(date))

cat("Loaded:", nrow(df), "observations,",
    format(min(df$date)), "to", format(max(df$date)), "\n\n")

# ---- Convert to ts objects ---------------------------------
start_yr  <- year(min(df$date))
start_mo  <- month(min(df$date))

y_ts    <- ts(df$y,    start = c(start_yr, start_mo), frequency = 12)
gdp_ts  <- ts(df$gdp,  start = c(start_yr, start_mo), frequency = 12)
debt_ts <- ts(df$debt, start = c(start_yr, start_mo), frequency = 12)

series_list <- list(
  "y    (log exchange rate return)"                        = y_ts,
  "gdp  (Δ monthly GDP growth, 1st-differenced interp.)"  = gdp_ts,
  "debt (Δ monthly debt service %, 1st-differenced interp.)" = debt_ts
)

# ---- Helper: run one suite of tests and return a data frame ----
run_tests <- function(x, series_name) {

  results <- list()

  # --- ADF (urca) ---
  # Specification: "none", "drift", "trend"; lag selection by AIC
  for (spec in c("none", "drift", "trend")) {
    adf <- tryCatch(
      ur.df(x, type = spec, lags = 12, selectlags = "AIC"),
      error = function(e) NULL
    )
    if (!is.null(adf)) {
      tau_stat <- adf@teststat[1]   # tau statistic (first row = tau)
      cv_1pct  <- adf@cval[1, 1]
      cv_5pct  <- adf@cval[1, 2]
      cv_10pct <- adf@cval[1, 3]
      decision <- ifelse(tau_stat < cv_5pct, "Reject H0 (stationary)", "Fail to reject H0 (unit root)")
      results[[length(results) + 1]] <- data.frame(
        series   = series_name,
        test     = paste0("ADF (", spec, ")"),
        statistic = round(tau_stat, 4),
        cv_1pct   = round(cv_1pct,  4),
        cv_5pct   = round(cv_5pct,  4),
        cv_10pct  = round(cv_10pct, 4),
        decision  = decision,
        stringsAsFactors = FALSE
      )
    }
  }

  # --- PP (urca) ---
  for (model in c("constant", "trend")) {
    pp <- tryCatch(
      ur.pp(x, type = "Z-tau", model = model, lags = "long"),
      error = function(e) NULL
    )
    if (!is.null(pp)) {
      stat    <- pp@teststat[1]
      cv_5pct <- pp@cval[, 2]
      decision <- ifelse(stat < cv_5pct, "Reject H0 (stationary)", "Fail to reject H0 (unit root)")
      results[[length(results) + 1]] <- data.frame(
        series    = series_name,
        test      = paste0("PP (", model, ")"),
        statistic = round(stat,             4),
        cv_1pct   = round(pp@cval[, 1],    4),
        cv_5pct   = round(cv_5pct,          4),
        cv_10pct  = round(pp@cval[, 3],    4),
        decision  = decision,
        stringsAsFactors = FALSE
      )
    }
  }

  # --- KPSS (urca) ---
  # H0: series IS stationary; reject H0 means unit root
  for (trend_type in c("mu", "tau")) {
    kp <- tryCatch(
      ur.kpss(x, type = trend_type, lags = "long"),
      error = function(e) NULL
    )
    if (!is.null(kp)) {
      stat    <- kp@teststat[1]
      cv_5pct <- kp@cval[, 2]
      # KPSS: REJECT H0 (non-stationary) if stat > critical value
      decision <- ifelse(stat > cv_5pct,
                         "Reject H0 (unit root / non-stationary)",
                         "Fail to reject H0 (stationary)")
      results[[length(results) + 1]] <- data.frame(
        series    = series_name,
        test      = paste0("KPSS (", trend_type, ")"),
        statistic = round(stat,           4),
        cv_1pct   = round(kp@cval[, 1],  4),
        cv_5pct   = round(cv_5pct,         4),
        cv_10pct  = round(kp@cval[, 3],  4),
        decision  = decision,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, results)
}

# ---- Run all tests -----------------------------------------
cat("Running unit root tests...\n\n")
all_results <- do.call(rbind, lapply(
  names(series_list),
  function(nm) run_tests(series_list[[nm]], nm)
))

# ---- Pretty-print to console -------------------------------
cat(strrep("=", 80), "\n")
cat("UNIT ROOT TEST RESULTS — Post-1999 Monthly Series\n")
cat(strrep("=", 80), "\n\n")

for (nm in names(series_list)) {
  sub <- all_results[all_results$series == nm, ]
  cat(strrep("-", 70), "\n")
  cat("Series:", nm, "\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("%-22s  %10s  %8s  %8s  %8s   %s\n",
              "Test", "Statistic", "CV 1%", "CV 5%", "CV 10%", "Decision"))
  cat(strrep("-", 70), "\n")
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("%-22s  %10.4f  %8.4f  %8.4f  %8.4f   %s\n",
                sub$test[i], sub$statistic[i],
                sub$cv_1pct[i], sub$cv_5pct[i], sub$cv_10pct[i],
                sub$decision[i]))
  }
  cat("\n")
}

# ---- Cross-check with tseries ------------------------------
cat(strrep("=", 80), "\n")
cat("CROSS-CHECK — tseries package (ADF and KPSS)\n")
cat(strrep("=", 80), "\n\n")

for (nm in names(series_list)) {
  x <- series_list[[nm]]
  cat("Series:", nm, "\n")

  adf_ts <- tryCatch(
    adf.test(x, alternative = "stationary"),
    error = function(e) NULL
  )
  if (!is.null(adf_ts))
    cat(sprintf("  ADF  : statistic = %7.4f  p-value = %.4f  [%s]\n",
                adf_ts$statistic, adf_ts$p.value,
                ifelse(adf_ts$p.value < 0.05, "REJECT H0", "fail to reject")))

  kpss_ts <- tryCatch(
    kpss.test(x, null = "Level"),
    error = function(e) NULL
  )
  if (!is.null(kpss_ts))
    cat(sprintf("  KPSS : statistic = %7.4f  p-value = %.4f  [%s]\n",
                kpss_ts$statistic, kpss_ts$p.value,
                ifelse(kpss_ts$p.value < 0.05, "REJECT H0 (unit root)", "stationary")))
  cat("\n")
}

# ---- Save results to file ----------------------------------
out_path <- tryCatch({
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
  file.path(script_dir, "../Preliminary_Tests/stationarity_results.txt")
}, error = function(e) {
  args      <- commandArgs(trailingOnly = FALSE)
  file_flag <- grep("--file=", args, value = TRUE)
  if (length(file_flag) > 0) {
    file.path(dirname(normalizePath(sub("--file=", "", file_flag[1]))),
              "../Preliminary_Tests/stationarity_results.txt")
  } else {
    file.path(getwd(), "../Preliminary_Tests/stationarity_results.txt")
  }
})
if (!dir.exists(dirname(out_path))) dir.create(dirname(out_path), recursive = TRUE)

sink(out_path)
cat("UNIT ROOT TEST RESULTS\n")
cat("Project: Impacts of GDP and Debt Servicing on Exchange Rate Volatility\n")
cat("Sample : Post-1999 managed-float era (Feb 1999 – Dec 2020)\n")
cat("Date   :", format(Sys.Date()), "\n\n")

for (nm in names(series_list)) {
  sub <- all_results[all_results$series == nm, ]
  cat(strrep("-", 90), "\n")
  cat("Series:", nm, "\n")
  cat(strrep("-", 90), "\n")
  cat(sprintf("%-22s  %10s  %8s  %8s  %8s   %s\n",
              "Test", "Statistic", "CV 1%", "CV 5%", "CV 10%", "Decision"))
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("%-22s  %10.4f  %8.4f  %8.4f  %8.4f   %s\n",
                sub$test[i], sub$statistic[i],
                sub$cv_1pct[i], sub$cv_5pct[i], sub$cv_10pct[i],
                sub$decision[i]))
  }
  cat("\n")
}
sink()

cat("Results saved to:", out_path, "\n")
cat("\nInterpretation reminder:\n")
cat("  ADF / PP  : H0 = unit root. Reject H0 (p < 0.05) => stationary.\n")
cat("  KPSS      : H0 = stationary. Reject H0 (p < 0.05) => unit root.\n")
cat("  For GARCH : y (log return) must be stationary at level.\n")
cat("              Macro regressors (gdp, debt) used in variance equation;\n")
cat("              covariance-stationarity is less critical but should be noted.\n")
