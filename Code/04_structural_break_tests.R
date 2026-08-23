# ============================================================
# Script: structural_break_tests.R
# Purpose: Detect structural breaks in the post-1999 log
#          return series before Bayesian EGARCH-X estimation
#
# Tests applied:
#   Zivot-Andrews (1992) — unit root test allowing for a single
#     structural break in intercept and/or trend (urca::ur.za)
#   Bai-Perron (1998, 2003) — multiple structural break test
#     (strucchange::breakpoints)
#   Chow test at a priori dates — June 2016 CBN float
#     (strucchange::sctest)
#
# Input:  Analysis/Data/nigeria_monthly_stan_input_post1999.csv
# Output: Analysis/Data/structural_break_results.txt
#         Analysis/Data/structural_break_plot.pdf
# ============================================================

required <- c("urca", "strucchange", "dplyr", "lubridate", "zoo")
new_pkg  <- required[!(required %in% installed.packages()[, "Package"])]
if (length(new_pkg)) install.packages(new_pkg)

library(urca)
library(strucchange)
library(dplyr)
library(lubridate)
library(zoo)

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
       "\nRun data_fetch_monthly.R first to generate the post-1999 dataset.")
}

df <- read.csv(input_csv, stringsAsFactors = FALSE) %>%
  mutate(date = as.Date(date))

cat("Loaded:", nrow(df), "observations,",
    format(min(df$date)), "to", format(max(df$date)), "\n\n")

start_yr <- year(min(df$date))
start_mo <- month(min(df$date))
y_ts     <- ts(df$y, start = c(start_yr, start_mo), frequency = 12)

# ============================================================
# TEST 1: Zivot-Andrews (1992)
# H0: unit root with no structural break
# Allows break in: intercept ("intercept"), trend ("trend"), both ("both")
# ============================================================
cat(strrep("=", 70), "\n")
cat("TEST 1: Zivot-Andrews Unit Root Test (with structural break)\n")
cat(strrep("=", 70), "\n\n")

za_results <- list()
for (model_type in c("intercept", "trend", "both")) {
  za <- tryCatch(
    ur.za(y_ts, model = model_type, lag = 12),
    error = function(e) {
      cat("  ur.za(model='", model_type, "') failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (!is.null(za)) {
    stat    <- za@teststat
    cv_5pct <- za@cval[2]   # 5% critical value
    break_pt <- za@bpoint

    # Convert breakpoint index to date
    break_date <- df$date[break_pt]

    decision <- ifelse(stat < cv_5pct,
                       "Reject H0 (stationary WITH break)",
                       "Fail to reject H0 (unit root)")

    cat(sprintf("Model: %-12s | Statistic: %8.4f | CV 5%%: %7.4f | ",
                model_type, stat, cv_5pct))
    cat(sprintf("Break: %s | %s\n", format(break_date), decision))

    za_results[[model_type]] <- list(
      stat = stat, cv_5pct = cv_5pct,
      break_date = break_date, decision = decision
    )
  }
}

# ============================================================
# TEST 2: Bai-Perron (1998, 2003) — multiple break detection
# Tests for up to 5 structural breaks in the mean of y
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("TEST 2: Bai-Perron Multiple Structural Break Test\n")
cat(strrep("=", 70), "\n\n")

bp <- tryCatch(
  breakpoints(y_ts ~ 1, h = 0.10, breaks = 5),
  error = function(e) {
    cat("breakpoints() failed:", conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(bp)) {
  cat("BIC-selected number of breaks:", bp$breakpoints[which.min(BIC(bp))], "\n")
  cat("\nBIC values by number of breaks:\n")
  bic_vals <- BIC(bp)
  print(round(bic_vals, 2))

  best_breaks <- bp$breakpoints[which.min(bic_vals)]
  if (!is.na(best_breaks) && best_breaks > 0) {
    # Fit the selected model and extract break dates
    bp_fit  <- breakpoints(bp, breaks = best_breaks)
    bp_idx  <- bp_fit$breakpoints
    bp_dates <- df$date[bp_idx]
    cat("\nDetected break date(s):\n")
    for (d in format(bp_dates)) cat(" ", d, "\n")
  } else {
    cat("BIC selects 0 breaks (no structural break detected).\n")
  }
}

# ============================================================
# TEST 3: Chow test at a priori break date — June 2016
# The CBN floated the naira on June 20, 2016; the July 2016
# observation captures the first full month post-float.
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("TEST 3: Chow Test at a priori break — June/July 2016\n")
cat(strrep("=", 70), "\n\n")

# Index of the June 2016 observation in the post-1999 sample
break_2016_idx <- which(df$date == as.Date("2016-06-01"))
if (length(break_2016_idx) == 0) {
  break_2016_idx <- which(df$date == as.Date("2016-07-01"))
}

if (length(break_2016_idx) > 0) {
  chow <- tryCatch(
    sctest(y_ts ~ 1, type = "Chow", point = break_2016_idx),
    error = function(e) NULL
  )
  if (!is.null(chow)) {
    cat(sprintf("Chow F-statistic: %.4f   p-value: %.4f   [%s]\n",
                chow$statistic, chow$p.value,
                ifelse(chow$p.value < 0.05,
                       "BREAK CONFIRMED at 5%",
                       "No significant break at 5%")))
  }

  # CUSUM test (visual and formal)
  cusum <- tryCatch(efp(y_ts ~ 1, type = "OLS-CUSUM"), error = function(e) NULL)
  if (!is.null(cusum)) {
    cusum_test <- sctest(cusum)
    cat(sprintf("CUSUM test: statistic = %.4f   p-value = %.4f   [%s]\n",
                cusum_test$statistic, cusum_test$p.value,
                ifelse(cusum_test$p.value < 0.05,
                       "Parameter instability detected",
                       "No parameter instability")))
  }
} else {
  cat("Note: June 2016 not in sample — cannot run Chow test.\n")
}

# ============================================================
# Plot: series with detected break dates
# ============================================================
out_plot <- tryCatch({
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
  file.path(script_dir, "../Preliminary_Tests/structural_break_plot.pdf")
}, error = function(e) file.path(getwd(), "../Preliminary_Tests/structural_break_plot.pdf"))
if (!dir.exists(dirname(out_plot))) dir.create(dirname(out_plot), recursive = TRUE)

tryCatch({
  pdf(out_plot, width = 10, height = 5)

  plot(df$date, df$y,
       type = "l", col = "steelblue", lwd = 0.8,
       xlab = "Date", ylab = "Log return (y)",
       main = "Monthly log return — NGN/USD (Feb 1999–Dec 2020)\nwith detected structural breaks")

  # Mark the a priori 2016 break
  abline(v = as.Date("2016-06-01"), col = "red",  lty = 2, lwd = 1.5)
  text(as.Date("2016-06-01"), max(df$y) * 0.85, "Jun 2016\n(CBN float)",
       col = "red", cex = 0.75, pos = 4)

  # Mark Bai-Perron breaks (if any)
  if (!is.null(bp) && !is.na(best_breaks) && best_breaks > 0) {
    abline(v = bp_dates, col = "darkgreen", lty = 3, lwd = 1.2)
    for (d in bp_dates)
      text(d, max(df$y) * 0.70, format(d, "%Y-%m"),
           col = "darkgreen", cex = 0.7, pos = 4)
  }

  # Mark Zivot-Andrews breaks
  for (nm in names(za_results)) {
    abline(v = za_results[[nm]]$break_date, col = "orange", lty = 4, lwd = 1)
  }

  legend("topleft",
         legend = c("y (log return)", "A priori 2016 break",
                    "Bai-Perron breaks", "ZA break point"),
         col    = c("steelblue", "red", "darkgreen", "orange"),
         lty    = c(1, 2, 3, 4), cex = 0.75, bty = "n")

  dev.off()
  cat("\nPlot saved to:", out_plot, "\n")
}, error = function(e) {
  cat("\nPlot generation failed:", conditionMessage(e), "\n")
})

# ============================================================
# Save text results
# ============================================================
out_txt <- tryCatch({
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
  file.path(script_dir, "../Preliminary_Tests/structural_break_results.txt")
}, error = function(e) file.path(getwd(), "../Preliminary_Tests/structural_break_results.txt"))
if (!dir.exists(dirname(out_txt))) dir.create(dirname(out_txt), recursive = TRUE)

sink(out_txt)
cat("STRUCTURAL BREAK TEST RESULTS\n")
cat("Project: Impacts of GDP and Debt Servicing on Exchange Rate Volatility\n")
cat("Series : y (log return, NGN/USD official rate, post-1999)\n")
cat("Sample : Feb 1999 – Dec 2020 (", nrow(df), "obs)\n")
cat("Date   :", format(Sys.Date()), "\n\n")

cat(strrep("=", 70), "\n")
cat("Zivot-Andrews Results:\n")
cat(strrep("=", 70), "\n")
for (nm in names(za_results)) {
  cat(sprintf("  Model: %-12s | Stat: %8.4f | CV 5%%: %7.4f | Break: %s | %s\n",
              nm,
              za_results[[nm]]$stat,
              za_results[[nm]]$cv_5pct,
              format(za_results[[nm]]$break_date),
              za_results[[nm]]$decision))
}

cat("\n")
cat(strrep("=", 70), "\n")
cat("Bai-Perron Results:\n")
cat(strrep("=", 70), "\n")
if (!is.null(bp)) {
  cat("BIC-selected breaks:", best_breaks, "\n")
  if (!is.na(best_breaks) && best_breaks > 0)
    cat("Break dates:", paste(format(bp_dates), collapse = ", "), "\n")
}

cat("\n")
cat(strrep("=", 70), "\n")
cat("Chow Test at June 2016:\n")
cat(strrep("=", 70), "\n")
if (exists("chow") && !is.null(chow)) {
  cat(sprintf("F-statistic: %.4f  p-value: %.4f\n", chow$statistic, chow$p.value))
}
sink()

cat("Results saved to:", out_txt, "\n")
cat("\nKey question for the model:\n")
cat("  If June 2016 break is confirmed (Chow p < 0.05), add a dummy D_2016\n")
cat("  to the EGARCH variance equation: D_2016 = 1 from June 2016 onward.\n")
cat("  This is the EGARCH-X-break specification defended in the audit report.\n")
