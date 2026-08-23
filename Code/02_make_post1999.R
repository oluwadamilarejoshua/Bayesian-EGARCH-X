# ============================================================
# make_post1999.R
# Reads the already-generated full monthly dataset, filters to
# the post-1999 managed-float era, first-differences the macro
# regressors (gdp, debt) to achieve stationarity, and writes
# the final Stan-ready input files.
#
# Background: stationarity_tests.R confirmed that y is I(0)
# but the Denton-Cholette interpolated gdp and debt series are
# I(1) (non-stationary). Their first differences are I(0) and
# have clean economic interpretation:
#   d_gdp  = monthly change in GDP growth momentum
#   d_debt = monthly change in debt service burden (% exports)
#
# Run this before stationarity_tests.R and structural_break_tests.R.
# No internet connection required.
# ============================================================

library(dplyr)
library(urca)      # for ADF confirmation on differenced series

# ---- Locate the full merged file ---------------------------
full_csv <- tryCatch({
  d <- dirname(rstudioapi::getSourceEditorContext()$path)
  file.path(d, "../Data/processed/nigeria_monthly_merged.csv")
}, error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  flag <- grep("--file=", args, value = TRUE)
  if (length(flag) > 0)
    file.path(dirname(normalizePath(sub("--file=", "", flag[1]))),
              "../Data/processed/nigeria_monthly_merged.csv")
  else
    file.path(getwd(), "../Data/processed/nigeria_monthly_merged.csv")
})

if (!file.exists(full_csv))
  stop("Cannot find nigeria_monthly_merged.csv.\n",
       "Run data_fetch_monthly.R first to generate the full dataset.")

cat("Reading:", full_csv, "\n")
df <- read.csv(full_csv, stringsAsFactors = FALSE)
df$date <- as.Date(df$date)
cat("Full sample:", nrow(df), "obs (",
    format(min(df$date)), "to", format(max(df$date)), ")\n\n")

# ---- Filter to post-1999 managed-float era -----------------
df_post <- df %>%
  filter(date >= as.Date("1999-02-01")) %>%
  arrange(date)

n     <- nrow(df_post)
y     <- df_post$log_exr_return
zeros <- sum(y == 0)

cat("Post-1999 sample:", n, "obs (",
    format(min(df_post$date)), "to", format(max(df_post$date)), ")\n")
cat("Exact zeros in y:", zeros, "of", n,
    paste0("(", round(100 * zeros / n, 1), "%)\n\n"))

# ---- First-difference the macro regressors -----------------
# diff() returns n-1 values; pad the first row with NA then drop it.
df_post <- df_post %>%
  mutate(
    d_gdp  = c(NA, diff(gdp_growth_m)),
    d_debt = c(NA, diff(debt_service_m))
  ) %>%
  filter(!is.na(d_gdp))         # drops the first obs (Feb 1999)

cat("After differencing, final sample:", nrow(df_post), "obs (",
    format(min(df_post$date)), "to", format(max(df_post$date)), ")\n\n")

# ---- Confirm differenced series are I(0) -------------------
cat("ADF confirmation on first-differenced macro series:\n")

for (nm in c("d_gdp", "d_debt")) {
  x   <- df_post[[nm]]
  adf <- tryCatch(
    ur.df(ts(x, frequency = 12), type = "drift", lags = 12, selectlags = "AIC"),
    error = function(e) NULL
  )
  if (!is.null(adf)) {
    stat    <- round(adf@teststat[1], 4)
    cv_5pct <- adf@cval[1, 2]
    verdict <- ifelse(stat < cv_5pct, "STATIONARY (I(0))", "Still non-stationary — check manually")
    cat(sprintf("  %-8s: ADF tau = %7.4f  CV 5%% = %6.4f  => %s\n",
                nm, stat, cv_5pct, verdict))
  }
}

# ---- Save outputs ------------------------------------------
out_dir <- dirname(full_csv)

# Full merged file with differenced columns appended
merged_out <- file.path(out_dir, "nigeria_monthly_merged_post1999.csv")
write.csv(df_post, merged_out, row.names = FALSE)
cat("\nSaved:", merged_out, "\n")

# Stan-ready file: y, d_gdp, d_debt
stan_out <- file.path(out_dir, "nigeria_monthly_stan_input_post1999.csv")
df_stan  <- df_post %>%
  select(date,
         y    = log_exr_return,
         gdp  = d_gdp,
         debt = d_debt)

write.csv(df_stan, stan_out, row.names = FALSE)
cat("Saved:", stan_out, "\n")

# Summary of final Stan input
cat("\nFinal Stan input summary:\n")
cat(sprintf("  Observations : %d\n", nrow(df_stan)))
cat(sprintf("  Date range   : %s to %s\n",
            format(min(df_stan$date)), format(max(df_stan$date))))
cat(sprintf("  y   — mean %.5f  sd %.5f\n", mean(df_stan$y),    sd(df_stan$y)))
cat(sprintf("  gdp — mean %.5f  sd %.5f\n", mean(df_stan$gdp),  sd(df_stan$gdp)))
cat(sprintf("  debt— mean %.5f  sd %.5f\n", mean(df_stan$debt), sd(df_stan$debt)))

cat("\nDone. Re-run stationarity_tests.R to confirm all three series are now I(0).\n")
cat("Then proceed to structural_break_tests.R and Stan estimation.\n")
