# ============================================================
# Script: Fetch complete monthly data for improved analysis
# Purpose: Retrieve monthly official exchange rate (IMF IFS source)
#          and annual macro data from WDI, then prepare the
#          combined dataset for the Bayesian EGARCH-X model.
#
# WHAT CHANGED FROM v1:
#   The original script used FRED series "EXNAUS" which does not
#   exist — that is why you got the 404 error. The correct
#   approach is Method A below (rdbnomics / IMF IFS), with
#   Method B (FRED with corrected series IDs) as a fallback.
#
# Data sources:
#   Exchange rate (monthly, official):
#     METHOD A [RECOMMENDED]:
#       rdbnomics package → IMF IFS series M.NG.ENDA_XDC_USD_RATE
#       Definition: Nigeria Domestic Currency Units per USD,
#                   Period Average (official rate), Monthly
#       Coverage:   January 1957 – present
#       Citation:   International Monetary Fund (2024). International
#                   Financial Statistics. IMF Data, Washington D.C.
#                   Retrieved via DBnomics: https://db.nomics.world/IMF/IFS/M.NG.ENDA_XDC_USD_RATE
#
#     METHOD B [FRED fallback]:
#       FRED series FXRATENGA618NUPN  (IMF IFS source, 1950-2010)
#       FRED series XRNCUSNGA618NRUG  (market+estimated, 1950-2023)
#       quantmod::getSymbols(src = "FRED")
#       Citation:   Federal Reserve Bank of St. Louis (2024). FRED.
#                   https://fred.stlouisfed.org
#
#   GDP growth (annual %):
#       World Bank WDI indicator NY.GDP.MKTP.KD.ZG
#
#   Debt service (% of exports):
#       World Bank WDI indicator DT.TDS.DECT.EX.ZS
#
# Interpolation:
#   GDP and debt servicing are annual. They are disaggregated to
#   monthly using the Denton-Cholette method (tempdisagg package).
#   Reference: Sax & Steiner (2013), The R Journal, 5(2), 80-87.
# ============================================================

# --- Install required packages ------------------------------
required_packages <- c("rdbnomics", "quantmod", "WDI", "tempdisagg",
                       "dplyr", "lubridate", "zoo")
new_packages <- required_packages[!(required_packages %in%
                                      installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

library(dplyr)
library(lubridate)
library(zoo)
library(tempdisagg)
library(WDI)

# ============================================================
# STEP 1: Monthly exchange rate
# Try Method A (rdbnomics/IMF IFS) first; fall back to Method B (FRED)
# ============================================================

exr_monthly <- NULL

# --- Method A: rdbnomics (IMF IFS) --------------------------
tryCatch({
  library(rdbnomics)
  cat("Trying Method A: rdbnomics / IMF IFS...\n")

  raw_rdb <- rdb("IMF", "IFS", mask = "M.NG.ENDA_XDC_USD_RATE")

  exr_monthly <- raw_rdb %>%
    filter(!is.na(value)) %>%
    mutate(date = as.Date(paste0(period, "-01"))) %>%
    filter(date >= as.Date("1980-01-01"),
           date <= as.Date("2020-12-31")) %>%
    select(date, ngn_per_usd = value) %>%
    arrange(date)

  cat("Method A succeeded:", nrow(exr_monthly), "monthly observations.\n")

}, error = function(e) {
  cat("Method A failed:", conditionMessage(e), "\n")
  cat("Falling back to Method B (FRED)...\n")
})

# --- Method B: FRED (quantmod) fallback ---------------------
if (is.null(exr_monthly) || nrow(exr_monthly) == 0) {
  tryCatch({
    library(quantmod)

    # Try the primary FRED series (IMF IFS official rate, 1950-2010)
    # If this series is annual, switch to XRNCUSNGA618NRUG
    tryCatch({
      cat("Trying FRED series FXRATENGA618NUPN...\n")
      getSymbols("FXRATENGA618NUPN", src = "FRED", auto.assign = TRUE)
      raw_fred <- get("FXRATENGA618NUPN")
      freq_check <- nrow(raw_fred)
      cat("  Observations fetched:", freq_check, "\n")

      exr_monthly <<- data.frame(
        date        = as.Date(index(raw_fred)),
        ngn_per_usd = as.numeric(raw_fred)
      ) %>%
        filter(date >= as.Date("1980-01-01"),
               date <= as.Date("2020-12-31"),
               !is.na(ngn_per_usd)) %>%
        arrange(date)

      cat("FXRATENGA618NUPN loaded:", nrow(exr_monthly), "observations.\n")

    }, error = function(e) {
      cat("FXRATENGA618NUPN failed:", conditionMessage(e), "\n")

      # Second FRED fallback: market+estimated rate 1950-2023
      cat("Trying FRED series XRNCUSNGA618NRUG...\n")
      getSymbols("XRNCUSNGA618NRUG", src = "FRED", auto.assign = TRUE)
      raw_fred2 <- get("XRNCUSNGA618NRUG")

      exr_monthly <<- data.frame(
        date        = as.Date(index(raw_fred2)),
        ngn_per_usd = as.numeric(raw_fred2)
      ) %>%
        filter(date >= as.Date("1980-01-01"),
               date <= as.Date("2020-12-31"),
               !is.na(ngn_per_usd)) %>%
        arrange(date)

      cat("XRNCUSNGA618NRUG loaded:", nrow(exr_monthly), "observations.\n")
    })

  }, error = function(e) {
    cat("Method B also failed:", conditionMessage(e), "\n")
    cat("STOPPING: Could not obtain monthly exchange rate data.\n")
    cat("Please ensure you have internet access and retry.\n")
    stop("Exchange rate data unavailable. Cannot proceed.")
  })
}

# Confirm frequency — warn if likely annual rather than monthly
n_obs  <- nrow(exr_monthly)
n_years <- as.numeric(difftime(max(exr_monthly$date),
                                min(exr_monthly$date), units = "days")) / 365.25
if (n_obs / n_years < 6) {
  warning(paste(
    "Only", round(n_obs / n_years, 1), "observations per year detected.",
    "This series may be annual, not monthly.",
    "Monthly data requires ~12 obs/year.",
    "Check the series frequency before proceeding."
  ))
} else {
  cat("Frequency check passed:", round(n_obs / n_years, 1), "obs/year.\n")
}

cat("Exchange rate data range:", format(min(exr_monthly$date)),
    "to", format(max(exr_monthly$date)), "\n")

# ============================================================
# STEP 1b: Monthly indicator for GDP disaggregation — WTI oil price
# Nigeria's GDP and export earnings are heavily oil-driven, so the
# monthly WTI spot price is used as the within-year movement proxy
# for the Denton-Cholette disaggregation of annual GDP growth.
# Source: FRED series WTISPLC (Spot Crude Oil Price: WTI, dollars
# per barrel), monthly, Jan 1946-present, no gaps.
# Citation: Federal Reserve Bank of St. Louis (2024). FRED.
#           https://fred.stlouisfed.org/series/WTISPLC
# ============================================================
cat("\nFetching GDP indicator: WTI oil price (FRED WTISPLC)...\n")
library(quantmod)
getSymbols("WTISPLC", src = "FRED", auto.assign = TRUE)
oil_raw <- data.frame(
  date = as.Date(index(WTISPLC)),
  oil  = as.numeric(WTISPLC)
) %>%
  filter(date >= as.Date("1980-01-01"), date <= as.Date("2020-12-01"),
         !is.na(oil)) %>%
  arrange(date)
cat("Oil price indicator:", nrow(oil_raw), "monthly observations,",
    format(min(oil_raw$date)), "to", format(max(oil_raw$date)), "\n")

# ============================================================
# STEP 1c: Monthly indicator for debt-service disaggregation —
# US 10-Year Treasury Constant Maturity Rate
#
# REVISION NOTE: an earlier version of this script used Nigeria's
# own FX reserves (IMF IFS M.NG.RAXG_USD) as the debt indicator.
# That was dropped: reserves are the exact instrument the CBN uses
# to defend the naira band, so a reserves-shaped debt regressor
# shares a common driver with the dependent variable (exchange-rate
# returns) — a real circularity risk that could spuriously inflate
# betaDebt. The US 10-year Treasury yield has no such problem: it
# is set entirely by US monetary policy and global bond markets,
# with zero mechanical link to CBN operations or the naira market.
# It is still economically motivated, not an arbitrary substitute:
# most of Nigeria's external debt is dollar-denominated, and the
# prevailing global rate directly drives the cost of servicing and
# refinancing that debt — a higher US rate raises Nigeria's debt-
# service burden independent of anything happening in Nigeria.
# Source: FRED series GS10 (Market Yield on U.S. Treasury Securities
# at 10-Year Constant Maturity), monthly, percent, Apr 1953-present,
# no gaps in the 1980-2020 window.
# Citation: Federal Reserve Bank of St. Louis (2024). FRED.
#           https://fred.stlouisfed.org/series/GS10
# ============================================================
cat("\nFetching debt-service indicator: US 10Y Treasury yield (FRED GS10)...\n")
library(quantmod)
getSymbols("GS10", src = "FRED", auto.assign = TRUE)
gs10_raw <- data.frame(
  date  = as.Date(index(GS10)),
  gs10  = as.numeric(GS10)
) %>%
  filter(date >= as.Date("1980-01-01"), date <= as.Date("2020-12-01"),
         !is.na(gs10)) %>%
  arrange(date)
cat("US 10Y Treasury yield indicator:", nrow(gs10_raw), "monthly observations,",
    format(min(gs10_raw$date)), "to", format(max(gs10_raw$date)), "\n")

# ============================================================
# STEP 2 & 3: Annual GDP growth and debt service
# Primary source: local CSV (already downloaded — no API call needed)
# Fallback: WDI API (used only if the local file is missing)
# ============================================================

# Locate the pre-downloaded annual file — try multiple strategies in order:
#   1. RStudio active editor path (interactive session)
#   2. Path of the currently sourced file (source() call)
#   3. Rscript --file argument (command-line execution)
#   4. Current working directory (last resort)
annual_csv <- tryCatch({
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
  file.path(script_dir, "../Data/raw/nigeria_annual_WDI.csv")
}, error = function(e) {
  # Try sys.frame path (works when script is called via source())
  src_path <- tryCatch(
    normalizePath(sys.frame(1)$ofile),
    error = function(e2) NULL
  )
  if (!is.null(src_path)) {
    file.path(dirname(src_path), "../Data/raw/nigeria_annual_WDI.csv")
  } else {
    # Try Rscript --file argument
    args <- commandArgs(trailingOnly = FALSE)
    file_flag <- grep("--file=", args, value = TRUE)
    if (length(file_flag) > 0) {
      script_path <- normalizePath(sub("--file=", "", file_flag[1]))
      file.path(dirname(script_path), "../Data/raw/nigeria_annual_WDI.csv")
    } else {
      # Fall back to working directory
      file.path(getwd(), "../Data/raw/nigeria_annual_WDI.csv")
    }
  }
})

if (file.exists(annual_csv)) {
  cat("\nReading annual macro data from local file:", annual_csv, "\n")
  annual_raw <- read.csv(annual_csv, stringsAsFactors = FALSE)

  gdp_raw <- annual_raw %>%
    select(year, gdp_growth = gdp_growth_pct) %>%
    filter(year >= 1979, year <= 2020) %>%
    arrange(year)

  debt_raw <- annual_raw %>%
    select(year, debt_service_pct_exports) %>%
    filter(year >= 1979, year <= 2020) %>%
    arrange(year)

  cat("GDP observations loaded from file:", nrow(gdp_raw), "\n")
  cat("Debt service observations loaded from file:", nrow(debt_raw), "\n")

} else {
  # --- Fallback: download from WDI API ----------------------
  cat("\nLocal file not found. Downloading from World Bank WDI API...\n")
  cat("(This may time out — if it does, ensure nigeria_annual_WDI.csv is in the Data folder)\n")

  library(WDI)

  gdp_raw <- tryCatch(
    WDI(indicator = "NY.GDP.MKTP.KD.ZG", country = "NG",
        start = 1979, end = 2020, extra = FALSE) %>%
      select(year, gdp_growth = NY.GDP.MKTP.KD.ZG) %>%
      arrange(year),
    error = function(e) {
      stop("WDI download failed for GDP. Place nigeria_annual_WDI.csv in the Data folder and re-run.\n",
           "Original error: ", conditionMessage(e))
    }
  )

  debt_raw <- tryCatch(
    WDI(indicator = "DT.TDS.DECT.EX.ZS", country = "NG",
        start = 1979, end = 2020, extra = FALSE) %>%
      select(year, debt_service_pct_exports = DT.TDS.DECT.EX.ZS) %>%
      arrange(year),
    error = function(e) {
      stop("WDI download failed for debt service. Place nigeria_annual_WDI.csv in the Data folder and re-run.\n",
           "Original error: ", conditionMessage(e))
    }
  )

  cat("Annual GDP observations:", nrow(gdp_raw), "\n")
  cat("Annual debt service observations:", nrow(debt_raw), "\n")
}

# ============================================================
# STEP 4: Interpolate annual series to monthly (Denton-Cholette),
# WITH monthly indicators (oil price for GDP; US 10Y Treasury yield
# for debt service).
# Previously this used formula `~ 1` (no indicator), which forces
# the disaggregation to fall back on pure curvature-minimisation —
# the smoothest mathematically possible path consistent with the
# annual totals, with no real monthly economic signal shaping the
# within-year distribution. Supplying a genuine related indicator
# lets Denton-Cholette preserve the indicator's actual month-to-
# month movement instead, which is the method's intended use case.
# NOTE: Denton methods require `0 +` in the formula to suppress the
# intercept — without it, the model matrix gets 2 columns (intercept
# + indicator) and td() errors with "only one series allowed".
# ============================================================
cat("\nInterpolating annual macro series to monthly (indicator-based)...\n")

months_all <- seq.Date(as.Date("1980-01-01"), as.Date("2020-12-01"), by = "month")

# Align indicators to the full monthly grid and fail loudly on any gap
oil_aligned <- data.frame(date = months_all) %>%
  left_join(oil_raw, by = "date")
gs10_aligned <- data.frame(date = months_all) %>%
  left_join(gs10_raw, by = "date")

if (any(is.na(oil_aligned$oil)))
  stop("Oil indicator has missing months in the 1980-2020 window — check the fetch.")
if (any(is.na(gs10_aligned$gs10)))
  stop("GS10 indicator has missing months in the 1980-2020 window — check the fetch.")

oil_monthly_ts  <- ts(oil_aligned$oil,   start = c(1980, 1), frequency = 12)
gs10_monthly_ts <- ts(gs10_aligned$gs10, start = c(1980, 1), frequency = 12)

# GDP — indicator: WTI oil price
gdp_annual_ts <- ts(
  gdp_raw %>% filter(year >= 1980, year <= 2020) %>% pull(gdp_growth),
  start = 1980, frequency = 1
)
gdp_td <- td(gdp_annual_ts ~ 0 + oil_monthly_ts, to = 12, method = "denton-cholette")
gdp_monthly_ts <- predict(gdp_td)
cat("\nGDP disaggregation (indicator = WTI oil price):\n")
print(summary(gdp_td))

# Debt service — indicator: US 10Y Treasury yield (exogenous to Nigeria's FX market)
debt_annual_ts <- ts(
  debt_raw %>% filter(year >= 1980, year <= 2020) %>% pull(debt_service_pct_exports),
  start = 1980, frequency = 1
)
debt_td <- td(debt_annual_ts ~ 0 + gs10_monthly_ts, to = 12, method = "denton-cholette")
debt_monthly_ts <- predict(debt_td)
cat("\nDebt-service disaggregation (indicator = US 10Y Treasury yield):\n")
print(summary(debt_td))

macro_monthly <- data.frame(
  date                        = months_all,
  gdp_growth_monthly_interp   = as.numeric(gdp_monthly_ts),
  debt_service_monthly_interp = as.numeric(debt_monthly_ts)
)

cat("\nInterpolated monthly macro observations:", nrow(macro_monthly), "\n")

# ============================================================
# STEP 5: Merge all series
# ============================================================
cat("\nMerging series...\n")
df_monthly <- exr_monthly %>%
  left_join(macro_monthly, by = "date") %>%
  arrange(date) %>%
  filter(complete.cases(.))

cat("Merged dataset:", nrow(df_monthly), "observations,",
    format(min(df_monthly$date)), "to", format(max(df_monthly$date)), "\n")

# ============================================================
# STEP 6: Compute log returns and differences for modelling
# ============================================================
df_monthly <- df_monthly %>%
  arrange(date) %>%
  mutate(
    log_exr          = log(ngn_per_usd),
    log_exr_return   = c(NA, diff(log_exr)),
    gdp_growth_m     = gdp_growth_monthly_interp,
    debt_service_m   = debt_service_monthly_interp
  ) %>%
  filter(!is.na(log_exr_return))

cat("Final modelling dataset:", nrow(df_monthly), "observations\n")

# ============================================================
# STEP 7: Save outputs
# ============================================================
# All outputs are written to ../Data/processed/, relative to this
# script's own location (Code/), regardless of the caller's working
# directory. out_dir uses the same lookup strategy as annual_csv above.
out_dir <- tryCatch({
  file.path(dirname(rstudioapi::getSourceEditorContext()$path), "../Data/processed")
}, error = function(e) {
  src_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e2) NULL)
  if (!is.null(src_path)) {
    file.path(dirname(src_path), "../Data/processed")
  } else {
    args <- commandArgs(trailingOnly = FALSE)
    flag <- grep("--file=", args, value = TRUE)
    if (length(flag) > 0)
      file.path(dirname(normalizePath(sub("--file=", "", flag[1]))), "../Data/processed")
    else
      file.path(getwd(), "../Data/processed")
  }
})
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

write.csv(df_monthly, file.path(out_dir, "nigeria_monthly_merged.csv"), row.names = FALSE)
cat("Saved: nigeria_monthly_merged.csv\n")

df_stan_input <- df_monthly %>%
  select(date,
         y    = log_exr_return,
         gdp  = gdp_growth_m,
         debt = debt_service_m)

write.csv(df_stan_input, file.path(out_dir, "nigeria_monthly_stan_input.csv"), row.names = FALSE)
cat("Saved: nigeria_monthly_stan_input.csv\n")

cat("\nDone. Check the frequency warning above if it appeared.\n")
cat("If the series turned out to be annual (< 6 obs/year), switch to\n")
cat("rdbnomics Method A or contact the maintainer.\n")

# ============================================================
# STEP 8: Restricted sample — post-1999 managed-float era
# ============================================================
# The full sample (1980-2020) cannot be fed directly into EGARCH
# because the Abacha peg (1993-06 to 1998-12) produced 67
# consecutive months with y=0, and the 1986 SAP devaluation
# (y=+1.050) and end-of-peg jump (1999-01, y=+1.364) are
# discrete policy shocks 20-30 SD from the floating-era mean.
# Starting from February 1999 gives the clean managed-float era.
cat("\n--- STEP 8: Restricted sample (post-1999 managed-float era) ---\n")

df_post1999 <- df_monthly %>%
  filter(date >= as.Date("1999-02-01"))

n_post    <- nrow(df_post1999)
y_post    <- df_post1999$log_exr_return
n_zeros   <- sum(y_post == 0)
zero_pct  <- round(100 * n_zeros / n_post, 1)

cat("Sample:", n_post, "observations,",
    format(min(df_post1999$date)), "to", format(max(df_post1999$date)), "\n")
cat("\nDescriptive stats — y (log return, post-1999):\n")
cat("  Mean      :", round(mean(y_post), 6), "\n")
cat("  Std dev   :", round(sd(y_post),   6), "\n")
cat("  Min       :", round(min(y_post),  6), " on",
    format(df_post1999$date[which.min(y_post)]), "\n")
cat("  Max       :", round(max(y_post),  6), " on",
    format(df_post1999$date[which.max(y_post)]), "\n")
cat("  Exact zeros:", n_zeros, "of", n_post,
    paste0("(", zero_pct, "%)"), "\n")

if (zero_pct > 5)
  warning(paste0("Post-1999 sample still has ", zero_pct,
                 "% exact zeros — check for remaining peg sub-periods."))

df_stan_post1999 <- df_post1999 %>%
  select(date,
         y    = log_exr_return,
         gdp  = gdp_growth_m,
         debt = debt_service_m)

write.csv(df_post1999,       file.path(out_dir, "nigeria_monthly_merged_post1999.csv"),    row.names = FALSE)
write.csv(df_stan_post1999,  file.path(out_dir, "nigeria_monthly_stan_input_post1999.csv"), row.names = FALSE)

cat("\nSaved: nigeria_monthly_merged_post1999.csv\n")
cat("Saved: nigeria_monthly_stan_input_post1999.csv\n")
cat("\nNext step: run stationarity_tests.R on the post-1999 dataset.\n")
