# ============================================================
# smoothed_gs10_pipeline.R
# Rebuilds the full data pipeline (fetch -> interpolate -> post-1999
# restrict -> difference) using a SMOOTHED (centred 3-month moving
# average) version of GS10 as the debt-service Denton-Cholette
# indicator, instead of raw GS10. Goal: dilute the sharp June/July
# 2016 dip (a US/global-rate coincidence, unrelated to Nigeria) while
# keeping the genuine "global rates matter" signal GS10 provides.
# GDP's indicator (raw WTI oil price) is unchanged -- oil was never
# implicated in the fragility diagnosis.
# Output goes to a SEPARATE folder, not overwriting the current
# pipeline's files, since we haven't yet decided which specification
# to adopt.
# ============================================================

if (.Platform$OS.type == "windows" && nchar(Sys.which("make")) == 0) {
  Sys.setenv(PATH = paste("C:/rtools45/usr/bin", Sys.getenv("PATH"), sep = ";"))
}

suppressMessages({
  library(dplyr); library(zoo); library(tempdisagg); library(quantmod); library(rdbnomics); library(urca)
})
options(timeout = 120)

OUT_DIR <- "c:/Users/joshua.ogundairo/Downloads/My Research Works/Impacts_of_GDP_and_Debt_Servicing_on_Exchange_Rate_Volatility_of_the_Nigerian_Naira__A_Bayesian_EGARCH_Approach/Analysis/Data/smoothed_gs10_test"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

ANNUAL_CSV <- "c:/Users/joshua.ogundairo/Downloads/My Research Works/Impacts_of_GDP_and_Debt_Servicing_on_Exchange_Rate_Volatility_of_the_Nigerian_Naira__A_Bayesian_EGARCH_Approach/Analysis/Data/nigeria_annual_WDI.csv"

# ---- Step 1: exchange rate (same as main pipeline) ----------
cat("Fetching exchange rate...\n")
exr_monthly <- NULL
tryCatch({
  raw_rdb <- rdb("IMF", "IFS", mask = "M.NG.ENDA_XDC_USD_RATE")
  exr_monthly <- raw_rdb %>% filter(!is.na(value)) %>%
    mutate(date = as.Date(paste0(period, "-01"))) %>%
    filter(date >= as.Date("1980-01-01"), date <= as.Date("2020-12-31")) %>%
    select(date, ngn_per_usd = value) %>% arrange(date)
}, error = function(e) cat("Method A failed:", conditionMessage(e), "\n"))
cat("Exchange rate obs:", nrow(exr_monthly), "\n")

# ---- Step 2: GDP indicator (raw oil price, unchanged) --------
cat("Fetching oil price indicator...\n")
getSymbols("WTISPLC", src = "FRED", auto.assign = TRUE)
oil_raw <- data.frame(date = as.Date(index(WTISPLC)), oil = as.numeric(WTISPLC)) %>%
  filter(date >= as.Date("1980-01-01"), date <= as.Date("2020-12-01"), !is.na(oil)) %>%
  arrange(date)

# ---- Step 3: debt indicator (GS10, SMOOTHED) ------------------
cat("Fetching GS10 and applying centred 3-month moving average...\n")
getSymbols("GS10", src = "FRED", auto.assign = TRUE)
gs10_raw <- data.frame(date = as.Date(index(GS10)), gs10 = as.numeric(GS10)) %>%
  filter(date >= as.Date("1979-11-01"), date <= as.Date("2021-02-01"), !is.na(gs10)) %>%
  arrange(date)   # fetch a bit wider so the centred MA has neighbours at the 1980/2020 edges
gs10_raw$gs10_smooth <- rollapply(gs10_raw$gs10, width = 3, FUN = mean, fill = NA, align = "center")
# fall back to the raw value at the very edges where the 3-window has no neighbour
gs10_raw$gs10_smooth[is.na(gs10_raw$gs10_smooth)] <- gs10_raw$gs10[is.na(gs10_raw$gs10_smooth)]
gs10_raw <- gs10_raw %>% filter(date >= as.Date("1980-01-01"), date <= as.Date("2020-12-01"))

cat("Smoothing check around 2016 float:\n")
print(gs10_raw[gs10_raw$date >= as.Date("2016-04-01") & gs10_raw$date <= as.Date("2016-09-01"), ])

# ---- Step 4: annual GDP/debt ----------------------------------
annual_raw <- read.csv(ANNUAL_CSV, stringsAsFactors = FALSE)
gdp_raw  <- annual_raw %>% select(year, gdp_growth = gdp_growth_pct) %>% filter(year >= 1979, year <= 2020) %>% arrange(year)
debt_raw <- annual_raw %>% select(year, debt_service_pct_exports) %>% filter(year >= 1979, year <= 2020) %>% arrange(year)

# ---- Step 5: Denton-Cholette interpolation --------------------
months_all <- seq.Date(as.Date("1980-01-01"), as.Date("2020-12-01"), by = "month")
oil_aligned  <- data.frame(date = months_all) %>% left_join(oil_raw, by = "date")
gs10_aligned <- data.frame(date = months_all) %>% left_join(gs10_raw[, c("date", "gs10_smooth")], by = "date")
stopifnot(!any(is.na(oil_aligned$oil)), !any(is.na(gs10_aligned$gs10_smooth)))

oil_monthly_ts  <- ts(oil_aligned$oil,           start = c(1980, 1), frequency = 12)
gs10_monthly_ts <- ts(gs10_aligned$gs10_smooth,  start = c(1980, 1), frequency = 12)

gdp_annual_ts  <- ts(gdp_raw  %>% filter(year >= 1980, year <= 2020) %>% pull(gdp_growth), start = 1980, frequency = 1)
debt_annual_ts <- ts(debt_raw %>% filter(year >= 1980, year <= 2020) %>% pull(debt_service_pct_exports), start = 1980, frequency = 1)

gdp_td  <- td(gdp_annual_ts  ~ 0 + oil_monthly_ts,  to = 12, method = "denton-cholette")
debt_td <- td(debt_annual_ts ~ 0 + gs10_monthly_ts, to = 12, method = "denton-cholette")

macro_monthly <- data.frame(
  date = months_all,
  gdp_growth_monthly_interp   = as.numeric(predict(gdp_td)),
  debt_service_monthly_interp = as.numeric(predict(debt_td))
)

# ---- Step 6: merge, log returns, post-1999 restriction, differencing ----
df_monthly <- exr_monthly %>% left_join(macro_monthly, by = "date") %>% arrange(date) %>% filter(complete.cases(.)) %>%
  mutate(log_exr = log(ngn_per_usd), log_exr_return = c(NA, diff(log_exr))) %>% filter(!is.na(log_exr_return))

write.csv(df_monthly, file.path(OUT_DIR, "nigeria_monthly_merged_smoothed.csv"), row.names = FALSE)

df_post <- df_monthly %>% filter(date >= as.Date("1999-02-01")) %>%
  mutate(d_gdp = c(NA, diff(gdp_growth_monthly_interp)), d_debt = c(NA, diff(debt_service_monthly_interp))) %>%
  filter(!is.na(d_gdp))

df_stan <- df_post %>% select(date, y = log_exr_return, gdp = d_gdp, debt = d_debt)
write.csv(df_stan, file.path(OUT_DIR, "nigeria_monthly_stan_input_post1999_smoothed.csv"), row.names = FALSE)
cat("\nSaved:", file.path(OUT_DIR, "nigeria_monthly_stan_input_post1999_smoothed.csv"), "\n")
cat("Final sample:", nrow(df_stan), "obs,", format(min(df_stan$date)), "to", format(max(df_stan$date)), "\n")

# ---- Step 7: quick ADF confirmation on smoothed d_debt --------
cat("\nADF confirmation on smoothed-indicator series:\n")
for (nm in c("gdp", "debt")) {
  x <- df_stan[[nm]]
  adf <- ur.df(ts(x, frequency = 12), type = "drift", lags = 12, selectlags = "AIC")
  cat(sprintf("  d_%-4s: ADF tau = %7.4f  CV 5%% = %6.4f  => %s\n", nm,
              adf@teststat[1], adf@cval[1, 2],
              ifelse(adf@teststat[1] < adf@cval[1, 2], "STATIONARY (I(0))", "check manually")))
}

# ---- Step 8: side-by-side comparison with the RAW (unsmoothed) series around 2016 ----
raw_stan <- read.csv("c:/Users/joshua.ogundairo/Downloads/My Research Works/Impacts_of_GDP_and_Debt_Servicing_on_Exchange_Rate_Volatility_of_the_Nigerian_Naira__A_Bayesian_EGARCH_Approach/Analysis/Data/nigeria_monthly_stan_input_post1999.csv")
cmp <- merge(df_stan[df_stan$date >= "2016-04-01" & df_stan$date <= "2016-09-01", c("date","debt")],
             raw_stan[raw_stan$date >= "2016-04-01" & raw_stan$date <= "2016-09-01", c("date","debt")],
             by = "date", suffixes = c("_smoothed", "_raw"))
cat("\n2016 window comparison (smoothed vs raw d_debt):\n")
print(cmp)

cat("\nDone: smoothed_gs10_pipeline\n")
