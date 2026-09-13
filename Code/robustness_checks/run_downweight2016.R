# ============================================================
# run_downweight2016.R
# Robustness refit: baseline priors, but June & July 2016 excluded
# from the likelihood (weight = 0), to test whether betaDebt survives
# without those two high-leverage observations.
# ============================================================

if (.Platform$OS.type == "windows" && nchar(Sys.which("make")) == 0) {
  Sys.setenv(PATH = paste("C:/rtools45/usr/bin", Sys.getenv("PATH"), sep = ";"))
}

suppressMessages(library(rstan))
rstan_options(auto_write = TRUE)
options(mc.cores = 1)

data_path <- "c:/Users/joshua.ogundairo/Downloads/My Research Works/Impacts_of_GDP_and_Debt_Servicing_on_Exchange_Rate_Volatility_of_the_Nigerian_Naira__A_Bayesian_EGARCH_Approach/Analysis/Data/nigeria_monthly_stan_input_post1999.csv"
stan_file <- "C:/Users/JOSHUA~1.OGU/AppData/Local/Temp/claude/c--Users-joshua-ogundairo-Downloads-My-Research-Works-Impacts-of-GDP-and-Debt-Servicing-on-Exchange-Rate-Volatility-of-the-Nigerian-Naira--A-Bayesian-EGARCH-Approach/89ca73a5-75a5-4283-8ef8-9cb30222e7e7/scratchpad/sensitivity/egarch_x_downweight2016.stan"
out_csv   <- "C:/Users/JOSHUA~1.OGU/AppData/Local/Temp/claude/c--Users-joshua-ogundairo-Downloads-My-Research-Works-Impacts-of-GDP-and-Debt-Servicing-on-Exchange-Rate-Volatility-of-the-Nigerian-Naira--A-Bayesian-EGARCH-Approach/89ca73a5-75a5-4283-8ef8-9cb30222e7e7/scratchpad/sensitivity/downweight2016_summary.csv"

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

dates_for_y <- df$date[2:T_full]   # aligns 1:1 with y, gdp, debt, t index
weight <- rep(1, T_obs)
excl_idx <- which(dates_for_y %in% as.Date(c("2016-06-01", "2016-07-01")))
cat("Excluding indices:", excl_idx, "-> dates:", format(dates_for_y[excl_idx]), "\n")
weight[excl_idx] <- 0

stan_data <- list(T = T_obs, y = y, y_lag = y_lag, gdp = gdp, debt = debt,
                   h1_init = h1_init, weight = weight)

init_list <- lapply(1:4, function(i) list(
  beta0 = 0.006, betaExR = 0.20, omega = -0.50, alpha = 0.15,
  gamma = 0.10, phi_raw = 0.00, betaGDP = -0.30, betaDebt = 0.30
))

cat("Fitting downweight-2016 variant...\n")
t0 <- Sys.time()
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
cat("Fit time:", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")

rhat_vals <- summary(fit)$summary[, "Rhat"]
rhat_bad  <- rhat_vals[!is.na(rhat_vals) & rhat_vals > 1.01]
cat("R-hat > 1.01 count:", length(rhat_bad), "\n")

pars_all <- c("beta0", "betaExR", "omega", "alpha", "gamma", "phi", "betaGDP", "betaDebt")
post_summary <- as.data.frame(summary(fit, pars = pars_all)$summary)
post_summary$parameter <- rownames(post_summary)
post_summary$variant   <- "downweight2016"
post_summary <- post_summary[, c("variant", "parameter", "mean", "sd",
                                  "2.5%", "50%", "97.5%", "n_eff", "Rhat")]
colnames(post_summary) <- c("variant", "parameter", "mean", "sd",
                             "ci_2.5", "ci_50", "ci_97.5", "n_eff", "rhat")

write.csv(post_summary, out_csv, row.names = FALSE)
cat("Saved:", out_csv, "\n\n")
print(post_summary, row.names = FALSE)
cat("\nDone: downweight2016\n")
