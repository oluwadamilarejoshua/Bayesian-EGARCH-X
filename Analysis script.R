
# Load libraries
library(rstan)
library(bayesplot)
library(tidyverse)

# Improve speed of Stan models
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

data <- read.csv("../GARCH for Conference/df_diff.csv")

# Create lagged variables
y <- data[c(2:40), 1]
y_lag <- data[c(1:39), 1]
gdp_lag <- data[c(1:39), 2]
debt_lag <- data[c(1:39), 3]
T <- length(y)

# Data list for Stan
stan_data <- list(
  T = T,
  y = y,
  y_lag = y_lag,
  gdp_lag = gdp_lag,
  debt_lag = debt_lag
)

# fit <- stan(
#   file = "egarch_model.stan",
#   data = stan_data,
#   iter = 2000,
#   chains = 4,
#   seed = 123
# )

fit <- stan(
  file = "egarch_model.stan",
  data = stan_data,
  iter = 10000, 
  warmup = 2000,
  chains = 4,
  seed = 123,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)

# Posterior summaries
print(fit, pars = c("beta0", "betaExR", "betaGDP", "betaDebt", "omega", "alpha", "gamma", "phi"))


# Set bayesplot theme to black-and-white
color_scheme_set("gray")

# Traceplots in black & white
stan_trace(fit, pars = c("beta0", "betaExR", "betaGDP", "betaDebt",
                         "omega", "alpha", "gamma", "phi"))

# Posterior density plots in black & white
stan_dens(fit, pars = c("beta0", "betaExR", "betaGDP", "betaDebt", 
                        "omega", "alpha", "gamma", "phi"))



print(fit)
sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
acceptance_rates <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
mean_acceptance_rate <- mean(acceptance_rates)
mean_acceptance_rate
