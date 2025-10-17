

data {
  int<lower=1> T;
  vector[T] y;         // Exchange rate
  vector[T] y_lag;     // Lag of exchange rate
  vector[T] gdp_lag;   // Lag of GDP
  vector[T] debt_lag;  // Lag of debt servicing
}

parameters {
  real beta0;
  real betaExR;
  real betaGDP;
  real betaDebt;

  real omega;
  real alpha;
  real gamma;
  real phi;
}

transformed parameters {
  vector[T] epsilon;
  vector[T] log_h;
  vector[T] h;

  log_h[1] = log(0.01); // Initialize log volatility
  epsilon[1] = y[1] - (beta0 + betaExR * y_lag[1] + betaGDP * gdp_lag[1] + betaDebt * debt_lag[1]);

  for (t in 2:T) {
    log_h[t] = omega + alpha * fabs(epsilon[t-1] / sqrt(exp(log_h[t-1]))) +
                      gamma * (epsilon[t-1] / sqrt(exp(log_h[t-1]))) +
                      phi * log_h[t-1];

    epsilon[t] = y[t] - (beta0 + betaExR * y_lag[t] + betaGDP * gdp_lag[t] + betaDebt * debt_lag[t]);
  }

  for (t in 1:T) {
    h[t] = exp(log_h[t]);
  }
}

model {
  // Priors from frequentist estimates
  beta0 ~ normal(0, 1);
  betaExR ~ normal(0.4, 0.2);      // Lagged exchange rate
  betaGDP ~ normal(-0.25, 0.1);    // GDP
  betaDebt ~ normal(-0.06, 0.01);   // Debt servicing

  omega ~ normal(0, 1);
  alpha ~ normal(0, 1);
  gamma ~ normal(0, 1);
  phi   ~ normal(0.5, 0.2);

  for (t in 1:T) {
    target += normal_lpdf(epsilon[t] | 0, sqrt(h[t]));
  }
}

