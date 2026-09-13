// Prior sensitivity variant A: WIDENED
// Same prior means as the baseline model (betaGDP=-0.50, betaDebt=+0.50)
// but SD doubled (0.40 -> 0.80), i.e. much weaker confidence in the
// magnitude while keeping the theory-motivated sign.

data {
  int<lower=1> T;
  vector[T] y;
  vector[T] y_lag;
  vector[T] gdp;
  vector[T] debt;
  real<lower=0> h1_init;
}

parameters {
  real beta0;
  real betaExR;
  real omega;
  real alpha;
  real gamma;
  real phi_raw;
  real betaGDP;
  real betaDebt;
}

transformed parameters {
  real phi;
  vector[T] epsilon;
  vector[T] log_h;
  vector[T] h;

  phi = inv_logit(phi_raw);

  log_h[1]   = log(h1_init);
  epsilon[1] = y[1] - (beta0 + betaExR * y_lag[1]);

  for (t in 2:T) {
    real h_lag = fmax(exp(log_h[t-1]), 1e-10);
    real z_lag = epsilon[t-1] / sqrt(h_lag);

    real raw_log_h = omega
                   + alpha  * (fabs(z_lag) - sqrt(2.0 / pi()))
                   + gamma  * z_lag
                   + phi    * log_h[t-1]
                   + betaGDP  * gdp[t]
                   + betaDebt * debt[t];
    log_h[t] = fmin(raw_log_h, 5.0);

    epsilon[t] = y[t] - (beta0 + betaExR * y_lag[t]);
  }

  for (t in 1:T)
    h[t] = fmax(exp(log_h[t]), 1e-10);
}

model {
  beta0   ~ normal(0,     0.10);
  omega   ~ normal(0,     2.00);
  alpha   ~ normal(0,     0.50);
  gamma   ~ normal(0,     0.50);
  phi_raw ~ normal(0,     1.00);

  betaExR  ~ normal( 0.20, 0.15);
  betaGDP  ~ normal(-0.50, 0.80);   // <-- SD widened from 0.40
  betaDebt ~ normal( 0.50, 0.80);   // <-- SD widened from 0.40

  for (t in 1:T)
    target += normal_lpdf(epsilon[t] | 0, sqrt(h[t]));
}

generated quantities {
  vector[T] log_lik;
  for (t in 1:T)
    log_lik[t] = normal_lpdf(epsilon[t] | 0, sqrt(h[t]));
}
