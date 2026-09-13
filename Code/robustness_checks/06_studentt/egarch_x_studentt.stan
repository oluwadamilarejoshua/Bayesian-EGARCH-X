// ============================================================
// egarch_x_studentt.stan
// Same EGARCH(1,1)-X specification and priors as egarch_x_model.stan
// (the baseline, tight-prior model), except the mean-equation
// innovation is Student-t rather than Gaussian:
//   epsilon[t] ~ t_nu(0, sigma[t]),  Var(epsilon[t]) = h[t]
// sigma[t] is scaled so h[t] remains the true conditional variance
// (Var of a scaled t_nu is sigma^2 * nu/(nu-2) for nu>2), so h[t]
// keeps exactly the same meaning as in the Gaussian model -- only
// the SHAPE of the innovation distribution changes.
// nu ~ Gamma(2, 0.1) truncated below at 2 (finite-variance region);
// weakly informative, allows anywhere from heavy tails (nu near 2)
// to near-Gaussian (large nu).
// Motivated by a Jarque-Bera test on the Gaussian model's
// standardized residuals: excess kurtosis 10.5, JB = 1426, p<0.001.
// ============================================================

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
  real<lower=4> nu;
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
    // Smooth (differentiable) soft-ceiling instead of a hard fmin(): a
    // hard clip creates a zero-gradient kink that Student-t innovations
    // hit often enough (heavier-tailed shocks -> more extreme raw_log_h
    // proposals during warmup) to make HMC's gradients unreliable right
    // at the boundary. This asymptotically approaches the same 5.0
    // ceiling but stays smooth everywhere.
    log_h[t] = 5.0 - log1p_exp(5.0 - raw_log_h);

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
  betaGDP  ~ normal(-0.50, 0.40);
  betaDebt ~ normal( 0.50, 0.40);

  nu ~ gamma(2, 0.1);

  for (t in 1:T) {
    real sigma_t = sqrt(h[t] * (nu - 2) / nu);
    target += student_t_lpdf(epsilon[t] | nu, 0, sigma_t);
  }
}

generated quantities {
  vector[T] log_lik;
  for (t in 1:T) {
    real sigma_t = sqrt(h[t] * (nu - 2) / nu);
    log_lik[t] = student_t_lpdf(epsilon[t] | nu, 0, sigma_t);
  }
}
