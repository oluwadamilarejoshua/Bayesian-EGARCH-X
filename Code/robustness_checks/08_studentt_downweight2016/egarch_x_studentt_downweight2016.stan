// ============================================================
// egarch_x_studentt_downweight2016.stan
// Robustness variant of the ADOPTED baseline (egarch_x_studentt.stan:
// AR(1) mean eq, Student-t innovations, smooth soft-ceiling, smoothed-
// GS10 data) that EXCLUDES June & July 2016 from the likelihood.
// The recursive variance path (epsilon, log_h, h) is computed exactly
// as in the baseline, using the true data at every t -- so the
// time-series structure (each month's h depends on the prior month's
// actual shock) is preserved intact. Only the LIKELIHOOD contribution
// of the two flagged months is zeroed via a per-observation `weight`,
// so those two extreme observations no longer influence the posterior
// for theta, while still propagating forward into subsequent h[t].
// Purpose: test whether betaDebt (and the other structural parameters)
// survive without the two months containing the Naira float shock --
// i.e. that the debt-servicing result is not an artifact of GS10's
// value happening to coincide with the sample's largest FX shock.
// ============================================================

data {
  int<lower=1> T;
  vector[T] y;
  vector[T] y_lag;
  vector[T] gdp;
  vector[T] debt;
  real<lower=0> h1_init;
  vector[T] weight;        // 1 = included, 0 = excluded from target
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
    // Same smooth soft-ceiling as the adopted baseline.
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
    target += weight[t] * student_t_lpdf(epsilon[t] | nu, 0, sigma_t);
  }
}

generated quantities {
  vector[T] log_lik;
  for (t in 1:T) {
    real sigma_t = sqrt(h[t] * (nu - 2) / nu);
    log_lik[t] = student_t_lpdf(epsilon[t] | nu, 0, sigma_t);
  }
}
