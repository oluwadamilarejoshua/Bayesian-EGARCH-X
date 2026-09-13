// ============================================================
// egarch_x_ar2_studentt.stan
// Combined revision of egarch_x_model.stan, prompted by residual
// diagnostics on the AR(1)+Gaussian fit:
//   (1) Mean equation extended from AR(1) to AR(2). A Ljung-Box
//       test on that model's standardized residuals rejected at
//       every lag checked (6, 12, 24); the PACF of y_t itself has
//       a clear, significant spike at lag 2 (-0.196, exceeding the
//       ~0.121 significance bound for n=262), and a BIC scan across
//       AR(0)-AR(5) is minimised at AR(2). Both point the same way.
//   (2) Mean-equation innovation changed from Gaussian to Student-t.
//       A Jarque-Bera test on the same residuals rejected normality
//       overwhelmingly (excess kurtosis 10.5, JB=1426, p<0.001).
// Variance equation, and all its priors, are unchanged from
// egarch_x_model.stan. beta2's prior is weakly informative
// (N(0,0.15)) and does not presuppose a sign, unlike beta1's
// theory/OLS-motivated prior -- this lag was added in response to
// diagnostics, not theory, so we let the data speak for it.
// ============================================================

data {
  int<lower=1> T;
  vector[T] y;
  vector[T] y_lag;    // y[t-1]
  vector[T] y_lag2;   // y[t-2]
  vector[T] gdp;
  vector[T] debt;
  real<lower=0> h1_init;
}

parameters {
  real beta0;
  real betaExR;   // coefficient on y[t-1]
  real beta2;     // coefficient on y[t-2] (new)
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
  epsilon[1] = y[1] - (beta0 + betaExR * y_lag[1] + beta2 * y_lag2[1]);

  for (t in 2:T) {
    real h_lag = fmax(exp(log_h[t-1]), 1e-10);
    real z_lag = epsilon[t-1] / sqrt(h_lag);

    real raw_log_h = omega
                   + alpha  * (fabs(z_lag) - sqrt(2.0 / pi()))
                   + gamma  * z_lag
                   + phi    * log_h[t-1]
                   + betaGDP  * gdp[t]
                   + betaDebt * debt[t];
    // Smooth (differentiable) soft-ceiling instead of a hard fmin() --
    // see egarch_x_studentt.stan for why this matters under Student-t
    // innovations specifically.
    log_h[t] = 5.0 - log1p_exp(5.0 - raw_log_h);

    epsilon[t] = y[t] - (beta0 + betaExR * y_lag[t] + beta2 * y_lag2[t]);
  }

  for (t in 1:T)
    h[t] = fmax(exp(log_h[t]), 1e-10);
}

model {
  beta0   ~ normal(0,     0.10);
  beta2   ~ normal(0,     0.15);   // weakly informative, no presumed sign
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
