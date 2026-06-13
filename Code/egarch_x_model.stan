// ============================================================
// egarch_x_model.stan  (v4 — numerical ceiling on log_h, phi prior widened)
// Bayesian EGARCH(1,1)-X: Nelson (1991) specification.
//
// Mean equation (AR(1) in log returns):
//   y[t] = beta0 + betaExR * y_lag[t] + epsilon[t],  epsilon[t] ~ N(0, h[t])
//
// Variance equation (EGARCH-X):
//   log h[t] = omega
//            + alpha  * (|z[t-1]| - sqrt(2/pi))   // ARCH (Nelson-centred)
//            + gamma  * z[t-1]                      // asymmetry / leverage
//            + phi    * log h[t-1]                  // GARCH persistence
//            + betaGDP  * gdp[t]                    // GDP growth channel
//            + betaDebt * debt[t]                   // debt service channel
//   z[t] = epsilon[t] / sqrt(h[t])
//
// gdp and debt are standardised (mean 0, SD 1) in the R script before
// being passed here, so betaGDP and betaDebt are in units of
// "change in log conditional variance per 1-SD change in the regressor."
//
// phi is reparametrised via the logistic transform:
//   phi_raw ~ unconstrained real; phi = inv_logit(phi_raw) in (0,1).
// This removes the hard upper boundary at 1 that causes HMC gradient
// explosion when posterior mass concentrates near high persistence.
//
// Prior specification:
//   WEAKLY INFORMATIVE (system parameters — let the data decide):
//   beta0   : normal(0, 0.10) — monthly log return near zero under managed float
//   omega   : normal(0, 2.00) — log-variance intercept; covers plausible range
//   alpha   : normal(0, 0.50) — symmetric shock response near zero a priori
//   gamma   : normal(0, 0.50) — leverage can be positive or negative
//   phi_raw : normal(0, 1.00) — inv_logit(0)=0.50; 95% prior for phi ≈ [0.12,0.88]
//
//   INFORMATIVE (model variables of interest — theory-driven):
//   betaExR : normal(0.20, 0.15) — mild positive AR(1) persistence under CBN
//             managed float; OLS estimates for NGN/USD post-1999 typically 0.1–0.4
//   betaGDP : normal(-0.50, 0.40) — GDP growth reduces volatility (Aliyu 2009)
//   betaDebt: normal( 0.50, 0.40) — debt service raises volatility (Reinhart 2003)
// ============================================================

data {
  int<lower=1> T;
  vector[T] y;             // log return of NGN/USD (monthly)
  vector[T] y_lag;         // one-period lagged log return
  vector[T] gdp;           // standardised delta monthly GDP growth
  vector[T] debt;          // standardised delta monthly debt service %
  real<lower=0> h1_init;   // initial conditional variance (sample var of y)
}

parameters {
  // Mean equation
  real beta0;
  real betaExR;

  // Variance equation
  real omega;
  real alpha;
  real gamma;
  real phi_raw;   // unconstrained; phi = inv_logit(phi_raw) in (0,1)
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
    // Ceiling at 5.0: prevents h > exp(5) ≈ 148, i.e. cond. SD > 1200%/month.
    // Without this, high-phi draws cause the recursion to diverge numerically.
    log_h[t] = fmin(raw_log_h, 5.0);

    epsilon[t] = y[t] - (beta0 + betaExR * y_lag[t]);
  }

  for (t in 1:T)
    h[t] = fmax(exp(log_h[t]), 1e-10);
}

model {
  // ---- Weakly informative priors (system / GARCH parameters) ------
  beta0   ~ normal(0,     0.10);   // mean log return near zero
  omega   ~ normal(0,     2.00);   // log-variance intercept
  alpha   ~ normal(0,     0.50);   // symmetric ARCH effect
  gamma   ~ normal(0,     0.50);   // asymmetry / leverage
  phi_raw ~ normal(0,     1.00);   // persistence: inv_logit(0)=0.50 centre

  // ---- Informative priors (model variables of interest) -----------
  // Mild positive AR(1) persistence under CBN managed float (OLS 0.1–0.4)
  betaExR  ~ normal( 0.20, 0.15);
  // GDP growth reduces exchange rate uncertainty (Aliyu 2009 for Nigeria)
  betaGDP  ~ normal(-0.50, 0.40);
  // Rising debt service signals fiscal stress, increasing currency risk
  // (Reinhart, Rogoff, and Savastano 2003)
  betaDebt ~ normal( 0.50, 0.40);

  // ---- Likelihood ------------------------------------------
  for (t in 1:T)
    target += normal_lpdf(epsilon[t] | 0, sqrt(h[t]));
}

generated quantities {
  // Pointwise log-likelihood for WAIC / approximate LOO.
  // Note: standard LOO is not fully valid for time series
  // (non-exchangeable observations); treat as approximate.
  vector[T] log_lik;
  for (t in 1:T)
    log_lik[t] = normal_lpdf(epsilon[t] | 0, sqrt(h[t]));
}
