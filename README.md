# Bayesian EGARCH-X: GDP and Debt Servicing Effects on Nigerian Naira Volatility

Replication materials for:

> Ogundairo, J. O., Adewale, A., & Ojo, O. O. (2025). *Impacts of GDP and Debt Servicing on Exchange Rate Volatility of the Nigerian Naira: A Bayesian EGARCH-X Approach.*

## Repository Structure

```
Code/
  analysis_egarchx.R          # Main estimation script (RStan)
  egarch_x_model.stan         # Bayesian EGARCH(1,1)-X Stan model definition
Data/
  nigeria_monthly_stan_input_post1999.csv   # Estimation sample (Mar 1999–Dec 2020, T=262)
  nigeria_monthly_merged_post1999.csv       # Merged raw monthly series
Results/
  posterior_summary.csv       # Posterior means, 95% CIs, n_eff, R-hat (all 8 parameters)
  conditional_volatility.csv  # Posterior median and 95% credible bands for h_t
  loo_cv_metrics.csv          # LOO-CV: elpd_loo, p_loo, LOOIC
Figures/
  conditional_volatility_plot.pdf   # Figure 1: Posterior median conditional volatility
  trace_variance_eq.pdf             # Figure 2: MCMC trace plots (variance-equation parameters)
  posterior_densities.pdf           # Figure 3: Posterior marginal density estimates
```

## Model

The study estimates a Bayesian Exponential GARCH model augmented with macro-fiscal regressors in the **variance equation** (EGARCH-X). GDP growth and external debt service enter as structural determinants of exchange rate risk — not the conditional mean — using monthly NGN/USD log returns over the post-managed-float era.

**Variance equation:**

$$\ln h_t = \omega + \alpha\left(\frac{|\varepsilon_{t-1}|}{\sqrt{h_{t-1}}} - \sqrt{\frac{2}{\pi}}\right) + \gamma\frac{\varepsilon_{t-1}}{\sqrt{h_{t-1}}} + \varphi\ln h_{t-1} + \beta_g\,\widetilde{\text{GDP}}_t + \beta_d\,\widetilde{\text{Debt}}_t$$

Posterior inference is conducted via the No-U-Turn Sampler (NUTS) implemented in RStan (4 chains × 3,000 post-warmup draws = 12,000 total samples).

## Requirements

- R (≥ 4.0) with packages: `rstan`, `bayesplot`, `loo`, `dplyr`, `ggplot2`, `tempdisagg`
- [Rtools](https://cran.r-project.org/bin/windows/Rtools/) (Windows only — required for Stan C++ compilation)

## Data Sources

| Series | WDI Code |
|--------|----------|
| NGN/USD exchange rate (official, period average) | `PA.NUS.FCRF` |
| GDP growth (annual %) | `NY.GDP.MKTP.KD.ZG` |
| External debt service (% of GNI) | `DT.TDS.DECT.GN.ZS` |

Source: [World Development Indicators](https://databank.worldbank.org/source/world-development-indicators), World Bank.

Annual GDP and debt series are disaggregated to monthly frequency via the Denton–Cholette quadratic minimisation method (`tempdisagg` R package), then first-differenced and standardised before entering the model.

## Sample

Post-managed-float era: **March 1999 – December 2020**  
T = 262 raw observations; T = 261 used in estimation (one observation lost to differencing and lagging).

## How to Reproduce

1. Clone this repository.
2. Open `Code/analysis_egarchx.R` in RStudio.
3. Verify that `data_path` resolves to `Data/nigeria_monthly_stan_input_post1999.csv` relative to the script location (the script auto-detects this).
4. Run the script. The Stan model in `Code/egarch_x_model.stan` compiles automatically on first run.
5. Outputs are saved to the working directory (match `Results/` and `Figures/` contents in this repo).

## Key Results

| Parameter | Mean | 95% CI |
|-----------|------|--------|
| β₀ (intercept) | 0.001 | [−0.001, 0.002] |
| β₁ (AR(1)) | 0.357 | [0.178, 0.535] |
| ω (log-var intercept) | −1.711 | [−2.433, −1.085] |
| α (ARCH) | 0.450 | [0.280, 0.620] |
| γ (leverage) | 0.221 | [0.055, 0.391] |
| φ (persistence) | 0.780 | [0.696, 0.851] |
| βg (GDP) | −0.145 | [−0.231, −0.059] |
| βd (Debt service) | 0.275 | [0.163, 0.414] |

LOO-CV: elpd_loo = 680.2 (SE = 61.6); LOOIC = −1360.3. All R̂ ≤ 1.01.
