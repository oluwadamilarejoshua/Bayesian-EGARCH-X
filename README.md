# Bayesian EGARCH-X: Nigerian Naira Exchange Rate Volatility

Execution scripts for *"Impacts of GDP and Debt Servicing on Exchange Rate Volatility of the Nigerian Naira: A Bayesian EGARCH-X Approach"* (Ogundairo, Adewale & Ojo).

This repository contains **code only** — no bundled data, no result files. Scripts in `Code/` are numbered `01`–`05` in the order they must be run; each one fetches or derives its own inputs and writes its own outputs to local folders (`Data/`, `Preliminary_Tests/`, `Figures/`, `Results/`) that are created on first run and are not tracked in version control.

## Repository structure

```
Bayesian-EGARCH-X/
└── Code/
    ├── 01_data_fetch_monthly.R      Fetch monthly FX (IMF IFS/FRED) + annual WDI
    │                                 GDP/debt; Denton-Cholette interpolation to monthly
    ├── 02_make_post1999.R           Restrict to post-1999 managed-float era;
    │                                 first-difference GDP/debt to achieve stationarity
    ├── 03_stationarity_tests.R      ADF, Phillips-Perron, KPSS unit-root tests
    ├── 04_structural_break_tests.R  Zivot-Andrews, Bai-Perron, Chow, CUSUM tests
    ├── 05_analysis_egarchx.R        Fit the Bayesian EGARCH(1,1)-X model in Stan
    └── egarch_x_model.stan          Stan model definition (used by 05)
```

## Pipeline, in order

1. **`01_data_fetch_monthly.R`** — Downloads the monthly official NGN/USD exchange rate (IMF IFS via `rdbnomics`, with a FRED fallback) and the annual WDI GDP growth and debt-service series (via the `WDI` API). Disaggregates the annual GDP/debt series to monthly frequency using the **Denton-Cholette** temporal disaggregation method (`tempdisagg`, Denton 1971; Sax & Steiner 2013), which distributes annual totals to monthly values while minimizing period-to-period movement, so no spurious volatility is introduced by the interpolation itself. Writes the full 1980-2020 merged series and a provisional post-1999 subset (superseded by step 2) to `Data/processed/` (created automatically).

2. **`02_make_post1999.R`** — The full 1980-2020 sample cannot be used for EGARCH estimation directly: it contains 67 consecutive zero-return months from the 1993-1998 Abacha fixed-rate peg (degenerate for variance estimation) and two 20-30 SD outlier months (the 1986 SAP devaluation and the 1999 end-of-peg jump) that would dominate the variance equation. This script restricts the sample to **March 1999 - December 2020** (the continuous managed-float era) and first-differences the interpolated GDP and debt series, since their levels test as I(1) — an artifact of the smooth Denton-Cholette interpolation. The differenced series (`d_gdp`, `d_debt`) have a clean interpretation as *growth momentum* and *debt-service momentum*. Output: `nigeria_monthly_stan_input_post1999.csv`, the file the model actually consumes.

3. **`03_stationarity_tests.R`** — Runs ADF (`none`/`drift`/`trend` specifications), Phillips-Perron (`constant`/`trend`), and KPSS (`mu`/`tau`) on all three final series (`y`, `gdp`, `debt`), cross-checked against a second implementation (`tseries`). ADF and PP test H0: unit root; KPSS tests the opposite H0: stationarity — using both lets one test's blind spot get caught by the other. Result: `y` is unambiguously I(0) on all 7 tests; the differenced GDP/debt series pass PP and KPSS but show some ADF ambiguity attributable to autocorrelation induced by the Denton-Cholette interpolation. All three series are treated as I(0) in estimation.

4. **`04_structural_break_tests.R`** — Runs Zivot-Andrews (endogenous single-break unit-root test), Bai-Perron (BIC-selected multiple-break test), a Chow test at the a priori June 2016 CBN-float date, and an OLS-CUSUM stability test, all on `y`. Zivot-Andrews rejects the unit root even under its worst-case endogenous break (found at October 2014, the oil-price collapse); Bai-Perron selects zero breaks in the mean; the Chow test at June 2016 is not significant (p = 0.25). Conclusion: no break dummy is needed in the conditional **mean** equation — the 2016 devaluation is expected to appear instead as a spike in conditional **variance**, which the EGARCH leverage term is designed to capture. This is the empirical basis for placing GDP and debt in the variance equation rather than the mean.

5. **`05_analysis_egarchx.R`** — Standardizes the GDP/debt regressors (mean 0, SD 1) and fits the Bayesian EGARCH(1,1)-X model (`egarch_x_model.stan`) via NUTS in RStan: 4 chains, 1,000 warmup + 3,000 post-warmup draws each, `adapt_delta = 0.92`, `max_treedepth = 15`. Reports posterior summaries, R-hat/effective-sample-size diagnostics, and LOO-CV (`loo` package); saves trace/density/interval plots to `Figures/` and posterior summaries, raw draws, and conditional-volatility estimates to `Results/`.

To reproduce end to end, run the five `Code/` scripts in numeric order from an R session (RStudio recommended — the scripts locate their own paths via `rstudioapi`, with a fallback for `Rscript` command-line execution). Each script creates the local output folders it needs; none of those folders are committed to this repository.

## Data sources

| Source | Series | Frequency | Coverage |
|---|---|---|---|
| World Bank WDI | `PA.NUS.FCRF` (official exchange rate) | Annual | 1980-2020 |
| World Bank WDI | `NY.GDP.MKTP.KD.ZG` (GDP growth) | Annual | 1980-2020 |
| World Bank WDI | `DT.TDS.DECT.EX.ZS` (debt service, % exports) | Annual | 1980-2020 |
| IMF IFS via DBnomics (primary) | `M.NG.ENDA_XDC_USD_RATE` | Monthly | 1957-present |
| FRED (fallback) | `FXRATENGA618NUPN` / `XRNCUSNGA618NRUG` | Monthly | 1950-2010/2023 |

`01_data_fetch_monthly.R` pulls all of the above live via API on each run (`rdbnomics`, `quantmod`/FRED, `WDI`) — no data files are bundled with this repository.

## Priors

Weakly informative priors are used for the EGARCH system parameters (β₀, ω, α, γ, φ). Informative priors on the three parameters of substantive interest are grounded in the literature:

- **β_ExR** (AR(1) persistence) ~ N(0.20, 0.15²) — consistent with prior OLS estimates (0.1-0.4) for post-1999 NGN/USD returns.
- **β_GDP** (GDP growth's effect on log conditional variance) ~ N(-0.50, 0.40²) — Aliyu, S.U.R. (2009), "Impact of Oil Price Shock and Exchange Rate Volatility on Economic Growth in Nigeria," *Research Journal of International Studies*, 11, 4-15.
- **β_Debt** (debt service's effect on log conditional variance) ~ N(0.50, 0.40²) — Reinhart, C.M., Rogoff, K.S., & Savastano, M.A. (2003), "Debt Intolerance," *Brookings Papers on Economic Activity*, 2003(1).

The EGARCH functional form follows Nelson, D.B. (1991), "Conditional Heteroskedasticity in Asset Returns: A New Approach," *Econometrica*, 59(2), 347-370. Posterior sampling uses the No-U-Turn Sampler (Hoffman & Gelman, 2014, *JMLR*, 15, 1593-1623) via RStan.

## Requirements

R packages: `rstan`, `bayesplot`, `loo`, `dplyr`, `ggplot2`, `urca`, `tseries`, `strucchange`, `tempdisagg`, `WDI`, `rdbnomics`, `quantmod`, `lubridate`, `zoo`. On Windows, Stan requires Rtools (see the PATH guard at the top of `05_analysis_egarchx.R`).
