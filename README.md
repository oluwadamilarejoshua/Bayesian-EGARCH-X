# Bayesian EGARCH-X: Nigerian Naira Exchange Rate Volatility

Execution scripts for *"Impacts of GDP and Debt Servicing on Exchange Rate Volatility of the Nigerian Naira: A Bayesian EGARCH-X Approach"* (Ogundairo, Adewale & Ojo).

This repository contains **code only** for the main pipeline — no bundled data, no result files. Scripts in `Code/` are numbered `01`–`05` in the order they must be run; each one fetches or derives its own inputs and writes its own outputs to local folders (`Data/`, `Preliminary_Tests/`, `Figures/`, `Results/`) that are created on first run and are not tracked in version control. The one deliberate exception is `Code/robustness_checks/` — small summary CSVs, per-chain diagnostics, and plots documenting a specification-robustness investigation are kept in version control there, since they are the evidence trail behind a real methodological decision (see below), not easily-regenerable intermediate data.

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
    ├── egarch_x_model.stan          Stan model definition (used by 05)
    └── robustness_checks/           Specification-robustness investigation (see below)
```

## Pipeline, in order

1. **`01_data_fetch_monthly.R`** — Downloads the monthly official NGN/USD exchange rate (IMF IFS via `rdbnomics`, with a FRED fallback) and the annual WDI GDP growth and debt-service series (via the `WDI` API). Disaggregates the annual GDP/debt series to monthly frequency using the **Denton-Cholette** temporal disaggregation method (`tempdisagg`, Denton 1971; Sax & Steiner 2013) *with a genuine monthly indicator series for each* — the WTI oil price for GDP (Nigeria's output is heavily oil-driven) and a **smoothed (centred 3-month moving average) US 10-Year Treasury yield** for debt service (most of Nigeria's external debt is dollar-denominated, so global rates drive its servicing cost, and the series is fully exogenous to the naira market). Three prior versions of the debt indicator were tried and superseded, each for a documented reason — see `Code/robustness_checks/README.md` for the full investigation: (i) no indicator at all (`~ 1`), which produced pathologically smooth series that confused ADF unit-root tests; (ii) Nigeria's own FX reserves, dropped because reserves are the CBN's own instrument for defending the naira band (circularity risk); (iii) raw (unsmoothed) GS10, which worked well but turned out to have a sharp coincidental dip in June-July 2016 (the same two months as the naira's largest shock in the sample) that made the debt-service coefficient's sign sensitive to those two observations. The current smoothed-GS10 version dilutes that coincidence and was confirmed to *strengthen*, not weaken, the debt-service finding. Writes the full 1980-2020 merged series and a provisional post-1999 subset (superseded by step 2) to `Data/processed/` (created automatically).

2. **`02_make_post1999.R`** — The full 1980-2020 sample cannot be used for EGARCH estimation directly: it contains 67 consecutive zero-return months from the 1993-1998 Abacha fixed-rate peg (degenerate for variance estimation) and two 20-30 SD outlier months (the 1986 SAP devaluation and the 1999 end-of-peg jump) that would dominate the variance equation. This script restricts the sample to **March 1999 - December 2020** (the continuous managed-float era) and first-differences the interpolated GDP and debt series, since their levels test as I(1) — an artifact of the smooth Denton-Cholette interpolation. The differenced series (`d_gdp`, `d_debt`) have a clean interpretation as *growth momentum* and *debt-service momentum*. Output: `nigeria_monthly_stan_input_post1999.csv`, the file the model actually consumes.

3. **`03_stationarity_tests.R`** — Runs ADF (`none`/`drift`/`trend` specifications), Phillips-Perron (`constant`/`trend`), and KPSS (`mu`/`tau`) on all three final series (`y`, `gdp`, `debt`), cross-checked against a second implementation (`tseries`). ADF and PP test H0: unit root; KPSS tests the opposite H0: stationarity — using both lets one test's blind spot get caught by the other. Result: all three series are unambiguously I(0) on all 7 tests. (With the earlier no-indicator interpolation, ADF specifically had trouble with the GDP/debt series — PP and KPSS still passed, but ADF's parametric lag-augmentation couldn't fully whiten the artificial smoothness that a no-indicator Denton-Cholette run introduces. Adding real monthly indicators in step 1 resolved this cleanly.)

4. **`04_structural_break_tests.R`** — Runs Zivot-Andrews (endogenous single-break unit-root test), Bai-Perron (BIC-selected multiple-break test), a Chow test at the a priori June 2016 CBN-float date, and an OLS-CUSUM stability test, all on `y`. Zivot-Andrews rejects the unit root even under its worst-case endogenous break (found at October 2014, the oil-price collapse); Bai-Perron selects zero breaks in the mean; the Chow test at June 2016 is not significant (p = 0.25). Conclusion: no break dummy is needed in the conditional **mean** equation — the 2016 devaluation is expected to appear instead as a spike in conditional **variance**, which the EGARCH leverage term is designed to capture. This is the empirical basis for placing GDP and debt in the variance equation rather than the mean.

5. **`05_analysis_egarchx.R`** — Standardizes the GDP/debt regressors (mean 0, SD 1) and fits the Bayesian EGARCH(1,1)-X model (`egarch_x_model.stan`) via NUTS in RStan: 4 chains, 1,000 warmup + 3,000 post-warmup draws each, `adapt_delta = 0.92`, `max_treedepth = 15`. Reports posterior summaries, R-hat/effective-sample-size diagnostics, and LOO-CV (`loo` package); saves trace/density/interval plots to `Figures/` and posterior summaries, raw draws, and conditional-volatility estimates to `Results/`.

To reproduce end to end, run the five `Code/` scripts in numeric order from an R session (RStudio recommended — the scripts locate their own paths via `rstudioapi`, with a fallback for `Rscript` command-line execution). Each script creates the local output folders it needs; none of those folders are committed to this repository.

## Robustness checks

`Code/robustness_checks/` documents a specification-robustness investigation carried out after
switching to indicator-based interpolation, prompted by a simple question: is the debt-service
finding stable, or an artifact of specific modelling choices? In order:

1. Widened and fully flat priors on `betaGDP`/`betaDebt` were tried. The flat-prior run,
   inspected **per chain** rather than pooled, revealed the model's likelihood is genuinely
   **bimodal** — one mode matching theory (high persistence, positive debt effect), another
   with lower persistence and a negative debt effect. Both are real, internally-consistent
   explanations of the data.
2. A dispersed-initialization check (chains deliberately seeded at each mode's values, under
   the model's actual tight priors) showed the theory-consistent mode is a genuine attractor
   on the full dataset — not an artifact of always starting the sampler from the same point.
3. Excluding just the two largest-shock months (June-July 2016) from the likelihood flipped
   the debt coefficient's sign entirely, raising a real concern about high-leverage sensitivity.
4. Smoothing the debt indicator (see above) rather than deleting those months tested this more
   gently — and the debt-service finding came back *stronger*, with a tighter, more credible
   interval and better out-of-sample fit (LOO-CV). This is why the smoothed indicator was
   adopted as the final specification.

See `Code/robustness_checks/README.md` for the full narrative, every intermediate result, and
the exact numbers at each step.

## Data sources

| Source | Series | Frequency | Coverage |
|---|---|---|---|
| World Bank WDI | `PA.NUS.FCRF` (official exchange rate) | Annual | 1980-2020 |
| World Bank WDI | `NY.GDP.MKTP.KD.ZG` (GDP growth) | Annual | 1980-2020 |
| World Bank WDI | `DT.TDS.DECT.EX.ZS` (debt service, % exports) | Annual | 1980-2020 |
| IMF IFS via DBnomics (primary) | `M.NG.ENDA_XDC_USD_RATE` | Monthly | 1957-present |
| FRED (fallback) | `FXRATENGA618NUPN` / `XRNCUSNGA618NRUG` | Monthly | 1950-2010/2023 |
| FRED | `WTISPLC` (WTI crude oil price) — Denton-Cholette indicator for GDP | Monthly | 1946-present |
| FRED | `GS10` (US 10-Year Treasury yield) — Denton-Cholette indicator for debt service | Monthly | 1953-present |

`01_data_fetch_monthly.R` pulls all of the above live via API on each run (`rdbnomics`, `quantmod`/FRED, `WDI`) — no data files are bundled with this repository.

## Priors

Weakly informative priors are used for the EGARCH system parameters (β₀, ω, α, γ, φ). Informative priors on the three parameters of substantive interest are grounded in the literature:

- **β_ExR** (AR(1) persistence) ~ N(0.20, 0.15²) — consistent with prior OLS estimates (0.1-0.4) for post-1999 NGN/USD returns.
- **β_GDP** (GDP growth's effect on log conditional variance) ~ N(-0.50, 0.40²) — direction motivated by Aliyu, S.U.R. (2009), "Impact of Oil Price Shock and Exchange Rate Volatility on Economic Growth in Nigeria," *Research Journal of International Studies*, 11, 4-15.
- **β_Debt** (debt service's effect on log conditional variance) ~ N(0.50, 0.40²) — direction motivated by Reinhart, C.M., Rogoff, K.S., & Savastano, M.A. (2003), "Debt Intolerance," *Brookings Papers on Economic Activity*, 2003(1).

Both citations motivate the *sign* of these priors; neither estimates this exact coefficient, so the numeric magnitude (mean and SD) is a subjective, theory-motivated choice rather than an empirical estimate transplanted from either source. This matters because the likelihood these priors are combined with is genuinely bimodal (see Robustness checks below) — the tight prior is doing real, disclosed work selecting between two statistically competing explanations of the data, not merely refining an already-clear answer.

The EGARCH functional form follows Nelson, D.B. (1991), "Conditional Heteroskedasticity in Asset Returns: A New Approach," *Econometrica*, 59(2), 347-370. Posterior sampling uses the No-U-Turn Sampler (Hoffman & Gelman, 2014, *JMLR*, 15, 1593-1623) via RStan.

## Requirements

R packages: `rstan`, `bayesplot`, `loo`, `dplyr`, `ggplot2`, `urca`, `tseries`, `strucchange`, `tempdisagg`, `WDI`, `rdbnomics`, `quantmod`, `lubridate`, `zoo`. On Windows, Stan requires Rtools (see the PATH guard at the top of `05_analysis_egarchx.R`).
