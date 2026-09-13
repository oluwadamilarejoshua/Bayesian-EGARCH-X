```
Bayesian-EGARCH-X/
├── Code/
│   ├── analysis_egarchx.R                          # original baseline (Gaussian) model
│   ├── egarch_x_model.stan
│   └── robustness_checks/
│       ├── 06_studentt/                            # adopted model: AR(1) mean eq, Student-t innovations
│       │   ├── run_studentt.R
│       │   ├── egarch_x_studentt.stan
│       │   └── residual_diagnostics_studentt.R      # ARCH-LM + Ljung-Box verification (raw and winsorized)
│       ├── 07_ar2_studentt/                        # robustness check: AR(2) mean equation
│       │   ├── run_ar2_studentt.R
│       │   └── egarch_x_ar2_studentt.stan
│       └── 08_studentt_downweight2016/             # robustness check: 2016 float excluded from likelihood
│           ├── run_studentt_downweight2016.R
│           └── egarch_x_studentt_downweight2016.stan
├── Data/
│   ├── nigeria_monthly_merged_post1999.csv
│   └── nigeria_monthly_stan_input_post1999.csv
├── Figures/
│   ├── conditional_volatility_plot.pdf              # original baseline model
│   ├── posterior_densities.pdf
│   ├── trace_variance_eq.pdf
│   └── robustness_checks/
│       └── 06_studentt/                            # figures for the adopted model, as used in the manuscript
│           ├── conditional_volatility_plot_studentt.pdf
│           ├── trace_variance_eq_studentt.pdf
│           └── posterior_densities_studentt.pdf
└── Results/
    ├── conditional_volatility.csv                   # original baseline model
    ├── loo_cv_metrics.csv
    ├── posterior_summary.csv
    └── robustness_checks/
        ├── 06_studentt/                             # adopted model: posterior summary, LOO-CV, residuals
        ├── 07_ar2_studentt/                         # AR(2) comparison
        └── 08_studentt_downweight2016/              # 2016-exclusion comparison
```

## Adopted model

The model reported in the manuscript is `Code/robustness_checks/06_studentt/`: an EGARCH(1,1)-X
with an AR(1) mean equation and Student-t innovations, estimated on the smoothed-GS10 debt
indicator. `07_ar2_studentt/` and `08_studentt_downweight2016/` are the two robustness checks
reported alongside it (mean-equation lag length, and sensitivity to the June/July 2016
devaluation). `Code/analysis_egarchx.R` and `Code/egarch_x_model.stan` are an earlier baseline
specification (Gaussian innovations) kept for reference; it is not the model reported in the
manuscript.

Each `run_*.R` script is self-contained: it reads the data, fits the model in RStan, prints
posterior summaries and diagnostics, and writes its outputs to the matching `Results/` and
`Figures/` subfolder. The very large (~48MB) fitted-model `.rds` object used to regenerate the
trace and density plots is not included in this repository; re-running `run_studentt.R`
regenerates it locally.
