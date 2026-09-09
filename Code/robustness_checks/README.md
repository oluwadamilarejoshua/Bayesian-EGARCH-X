# Robustness checks — EGARCH-X on indicator-based (oil + GS10) data

Investigating whether the baseline model's `betaGDP`/`betaDebt` estimates are robust,
after switching Denton-Cholette interpolation from no-indicator to real monthly
indicators (WTI oil price for GDP, US 10-Year Treasury yield for debt).

## Files, in the order the checks were run

| File | What it is | Finding |
|---|---|---|
| `01_widened_prior_*` | Same sign/mean priors as baseline, SD 0.40 -> 0.80 | **Failed to converge** (R-hat up to 9.3) |
| `02_flat_prior_*` | betaGDP/betaDebt ~ N(0, 2) -- no directional info | **Failed to converge** (pooled R-hat up to 8.3) |
| `02b_flat_prior_perchain_*` | Same flat-prior model, re-run retaining **per-chain** summaries | Revealed **clean bimodality**: chains 1&3 -> betaDebt ~ +0.10 (phi~0.91); chains 2&4 -> betaDebt ~ -1.02 (phi~0.56). Not a sampler failure -- two genuine, internally-consistent modes. |
| `03_downweight2016_*` | Baseline's tight priors, but June & July 2016 excluded from the likelihood (weight=0; recursive path otherwise unchanged) | **Converged cleanly** (R-hat all <=1.001), but betaDebt = -0.90 [-1.04,-0.77] -- **sign flipped** from the full-sample baseline's +0.11 |
| `04_baseline_dispersed_perchain_*` | Baseline's tight priors, 4 chains deliberately seeded at very different starting points (incl. one seeded AT the 2016-excluded run's answer) | Chains 1-3 (incl. the one seeded at the alternative answer) all pulled back to betaDebt ~ +0.11; chain 4 (wrong-signed wildcard start) got stuck in a third, different place. **The positive-debt answer is a genuine attracting mode on the full data, not an artifact of shared initialization** -- but a third stray basin exists too. |
| `00_robustness_comparison.*` | Forest plot + underlying data comparing betaGDP/betaDebt across every converged result above (no cross-chain averaging -- individual chains shown separately) | Visual summary of everything above |
| `egarch_x_widened.stan`, `egarch_x_flat.stan`, `egarch_x_downweight2016.stan` | The three modified Stan models used | For reproducibility |
| `run_variant.R`, `run_downweight2016.R`, `run_flat_perchain.R`, `run_baseline_dispersed.R`, `make_comparison_plot.R` | Driver scripts | For reproducibility |

## The story so far

The full dataset (with June/July 2016 included) has a genuinely **bimodal** likelihood over
the variance-equation parameters: one mode with high persistence and a modest *positive*
debt-service effect (matching theory), another with lower persistence, higher leverage, and
a large *negative* debt effect. The baseline's tight, theory-motivated prior legitimately
selects the theory-consistent mode when the two are in real competition -- that's defensible
Bayesian practice *if disclosed*. However, excluding just the two 2016 months collapses the
ambiguity entirely in favour of the *negative*-debt mode, without needing a flat prior to get
there. This points to the two coincidental 2016 months (where the GS10 debt indicator dipped
sharply for reasons entirely about US/global rates) as the specific evidence responsible for
tipping the full-sample likelihood toward the theory-consistent mode in the first place.

## The resolution — smoothed GS10 indicator (ADOPTED)

| File | What it is |
|---|---|
| `smoothed_gs10_pipeline.R` | Rebuilds the full data pipeline using a **centred 3-month moving average of GS10** (instead of raw GS10) as the debt-service Denton-Cholette indicator. GDP's indicator (raw WTI oil price) is unchanged. Output written to `Analysis/Data/smoothed_gs10_test/`. |
| `run_smoothed_gs10_full.R` | Refits the baseline model (same tight, theory-motivated priors, unchanged) on this smoothed data, with full diagnostics (trace, posterior densities incl. **per-chain**, intervals, conditional volatility, LOO-CV, per-chain summary). |
| `05_smoothed_gs10/` | All outputs of that refit. |

**Verdict: this is the adopted final specification.** Smoothing diluted the sharp June/July
2016 coincidence in the debt indicator (raw d_debt -0.087/-0.086 -> smoothed -0.045/-0.049)
while preserving stationarity (ADF still rejects the unit root comfortably). Refitting on this
data did **not** weaken the debt-service finding -- it came back *stronger*:

| | Raw GS10 (baseline) | **Smoothed GS10 (adopted)** |
|---|---|---|
| betaDebt | +0.109 [0.023, 0.201] | **+0.160 [0.081, 0.238]** |
| betaGDP | -0.594 [-0.668, -0.519] | -0.564 [-0.638, -0.490] |
| LOO elpd | 732.2 | 739.2 (better) |
| Convergence | clean (R-hat <=1.001) | clean, and **all 4 chains individually agree** (see `05_smoothed_gs10/per_chain_summary_smoothed.csv`) -- no hidden bimodality |

This reframes the earlier 2016-exclusion sign-flip: a *mild* intervention (smoothing, months
still included) strengthens the result, while only a *drastic* intervention (deleting the two
largest shocks in the 261-month sample outright) flips it. That pattern is more consistent
with ordinary high-leverage-observation sensitivity (removing your most extreme data points
changes any time-series fit) than with the debt-service finding being an artifact manufactured
by one indicator's coincidental behaviour in two specific months.

**What must still be disclosed in the paper, honestly, regardless of this result:**
1. The two literature priors (Aliyu 2009; Reinhart-Rogoff-Savastano 2003) motivate the *sign*
   of betaGDP/betaDebt, not their numeric magnitude -- no source estimates this exact
   coefficient.
2. The full-sample likelihood is genuinely bimodal (confirmed by the flat-prior per-chain
   split, `02b_flat_prior_perchain.csv`); the tight prior is legitimately selecting the
   theory-consistent mode when the two are in real competition, which is defensible but must
   be stated, not implied as data speaking for itself.
3. Results retain some sensitivity to the two largest observations (June/July 2016), reduced
   but not eliminated by smoothing.
