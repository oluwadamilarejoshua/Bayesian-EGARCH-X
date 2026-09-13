# ============================================================
# make_comparison_plot.R
# Builds a transparent, no-aggregation forest plot of betaGDP and
# betaDebt across every robustness check run so far. Every row is
# either a full pooled fit (only shown when it actually converged)
# or one individual chain's per-chain summary (for the two runs
# where we deliberately inspected per-chain behaviour). Nothing is
# averaged across chains -- averaging two different modes together
# would misrepresent what was actually found.
# ============================================================

suppressMessages({ library(ggplot2); library(dplyr) })

BEG <- "c:/Users/joshua.ogundairo/Downloads/My Research Works/Impacts_of_GDP_and_Debt_Servicing_on_Exchange_Rate_Volatility_of_the_Nigerian_Naira__A_Bayesian_EGARCH_Approach/Analysis/Bayesian-EGARCH"
RC  <- file.path(BEG, "robustness_checks")

baseline   <- read.csv(file.path(BEG, "posterior_summary.csv"))
dw2016     <- read.csv(file.path(RC, "03_downweight2016_summary.csv"))
flat_pc    <- read.csv(file.path(RC, "02b_flat_prior_perchain.csv"))
disp_pc    <- read.csv(file.path(RC, "04_baseline_dispersed_perchain.csv"))

rows <- list()

add_row <- function(label, mean, lo, hi, group, converged = TRUE) {
  rows[[length(rows) + 1]] <<- data.frame(label = label, mean = mean, lo = lo, hi = hi,
                                           group = group, converged = converged)
}

for (par in c("betaGDP", "betaDebt")) {
  b <- baseline[baseline$parameter == par, ]
  add_row(paste0(par, " | baseline (tight prior, full data)"), b$mean, b$ci_2.5, b$ci_97.5, "baseline", TRUE)

  d <- dw2016[dw2016$parameter == par, ]
  add_row(paste0(par, " | 2016 excluded (tight prior)"), d$mean, d$ci_2.5, d$ci_97.5, "2016-excluded", TRUE)

  for (ch in 1:4) {
    r <- flat_pc[flat_pc$parameter == par & flat_pc$chain == ch, ]
    mode <- if (ch %in% c(1, 3)) "flat prior -- Mode A (theory-consistent)" else "flat prior -- Mode B (opposite sign)"
    add_row(paste0(par, " | flat prior, chain ", ch), r$mean, r$q2.5, r$q97.5, mode, TRUE)
  }

  for (ch in 1:4) {
    r <- disp_pc[disp_pc$parameter == par & disp_pc$chain == ch, ]
    grp <- if (ch %in% c(1, 2, 3)) "dispersed-init baseline -- Mode A" else "dispersed-init baseline -- stray chain"
    add_row(paste0(par, " | dispersed-init baseline, chain ", ch), r$mean, r$q2.5, r$q97.5, grp, TRUE)
  }
}

df <- do.call(rbind, rows)
df$param <- ifelse(grepl("^betaGDP", df$label), "betaGDP (GDP -> volatility)", "betaDebt (Debt -> volatility)")
df$label <- factor(df$label, levels = rev(df$label))

pal <- c(
  "baseline" = "#1b7837",
  "2016-excluded" = "#762a83",
  "flat prior -- Mode A (theory-consistent)" = "#2166ac",
  "flat prior -- Mode B (opposite sign)" = "#b2182b",
  "dispersed-init baseline -- Mode A" = "#2166ac",
  "dispersed-init baseline -- stray chain" = "#f4a582"
)

p <- ggplot(df, aes(x = mean, y = label, xmin = lo, xmax = hi, color = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(size = 0.5) +
  scale_color_manual(values = pal, name = "Run / mode") +
  facet_wrap(~param, scales = "free", ncol = 1) +
  labs(
    title = "Robustness of betaGDP and betaDebt across all checks run so far",
    subtitle = "No values are averaged across chains -- every row is either a converged pooled fit or one individual chain",
    x = "Posterior mean (95% CI)", y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", legend.direction = "vertical",
        plot.title = element_text(face = "bold"))

out_pdf <- file.path(RC, "00_robustness_comparison.pdf")
ggsave(out_pdf, p, width = 10, height = 9)
cat("Saved:", out_pdf, "\n")

write.csv(df, file.path(RC, "00_robustness_comparison_data.csv"), row.names = FALSE)
cat("Saved:", file.path(RC, "00_robustness_comparison_data.csv"), "\n")
