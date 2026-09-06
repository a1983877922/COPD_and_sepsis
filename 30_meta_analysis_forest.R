## =========================================================================
## 30_meta_analysis_forest.R
## Three-cohort COPD inflammatory signature -> 28-day mortality meta-analysis
##   (formal recalculation + forest plot)
##
## Purpose:
##   1) Recalculate the three-cohort fixed/random effects pooled HR with meta::metagen
##      (previously logHR/SE were back-derived manually in Python, giving pooled HR 2.68,
##      P=0.019; this script performs the formal check);
##   2) Produce the forest plot, echoing Fig 6d (prognostic forest plot);
##   3) Report heterogeneity tests (Q / I^2) and directional consistency as
##      standard statistics for submission.
##
## Input: univariable Cox HR and 95% CI for the three cohorts (cox_stats from path14/16/24)
## Output: path30_meta_stats.txt / path30_meta_forest.pdf
## =========================================================================

out_dir <- "."   # run from the script directory

## ---- Dependency installation (meta package) ----
if (!requireNamespace("meta", quietly = TRUE)) {
  install.packages("meta", repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(meta))

## ---- Three-cohort data (COPD inflammatory signature, univariable Cox) ----
## GSE65682:   HR 2.81 (0.99-7.97), n=468 (28-day follow-up, 103 deaths)
## E-MTAB-4451:HR 2.00 (0.47-8.55), n=106
## GSE95233:   HR 8.01 (0.27-240.4), n=51
study <- c("GSE65682", "E-MTAB-4451", "GSE95233")
hr  <- c(2.81, 2.00, 8.01)
lo  <- c(0.99, 0.47, 0.27)
hi  <- c(7.97, 8.55, 240.4)
n   <- c(468, 106, 51)

## ---- metagen (standard inverse-variance, fixed + random effects) ----
m <- metagen(
  TE        = log(hr),
  seTE      = (log(hi) - log(lo)) / (2 * qnorm(0.975)),
  studlab   = study,
  sm        = "HR",
  fixed     = TRUE,
  random    = TRUE,
  method.tau = "DL",       # DerSimonian-Laird
  n.e       = n
)

## ---- Write results to text ----
sink(file.path(out_dir, "path30_meta_stats.txt"))
cat("=========== Three-cohort COPD inflammatory signature meta-analysis ===========\n\n")
print(summary(m))
cat("\n\n----- Human-readable summary -----\n")
cat(sprintf("Fixed effects: HR = %.3f (95%% CI %.3f-%.3f), P = %.4f\n",
            exp(m$TE.fixed), exp(m$lower.fixed), exp(m$upper.fixed), m$pval.fixed))
cat(sprintf("Random effects: HR = %.3f (95%% CI %.3f-%.3f), P = %.4f\n",
            exp(m$TE.random), exp(m$lower.random), exp(m$upper.random), m$pval.random))
cat(sprintf("Heterogeneity: Q = %.3f, df = %d, P = %.4f, I^2 = %.1f%%\n",
            m$Q, m$df.Q, m$pval.Q, m$I2))
cat(sprintf("tau^2 = %.4f\n", m$tau2))
sink()

## ---- Forest plot ----
pdf(file.path(out_dir, "path30_meta_forest.pdf"), width = 8, height = 4.5)
forest(m,
       xlab = "Hazard ratio (log scale)",
       col.square = "#d62728",
       col.diamond = "#1f77b4",
       col.diamond.lines = "#1f77b4",
       overall = TRUE,
       overall.hetstat = TRUE,
       fontsize = 10)
dev.off()

cat("Done. Output: path30_meta_stats.txt / path30_meta_forest.pdf\n")
