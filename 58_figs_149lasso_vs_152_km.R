# =============================================================================
# 58_FigS_149Lasso_vs_152_KM.R   (v3: pure archival re-render, zero recomputation)
# Purpose: GSE65682 (n=468) two-panel KM supplementary figure -- no signature score is
#       recomputed; each panel reads the archive of its own validated analysis, so the
#       numbers match the manuscript digit for digit:
#   Panel A: 152 COPD signature  <- path15_GSE65682_train_scores.csv (original path15 analysis:
#                               score_geneset, whole-cohort median split, log-rank 0.020/HR 2.81)
#   Panel B: 149-Lasso signature <- path57_lasso_plotdata.rds (final version of 57, risk;
#                               lambda.min 12 genes, log-rank 0.206/HR 1.07)
# Consistency check: Panel B re-estimates Cox/HR of risk using time/status from the path15
#             archive; agreement with the cox_p/logrank_p archived by 57 (1.07/0.206)
#             confirms the samples are aligned.
# Output: FigS_Lasso149_vs_152_KM.pdf / path58_stats.txt / path58_plotdata.rds
# Run:   Rscript 58_FigS_149Lasso_vs_152_KM.R [OUT_DIR]   (only needs survival+ggplot2)
# =============================================================================
suppressPackageStartupMessages({
  library(survival); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
if (length(args) >= 1 && nzchar(args[1])) script_dir <- args[1]
OUT_DIR <- if (length(args) >= 2 && nzchar(args[2])) args[2] else script_dir
out_dir <- OUT_DIR
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("===== FigS: 149-Lasso vs 152 COPD signature (archival re-render) =====\n")
s15_file <- file.path(out_dir, "path15_GSE65682_train_scores.csv")
p57_file <- file.path(out_dir, "path57_lasso_plotdata.rds")
if (!file.exists(s15_file)) { cat("!! missing", s15_file, "\n"); quit(status = 1) }
if (!file.exists(p57_file)) { cat("!! missing", p57_file, "\n"); quit(status = 1) }

## Read the two archives
dfA <- read.csv(s15_file, stringsAsFactors = FALSE)          # time,status,copd,grp
p57 <- readRDS(p57_file)                                      # risk, cox_p, logrank_p...
risk <- as.numeric(p57$risk)
n <- length(risk)
cat("path15 scores rows:", nrow(dfA), "| path57 risk length:", n, "\n")
if (nrow(dfA) != n) { cat("!! length mismatch, aborting\n"); quit(status = 1) }

## Panel B re-estimates risk with time/status from path15 -> verify sample alignment
dfB <- data.frame(time = dfA$time, status = dfA$status,
                  risk = risk,
                  grp = ifelse(risk >= median(risk), "High", "Low"))
fitB <- coxph(Surv(time, status) ~ risk, data = dfB)
lrB  <- survdiff(Surv(time, status) ~ grp, data = dfB)
pB_cox <- summary(fitB)$coefficients[5]
pB_lr  <- 1 - pchisq(lrB$chisq, 1)
cat("alignment check: Panel B re-estimated Cox HR =", round(exp(coef(fitB)), 3),
    "| expected (archive of 57) HR = 1.071, cox_p =", format(p57$cox_p, digits = 3),
    "| re-estimated logrank P =", format(pB_lr, digits = 3),
    "| expected logrank =", format(p57$logrank_p, digits = 3), "\n")

## Panel A statistics (computed directly from the archived df; must match the original path15)
dfA <- dfA[!is.na(dfA$copd), ]
fitA <- coxph(Surv(time, status) ~ copd, data = dfA)
lrA  <- survdiff(Surv(time, status) ~ grp, data = dfA)
csA  <- summary(fitA)$conf.int
pA_lr <- 1 - pchisq(lrA$chisq, 1)
cat("Panel A (152) Cox HR =", round(csA[1], 3), "| log-rank P =",
    format(pA_lr, digits = 3), "(original path15: HR 2.81, P 0.020)\n")

## KM curves (extracted via summary)
km_curve <- function(d) {
  fit <- survfit(Surv(time, status) ~ grp, data = d)
  sf  <- summary(fit)
  g   <- sub("^.*=", "", as.character(sf$strata))
  data.frame(t = sf$time, s = sf$surv,
             grp = factor(g, levels = c("High", "Low")))
}
cA <- km_curve(dfA); cB <- km_curve(dfB)
cA$sig <- "A"; cB$sig <- "B"
curves <- rbind(cA, cB)

lab_df <- data.frame(
  sig = c("A", "B"),
  lab = c(sprintf("COPD signature (152 genes)\nlog-rank P = %.3f, HR = %.2f",
                  pA_lr, csA[1]),
          sprintf("Lasso signature from 149-gene pool\n12 genes (lambda.min)\nlog-rank P = %.3f, HR = %.2f",
                  pB_lr, exp(coef(fitB))))
)
curves$sig2 <- factor(ifelse(curves$sig == "A", lab_df$lab[1], lab_df$lab[2]),
                      levels = lab_df$lab)

p <- ggplot(curves, aes(t, s, colour = grp)) +
  geom_step(linewidth = 0.9) +
  facet_wrap(~ sig2, ncol = 2) +
  scale_colour_manual(values = c("High" = "#B2182B", "Low" = "#2166AC"),
                      breaks = c("High", "Low"), name = "Risk group",
                      drop = FALSE) +
  scale_x_continuous(limits = c(0, NA)) +
  labs(x = "Days", y = "Survival probability") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 9.5, lineheight = 0.95),
        legend.position = "bottom")
ggsave(file.path(out_dir, "FigS_Lasso149_vs_152_KM.pdf"), p,
       width = 10, height = 4.8)

stats_txt <- sprintf(
  "FigS: 149-Lasso vs 152 COPD signature, GSE65682 (n=%d)\n", n)
stats_txt <- paste0(stats_txt, sprintf(
  "A) 152 signature (archival re-render of path15): Cox HR=%.3f (%.3f-%.3f) | log-rank P=%.3e\n",
  csA[1], csA[3], csA[4], pA_lr))
stats_txt <- paste0(stats_txt, sprintf(
  "B) 149-Lasso 12 genes (archival re-render of 57 plotdata): Cox HR=%.3f (%.3f-%.3f) | log-rank P=%.3e\n",
  exp(coef(fitB)), summary(fitB)$conf.int[3], summary(fitB)$conf.int[4], pB_lr))
writeLines(stats_txt, file.path(out_dir, "path58_stats.txt"))

saveRDS(list(curves = curves, dfA = dfA, dfB = dfB,
             stA = list(HR = csA[1], lo = csA[3], hi = csA[4], p_lr = pA_lr),
             stB = list(HR = exp(coef(fitB)), p_cox = pB_cox, p_lr = pB_lr)),
        file.path(out_dir, "path58_plotdata.rds"))
cat("\noutput: FigS_Lasso149_vs_152_KM.pdf / path58_stats.txt / path58_plotdata.rds\n")
cat("===== script 58 done =====\n")
