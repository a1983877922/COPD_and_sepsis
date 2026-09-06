# =========================================================================
# 51_Fig6D_external cohort 28-day mortality.R
#
# Purpose: Fig6 panel D -- COPD inflammatory signature 28-day mortality risk
#   in independent external cohorts
#   E-MTAB-4451 (Davenport 2016; UK ICU sepsis, n=106)  -- 28-day binary outcome
#   GSE95233   (n=51)                                   -- 28-day binary outcome
#
# Display note: Both cohorts have 28-day binary follow-up (no daily death
#   times; sdrf only survivor/non-survivor); Kaplan-Meier yields only flat
#   curves (non-informative), so instead plot "28-day mortality of High vs
#   Low signature groups" as grouped bars (Fisher's exact test), and keep
#   univariate Cox HR in stats (consistent with meta analysis).
#
# Input (server out_dir, produced by scripts 16 / 24):
#   path16_EMTAB4451_scores.csv    (columns include time,status,copd)
#   path24_GSE95233_valid_scores.csv
# Grouping: copd >= within-cohort median -> High, else Low (consistent with scripts 15/16)
#
# Output:
#   Fig6D_external_28d_mortality.pdf/.png   (two cohorts side by side, patchwork; single plot if missing)
#   path51_fig6d_stats.txt / path51_fig6d_plotdata.rds
#
# Run: Rscript 51_Fig6D_external_cohort_28d_mortality.R
# =========================================================================

.this_file <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
if (.this_file == "" || !file.exists(.this_file)) {
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1])
}
SCRIPT_DIR <- dirname(.this_file)
config_file <- file.path(SCRIPT_DIR, "00_config.R")
if (!file.exists(config_file)) stop("Cannot find 00_config.R: ", config_file)
source(config_file)

suppressPackageStartupMessages({
  library(survival); library(ggplot2); library(dplyr)
})
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)

cat("\n==============================================================\n")
cat("Fig6D: COPD signature external cohort 28-day mortality comparison\n")
cat("==============================================================\n")
stats_lines <- c("===== Fig6D external cohort 28-day mortality =====")

mort <- function(df, cohort_tag) {
  df <- df[!is.na(df$time) & !is.na(df$status) & !is.na(df$copd), ]
  if (nrow(df) < 10) { cat("  [", cohort_tag, "] insufficient samples\n", sep=""); return(NULL) }
  med <- median(df$copd, na.rm = TRUE)
  df$grp <- factor(ifelse(df$copd >= med, "High", "Low"), levels = c("Low", "High"))
  tab <- df %>% group_by(grp) %>%
    summarise(n = n(), deaths = sum(status), .groups = "drop")
  tab$pct <- tab$deaths / tab$n * 100
  m <- matrix(c(tab$deaths[2], tab$n[2] - tab$deaths[2],
                tab$deaths[1], tab$n[1] - tab$deaths[1]), nrow = 2)
  ft <- tryCatch(fisher.test(m), error = function(e) NULL)
  p <- if (!is.null(ft)) ft$p.value else NA_real_
  cox <- tryCatch(summary(coxph(Surv(time, status) ~ copd, data = df))$conf.int[1, ],
                  error = function(e) NULL)
  if (!is.null(ft)) {
    l1 <- sprintf("[%s] 28d mortality High %d/%d (%.1f%%) vs Low %d/%d (%.1f%%); Fisher P=%.3e",
                  cohort_tag, tab$deaths[2], tab$n[2], tab$pct[2],
                  tab$deaths[1], tab$n[1], tab$pct[1], p)
  } else {
    l1 <- sprintf("[%s] 28d mortality High %d/%d (%.1f%%) vs Low %d/%d (%.1f%%)", cohort_tag,
                  tab$deaths[2], tab$n[2], tab$pct[2], tab$deaths[1], tab$n[1], tab$pct[1])
  }
  cat(" ", l1, "\n"); stats_lines <<- c(stats_lines, l1)
  if (!is.null(cox)) {
    l2 <- sprintf("[%s] COPD univariate Cox: HR=%.3f (%.3f-%.3f)", cohort_tag, cox[1], cox[3], cox[4])
    cat(" ", l2, "\n"); stats_lines <<- c(stats_lines, l2)
  }
  list(tab = tab, p = p, tag = cohort_tag, hr = if (!is.null(cox)) cox[1] else NA)
}

plot_mort <- function(r) {
  pv <- if (!is.na(r$p)) sprintf("Fisher P = %.3f", r$p) else "Fisher P = NA"
  ggplot(r$tab, aes(grp, pct)) +
    geom_col(aes(fill = grp), width = 0.6, alpha = 0.9) +
    geom_text(aes(label = sprintf("%d/%d\n(%.0f%%)", deaths, n, pct)),
              vjust = -0.4, size = 3.6) +
    scale_fill_manual(values = c("Low" = "#4C8BF5", "High" = "#D64545"), guide = "none") +
    labs(title = r$tag,
         subtitle = pv,
         x = "COPD signature group", y = "28-day mortality (%)") +
    coord_cartesian(ylim = c(0, max(r$tab$pct) * 1.35 + 5)) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 10, colour = "grey30"))
}

## ---- E-MTAB-4451 ----
cat("\n===== E-MTAB-4451 =====\n")
f1 <- file.path(out_dir, "path16_EMTAB4451_scores.csv")
if (!file.exists(f1)) {
  cat("!! Missing path16_EMTAB4451_scores.csv -- please run 16_EMTAB4451_external_validation.R first\n")
  r1 <- NULL
} else {
  s1 <- read.csv(f1, check.names = FALSE, stringsAsFactors = FALSE)
  cat("E-MTAB scores:", nrow(s1), "rows\n")
  r1 <- mort(s1, "E-MTAB-4451")
}

## ---- GSE95233 ----
cat("\n===== GSE95233 =====\n")
f2 <- file.path(out_dir, "path24_GSE95233_valid_scores.csv")
if (!file.exists(f2)) {
  cat("!! Missing path24_GSE95233_valid_scores.csv -- run 24_GSE95233_prognosis_validation.R first\n")
  r2 <- NULL
} else {
  s2 <- read.csv(f2, check.names = FALSE, stringsAsFactors = FALSE)
  cat("GSE95233 scores:", nrow(s2), "rows\n")
  r2 <- mort(s2, "GSE95233")
}

## ---- Output ----
pl <- list()
if (!is.null(r1)) pl[["E"]] <- plot_mort(r1)
if (!is.null(r2)) pl[["G"]] <- plot_mort(r2)
if (length(pl) == 0) {
  cat("!! Both cohorts unavailable, no output\n")
} else {
  if (length(pl) == 2 && has_patchwork) {
    suppressPackageStartupMessages(library(patchwork))
    combo <- pl[["E"]] + pl[["G"]] + plot_layout(ncol = 2)
    ggsave(file.path(out_dir, "Fig6D_external_28d_mortality.pdf"), combo,
           width = 9, height = 4.4)
    ggsave(file.path(out_dir, "Fig6D_external_28d_mortality.png"), combo,
           width = 9, height = 4.4, dpi = 300, device = grDevices::png)
    cat("Saved: Fig6D_external_28d_mortality.{pdf,png}\n")
  } else {
    ggsave(file.path(out_dir, "Fig6D_external_28d_mortality.pdf"),
           pl[[1]], width = 5, height = 4.4)
    ggsave(file.path(out_dir, "Fig6D_external_28d_mortality.png"),
           pl[[1]], width = 5, height = 4.4, dpi = 300, device = grDevices::png)
    cat("Saved: Fig6D_external_28d_mortality.{pdf,png} (single cohort)\n")
  }
}

writeLines(stats_lines, file.path(out_dir, "path51_fig6d_stats.txt"))
cat("\n----- Statistics summary -----\n"); cat(paste(stats_lines, collapse = "\n"), "\n")
plotdata <- list(
  emtab = if (!is.null(r1)) list(tab = r1$tab, fisher_p = r1$p, hr = r1$hr) else NULL,
  gse95233 = if (!is.null(r2)) list(tab = r2$tab, fisher_p = r2$p, hr = r2$hr) else NULL,
  stats = stats_lines
)
saveRDS(plotdata, file.path(out_dir, "path51_fig6d_plotdata.rds"))
cat("\nSaved path51_fig6d_plotdata.rds\n")
cat("===== Script 51 complete =====\n")
