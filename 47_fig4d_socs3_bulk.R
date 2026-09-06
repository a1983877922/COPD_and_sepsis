# =========================================================================
# 47 Fig4D SOCS3 bulk validation
#
# Purpose: Fig4 panel D — validate SOCS3 single-gene up-regulation in independent bulk cohorts
#   1) GSE57148 (lung tissue RNA-seq, FPKM): COPD (n≈98) vs smoking controls (n≈91)
#   2) GSE66099 (whole blood microarray GPL570): Sepsis (n≈199) vs Control (n≈47)
#
# Note: SOCS3 is a severity-grading regulator discovered on the blood side by single-cell
#   analysis. Here we test whether its single-gene expression is concordantly up-regulated in
#   independent bulk cohorts (Wilcoxon; the direction of the one-sided claim follows the data).
#
# Outputs (Fig4D_* prefix, PDF+PNG):
#   Fig4D_SOCS3_bulk_lung.pdf/.png     GSE57148 lung (COPD vs Control)
#   Fig4D_SOCS3_bulk_blood.pdf/.png    GSE66099 blood (Sepsis vs Control)
#   Fig4D_SOCS3_bulk.pdf/.png          combined plot (needs patchwork; otherwise only the
#                                      individual panels are produced)
#   Fig4D_SOCS3_bulk_stats.txt         statistical summary
#   Fig4D_SOCS3_bulk_long.csv          pre-plot long table (for local re-rendering)
#   path47_fig4d_plotdata.rds          pre-plot data (including both cohort statistics)
#
# Server: Rscript 47_fig4d_socs3_bulk.R
# Dependency: GEOquery (GSE66099 series matrix parsing)
# =========================================================================

.this_file <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
if (.this_file == "" || !file.exists(.this_file)) {
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1])
}
SCRIPT_DIR <- dirname(.this_file)
config_file <- file.path(SCRIPT_DIR, "00_config.R")
if (!file.exists(config_file)) stop("cannot find 00_config.R: ", config_file)
source(config_file)

suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)

cat("\n==============================================================\n")
cat("Fig4D: validation of SOCS3 as a single gene in independent bulk cohorts\n")
cat("==============================================================\n")

TARGET <- "SOCS3"
CTRL_COL <- "#3C5488"   # Control / Healthy blue
COPD_COL <- "#00A087"   # COPD green (consistent with the five-group palette of Fig2)
SEPS_COL <- "#E64B35"   # Sepsis red

## ---------------------------------------------------------------------------
## Shared: single-gene boxplot
## ---------------------------------------------------------------------------
socs3_box <- function(df, xvar, ycol, grp_ord, fills, ylab, subtitle, pval, title) {
  df <- df[!is.na(df[[ycol]]) & is.finite(df[[ycol]]), ]
  df[[xvar]] <- factor(df[[xvar]], levels = grp_ord)
  p <- ggplot(df, aes(x = .data[[xvar]], y = .data[[ycol]])) +
    geom_boxplot(aes(fill = .data[[xvar]]), outlier.shape = NA, width = 0.55,
                 alpha = 0.75, linewidth = 0.35) +
    geom_jitter(width = 0.12, size = 0.7, alpha = 0.35, colour = "grey25") +
    scale_fill_manual(values = fills, guide = "none") +
    labs(title = title, subtitle = subtitle, x = NULL, y = ylab) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 8.5, colour = "grey30"))
  if (!is.null(pval)) {
    ymax <- max(df[[ycol]], na.rm = TRUE)
    p <- p + annotate("text", x = 1.5, y = ymax * 1.12, hjust = 0.5,
                      label = pval, size = 3.6, fontface = "bold")
  }
  p
}

fmt_p <- function(p) {
  if (p < 2.2e-16) "P < 2.2e-16" else if (p < 0.001) sprintf("P = %.2e", p)
  else sprintf("P = %.4f", p)
}

## ===========================================================================
## Cohort 1: GSE57148 lung tissue RNA-seq (FPKM) — COPD vs smoking controls
## ===========================================================================
cat("\n########## GSE57148 (lung, COPD vs Control) ##########\n")

fpkm_file <- file.path(out_dir, "GSE57148_COPD_FPKM_Normalized.txt.gz")
matrix_file <- file.path(out_dir, "GSE57148_series_matrix.txt.gz")
if (!file.exists(fpkm_file)) stop("no local GSE57148 FPKM: ", fpkm_file,
                                  "\nplease run the GSE57148 lung-side validation script first, or download the file manually")
if (!file.exists(matrix_file)) stop("no local GSE57148 series matrix: ", matrix_file)

fpkm <- read.delim(gzfile(fpkm_file), row.names = 1, check.names = FALSE)
expr_lung <- log2(fpkm + 1)
cat("FPKM matrix:", nrow(expr_lung), "genes x", ncol(expr_lung), "samples\n")

gse_l <- getGEO(filename = matrix_file, getGPL = FALSE)
pd_l <- pData(gse_l)
title_raw <- as.character(pd_l$title)
sample_alias <- vapply(strsplit(title_raw, " "), function(x) x[length(x)], character(1))
sample_alias <- gsub("^\"|\"$", "", sample_alias)
grp_raw <- as.character(pd_l$title)
is_copd <- grepl("copd", grp_raw, ignore.case = TRUE)

common_l <- intersect(colnames(expr_lung), sample_alias)
cat("Samples aligned:", length(common_l), "\n")
if (length(common_l) == 0) stop("GSE57148 sample names cannot be aligned with the title aliases")
m_l <- match(common_l, sample_alias)
expr_lung <- expr_lung[, common_l, drop = FALSE]
grp_lung <- ifelse(is_copd[m_l], "COPD", "Control")
cat("Grouping: COPD =", sum(grp_lung == "COPD"), " / Control =", sum(grp_lung == "Control"), "\n")

if (!TARGET %in% rownames(expr_lung)) {
  cat("!! SOCS3 not found in GSE57148 (not included in the FPKM gene set?), skipping the lung cohort\n")
  lung_res <- NULL
} else {
  v <- as.numeric(expr_lung[TARGET, ])
  w_l <- wilcox.test(v[grp_lung == "COPD"], v[grp_lung == "Control"])
  med_copd <- median(v[grp_lung == "COPD"]); med_ctrl <- median(v[grp_lung == "Control"])
  cat(sprintf("SOCS3: COPD med=%.3f vs Control med=%.3f, Wilcoxon P = %.3e\n",
              med_copd, med_ctrl, w_l$p.value))
  df_l <- data.frame(sample = colnames(expr_lung), cohort = "GSE57148_lung",
                     group = grp_lung, socs3 = v, stringsAsFactors = FALSE)
  p_l <- socs3_box(df_l, "group", "socs3",
                   grp_ord = c("Control", "COPD"),
                   fills = c("Control" = CTRL_COL, "COPD" = COPD_COL),
                   ylab = paste0(TARGET, " expression\n(log2 FPKM + 1)"),
                   subtitle = sprintf("GSE57148 lung, COPD vs smoking controls\nn = %d vs %d",
                                      sum(grp_lung == "COPD"), sum(grp_lung == "Control")),
                   pval = fmt_p(w_l$p.value),
                   title = "SOCS3 in COPD lung tissue (bulk)")
  lung_res <- list(p = p_l, df = df_l, pval = w_l$p.value,
                   med_copd = med_copd, med_ctrl = med_ctrl)
  ggsave(file.path(out_dir, "Fig4D_SOCS3_bulk_lung.pdf"), p_l, width = 4.2, height = 4.4)
  ggsave(file.path(out_dir, "Fig4D_SOCS3_bulk_lung.png"), p_l, width = 4.2, height = 4.4,
         dpi = 300, device = grDevices::png)
  cat("Saved: Fig4D_SOCS3_bulk_lung.{pdf,png}\n")
}

## ===========================================================================
## Cohort 2: GSE66099 whole blood microarray (GPL570) — Sepsis vs Control
## ===========================================================================
cat("\n########## GSE66099 (whole blood, Sepsis vs Control) ##########\n")

matrix_file2 <- file.path(out_dir, "GSE66099_matrix.txt.gz")
if (!file.exists(matrix_file2)) {
  cat("no local series matrix, downloading...\n")
  gse_b <- getGEO("GSE66099", GSEMatrix = TRUE, getGPL = TRUE)[[1]]
} else {
  gse_b <- getGEO(filename = matrix_file2, getGPL = TRUE)
}
expr_b <- exprs(gse_b)
cat("Probe matrix:", nrow(expr_b), "x", ncol(expr_b), "\n")

pd_b <- pData(gse_b)
fd_b <- fData(gse_b)
sym_col <- NULL
for (cand in c("Gene Symbol", "GeneSymbol", "Symbol", "Gene.Symbol", "gene_symbol")) {
  if (cand %in% colnames(fd_b)) { sym_col <- cand; break }
}
if (is.null(sym_col)) {
  m <- grep("symbol", colnames(fd_b), ignore.case = TRUE)
  if (length(m) > 0) sym_col <- colnames(fd_b)[m[1]]
}
if (is.null(sym_col)) stop("no gene symbol column found in GSE66099 fData")
cat("symbol column:", sym_col, "\n")

sym <- as.character(fd_b[[sym_col]]); names(sym) <- rownames(fd_b)
keep <- sym != "" & !is.na(sym) & sym != "---"
expr_b <- expr_b[keep, , drop = FALSE]; sym <- sym[keep]
s3_probes <- names(sym)[sym == TARGET]
cat("SOCS3 probe count:", length(s3_probes), "\n")
if (length(s3_probes) == 0) {
  cat("!! no SOCS3 probe in GSE66099 (GPL570)? skipping the blood cohort after checking\n")
  blood_res <- NULL
} else {
  if (length(s3_probes) > 1) {
    v_b <- colMeans(expr_b[s3_probes, , drop = FALSE])
  } else {
    v_b <- as.numeric(expr_b[s3_probes, ])
  }
  if (!"disease:ch1" %in% colnames(pd_b)) stop("GSE66099 has no disease:ch1 field")
  grp_b <- as.character(pd_b[["disease:ch1"]]); names(grp_b) <- rownames(pd_b)
  grp_b[grp_b == "SepticShock"] <- "Sepsis"
  common_b <- intersect(names(v_b), names(grp_b))
  v_b <- v_b[common_b]; grp_b <- grp_b[common_b]
  keep2 <- grp_b %in% c("Sepsis", "Control")
  v_b <- v_b[keep2]; grp_b <- grp_b[keep2]
  cat("after alignment + filtering:\n"); print(table(grp_b))

  if (sum(grp_b == "Sepsis") < 3 || sum(grp_b == "Control") < 3) {
    cat("!! too few samples per group, skipping the blood cohort\n"); blood_res <- NULL
  } else {
    w_b <- wilcox.test(v_b[grp_b == "Sepsis"], v_b[grp_b == "Control"])
    med_sep <- median(v_b[grp_b == "Sepsis"]); med_ctrl <- median(v_b[grp_b == "Control"])
    cat(sprintf("SOCS3: Sepsis med=%.3f vs Control med=%.3f, Wilcoxon P = %.3e\n",
                med_sep, med_ctrl, w_b$p.value))
    df_b <- data.frame(sample = names(v_b), cohort = "GSE66099_blood",
                       group = grp_b, socs3 = v_b, stringsAsFactors = FALSE)
    p_b <- socs3_box(df_b, "group", "socs3",
                     grp_ord = c("Control", "Sepsis"),
                     fills = c("Control" = CTRL_COL, "Sepsis" = SEPS_COL),
                     ylab = paste0(TARGET, " expression\n(log2 microarray)"),
                     subtitle = sprintf("GSE66099 whole blood, sepsis vs control\nn = %d vs %d",
                                        sum(grp_b == "Sepsis"), sum(grp_b == "Control")),
                     pval = fmt_p(w_b$p.value),
                     title = "SOCS3 in septic whole blood (bulk)")
    blood_res <- list(p = p_b, df = df_b, pval = w_b$p.value,
                      med_sep = med_sep, med_ctrl = med_ctrl)
    ggsave(file.path(out_dir, "Fig4D_SOCS3_bulk_blood.pdf"), p_b, width = 4.2, height = 4.4)
    ggsave(file.path(out_dir, "Fig4D_SOCS3_bulk_blood.png"), p_b, width = 4.2, height = 4.4,
           dpi = 300, device = grDevices::png)
    cat("Saved: Fig4D_SOCS3_bulk_blood.{pdf,png}\n")
  }
}

## ===========================================================================
## Combined plot (patchwork optional)
## ===========================================================================
if (has_patchwork && !is.null(lung_res) && !is.null(blood_res)) {
  suppressPackageStartupMessages(library(patchwork))
  combo <- lung_res$p + blood_res$p + plot_layout(ncol = 2)
  ggsave(file.path(out_dir, "Fig4D_SOCS3_bulk.pdf"), combo, width = 9, height = 4.6)
  ggsave(file.path(out_dir, "Fig4D_SOCS3_bulk.png"), combo, width = 9, height = 4.6,
         dpi = 300, device = grDevices::png)
  cat("Saved: Fig4D_SOCS3_bulk.{pdf,png} (combined plot)\n")
} else {
  cat("combined plot skipped (patchwork missing or one cohort has no result)\n")
}

## ===========================================================================
## Statistical summary + pre-plot data
## ===========================================================================
cat("\n===== statistical summary =====\n")
stats_lines <- c("===== Fig4D: SOCS3 independent bulk validation (Wilcoxon) =====")
if (!is.null(lung_res)) {
  stats_lines <- c(stats_lines,
    sprintf("GSE57148 lung   : COPD med=%.3f (n=%d) vs Control med=%.3f (n=%d), P = %.3e",
            lung_res$med_copd, sum(lung_res$df$group == "COPD"),
            lung_res$med_ctrl, sum(lung_res$df$group == "Control"), lung_res$pval))
}
if (!is.null(blood_res)) {
  stats_lines <- c(stats_lines,
    sprintf("GSE66099 blood  : Sepsis med=%.3f (n=%d) vs Control med=%.3f (n=%d), P = %.3e",
            blood_res$med_sep, sum(blood_res$df$group == "Sepsis"),
            blood_res$med_ctrl, sum(blood_res$df$group == "Control"), blood_res$pval))
}
if (is.null(lung_res) && is.null(blood_res)) {
  stats_lines <- c(stats_lines, "no SOCS3 signal in either cohort; check the gene name / probe mapping")
}
writeLines(stats_lines, file.path(out_dir, "Fig4D_SOCS3_bulk_stats.txt"))
cat(paste(stats_lines, collapse = "\n"), "\n")

all_df <- rbind(if (!is.null(lung_res)) lung_res$df else NULL,
                if (!is.null(blood_res)) blood_res$df else NULL)
if (!is.null(all_df)) {
  write.csv(all_df, file.path(out_dir, "Fig4D_SOCS3_bulk_long.csv"), row.names = FALSE)
  plotdata <- list(
    target = TARGET,
    lung  = lung_res[c("df", "pval", "med_copd", "med_ctrl")],
    blood = blood_res[c("df", "pval", "med_sep", "med_ctrl")],
    stats = stats_lines,
    ggplots = list(lung = if (!is.null(lung_res)) lung_res$p else NULL,
                   blood = if (!is.null(blood_res)) blood_res$p else NULL),
    ctrl_col = CTRL_COL, copd_col = COPD_COL, seps_col = SEPS_COL
  )
  saveRDS(plotdata, file.path(out_dir, "path47_fig4d_plotdata.rds"))
  sz <- round(file.info(file.path(out_dir, "path47_fig4d_plotdata.rds"))$size / 1e6, 2)
  cat("Saved path47_fig4d_plotdata.rds (", sz, " MB ) — local re-render: readRDS then ggsave(plotdata$ggplots$lung, ...)\n")
}

cat("\n===== script 47 done =====\n")
