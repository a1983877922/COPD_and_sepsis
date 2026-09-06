# =========================================================================
# 17_socs3_virtual_knockout.R (SOCS3 virtual knockout)
# SOCS3 is a negative-feedback inhibitor of the IFN/JAK-STAT pathway. This
# script validates on pDC data the relationship between SOCS3 and downstream
# IFN signaling (an association-level 'virtual knockout'):
#   1) SOCS3 high/low expression groups -> compare ISG/IFN signaling scores
#   2) Correlation between SOCS3 expression and ISG score (Spearman)
#   3) Compare STAT1/STAT3 downstream target gene scores
#
# Expectation (two possibilities; let the data decide):
#   - Negative correlation: SOCS3 high -> ISG low = SOCS3 is suppressing IFN
#     signaling (the 'brake' of immunoparalysis is engaged)
#   - Positive correlation: SOCS3 and ISG co-upregulated = SOCS3 is the
#     accompanying negative feedback after IFN activation (snapshot not yet effective)
# Either way, this supports 'SOCS3 is a key regulator of pDC IFN signaling'.
#
# Note: this is an 'associational' validation. For a true GRN-level in silico
# KO (CellOracle perturbation), CellOracle + base GRN must be installed
# separately, and a separate script written.
#
# Dependencies: Seurat (already present); reads path11_pdc_subset.rds (output of script 11)
# =========================================================================

.this_file <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
if (.this_file == "" || !file.exists(.this_file)) {
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1])
}
SCRIPT_DIR  <- dirname(.this_file)
config_file <- file.path(SCRIPT_DIR, "00_config.R")
if (!file.exists(config_file)) stop("00_config.R not found: ", config_file)
source(config_file)

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))

cat("\n==============================================================\n")
cat("SOCS3 downstream signaling and virtual knockout validation (pDC)\n")
cat("==============================================================\n")

pdc_file <- file.path(out_dir, "path11_pdc_subset.rds")
if (!file.exists(pdc_file)) stop("path11_pdc_subset.rds not found, please run 11 first")
pdc <- readRDS(pdc_file)
cat("pDC cell count:", ncol(pdc), "\n")

## ISG / IFN signaling geneset (representative readout of type I interferon)
isg_genes <- c("ISG15","MX1","MX2","OAS1","OAS2","OAS3","IFIT1","IFIT2","IFIT3",
               "IFI44","IFI44L","IFI6","ISG20","IRF7","IRF9","STAT1","STAT2",
               "GBP1","BST2","IFIH1","USP18","RSAD2")
## STAT1/STAT3 downstream target genes (simplified, includes SOCS3 itself)
stat_targets <- c("IRF1","IRF7","IRF9","MX1","OAS1","ISG15","STAT1","BCL3","SOCS3","IL6ST")

# Hit check
isg_hit <- intersect(isg_genes, rownames(pdc))
stat_hit <- intersect(stat_targets, rownames(pdc))
cat("ISG genes hit:", length(isg_hit), "/", length(isg_genes), "\n")
cat("STAT target genes hit:", length(stat_hit), "/", length(stat_targets), "\n")

## =========================================================================
## Step 1: SOCS3 expression + ISG/STAT module scoring
## =========================================================================
cat("\n===== Step 1: Compute SOCS3 expression + signaling scores =====\n")
if (!("SOCS3" %in% rownames(pdc))) stop("No SOCS3 gene in pDC object")

expr_soc3 <- GetAssayData(pdc, assay = "RNA", layer = "data")["SOCS3", ]

# ISG module score (z-score average)
score_module <- function(obj, genes) {
  genes <- intersect(genes, rownames(obj))
  if (length(genes) < 2) return(rep(NA_real_, ncol(obj)))
  mat <- GetAssayData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE]
  z <- t(scale(t(as.matrix(mat))))
  z[is.na(z)] <- 0
  colMeans(z)
}
isg_score <- score_module(pdc, isg_genes)
stat_score <- score_module(pdc, stat_targets)

meta <- data.frame(sample = colnames(pdc), group = as.character(pdc$group),
                   SOCS3 = as.numeric(expr_soc3),
                   ISG_score = isg_score, STAT_score = stat_score,
                   stringsAsFactors = FALSE)

## =========================================================================
## Step 2: SOCS3 high/low grouping (top/bottom 25%) -> ISG/STAT score comparison
## =========================================================================
cat("\n===== Step 2: SOCS3 high/low group comparison =====\n")
q <- quantile(meta$SOCS3, probs = c(0.25, 0.75), na.rm = TRUE)
meta$SOCS3_grp <- ifelse(meta$SOCS3 >= q[2], "SOCS3_high",
                  ifelse(meta$SOCS3 <= q[1], "SOCS3_low", "mid"))
cat("SOCS3 grouping thresholds (Q25/Q75):", round(q, 3), "\n")
print(table(meta$SOCS3_grp))

res_lines <- character()
for (score_nm in c("ISG_score", "STAT_score")) {
  v_hi <- meta[[score_nm]][meta$SOCS3_grp == "SOCS3_high"]
  v_lo <- meta[[score_nm]][meta$SOCS3_grp == "SOCS3_low"]
  if (length(v_hi) >= 5 && length(v_lo) >= 5) {
    wt <- wilcox.test(v_hi, v_lo)
    line <- sprintf("%s: SOCS3_high mean=%.3f vs SOCS3_low mean=%.3f, wilcox p=%.3e",
                    score_nm, mean(v_hi, na.rm=TRUE), mean(v_lo, na.rm=TRUE), wt$p.value)
    cat(" ", line, "\n"); res_lines <- c(res_lines, line)
  }
}

## =========================================================================
## Step 3: SOCS3 vs ISG score correlation (Spearman)
## =========================================================================
cat("\n===== Step 3: Correlation (Spearman) =====\n")
for (score_nm in c("ISG_score", "STAT_score")) {
  ct <- cor.test(meta$SOCS3, meta[[score_nm]], method = "spearman")
  line <- sprintf("SOCS3 vs %s: rho=%.3f, p=%.3e", score_nm, ct$estimate, ct$p.value)
  cat(" ", line, "\n"); res_lines <- c(res_lines, line)
}

write.csv(meta, file.path(out_dir, "path17_SOCS3_signal_scores.csv"), row.names = FALSE)
writeLines(res_lines, file.path(out_dir, "path17_SOCS3_signal_stats.txt"))

## =========================================================================
## Step 4: Visualization
## =========================================================================
cat("\n===== Step 4: Visualization =====\n")

# 4.1 ISG score boxplot for SOCS3 high/low groups
p1 <- ggplot(meta[meta$SOCS3_grp != "mid", ], aes(x = SOCS3_grp, y = ISG_score, fill = SOCS3_grp)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.2, alpha = 0.3) +
  scale_fill_manual(values = c("SOCS3_high" = "#d64545", "SOCS3_low" = "#4c8bf5")) +
  labs(x = NULL, y = "ISG signature score", title = "SOCS3 high vs low pDC: ISG score") +
  theme_bw(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(out_dir, "path17_SOCS3_ISG_boxplot.pdf"), p1, width = 5.5, height = 5)

# 4.2 SOCS3 vs ISG scatter plot
p2 <- ggplot(meta, aes(x = SOCS3, y = ISG_score)) +
  geom_point(size = 0.3, alpha = 0.3, color = "#5b6472") +
  geom_smooth(method = "lm", se = TRUE, color = "#d64545") +
  labs(x = "SOCS3 expression", y = "ISG signature score",
       title = "SOCS3 vs ISG score (pDC)") +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "path17_SOCS3_ISG_scatter.pdf"), p2, width = 5.5, height = 5)

# 4.3 SOCS3 by group
p3 <- ggplot(meta, aes(x = group, y = SOCS3, fill = group)) +
  geom_violin(alpha = 0.8, scale = "width") +
  labs(x = NULL, y = "SOCS3 expression", title = "SOCS3 by group (pDC)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out_dir, "path17_SOCS3_by_group.pdf"), p3, width = 6, height = 5)

cat("\nSOCS3 downstream signaling analysis complete\n")
cat("Outputs: path17_SOCS3_signal_scores.csv / path17_SOCS3_signal_stats.txt /\n")
cat("      path17_SOCS3_ISG_boxplot.pdf / path17_SOCS3_ISG_scatter.pdf / path17_SOCS3_by_group.pdf\n")
cat("Interpretation: if the SOCS3_high group's ISG_score is significantly lower than SOCS3_low (negative correlation), then SOCS3 is suppressing IFN signaling\n")
cat("      (immunoparalysis brake engaged); if positive, it is the accompanying negative feedback after IFN activation.\n")
cat("Note: a true GRN-level virtual knockout requires CellOracle (in silico perturbation), needing separate package install + base GRN.\n")
