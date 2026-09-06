# =========================================================================
# 25_batch_sensitivity.R
#
# Purpose: answer the reviewer's core objection — "group and dataset are fully
#          confounded; between-group differences may be batch effects".
#
# Argumentation strategy (three lines of evidence):
#   1) Dataset × group cross-table: show it is NOT fully confounded
#      (multiple datasets contain multiple groups internally).
#   2) Within-dataset contrast (most critical): compare pDC proportions between
#      Healthy and disease groups within the same dataset. If pDC also drops
#      with group "within the same batch", the batch-effect explanation is ruled out.
#   3) Variance partition: compare how much variance in pDC labels group vs
#      dataset each explains (logistic regression pseudo-R² / likelihood ratio),
#      showing group explanatory power ≥ batch.
#
# Input: path1_sepsis_copd_integrated.rds (blood object, with cell_type + group + dataset)
# Output: path25_ series
# =========================================================================

.this_file <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
if (.this_file == "" || !file.exists(.this_file)) {
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1])
}
SCRIPT_DIR  <- dirname(.this_file)
config_file <- file.path(SCRIPT_DIR, "00_config.R")
if (!file.exists(config_file)) stop("cannot find 00_config.R: ", config_file)
source(config_file)

suppressPackageStartupMessages(library(ggplot2))

cat("\n==============================================================\n")
cat("Batch sensitivity analysis (group vs dataset confounding)\n")
cat("==============================================================\n")

dc_markers <- list(
  "cDC1" = c("CLEC9A", "XCR1", "BATF3", "CADM1"),
  "cDC2" = c("CD1C", "CLEC10A", "FCER1A", "ITGAX"),
  "pDC"  = c("LILRA4", "CLEC4C", "IL3RA", "TCF4", "IRF7")
)

## =========================================================================
## Step 1: read the blood object → DC subclustering → pDC labels
## =========================================================================
cat("\n===== Step 1: read blood object + DC subclustering =====\n")

blood_file <- file.path(out_dir, "path1_sepsis_copd_integrated.rds")
if (!file.exists(blood_file)) stop("cannot find ", blood_file, "; run script 01 first")
blood_seu <- readRDS(blood_file)
cat("blood object:", ncol(blood_seu), "cells\n")
cat("meta fields:", paste(colnames(blood_seu@meta.data), collapse = ", "), "\n")

# Check key fields
need <- c("group", "dataset", "cell_type")
miss <- need[!need %in% colnames(blood_seu@meta.data)]
if (length(miss) > 0) stop("missing fields: ", paste(miss, collapse = ", "))

blood_dc <- subset(blood_seu, subset = cell_type == "DC")
rm(blood_seu); gc()
dc_counts <- GetAssayData(blood_dc, assay = "RNA", layer = "counts")
blood_dc[["RNA"]] <- CreateAssay5Object(counts = dc_counts)
rm(dc_counts); gc()

old_res <- resolution
resolution <- 0.6
blood_dc <- integrate_and_annotate(blood_dc, marker_list = dc_markers,
                                   col_name = "dc_subtype")
resolution <- old_res
cat("DC subtype distribution:\n")
print(table(blood_dc$dc_subtype))

## =========================================================================
## Step 2: dataset × group cross-table
## =========================================================================
cat("\n===== Step 2: dataset × group cross-table (pDC cell counts) =====\n")

meta <- blood_dc@meta.data
meta$group <- as.character(meta$group)
meta$dataset <- as.character(meta$dataset)

pdc_cells <- meta[meta$dc_subtype == "pDC", , drop = FALSE]
cross <- table(pdc_cells$dataset, pdc_cells$group)
cat("pDC cell counts (dataset × group):\n")
print(cross)
write.csv(cross, file.path(out_dir, "path25_pdc_dataset_group_cross.csv"))

# Total DC cell cross-table (for computing proportions)
cross_total <- table(meta$dataset, meta$group)
cat("\ntotal DC cell counts (dataset × group):\n")
print(cross_total)

## =========================================================================
## Step 3: within-dataset contrast (most critical) — Healthy vs disease groups
##         within the same dataset
## =========================================================================
cat("\n===== Step 3: within-dataset contrast (Fisher exact test) =====\n")

res_lines <- c("=== Batch sensitivity analysis: within-dataset contrast ===", "",
               "Strategy: for each dataset that internally contains both Healthy and a disease group,",
               "      compare pDC proportions within the same dataset (same batch) to rule out batch effects.", "")

datasets <- sort(unique(meta$dataset))
for (ds in datasets) {
  sub <- meta[meta$dataset == ds, , drop = FALSE]
  grps <- names(table(sub$group))
  if (!"Healthy" %in% grps) {
    cat(sprintf("  [%s] no Healthy control, skipping within-dataset contrast\n", ds))
    next
  }
  for (g in setdiff(grps, "Healthy")) {
    a <- sum(sub$group == "Healthy" & sub$dc_subtype == "pDC")
    b <- sum(sub$group == "Healthy" & sub$dc_subtype != "pDC")
    c <- sum(sub$group == g & sub$dc_subtype == "pDC")
    d <- sum(sub$group == g & sub$dc_subtype != "pDC")
    if (min(c(a, b, c, d)) < 5) {
      cat(sprintf("  [%s] %s vs Healthy: too few cells, skipping\n", ds, g))
      next
    }
    m <- matrix(c(a, b, c, d), nrow = 2,
                dimnames = list(c("pDC", "non-pDC"), c("Healthy", g)))
    ft <- fisher.test(m)
    p_healthy <- a / (a + b)
    p_disease <- c / (c + d)
    line <- sprintf("  [%s] %s vs Healthy: pDC%% %.2f%% vs %.2f%%, OR=%.2f, Fisher p=%.3e",
                    ds, g, p_healthy * 100, p_disease * 100,
                    ft$estimate, ft$p.value)
    cat(line, "\n"); res_lines <- c(res_lines, line)
  }
}

writeLines(res_lines, file.path(out_dir, "path25_within_dataset_contrast.txt"))

## =========================================================================
## Step 4: variance partition — explanatory power of group vs dataset for pDC labels
## =========================================================================
cat("\n===== Step 4: variance partition (group vs dataset) =====\n")

# Compare via logistic regression: univariate group / univariate dataset / bivariate
# McFadden pseudo-R² = 1 - LL_model / LL_null
meta$is_pdc <- as.integer(meta$dc_subtype == "pDC")

# Subsample to avoid overly many cells slowing computation (at most 50k cells)
if (nrow(meta) > 50000) {
  set.seed(123)
  meta_sub <- meta[sample(nrow(meta), 50000), ]
} else {
  meta_sub <- meta
}

ll_null <- tryCatch(logLik(glm(is_pdc ~ 1, data = meta_sub, family = binomial())),
                    error = function(e) NA)
ll_grp   <- tryCatch(logLik(glm(is_pdc ~ group, data = meta_sub, family = binomial())),
                     error = function(e) NA)
ll_ds    <- tryCatch(logLik(glm(is_pdc ~ dataset, data = meta_sub, family = binomial())),
                     error = function(e) NA)
ll_both  <- tryCatch(logLik(glm(is_pdc ~ group + dataset, data = meta_sub, family = binomial())),
                     error = function(e) NA)

pseudo_r2 <- function(ll_m, ll_null) {
  if (is.na(ll_m[1]) || is.na(ll_null[1])) return(NA_real_)
  1 - as.numeric(ll_m) / as.numeric(ll_null)
}
r2_grp  <- pseudo_r2(ll_grp, ll_null)
r2_ds   <- pseudo_r2(ll_ds, ll_null)
r2_both <- pseudo_r2(ll_both, ll_null)

cat(sprintf("McFadden pseudo-R²: group=%.4f | dataset=%.4f | group+dataset=%.4f\n",
            r2_grp, r2_ds, r2_both))
cat(sprintf("→ group/batch explanatory power ratio = %.2f (larger favors a group effect)\n",
            r2_grp / r2_ds))

res_lines2 <- c("=== Variance partition (logistic McFadden pseudo-R²) ===",
                sprintf("group univariate pseudo-R² = %.4f", r2_grp),
                sprintf("dataset univariate pseudo-R² = %.4f", r2_ds),
                sprintf("group+dataset pseudo-R² = %.4f", r2_both),
                sprintf("group/batch explanatory power ratio = %.2f", r2_grp / r2_ds),
                "Interpretation: if group's pseudo-R² clearly exceeds dataset's, the pDC differences are mainly group-driven,",
                "      not batch; if the two are close or dataset is larger, be cautious — batch confounding cannot be ruled out.")
writeLines(res_lines2, file.path(out_dir, "path25_variance_partition.txt"))
cat("saved: path25_variance_partition.txt\n")

## =========================================================================
## Step 5: visualization — pDC proportion of Healthy vs disease groups per dataset
## =========================================================================
cat("\n===== Step 5: visualization =====\n")

# pDC proportion per (dataset, group)
summ <- aggregate(is_pdc ~ dataset + group, data = meta,
                  FUN = function(x) c(pdc_prop = mean(x), n = length(x)))
summ <- do.call(data.frame, summ)
colnames(summ) <- c("dataset", "group", "pdc_prop", "n")

p <- ggplot(summ, aes(x = group, y = pdc_prop, fill = group)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = paste0("n=", n)), vjust = -0.3, size = 2.5) +
  facet_wrap(~ dataset, scales = "free_x") +
  labs(x = NULL, y = "pDC proportion (of DC)", title = "pDC proportion by group within each dataset") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave(file.path(out_dir, "path25_pdc_prop_by_dataset_group.pdf"), p,
       width = 10, height = 6)
ggsave(file.path(out_dir, "path25_pdc_prop_by_dataset_group.png"), p,
       width = 10, height = 6, dpi = 300)

cat("\nBatch sensitivity analysis done\n")
cat("Output: path25_pdc_dataset_group_cross.csv / path25_within_dataset_contrast.txt\n")
cat("      path25_variance_partition.txt / path25_pdc_prop_by_dataset_group.pdf/.png\n")
cat("Interpretation: if the within-dataset contrast (Step 3) also shows pDC declining with group, and group\n")
cat("      explanatory power (Step 4) exceeds batch, then the \"healthy < COPD < sepsis\" gradient is not a batch\n")
cat("      effect and can be written into the paper's response to the objection.\n")
