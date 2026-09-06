# =========================================================================
# 26_pdc_annotation_audit.R
#
# Purpose: Diagnose pDC sub-clustering annotation issues. Sensitivity analysis
#   (script 25) found:
#   - SCP548 pDC accounts for 49-59% of DC (abnormal; normal should be <20%)
#   - Contradictory control direction within datasets (in GSE167363/SCP548,
#     Sepsis pDC is paradoxically higher)
#   - Variance decomposition shows dataset explains variance >> group
#   These all point to inconsistent pDC annotation across datasets, needs audit.
#
# Diagnostic approach: Use AddModuleScore scoring of pDC/cDC1/cDC2 marker genes
#   to verify whether dc_subtype annotation is consistent with marker scores;
#   also check whether the DC major class (cell_type=="DC") itself is reasonable
#   across datasets.
#
# Input: path1_sepsis_copd_integrated.rds
# Output: path26_ series
# =========================================================================

.this_file <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
if (.this_file == "" || !file.exists(.this_file)) {
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1])
}
SCRIPT_DIR  <- dirname(.this_file)
config_file <- file.path(SCRIPT_DIR, "00_config.R")
if (!file.exists(config_file)) stop("Cannot find 00_config.R: ", config_file)
source(config_file)

suppressPackageStartupMessages(library(ggplot2))

cat("\n==============================================================\n")
cat("pDC annotation audit (marker score validation)\n")
cat("==============================================================\n")

dc_markers <- list(
  "cDC1" = c("CLEC9A", "XCR1", "BATF3", "CADM1"),
  "cDC2" = c("CD1C", "CLEC10A", "FCER1A", "ITGAX"),
  "pDC"  = c("LILRA4", "CLEC4C", "IL3RA", "TCF4", "IRF7")
)

## =========================================================================
## Step 1: Load blood object + print cell_type distribution (is DC major class reasonable)
## =========================================================================
cat("\n===== Step 1: cell_type distribution (DC proportion per dataset) =====\n")

blood_file <- file.path(out_dir, "path1_sepsis_copd_integrated.rds")
if (!file.exists(blood_file)) stop("Cannot find ", blood_file, " please run script 01 first")
blood_seu <- readRDS(blood_file)
meta_all <- blood_seu@meta.data

# cell_type distribution per dataset
cat("cell_type distribution per dataset (top 8 types):\n")
for (ds in sort(unique(meta_all$dataset))) {
  sub <- meta_all[meta_all$dataset == ds, ]
  tt <- sort(table(sub$cell_type), decreasing = TRUE)
  cat(sprintf("  [%s] total cells %d | DC proportion %.2f%% | ", ds, nrow(sub),
              sum(sub$cell_type == "DC") / nrow(sub) * 100))
  cat(paste(names(tt)[1:min(6, length(tt))], collapse = ", "), "\n")
}

## =========================================================================
## Step 2: subset DC + sub-clustering (reuse script 19)
## =========================================================================
cat("\n===== Step 2: DC sub-clustering =====\n")

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

cat("DC subtype distribution (per dataset):\n")
tt <- table(blood_dc$dataset, blood_dc$dc_subtype)
print(tt)
write.csv(tt, file.path(out_dir, "path26_dc_subtype_by_dataset.csv"))

## =========================================================================
## Step 3: AddModuleScore scoring validation
## =========================================================================
cat("\n===== Step 3: marker score validation =====\n")

# Ensure data layer exists (integrate_and_annotate should have done it; do once more just in case)
if (!"data" %in% Layers(blood_dc)) {
  blood_dc <- NormalizeData(blood_dc, verbose = FALSE)
}

for (nm in names(dc_markers)) {
  genes <- intersect(dc_markers[[nm]], rownames(blood_dc))
  if (length(genes) < 2) next
  blood_dc <- AddModuleScore(blood_dc, features = list(genes),
                             name = paste0("score_", nm))
  colnames(blood_dc@meta.data)[colnames(blood_dc@meta.data) == paste0("score_", nm, "1")] <-
    paste0("score_", nm)
}

## =========================================================================
## Step 4: Diagnostic output —— consistency of dc_subtype vs marker score
## =========================================================================
cat("\n===== Step 4: dc_subtype vs marker score =====\n")

m <- blood_dc@meta.data
m <- m[, !duplicated(colnames(m)), drop = FALSE]  # deduplicate columns (AddModuleScore rename may leave duplicates)
# For each dc_subtype, look at the median of the corresponding marker score
res_lines <- c("=== pDC annotation audit: dc_subtype vs marker score ===", "",
               "If cells with dc_subtype=='pDC' have pDC score median clearly higher than cDC score,",
               "then annotation is reasonable; if SCP548's pDC score is abnormally low or cDC score abnormally high,",
               "then sub-clustering mislabeled cDC as pDC.", "")

for (ds in sort(unique(m$dataset))) {
  sub <- m[m$dataset == ds, ]
  line <- sprintf("[%s] marker score median per dc_subtype:", ds)
  cat(line, "\n"); res_lines <- c(res_lines, line)
  for (st in c("pDC", "cDC1", "cDC2")) {
    ss <- sub[sub$dc_subtype == st, , drop = FALSE]
    if (nrow(ss) < 10) next
    pdc_m <- median(ss$score_pDC, na.rm = TRUE)
    cdc1_m <- median(ss$score_cDC1, na.rm = TRUE)
    cdc2_m <- median(ss$score_cDC2, na.rm = TRUE)
    line2 <- sprintf("    %-5s (n=%d): pDC score=%.3f | cDC1 score=%.3f | cDC2 score=%.3f",
                     st, nrow(ss), pdc_m, cdc1_m, cdc2_m)
    cat(line2, "\n"); res_lines <- c(res_lines, line2)
  }
}

writeLines(res_lines, file.path(out_dir, "path26_pdc_annotation_check.txt"))

## =========================================================================
## Step 4b: Within-dataset marker score control (truly clean validation, removes batch)
##   Use pDC marker score (score_pDC) instead of sub-cluster labels, compare
##   Healthy vs disease group scores within a single dataset —
##   objective (does not depend on annotation) + single dataset (removes batch).
## =========================================================================
cat("\n===== Step 4b: within-dataset marker score control (objective + batch-removed) =====\n")

res2_lines <- c("=== within-dataset marker score control (objective + batch-removed) ===", "",
                "Use pDC marker score (score_pDC, independent of sub-cluster annotation),",
                "compare Healthy vs disease group scores within a single dataset.",
                "If within SCP548 Healthy score is also significantly higher than Sepsis → supports pDC signal depletion;",
                "If opposite → pDC depletion does not hold.", "")

for (ds in sort(unique(m$dataset))) {
  sub <- m[m$dataset == ds, , drop = FALSE]
  grps <- names(table(sub$group))
  if (!"Healthy" %in% grps) {
    cat(sprintf("  [%s] no Healthy control, skip\n", ds))
    next
  }
  for (g in setdiff(grps, "Healthy")) {
    h <- sub$score_pDC[sub$group == "Healthy"]
    d <- sub$score_pDC[sub$group == g]
    if (length(h) < 20 || length(d) < 20) {
      cat(sprintf("  [%s] %s vs Healthy: too few cells, skip\n", ds, g))
      next
    }
    wt <- wilcox.test(d, h)
    line <- sprintf("  [%s] %s vs Healthy: pDC score median %.3f vs %.3f, Wilcoxon p=%.3e",
                    ds, g, median(d), median(h), wt$p.value)
    cat(line, "\n"); res2_lines <- c(res2_lines, line)
  }
}

writeLines(res2_lines, file.path(out_dir, "path26_marker_within_dataset_contrast.txt"))
cat("Saved: path26_marker_within_dataset_contrast.txt\n")

## =========================================================================
## Step 5: Visualization —— pDC score distribution (by dataset × dc_subtype)
## =========================================================================
cat("\n===== Step 5: visualization =====\n")

# pDC marker score violin plot: compare pDC vs cDC score per dataset
m$dc_subtype <- factor(m$dc_subtype, levels = c("pDC", "cDC1", "cDC2"))
p <- ggplot(m, aes(x = dc_subtype, y = score_pDC, fill = dc_subtype)) +
  geom_violin(scale = "width", alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
  facet_wrap(~ dataset, scales = "free_y") +
  scale_fill_manual(values = c("pDC" = "#E53935", "cDC1" = "#43A047", "cDC2" = "#1E88E5")) +
  labs(x = NULL, y = "pDC marker score", title = "pDC marker score × dataset × dc_subtype") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "path26_pdc_score_by_dataset.pdf"), p, width = 12, height = 6)
ggsave(file.path(out_dir, "path26_pdc_score_by_dataset.png"), p, width = 12, height = 6, dpi = 300)

# SCP548 author annotation vs unified annotation cross-tab (if scp548_cell_type exists)
if ("scp548_cell_type" %in% colnames(m)) {
  scp <- m[m$dataset == "SCP548", ]
  cat("\nSCP548 author annotation vs unified cell_type cross-tab:\n")
  print(table(scp$scp548_cell_type, scp$cell_type))
  cat("SCP548 author annotation vs dc_subtype cross-tab:\n")
  print(table(scp$scp548_cell_type, scp$dc_subtype))
}

cat("\npDC annotation audit complete\n")
cat("Output: path26_dc_subtype_by_dataset.csv / path26_pdc_annotation_check.txt\n")
cat("      path26_pdc_score_by_dataset.pdf/.png\n")
