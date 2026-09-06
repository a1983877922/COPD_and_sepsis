# =========================================================================
# 07_dc_subclustering.R
# Verify, in the "blood-lung DC redistribution", which specific DC subtype
# is migrating/exhausted (cDC1 vs cDC2 vs pDC).
#
# Blood side: read path1 integrated object, subset DC -> subclustering -> annotate cDC1/cDC2/pDC
#       -> observe subtype proportion changes across five groups (Healthy/Infection/COPD/Sepsis/Pneumonia-source)
# Lung side: read path2 lung object (with built-in Manuscript_Identity annotation, including cDC1/cDC2/pDC
#       /DC_Mature/DC_Langerhans) -> compare Control vs COPD proportions
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

cat("\n==============================================================\n")
cat("DC subset deep analysis (blood-side subclustering + lung-side subtype proportions)\n")
cat("==============================================================\n")

## DC subtype marker genes (all verified present in the object; CLEC9A/FCER1A may be missing
## due to min.cells, and annotate_by_markers will automatically skip insufficient ones)
dc_markers <- list(
  "cDC1" = c("CLEC9A", "XCR1", "BATF3", "CADM1"),
  "cDC2" = c("CD1C", "CLEC10A", "FCER1A", "ITGAX", "HLA-DQA1"),
  "pDC"  = c("LILRA4", "CLEC4C", "IL3RA", "TCF4", "IRF7")
)

## =========================================================================
## Step 1: blood-side DC subclustering
## =========================================================================
cat("\n===== Step 1: blood-side DC subclustering =====\n")

blood_seu <- readRDS(file.path(out_dir, "path1_sepsis_copd_integrated.rds"))
cat("Blood object:", ncol(blood_seu), "cells\n")

blood_dc <- subset(blood_seu, subset = cell_type == "DC")
rm(blood_seu); gc()
cat("Blood-side DC count:", ncol(blood_dc), "\n")
print(table(blood_dc$group))

# Fix out-of-sync @cells mapping after subset (Seurat v5 bug)
dc_counts <- GetAssayData(blood_dc, assay = "RNA", layer = "counts")
blood_dc[["RNA"]] <- CreateAssay5Object(counts = dc_counts)
rm(dc_counts); gc()

if (ncol(blood_dc) < 100) {
  cat("Too few blood-side DCs, skip subclustering, score directly by known markers\n")
} else {
  # Subclustering (lower resolution, DC population is small)
  old_res <- resolution
  resolution <- 0.6
  blood_dc <- integrate_and_annotate(blood_dc, marker_list = dc_markers,
                                     col_name = "dc_subtype")
  resolution <- old_res
  cat("\nBlood-side DC subtypes:\n")
  print(table(blood_dc$dc_subtype))
}

## =========================================================================
## Step 2: blood-side DC subtype proportions across five groups
## =========================================================================
cat("\n===== Step 2: blood-side DC subtype proportions across five groups =====\n")

if ("dc_subtype" %in% colnames(blood_dc@meta.data)) {
  # Subtype x group contingency table
  dc_tab <- table(blood_dc$dc_subtype, blood_dc$group)
  print(dc_tab)
  write.csv(dc_tab, file.path(out_dir, "path7_blood_dc_subtype_by_group.csv"))

  # Proportion (subtype fraction within each group)
  dc_prop <- sweep(dc_tab, 2, colSums(dc_tab), "/")
  print(round(dc_prop, 4))

  # Ro/e (observed/expected)
  expected <- outer(rowSums(dc_tab), colSums(dc_tab)) / sum(dc_tab)
  roe <- dc_tab / expected
  print(round(roe, 3))
  write.csv(as.data.frame.matrix(round(roe, 3)),
            file.path(out_dir, "path7_blood_dc_subtype_roe.csv"))
  cat("\n(Interpretation: if pDC Ro/e in Sepsis is far <1, it means pDC is most severely exhausted)\n")
}

## =========================================================================
## Step 3: lung-side DC subtype proportions (Control vs COPD)
## =========================================================================
cat("\n===== Step 3: lung-side DC subtype proportions =====\n")

lung_seu <- readRDS(file.path(out_dir, "path2_copd_lung.rds"))
cat("Lung object:", ncol(lung_seu), "cells\n")

# Lung-side DC-related cell types (Manuscript_Identity)
dc_types <- c("cDC1", "cDC2", "pDC", "DC_Mature", "DC_Langerhans")
lung_dc <- subset(lung_seu, subset = cell_type %in% dc_types &
                    disease %in% c("Control", "COPD"))
rm(lung_seu); gc()
cat("Lung-side DC count:", ncol(lung_dc), "\n")

lung_dc_tab <- table(lung_dc$cell_type, lung_dc$disease)
print(lung_dc_tab)
write.csv(lung_dc_tab, file.path(out_dir, "path7_lung_dc_subtype_by_disease.csv"))

# Proportion (fraction of each subtype within Control/COPD)
lung_dc_prop <- sweep(lung_dc_tab, 2, colSums(lung_dc_tab), "/")
print(round(lung_dc_prop, 4))

# Fold change (COPD / Control)
fc <- lung_dc_prop[, "COPD"] / lung_dc_prop[, "Control"]
cat("\nLung-side DC subtype COPD/Control proportion fold change:\n")
print(round(fc, 2))

## =========================================================================
## Step 4: summarize conclusions
## =========================================================================
cat("\n===== Step 4: summary =====\n")
cat("Core question: which specific subtype of DC is migrating/exhausted?\n")
cat("- Blood side: see which subtype (cDC1/cDC2/pDC) has the lowest Ro/e in the Sepsis group (most severe exhaustion)\n")
cat("- Lung side: see which subtype has the highest enrichment fold in COPD lung\n")
cat("- If the same subtype on both sides (especially pDC), the migration direction is clear\n")

cat("\nDC subset deep analysis complete\n")
cat("Output: path7_blood_dc_subtype_by_group.csv / path7_blood_dc_subtype_roe.csv / path7_lung_dc_subtype_by_disease.csv\n")
