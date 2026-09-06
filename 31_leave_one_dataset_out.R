#!/usr/bin/env Rscript
# =============================================================================
# 31 leave-one-dataset-out sensitivity analysis
# Sepsis + COPD comorbidity analysis — leave-one-dataset-out sensitivity analysis
# =============================================================================
# Purpose: address the reviewer's concern about "group confounded with batch". Remove
#          one dataset at a time (leave-one-dataset-out) and recompute the three core
#          conclusions to check whether they are robust:
#        ① blood-side direction concordance of the 149 shared genes (donor-level edgeR,
#          number surviving with logFC>0)
#        ② donor-level up-regulation of SOCS3 (Sepsis vs Healthy logFC/FDR)
#        ③ blood monocyte inflammatory imprint gradient (healthy < COPD < sepsis,
#          donor-level Kruskal-Wallis)
#
# Criteria (same as 28b):
#   blood-side DEG: cell_type == "Monocyte", condition in (Sepsis, Healthy)
#   imprint gradient: cell_type == "Monocyte", group in (Healthy, COPD, Sepsis)
#   donor field: patient first, sample as fallback
#
# Outputs:
#   path31_leave_one_out_stats.txt     main result summary
#   path31_lo_149_consistency.csv      149 direction-concordance count + SOCS3 after
#                                      removing each dataset
#   path31_lo_imprint.csv              imprint medians + KW p after removing each dataset
# =============================================================================

## ---- Auto-install dependencies ----
if (!requireNamespace("edgeR", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("edgeR", update = FALSE, ask = FALSE)
}
suppressPackageStartupMessages({
  library(Seurat); library(edgeR); library(Matrix)
})

## ---- Paths ----
this_file <- commandArgs(trailingOnly = FALSE)
.f <- grep("--file=", this_file, value = TRUE)
if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1]) else .this_file <- "."
out_dir <- dirname(normalizePath(.this_file))
source(file.path(out_dir, "00_config.R"))

## ---- Donor field detection (patient first, sample as fallback) [same as 28b] ----
detect_donor <- function(seu) {
  meta <- seu@meta.data
  get <- function(f) if (f %in% colnames(meta)) as.character(meta[[f]]) else NULL
  patient <- get("patient"); sample <- get("sample")
  if (!is.null(patient)) {
    if (!is.null(sample)) {
      na_idx <- is.na(patient) | patient == ""
      patient[na_idx] <- sample[na_idx]
      return(list(field = "patient(+sample fallback)", vec = patient))
    }
    return(list(field = "patient", vec = patient))
  }
  for (cand in c("donor_id","donor","subject","subject_id","Subject_Identity","sample")) {
    v <- get(cand)
    if (!is.null(v)) return(list(field = cand, vec = v))
  }
  NULL
}

## ---- Donor-level pseudobulk + edgeR, returns the full DEG table [same as 28b] ----
pseudobulk_deg <- function(seu, keep_cell, donor_vec, group_vec, groups, tag) {
  meta <- seu@meta.data
  d <- donor_vec[keep_cell]; g <- group_vec[keep_cell]
  counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
  counts <- counts[, rownames(meta)[keep_cell], drop = FALSE]

  donors <- unique(d); donors <- donors[!is.na(donors) & donors != ""]
  cts <- do.call(cbind, lapply(donors, function(x)
    Matrix::rowSums(counts[, d == x, drop = FALSE])))
  colnames(cts) <- donors
  grp <- vapply(donors, function(x) unique(g[d == x])[1], character(1))
  cat("    [", tag, "] donors:", ncol(cts), " ", sep = "")
  print(table(grp))

  present <- groups[groups %in% unique(grp)]
  if (length(present) < 2 || min(table(grp)) < 2) {
    cat("    [", tag, "] too few donors, skipping\n", sep = ""); return(NULL)
  }
  grp_f <- factor(grp, levels = present)
  y <- DGEList(counts = cts, group = grp_f)
  keep <- filterByExpr(y, min.count = 10, min.total.count = 15)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~ grp_f)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = paste0("grp_f", present[2]))
  tab <- topTags(qlf, n = Inf)$table
  tab$gene <- rownames(tab)
  tab
}

## ---- Donor-level imprint scoring (mean z-score of the 149 genes) ----
imprint_donor <- function(seu, keep_cell, donor_vec, genes) {
  meta <- seu@meta.data
  data_mat <- GetAssayData(seu, assay = "RNA", layer = "data")
  genes <- genes[genes %in% rownames(data_mat)]
  if (length(genes) < 10) return(NULL)
  data_mat <- data_mat[genes, , drop = FALSE]
  # per-gene z-score across cells (data layer is already log-normalized)
  data_z <- t(scale(t(as.matrix(data_mat))))
  cells <- rownames(meta)[keep_cell]
  data_z <- data_z[, cells, drop = FALSE]
  d <- donor_vec[keep_cell]
  donors <- unique(d); donors <- donors[!is.na(donors) & donors != ""]
  scores <- vapply(donors, function(x)
    mean(colMeans(data_z[, d == x, drop = FALSE], na.rm = TRUE)), numeric(1))
  names(scores) <- donors
  scores
}

## ============================================================================
## Read the blood object + the original 149
## ============================================================================
cat("\n########## reading blood object + original 149 ##########\n")
blood_file <- file.path(out_dir, "path1_sepsis_copd_integrated.rds")
if (!file.exists(blood_file)) stop("cannot find path1_sepsis_copd_integrated.rds")
blood <- readRDS(blood_file)
blood_donor <- detect_donor(blood)
cat("Blood object donor field:", blood_donor$field, "\n")

ds_field <- if ("dataset" %in% colnames(blood@meta.data)) "dataset" else "orig.ident"
cat("Dataset field:", ds_field, "\n")
datasets <- sort(unique(as.character(blood@meta.data[[ds_field]])))
datasets <- datasets[!is.na(datasets) & datasets != ""]
cat("Blood object datasets (", length(datasets), "): ", paste(datasets, collapse=", "), "\n", sep="")

orig_file <- file.path(out_dir, "path2_shared_myeloid_genes.txt")
if (!file.exists(orig_file)) orig_file <- file.path(out_dir, "path21_shared_up_genes.txt")
orig_149 <- if (file.exists(orig_file)) readLines(orig_file) else character(0)
orig_149 <- orig_149[orig_149 != ""]
cat("Original 149 gene count:", length(orig_149), "\n")

## ============================================================================
## Baseline (full data, no dataset removed)
## ============================================================================
cat("\n########## baseline (full data) ##########\n")
base_keep <- !is.na(blood$cell_type) & blood$cell_type == "Monocyte" &
             !is.na(blood$condition) & blood$condition %in% c("Sepsis", "Healthy")
base_deg <- pseudobulk_deg(blood, base_keep, blood_donor$vec, blood$condition,
                           c("Healthy", "Sepsis"), "baseline blood monocyte")
if (is.null(base_deg)) stop("baseline pseudobulk failed")
base_149_ok <- sum(orig_149 %in% base_deg$gene[base_deg$logFC > 0])
base_socs3 <- base_deg[base_deg$gene == "SOCS3", ]
cat("Baseline: original 149 blood-side direction concordant =", base_149_ok, "/", length(orig_149), "\n")

## ============================================================================
## Part 1: leave-one-dataset-out (149 direction concordance + SOCS3)
## ============================================================================
cat("\n########## Part 1: leave-one-dataset-out (149 + SOCS3) ##########\n")
lo_rows <- list()
for (ds in datasets) {
  cat("\n--- removing", ds, "---\n")
  keep <- base_keep & blood@meta.data[[ds_field]] != ds
  deg <- pseudobulk_deg(blood, keep, blood_donor$vec, blood$condition,
                        c("Healthy", "Sepsis"), paste0("remove ", ds))
  if (is.null(deg)) {
    lo_rows[[ds]] <- data.frame(dataset_removed = ds, n_149_ok = NA,
                                socs3_logFC = NA, socs3_FDR = NA)
    next
  }
  n_ok <- sum(orig_149 %in% deg$gene[deg$logFC > 0])
  s3 <- deg[deg$gene == "SOCS3", ]
  lo_rows[[ds]] <- data.frame(
    dataset_removed = ds,
    n_149_ok = n_ok,
    socs3_logFC = if (nrow(s3) > 0) s3$logFC[1] else NA,
    socs3_FDR  = if (nrow(s3) > 0) s3$FDR[1] else NA)
  cat("  original 149 direction concordant:", n_ok, "/", length(orig_149),
      " | SOCS3 logFC=", if (nrow(s3)>0) round(s3$logFC[1],3) else NA, "\n", sep="")
}
lo_df <- do.call(rbind, lo_rows)

## ============================================================================
## Part 2: leave-one-dataset-out (imprint gradient)
## ============================================================================
cat("\n########## Part 2: leave-one-dataset-out (imprint gradient) ##########\n")
imp_keep <- !is.na(blood$cell_type) & blood$cell_type == "Monocyte" &
            !is.na(blood$group) & blood$group %in% c("Healthy", "COPD", "Sepsis")
imp_rows <- list()
for (ds in c(datasets, "NONE")) {
  if (ds == "NONE") {
    tag <- "baseline (full data)"
    keep <- imp_keep
  } else {
    tag <- paste0("remove ", ds)
    keep <- imp_keep & blood@meta.data[[ds_field]] != ds
  }
  cat("\n--- ", tag, " ---\n", sep="")
  grp <- blood$group[keep]
  if (length(unique(grp)) < 2) { cat("  fewer than 2 groups, skipping\n"); next }
  sc <- imprint_donor(blood, keep, blood_donor$vec, orig_149)
  if (is.null(sc)) { cat("  imprint failed, skipping\n"); next }
  donor_grp <- vapply(names(sc), function(x) {
    g <- unique(blood$group[blood_donor$vec == x & keep])
    if (length(g) > 0) g[1] else NA_character_
  }, character(1))
  med <- tapply(sc, donor_grp, median, na.rm = TRUE)
  kw <- tryCatch(kruskal.test(sc ~ as.factor(donor_grp))$p.value,
                 error = function(e) NA)
  imp_rows[[ds]] <- data.frame(
    dataset_removed = ds,
    Healthy_med = if ("Healthy" %in% names(med)) round(med["Healthy"],4) else NA,
    COPD_med    = if ("COPD" %in% names(med)) round(med["COPD"],4) else NA,
    Sepsis_med  = if ("Sepsis" %in% names(med)) round(med["Sepsis"],4) else NA,
    KW_p = kw)
  cat("  medians:", paste(names(med), round(med,4), sep="=", collapse=" "),
      " | KW p=", if (is.na(kw)) "NA" else format(kw, scientific=TRUE, digits=3), "\n", sep="")
}
imp_df <- do.call(rbind, imp_rows)

## ============================================================================
## Output
## ============================================================================
write.csv(lo_df, file.path(out_dir, "path31_lo_149_consistency.csv"), row.names = FALSE)
write.csv(imp_df, file.path(out_dir, "path31_lo_imprint.csv"), row.names = FALSE)

lines <- c(
  "===== leave-one-dataset-out sensitivity analysis =====",
  paste0("Blood object datasets (", length(datasets), "): ", paste(datasets, collapse=", ")),
  paste0("Original 149 gene count: ", length(orig_149)),
  "",
  "[Baseline (full data)]",
  paste0("  original 149 blood-side direction concordant: ", base_149_ok, "/", length(orig_149)),
  paste0("  SOCS3 logFC=", round(base_socs3$logFC[1],3),
         " FDR=", format(base_socs3$FDR[1], scientific=TRUE, digits=3)),
  "",
  "[Part 1: 149 direction concordance + SOCS3 leave-one-out]",
  "  (the original-149 blood-side concordance count should stay near baseline; SOCS3 logFC should stay >0 with significant FDR)",
  paste0("  see path31_lo_149_consistency.csv (", nrow(lo_df), " datasets)"),
  "",
  "[Part 2: imprint gradient leave-one-out]",
  "  (the Healthy < COPD < Sepsis median ordering should hold; KW p should stay significant)",
  paste0("  see path31_lo_imprint.csv (", nrow(imp_df), " rows)")
)
writeLines(lines, file.path(out_dir, "path31_leave_one_out_stats.txt"))
cat("\n", paste(lines, collapse = "\n"), "\n", sep = "")

cat("\n31 leave-one-dataset-out sensitivity analysis done\n")
