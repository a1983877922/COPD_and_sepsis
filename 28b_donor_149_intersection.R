#!/usr/bin/env Rscript
# =============================================================================
# 28b - Donor-level re-derivation of the 149-gene intersection
# Sepsis + COPD comorbidity analysis - donor-level re-derivation of the "149 shared myeloid gene" intersection
# =============================================================================
# Purpose: the original 149 shared genes were obtained by "per-cell Wilcoxon, avg_log2FC>0
#       intersection", which is affected by pseudo-replication. This script re-derives the same
#       intersection using donor-level pseudobulk + edgeR, answering the reviewers' inevitable
#       question: "how much of the 149 shared program survives at donor level?"
#
# Definition (identical to the original cross-tissue script 02):
#   Lung side: cell_category == "Myeloid", COPD vs Control
#   Blood side: cell_type == "Monocyte", Sepsis vs Healthy
#   Intersection: logFC > 0 on both sides (consistent direction)
#
# Output:
#   path28b_blood_mono_pseudobulk_deg.csv    blood monocyte donor-level DEG (all genes)
#   path28b_lung_myeloid_pseudobulk_deg.csv  lung myeloid donor-level DEG (all genes)
#   path28b_intersection_stats.txt           intersection statistics + overlap with the original 149
#   path28b_shared_donor_genes.txt           donor-level shared upregulated genes
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

## ---- Donor field detection (patient first, sample as fallback where missing) ----
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

## ---- Donor-level pseudobulk + edgeR, returns the full DEG table ----
# keep_cell: logical vector (whole object); donor_vec/group_vec: vectors aligned with object cells
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
  cat("\n[", tag, "] donor count:", ncol(cts), "\n", sep = "")
  print(table(grp))

  present <- groups[groups %in% unique(grp)]
  if (length(present) < 2 || min(table(grp)) < 2) {
    cat("[", tag, "] insufficient donors, skipping\n", sep = ""); return(NULL)
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
  list(tab = tab, grp = grp, present = present)
}

## ============================================================================
## Part 1: blood monocytes Sepsis vs Healthy (donor level)
## ============================================================================
cat("\n########## Part 1: blood monocytes Sepsis vs Healthy (donor level) ##########\n")
blood_file <- file.path(out_dir, "path1_sepsis_copd_integrated.rds")
if (!file.exists(blood_file)) stop("Cannot find path1_sepsis_copd_integrated.rds")
blood <- readRDS(blood_file)
blood_donor <- detect_donor(blood)
cat("Blood object donor field:", blood_donor$field, "\n")
# Use the condition field (consistent with the original script 02), print values for checking
cat("condition values:\n"); print(table(blood$condition))

blood_keep <- !is.na(blood$cell_type) & blood$cell_type == "Monocyte" &
              !is.na(blood$condition) & blood$condition %in% c("Sepsis", "Healthy")
blood_res <- pseudobulk_deg(blood, blood_keep, blood_donor$vec, blood$condition,
                            c("Healthy", "Sepsis"), "blood monocytes")
if (!is.null(blood_res)) {
  write.csv(blood_res$tab, file.path(out_dir, "path28b_blood_mono_pseudobulk_deg.csv"),
            row.names = FALSE)
}

## ============================================================================
## Part 2: lung myeloid COPD vs Control (donor level)
## ============================================================================
cat("\n########## Part 2: lung myeloid COPD vs Control (donor level) ##########\n")
lung_seu <- tryCatch(
  read_gse136831(file.path(copd_dir, "GSE136831"), gse_id = "GSE136831"),
  error = function(e) { cat("Failed to read the lung object:", conditionMessage(e), "\n"); NULL })
lung_res <- NULL
if (!is.null(lung_seu)) {
  lung_donor <- detect_donor(lung_seu)
  cat("Lung object donor field:", if (is.null(lung_donor)) "none" else lung_donor$field, "\n")
  if (!"cell_category" %in% colnames(lung_seu@meta.data)) {
    cat("Lung object has no cell_category field, defining it from myeloid terms in cell_type\n")
    lung_seu$cell_category <- ifelse(
      grepl("Macro|Mono|DC|Myeloid|mast|mast", lung_seu$cell_type, ignore.case = TRUE),
      "Myeloid", "Other")
  }
  cat("disease values:\n"); print(table(lung_seu$disease))
  lung_keep <- !is.na(lung_seu$cell_category) & lung_seu$cell_category == "Myeloid" &
               !is.na(lung_seu$disease) & lung_seu$disease %in% c("COPD", "Control")
  lung_res <- pseudobulk_deg(lung_seu, lung_keep, lung_donor$vec, lung_seu$disease,
                             c("Control", "COPD"), "lung myeloid")
  if (!is.null(lung_res)) {
    write.csv(lung_res$tab, file.path(out_dir, "path28b_lung_myeloid_pseudobulk_deg.csv"),
              row.names = FALSE)
  }
}

## ============================================================================
## Part 3: intersection (donor level) + overlap with the original 149
## ============================================================================
cat("\n########## Part 3: Intersection + overlap with the original 149 ##########\n")
if (is.null(blood_res) || is.null(lung_res)) {
  stop("Blood or lung pseudobulk failed, cannot compute the intersection")
}

blood_up  <- blood_res$tab$gene[blood_res$tab$logFC > 0]
lung_up   <- lung_res$tab$gene[lung_res$tab$logFC > 0]
shared_donor <- intersect(blood_up, lung_up)

blood_up_sig <- blood_res$tab$gene[blood_res$tab$logFC > 0 & blood_res$tab$FDR < 0.05]
lung_up_sig  <- lung_res$tab$gene[lung_res$tab$logFC > 0 & lung_res$tab$FDR < 0.05]
shared_donor_sig <- intersect(blood_up_sig, lung_up_sig)

# Read the original 149
orig_file <- file.path(out_dir, "path2_shared_myeloid_genes.txt")
if (!file.exists(orig_file)) orig_file <- file.path(out_dir, "path21_shared_up_genes.txt")
orig_149 <- if (file.exists(orig_file)) readLines(orig_file) else character(0)
orig_149 <- orig_149[orig_149 != ""]

overlap <- intersect(orig_149, shared_donor)
dropped <- setdiff(orig_149, shared_donor)     # within the original 149, no longer concordant at donor level
new     <- setdiff(shared_donor, orig_149)     # newly shared at donor level

# Direction distribution of the original 149 at donor level on both sides (how many have logFC>0 in blood/lung)
in_blood_up <- orig_149[orig_149 %in% blood_up]
in_lung_up  <- orig_149[orig_149 %in% lung_up]
both_ok     <- orig_149[orig_149 %in% shared_donor]

lines <- c(
  "===== Donor-level re-derivation of the 149-gene intersection =====",
  paste0("Blood monocyte donor-level upregulated genes (logFC>0): ", length(blood_up),
         " (significant FDR<0.05: ", length(blood_up_sig), ")"),
  paste0("Lung myeloid donor-level upregulated genes (logFC>0): ", length(lung_up),
         " (significant FDR<0.05: ", length(lung_up_sig), ")"),
  "",
  paste0("[Donor-level shared upregulation] (blood ∩ lung, logFC>0): ", length(shared_donor)),
  paste0("[Donor-level shared upregulation + significant] (blood ∩ lung, logFC>0 & FDR<0.05): ", length(shared_donor_sig)),
  "",
  paste0("Original 149 gene count: ", length(orig_149)),
  paste0("Of the original 149, still upregulated in blood at donor level: ", length(in_blood_up), " / ", length(orig_149)),
  paste0("Of the original 149, still upregulated in lung at donor level: ", length(in_lung_up), " / ", length(orig_149)),
  paste0("Of the original 149, still upregulated on both sides at donor level (surviving): ", length(both_ok), " / ", length(orig_149)),
  paste0("Of the original 149, no longer shared at donor level (dropped): ", length(dropped)),
  paste0("Newly shared genes at donor level: ", length(new)),
  "",
  paste0("Surviving genes (top 60): ", paste(head(both_ok, 60), collapse = ", ")),
  paste0("Dropped genes (in the original 149 but not shared at donor level): ", paste(head(dropped, 60), collapse = ", "))
)
writeLines(lines, file.path(out_dir, "path28b_intersection_stats.txt"))
cat(paste(lines, collapse = "\n"), "\n")

if (length(shared_donor) > 0)
  writeLines(shared_donor, file.path(out_dir, "path28b_shared_donor_genes.txt"))

cat("\n28b donor-level re-derivation of the 149-gene intersection finished\n")
