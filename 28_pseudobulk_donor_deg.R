#!/usr/bin/env Rscript
# =============================================================================
# 28_pseudobulk_donor_deg.R
# Sepsis + COPD comorbidity analysis — donor-level pseudobulk differential expression
# =============================================================================
# Purpose: Switch key DEGs from "per-cell Wilcoxon" to "donor-level pseudobulk" to
#       eliminate pseudoreplication. Per-cell Wilcoxon on ~960k cells yields P < 1e-300
#       which carries almost no information and is always checked by reviewers. Standard
#       practice: aggregate raw counts per donor → edgeR/DESeq2, and report donor counts.
#
# Analysis:
#   Part 1: Blood monocyte pseudobulk — donor-level DEG for Sepsis/COPD/IC vs Healthy (edgeR)
#   Part 2: Lung myeloid pseudobulk — donor-level DEG for COPD vs Control (edgeR)
#   Part 3: Key gene check — whether SOCS3 / ISG / 149 shared genes remain significant at donor level
#   Part 4: 149-gene imprint — donor-level group comparison (replaces per-cell Wilcoxon)
#
# Dependencies: Seurat, edgeR (Bioconductor)
# =============================================================================

## ---- Auto-install dependencies (edgeR via Bioconductor) ----
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

## ---- Donor field detection (returns donor vector aligned with meta row names) ----
# Blood object: `patient` is the true donor, but some datasets (e.g. GSE151263/GSE175453)
#         only set `sample` (library-level), in which case patient is NA. So patient is
#         preferred, with sample as fallback where missing.
# Lung object: no patient, use sample (= Subject_Identity).
detect_donor <- function(seu) {
  meta <- seu@meta.data
  get <- function(f) if (f %in% colnames(meta)) as.character(meta[[f]]) else NULL
  patient <- get("patient"); sample <- get("sample")
  if (!is.null(patient)) {
    if (!is.null(sample)) {
      na_idx <- is.na(patient) | patient == ""
      patient[na_idx] <- sample[na_idx]
      return(list(field = "patient (sample as fallback)", vec = patient))
    }
    return(list(field = "patient", vec = patient))
  }
  for (cand in c("donor_id","donor","subject","subject_id","Subject_Identity","sample")) {
    v <- get(cand)
    if (!is.null(v)) return(list(field = cand, vec = v))
  }
  NULL
}

## ---- Donor-level aggregation + edgeR ----
run_pseudobulk <- function(seu, celltype, groups, donor_info, out_tag) {
  meta <- seu@meta.data
  donor_vec <- donor_info$vec
  keep_cell <- !is.na(meta$cell_type) & meta$cell_type == celltype &
               !is.na(meta$group) & meta$group %in% groups
  if (sum(keep_cell, na.rm = TRUE) == 0) {
    cat("[", out_tag, "] no ", celltype, " cells, skip\n", sep = ""); return(NULL)
  }
  meta <- meta[keep_cell, , drop = FALSE]
  donor_vec <- donor_vec[keep_cell]
  counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
  counts <- counts[, rownames(meta), drop = FALSE]

  donors <- unique(donor_vec)
  donors <- donors[!is.na(donors) & donors != ""]
  cts_list <- lapply(donors, function(d) {
    cells <- rownames(meta)[donor_vec == d]
    Matrix::rowSums(counts[, cells, drop = FALSE])
  })
  cts <- do.call(cbind, cts_list)
  colnames(cts) <- donors
  grp <- vapply(donors, function(d) unique(meta$group[donor_vec == d])[1],
                character(1))

  cat("\n[", out_tag, "] ", celltype, " number of donors:", ncol(cts), "\n", sep = "")
  print(table(grp))

  # Need at least 2 groups and >= 2 donors per group to run edgeR
  if (length(unique(grp)) < 2 || min(table(grp)) < 2) {
    cat("[", out_tag, "] insufficient donors, skip edgeR\n", sep = ""); return(NULL)
  }

  present_groups <- groups[groups %in% unique(grp)]
  grp_f <- factor(grp, levels = present_groups)
  y <- DGEList(counts = cts, group = grp_f)
  keep <- filterByExpr(y, min.count = 10, min.total.count = 15)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)                 # TMM
  design <- model.matrix(~ grp_f)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)

  res_all <- list()
  for (g in present_groups[-1]) {
    coef_name <- paste0("grp_f", g)
    if (!coef_name %in% colnames(design)) next
    qlf <- glmQLFTest(fit, coef = coef_name)
    tab <- topTags(qlf, n = Inf)$table
    tab$gene <- rownames(tab)
    res_all[[g]] <- tab
    write.csv(tab, file.path(out_dir,
              paste0("path28_", out_tag, "_pseudobulk_", g, "_vs_", present_groups[1], ".csv")),
              row.names = FALSE)
    cat("[", out_tag, "] ", g, " vs ", present_groups[1],
        ": significant genes (FDR<0.05) ", sum(tab$FDR < 0.05), " / ", nrow(tab),
        " (donors ", sum(grp_f == g), " vs ", sum(grp_f == present_groups[1]), ")\n", sep = "")
  }
  list(counts = cts, grp = grp, res = res_all)
}

## ============================================================================
## Part 1: Blood monocyte pseudobulk
## ============================================================================
cat("\n########## Part 1: blood monocyte pseudobulk ##########\n")
blood_file <- file.path(out_dir, "path1_sepsis_copd_integrated.rds")
if (!file.exists(blood_file)) {
  stop("Cannot find blood object path1_sepsis_copd_integrated.rds, please confirm script 01 has generated it")
}
blood <- readRDS(blood_file)
cat("blood object meta columns:\n  ", paste(colnames(blood@meta.data), collapse = "\n   "), "\n")
blood_donor <- detect_donor(blood)
if (is.null(blood_donor)) stop("blood object has no donor-level field (patient/sample), cannot run pseudobulk")
cat("blood object donor field:", blood_donor$field, "\n")

blood_groups <- c("Healthy", "Infection_Control", "COPD", "Sepsis", "Sepsis_Pneumonia")
pb_blood <- run_pseudobulk(blood, "Monocyte", blood_groups, blood_donor, "blood_mono")

## ============================================================================
## Part 2: Lung myeloid pseudobulk
## ============================================================================
cat("\n########## Part 2: lung myeloid pseudobulk ##########\n")
lung_seu <- tryCatch(
  read_gse136831(file.path(copd_dir, "GSE136831"), gse_id = "GSE136831"),
  error = function(e) { cat("Failed to read lung object:", conditionMessage(e), "\n"); NULL })
if (!is.null(lung_seu)) {
  cat("lung object meta columns:\n  ", paste(colnames(lung_seu@meta.data), collapse = "\n   "), "\n")
  lung_donor <- detect_donor(lung_seu)
  if (!"cell_type" %in% colnames(lung_seu@meta.data)) {
    cat("lung object has no cell_type field, define myeloid by myeloid marker score\n")
    lung_seu <- AddModuleScore(lung_seu,
      features = list(myeloid = c("LYZ","CD68","CD14","FCGR3A","CSF1R","ITGAX","CLEC4C")),
      name = "myeloid_score")
    lung_seu$cell_type <- ifelse(lung_seu$myeloid_score1 > 0, "Myeloid", "Other")
  }
  if ("disease" %in% colnames(lung_seu@meta.data) && !"group" %in% colnames(lung_seu@meta.data)) {
    lung_seu$group <- lung_seu$disease
  }
  if (is.null(lung_donor)) {
    cat("lung object has no donor-level field, skip lung pseudobulk\n")
  } else {
    cat("lung object donor field:", lung_donor$field, "\n")
    lung_groups <- sort(unique(lung_seu$group))
    pb_lung <- run_pseudobulk(lung_seu, "Myeloid", lung_groups, lung_donor, "lung_myeloid")
  }
}

## ============================================================================
## Part 3: Key gene check (are SOCS3 / ISG / 149 shared genes still significant at donor level)
## ============================================================================
cat("\n########## Part 3: key gene check ##########\n")
key_genes <- c("SOCS3", "SOCS1", "IRF7", "TLR9", "ISG15", "MX1", "OAS2",
               "IFI44L", "ISG20", "STAT1", "STAT2", "IFNAR2", "JAK1")

check_genes <- function(res_list, genes, tag) {
  lines <- c(paste0("===== ", tag, ": key genes in donor-level pseudobulk ====="))
  for (g in names(res_list)) {
    tab <- res_list[[g]]
    lines <- c(lines, paste0("--- ", tag, " ", g, " (donor-level edgeR) ---"))
    for (gn in genes) {
      if (gn %in% tab$gene) {
        r <- tab[tab$gene == gn, ]
        lines <- c(lines, sprintf("  %-10s logFC=%+.3f  FDR=%s",
                                  gn, r$logFC, format(r$FDR, scientific = TRUE, digits = 3)))
      } else {
        lines <- c(lines, sprintf("  %-10s not detected", gn))
      }
    }
  }
  lines
}

check_lines <- character()
if (!is.null(pb_blood)) check_lines <- c(check_lines,
  check_genes(pb_blood$res, key_genes, "blood_mono"))
if (exists("pb_lung") && !is.null(pb_lung)) check_lines <- c(check_lines,
  check_genes(pb_lung$res, key_genes, "lung_myeloid"))
writeLines(check_lines, file.path(out_dir, "path28_key_genes_pseudobulk.txt"))
cat(paste(check_lines, collapse = "\n"), "\n")

## ============================================================================
## Part 4: 149-gene imprint — donor-level group comparison
## ============================================================================
cat("\n########## Part 4: 149-gene imprint (donor-level) ##########\n")
sig_file <- file.path(out_dir, "path21_shared_up_genes.txt")
if (!file.exists(sig_file)) sig_file <- file.path(out_dir, "path2_shared_myeloid_genes.txt")
if (file.exists(sig_file)) {
  sig_genes <- readLines(sig_file)
  sig_genes <- sig_genes[sig_genes %in% rownames(blood)]
  cat("149 shared genes matched in the blood object:", length(sig_genes), "\n")

  data_mat <- GetAssayData(blood, assay = "RNA", layer = "data")
  mono_keep <- !is.na(blood$cell_type) & blood$cell_type == "Monocyte" &
               !is.na(blood$group) & blood$group %in% blood_groups
  mono_meta <- blood@meta.data[mono_keep, , drop = FALSE]
  donor_vec <- blood_donor$vec[mono_keep]
  data_mono <- data_mat[sig_genes, rownames(mono_meta), drop = FALSE]

  donors <- unique(donor_vec)
  donors <- donors[!is.na(donors) & donors != ""]
  imp <- sapply(donors, function(d) {
    cells <- rownames(mono_meta)[donor_vec == d]
    mean(colMeans(data_mono[, cells, drop = FALSE]))
  })
  imp_grp <- vapply(donors, function(d)
    unique(mono_meta$group[donor_vec == d])[1], character(1))
  imp_df <- data.frame(donor = donors, imprint = as.numeric(imp), group = imp_grp,
                       stringsAsFactors = FALSE)

  kw <- kruskal.test(imprint ~ group, data = imp_df)
  wt <- lapply(blood_groups[-1], function(g) {
    x <- imp_df$imprint[imp_df$group == g]; y <- imp_df$imprint[imp_df$group == "Healthy"]
    if (length(x) >= 2 && length(y) >= 2) {
      w <- wilcox.test(x, y)
      sprintf("  %s vs Healthy: donors %d vs %d, median %.4f vs %.4f, p=%s",
              g, length(x), length(y), median(x), median(y),
              format(w$p.value, scientific = TRUE, digits = 3))
    } else sprintf("  %s vs Healthy: not enough donors (<2)", g)
  })
  imp_lines <- c("===== 149-gene imprint (donor-level) =====",
                 paste0("total donors: ", nrow(imp_df)),
                 paste0("donors per group: ", paste(names(table(imp_grp)), table(imp_grp),
                                                sep = "=", collapse = ", ")),
                 paste0("Kruskal-Wallis p = ", format(kw$p.value, scientific = TRUE, digits = 3)),
                 unlist(wt))
  writeLines(imp_lines, file.path(out_dir, "path28_imprint_donor_stats.txt"))
  cat(paste(imp_lines, collapse = "\n"), "\n")
  write.csv(imp_df, file.path(out_dir, "path28_imprint_donor.csv"), row.names = FALSE)
} else {
  cat("149 shared gene file not found, skipping Part 4\n")
}

cat("\npseudobulk donor-level DEG analysis finished\n")
cat("output: path28_*_pseudobulk_*.csv / path28_key_genes_pseudobulk.txt / path28_imprint_donor_stats.txt\n")
