#!/usr/bin/env Rscript
# =============================================================================
# 69b_Intersection_Gene_Table.R   (standalone Part C of 69 + dimension diagnostics)
# Purpose: output the 95 up + 55 down intersection gene table (with blood/lung
#       logFC and FDR) under the "true labels" of the 68 cached pipeline, for
#       rebuilding ST4. When 69 ran as a whole, Part C failed at estimateDisp with
#       "nrow(design) disagrees with ncol(y)". This script prints the dimensions of
#       every step to locate which side / which object goes wrong.
# Run:   Rscript 69b_Intersection_Gene_Table.R
# Deps:  Seurat, edgeR, Matrix (~5-10 minutes)
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(edgeR); library(Matrix) })

this_file <- commandArgs(trailingOnly = FALSE)
.f <- grep("--file=", this_file, value = TRUE)
if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1]) else .this_file <- "."
out_dir <- dirname(normalizePath(.this_file))
source(file.path(out_dir, "00_config.R"))

LOGFCT <- 0.25; FDR_CUT <- 0.05

detect_donor <- function(seu) {
  meta <- seu@meta.data
  get <- function(f) if (f %in% colnames(meta)) as.character(meta[[f]]) else NULL
  patient <- get("patient"); sample <- get("sample")
  if (!is.null(patient)) {
    if (!is.null(sample)) { na_idx <- is.na(patient) | patient == ""; patient[na_idx] <- sample[na_idx] }
    return(patient)
  }
  for (cand in c("donor_id","donor","subject","subject_id","Subject_Identity","sample")) {
    v <- get(cand); if (!is.null(v)) return(v)
  }
  NULL
}
build_pb <- function(seu, keep, tag = "") {
  meta <- seu@meta.data; donor_vec <- detect_donor(seu)
  if (is.null(donor_vec)) stop(sprintf("[%s] no donor field found", tag))
  keep[is.na(keep)] <- FALSE
  counts_all <- GetAssayData(seu, assay = "RNA", layer = "counts")
  donor_all <- as.character(donor_vec); names(donor_all) <- rownames(meta)
  cells <- rownames(meta)[keep]; cells <- cells[cells %in% colnames(counts_all)]
  donor_vec <- donor_all[cells]
  meta <- meta[cells, , drop = FALSE]
  counts <- counts_all[, cells, drop = FALSE]
  donors <- unique(donor_vec); donors <- donors[!is.na(donors) & donors != ""]
  cts_list <- lapply(donors, function(d) Matrix::rowSums(counts[, names(donor_vec)[donor_vec == d], drop = FALSE]))
  cts <- do.call(cbind, cts_list); colnames(cts) <- donors
  grp <- vapply(donors, function(d) unique(meta$group[donor_vec == d])[1], character(1))
  cat(sprintf("[%s] donors=%d | counts %d x %d | group table: ", tag, length(donors), nrow(cts), ncol(cts)))
  print(table(grp, useNA = "ifany"))
  list(counts = cts, grp = grp)
}
prep_cache <- function(pb, ref, target, tag = "") {
  grp_f <- factor(pb$grp, levels = c(ref, target))
  y0 <- DGEList(counts = pb$counts, group = grp_f)
  keep <- filterByExpr(y0, min.count = 10, min.total.count = 15)
  y0 <- y0[keep, , keep.lib.sizes = FALSE]
  y0 <- calcNormFactors(y0)
  cat(sprintf("[%s] cache: kept genes %d / %d, nf length %d\n", tag, sum(keep), nrow(pb$counts), length(y0$samples$norm.factors)))
  list(keep = keep, nf = y0$samples$norm.factors, levels = c(ref, target))
}
deg_arms <- function(pb, cache, grp_vec, tag = "") {
  grp_f <- factor(grp_vec, levels = cache$levels)
  cat(sprintf("[%s] deg_arms: grp length %d | counts of two levels %s\n", tag, length(grp_f),
              paste(table(grp_f, useNA = "ifany"), collapse = "/")))
  if (length(unique(grp_f)) < 2 || min(table(grp_f)) < 2) { cat("  !! fewer than two groups, returning NULL\n"); return(NULL) }
  counts_sub <- pb$counts[cache$keep, , drop = FALSE]
  cat(sprintf("[%s] counts_sub %d x %d | nf length %d\n", tag, nrow(counts_sub), ncol(counts_sub), length(cache$nf)))
  y <- DGEList(counts = counts_sub, group = grp_f, norm.factors = cache$nf)
  cat(sprintf("[%s] DGEList ncol=%d lib.size=%d items, duplicated sample names=%d\n", tag, ncol(y), length(y$samples$lib.size),
              sum(duplicated(colnames(y)))))
  design <- model.matrix(~ grp_f)
  cat(sprintf("[%s] design %d x %d | y ncol %d\n", tag, nrow(design), ncol(design), ncol(y)))
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  tab <- topTags(glmQLFTest(fit, coef = 2), n = Inf)$table
  list(tab = tab, up = rownames(tab)[tab$logFC > LOGFCT & tab$FDR < FDR_CUT],
       down = rownames(tab)[tab$logFC < -LOGFCT & tab$FDR < FDR_CUT])
}

cat("===== 69b: intersection gene table (with diagnostics) =====\n")

blood <- readRDS(file.path(out_dir, "path1_sepsis_copd_integrated.rds"))
b_keep <- blood$cell_type == "Monocyte" & blood$group %in% c("Healthy", "Sepsis")
pb_blood <- build_pb(blood, b_keep, tag = "blood_mono")
cache_b <- prep_cache(pb_blood, "Healthy", "Sepsis", tag = "blood_mono")
b <- tryCatch(deg_arms(pb_blood, cache_b, pb_blood$grp, tag = "blood_mono"), error = function(e) {
  cat("\n!! blood-side deg_arms error:", conditionMessage(e), "\n")
  NULL
})
cat("blood-side up=", if (is.null(b)) NA else length(b$up), "\n")

lung <- readRDS(file.path(out_dir, "path2_copd_lung.rds"))
if (!"group" %in% colnames(lung@meta.data)) lung$group <- lung$disease
ct_all <- unique(as.character(lung$cell_type))
myeloid_types <- grep("Macrophage|Monocyte|monocyte|macrophage|Myeloid|myeloid|DC|dendritic", ct_all, value = TRUE)
if (length(myeloid_types) == 0) myeloid_types <- intersect(c("Myeloid", "Myeloid cells"), ct_all)
cat("lung myeloid cell_type candidates:", paste(myeloid_types, collapse = ", "), "\n")
l_keep <- lung$cell_type %in% myeloid_types & lung$group %in% c("Control", "COPD")
pb_lung <- build_pb(lung, l_keep, tag = "lung_myeloid")
cache_l <- prep_cache(pb_lung, "Control", "COPD", tag = "lung_myeloid")
l <- tryCatch(deg_arms(pb_lung, cache_l, pb_lung$grp, tag = "lung_myeloid"), error = function(e) {
  cat("\n!! lung-side deg_arms error:", conditionMessage(e), "\n")
  NULL
})
cat("lung-side up=", if (is.null(l)) NA else length(l$up), "\n")

if (is.null(b) || is.null(l)) {
  cat("\nat least one side failed -> intersection table not written. Please paste back all [blood_mono]/[lung_myeloid] printed lines above.\n")
  quit(status = 1)
}

up <- intersect(b$up, l$up); dn <- intersect(b$down, l$down)
tb <- b$tab; tl <- l$tab
out_up <- data.frame(arm = "up_both", gene = up,
                     blood_logFC = tb[up, "logFC"], blood_FDR = tb[up, "FDR"],
                     lung_logFC  = tl[up, "logFC"], lung_FDR  = tl[up, "FDR"],
                     stringsAsFactors = FALSE)
out_dn <- data.frame(arm = "down_both", gene = dn,
                     blood_logFC = tb[dn, "logFC"], blood_FDR = tb[dn, "FDR"],
                     lung_logFC  = tl[dn, "logFC"], lung_FDR  = tl[dn, "FDR"],
                     stringsAsFactors = FALSE)
out_g <- rbind(out_up, out_dn)
write.csv(out_g, file.path(out_dir, "69_observed_intersection_genes.csv"), row.names = FALSE)
cat(sprintf("\nintersection genes: up=%d down=%d (%d rows total) -> 69_observed_intersection_genes.csv\n",
            length(up), length(dn), nrow(out_g)))
cat("===== 69b done =====\n")
