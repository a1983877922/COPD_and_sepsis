#!/usr/bin/env Rscript
# =============================================================================
# 33_Imprint scoring method comparison (z-score vs ssGSEA)
# Sepsis + COPD comorbidity analysis - imprint scoring method comparison (z-score vs ssGSEA)
# =============================================================================
# Purpose: the original imprint score used the "149-gene z-score mean", which is sensitive
#       to dataset baseline differences (script 32 showed within-group dataset differences
#       are comparable to the between-group gradient). This script re-scores with rank-based
#       ssGSEA and compares the two methods on:
#         (1) between-group gradient (Healthy < COPD < Sepsis medians)
#         (2) within-group dataset differences (imprint range across datasets, reflecting residual batch effects)
#       If the within-group dataset difference clearly shrinks with ssGSEA, the z-score gradient
#       was inflated by baseline differences, making ssGSEA the more credible version.
#
# ssGSEA implementation: prefer GSVA::gsva(method="ssgsea"); if GSVA is unavailable, fall
#       back to a hand-written rank-based score (mean normalized rank of the 149 genes within
#       each donor, equivalent to a simplified ssGSEA).
#
# Output:
#   path33_imprint_two_methods.csv     long table (donor/imprint/group/dataset/method)
#   path33_method_comparison.txt       between-group gradient + within-group dataset range for both methods
#   path33_imprint_ssgsea.pdf          ssGSEA facet plot (facet by group, x=dataset)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(ggplot2)
})

## ---- Paths ----
this_file <- commandArgs(trailingOnly = FALSE)
.f <- grep("--file=", this_file, value = TRUE)
if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1]) else .this_file <- "."
out_dir <- dirname(normalizePath(.this_file))
source(file.path(out_dir, "00_config.R"))

## ---- Donor field detection (patient first, sample as fallback) ----
detect_donor <- function(seu) {
  meta <- seu@meta.data
  get <- function(f) if (f %in% colnames(meta)) as.character(meta[[f]]) else NULL
  patient <- get("patient"); sample <- get("sample")
  if (!is.null(patient)) {
    if (!is.null(sample)) {
      na_idx <- is.na(patient) | patient == ""
      patient[na_idx] <- sample[na_idx]
      return(patient)
    }
    return(patient)
  }
  for (cand in c("donor_id","donor","subject","subject_id","Subject_Identity","sample")) {
    v <- get(cand); if (!is.null(v)) return(v)
  }
  NULL
}

## ---- Hand-written ssGSEA (rank-based, donor level) fallback ----
# For each donor expression vector, rank genes by expression; imprint = mean normalized
# rank (0-1) of the genes in the set, centered by subtracting 0.5. Equivalent to a
# simplified ssGSEA (mean within-set rank only).
rank_based_es <- function(expr_sym, genes) {
  # expr_sym: gene x donor matrix (logCPM)
  N <- nrow(expr_sym)
  genes <- genes[genes %in% rownames(expr_sym)]
  apply(expr_sym, 2, function(sample_expr) {
    r <- rank(sample_expr, ties.method = "average")
    names(r) <- rownames(expr_sym)
    mean(r[genes]) / N - 0.5
  })
}

## ============================================================================
## Load blood object + original 149 + donor-level pseudobulk
## ============================================================================
cat("\n########## Read data + donor-level pseudobulk ##########\n")
blood_file <- file.path(out_dir, "path1_sepsis_copd_integrated.rds")
if (!file.exists(blood_file)) stop("Cannot find path1_sepsis_copd_integrated.rds")
blood <- readRDS(blood_file)
donor_vec <- detect_donor(blood)
ds_field <- if ("dataset" %in% colnames(blood@meta.data)) "dataset" else "orig.ident"

orig_file <- file.path(out_dir, "path2_shared_myeloid_genes.txt")
if (!file.exists(orig_file)) orig_file <- file.path(out_dir, "path21_shared_up_genes.txt")
orig_149 <- if (file.exists(orig_file)) readLines(orig_file) else character(0)
orig_149 <- orig_149[orig_149 != ""]

keep <- !is.na(blood$cell_type) & blood$cell_type == "Monocyte" &
        !is.na(blood$group) & blood$group %in% c("Healthy", "COPD", "Sepsis")

# Donor-level pseudobulk counts
counts <- GetAssayData(blood, assay = "RNA", layer = "counts")
counts <- counts[, rownames(blood@meta.data)[keep], drop = FALSE]
d <- donor_vec[keep]
donors <- unique(d); donors <- donors[!is.na(donors) & donors != ""]
pb <- do.call(cbind, lapply(donors, function(x)
  Matrix::rowSums(counts[, d == x, drop = FALSE])))
colnames(pb) <- donors

# logCPM (donor level)
lib <- colSums(pb)
cpm <- sweep(pb, 2, lib, "/") * 1e6
logcpm <- log2(cpm + 1)
cat("Donor-level pseudobulk:", nrow(logcpm), "genes x", ncol(logcpm), "donors\n")

# Donor group / dataset
grp <- vapply(donors, function(x) {
  g <- unique(blood$group[donor_vec == x & keep]); if (length(g) > 0) g[1] else NA_character_
}, character(1))
ds <- vapply(donors, function(x) {
  v <- unique(blood@meta.data[[ds_field]][donor_vec == x & keep]); if (length(v) > 0) v[1] else NA_character_
}, character(1))

## ============================================================================
## Method 1: z-score mean (original method, cell level)
## ============================================================================
cat("\n########## Method 1: z-score mean (original method) ##########\n")
data_mat <- GetAssayData(blood, assay = "RNA", layer = "data")
genes <- orig_149[orig_149 %in% rownames(data_mat)]
data_mat <- data_mat[genes, , drop = FALSE]
data_z <- t(scale(t(as.matrix(data_mat))))
cells <- rownames(blood@meta.data)[keep]
data_z <- data_z[, cells, drop = FALSE]
d2 <- donor_vec[keep]
zscore_imp <- vapply(donors, function(x)
  mean(colMeans(data_z[, d2 == x, drop = FALSE], na.rm = TRUE)), numeric(1))

## ============================================================================
## Method 2: ssGSEA (rank-based)
## ============================================================================
cat("\n########## Method 2: ssGSEA (rank-based) ##########\n")
if (requireNamespace("GSVA", quietly = TRUE)) {
  suppressPackageStartupMessages(library(GSVA))
  gs <- list(shared_myeloid = orig_149)
  ssgsea_imp <- tryCatch({
    if (exists("ssgseaParam")) {
      cat("Using GSVA new ssgseaParam API\n")
      param <- ssgseaParam(logcpm, gs)
      ss <- gsva(param)
      as.numeric(ss["shared_myeloid", donors, drop = TRUE])
    } else {
      cat("Using GSVA legacy gsva(matrix, ...) API\n")
      ss <- gsva(logcpm, gs, method = "ssgsea", min.sz = 10, max.sz = 10000,
                 verbose = FALSE)
      as.numeric(ss["shared_myeloid", donors])
    }
  }, error = function(e) {
    cat("GSVA call failed, falling back to hand-written rank-based:", conditionMessage(e), "\n")
    rank_based_es(logcpm, orig_149)
  })
} else {
  cat("GSVA not installed, using hand-written rank-based scoring\n")
  ssgsea_imp <- rank_based_es(logcpm, orig_149)
}
names(ssgsea_imp) <- donors

## ============================================================================
## Assemble comparison long table + statistics
## ============================================================================
df_z <- data.frame(donor = donors, imprint = zscore_imp, group = grp, dataset = ds,
                   method = "z-score", stringsAsFactors = FALSE)
df_s <- data.frame(donor = donors, imprint = ssgsea_imp, group = grp, dataset = ds,
                   method = "ssGSEA", stringsAsFactors = FALSE)
df_all <- rbind(df_z, df_s)
df_all <- df_all[!is.na(df_all$group) & !is.na(df_all$dataset), ]
df_all$group <- factor(df_all$group, levels = c("Healthy", "COPD", "Sepsis"))
write.csv(df_all, file.path(out_dir, "path33_imprint_two_methods.csv"), row.names = FALSE)

## ---- Comparison: between-group gradient + within-group dataset range ----
summ <- function(df, m) {
  s <- df[df$method == m, ]
  med_group <- tapply(s$imprint, s$group, median, na.rm = TRUE)
  # Within-group dataset range (range across datasets for Healthy and Sepsis, take the max)
  range_ds <- tapply(s$imprint, list(s$group, s$dataset),
                     function(x) median(x, na.rm = TRUE))
  range_by_group <- sapply(rownames(range_ds), function(g) {
    v <- range_ds[g, !is.na(range_ds[g, ])]
    if (length(v) >= 2) diff(range(v)) else NA_real_
  })
  c(Healthy_med = med_group["Healthy"], COPD_med = med_group["COPD"],
    Sepsis_med = med_group["Sepsis"],
    Healthy_ds_range = range_by_group["Healthy"],
    Sepsis_ds_range = range_by_group["Sepsis"])
}
sum_z <- summ(df_all, "z-score")
sum_s <- summ(df_all, "ssGSEA")

lines <- c(
  "===== Imprint scoring method comparison (z-score vs ssGSEA) =====",
  "",
  "[z-score mean (original method)]",
  paste0("  Between-group gradient: Healthy ", round(sum_z["Healthy_med"],3),
         " < COPD ", round(sum_z["COPD_med"],3),
         " < Sepsis ", round(sum_z["Sepsis_med"],3)),
  paste0("  Within-group dataset range: Healthy ", round(sum_z["Healthy_ds_range"],3),
         " | Sepsis ", round(sum_z["Sepsis_ds_range"],3)),
  "",
  "[ssGSEA (rank-based)]",
  paste0("  Between-group gradient: Healthy ", round(sum_s["Healthy_med"],3),
         " < COPD ", round(sum_s["COPD_med"],3),
         " < Sepsis ", round(sum_s["Sepsis_med"],3)),
  paste0("  Within-group dataset range: Healthy ", round(sum_s["Healthy_ds_range"],3),
         " | Sepsis ", round(sum_s["Sepsis_ds_range"],3)),
  "",
  "Interpretation:",
  "  The within-group dataset range reflects residual batch effects. If the ssGSEA within-group range is clearly < the z-score one,",
  "  rank-based scoring has removed part of the dataset baseline differences, so ssGSEA is more credible;",
  "  if the two within-group ranges are similar, batch confounding is inherent to the design and cannot be fixed by any scoring method."
)
writeLines(lines, file.path(out_dir, "path33_method_comparison.txt"))
cat("\n", paste(lines, collapse = "\n"), "\n", sep = "")

## ---- ssGSEA facet plot ----
ds_levels <- sort(unique(df_s$dataset))
ds_cols <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02")
if (length(ds_levels) > 6) ds_cols <- colorRampPalette(ds_cols)(length(ds_levels))
p <- ggplot(df_s, aes(x = dataset, y = imprint, color = dataset)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  facet_wrap(~ group, scales = "free_x") +
  scale_color_manual(values = ds_cols) +
  labs(x = NULL, y = "ssGSEA enrichment score",
       title = "ssGSEA imprint by dataset within each group") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(size = 11, face = "bold"))
ggsave(file.path(out_dir, "path33_imprint_ssgsea.pdf"), p, width = 8, height = 4)

cat("\nDone. Output:\n")
cat("  path33_imprint_two_methods.csv   (long table for both methods)\n")
cat("  path33_method_comparison.txt     (gradient + within-group range comparison for both methods)\n")
cat("  path33_imprint_ssgsea.pdf        (ssGSEA facet plot)\n")
