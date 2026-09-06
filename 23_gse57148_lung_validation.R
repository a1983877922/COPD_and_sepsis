# =========================================================================
# 23 - GSE57148 lung-side validation
#
# Purpose: use independent bulk lung tissue RNA-seq (GSE57148) to validate the single-cell
#       finding of "upregulated myeloid / IFN / pDC signals on the COPD lung side".
#       -- On the lung side we currently have only the single-cell GSE136831 source,
#       lacking independent bulk validation.
#
# Data: GSE57148 (lung tissue RNA-seq, GPL11154/Illumina HiSeq)
#   98 COPD + 91 smoking controls with normal lung function (all male smokers, surgical lung resection specimens)
#   Preprocessed FPKM matrix: GSE57148_COPD_FPKM_Normalized.txt.gz (16739 genes x 189 samples)
#   Clinical data include FEV1% / FEV1-FVC / DLCO / pack-years (correlated if present in the series matrix)
#
# Output: path23_* series
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

suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(ggplot2))

cat("\n==============================================================\n")
cat("GSE57148 lung-side bulk validation (COPD vs smoking controls)\n")
cat("==============================================================\n")

## Gene signatures
# 149 shared myeloid genes (read the output of script 21; fall back to path2 from script 02 if absent)
shared_file <- file.path(out_dir, "path21_shared_up_genes.txt")
if (!file.exists(shared_file)) shared_file <- file.path(out_dir, "path2_shared_myeloid_genes.txt")
shared149 <- if (file.exists(shared_file)) {
  trimws(readLines(shared_file))
} else {
  cat("!! Cannot find the 149 shared gene file, using an empty set (skip the 149 signature)\n")
  character(0)
}
cat("149 shared gene count:", length(shared149), "\n")

ifn_genes <- c("ISG15","MX1","MX2","OAS1","OAS2","OAS3","IFI44L","IFI6",
               "IFIT1","IFIT2","IFIT3","GBP1","STAT1","IRF7","ISG20","USP18")
pdc_markers <- c("LILRA4","CLEC4C","IL3RA","TCF4","IRF7","IRF8","GZMB","SPIB")

score_geneset <- function(expr_mat, genes) {
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) < 2) return(rep(NA_real_, ncol(expr_mat)))
  sub <- expr_mat[genes, , drop = FALSE]
  z   <- t(scale(t(sub)))
  z[is.na(z)] <- 0
  colMeans(z)
}

## =========================================================================
## Step 1: read the FPKM matrix (local first, curl download as fallback)
## =========================================================================
cat("\n===== Step 1: Read FPKM matrix =====\n")

fpkm_file <- file.path(out_dir, "GSE57148_COPD_FPKM_Normalized.txt.gz")
fpkm_url  <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE57nnn/GSE57148/suppl/GSE57148_COPD_FPKM_Normalized.txt.gz"
if (!file.exists(fpkm_file)) {
  cat("No local FPKM file, downloading...\n")
  ok <- system(paste0("curl -sk --max-time 600 -o \"", fpkm_file, "\" \"", fpkm_url, "\""))
  if (ok != 0 || !file.exists(fpkm_file)) stop("Failed to download FPKM, please download manually to ", fpkm_file)
}
fpkm <- read.delim(gzfile(fpkm_file), row.names = 1, check.names = FALSE)
cat("FPKM matrix:", nrow(fpkm), "genes x", ncol(fpkm), "samples\n")

# log2(FPKM + 1)
expr <- log2(fpkm + 1)
rm(fpkm); gc()

## =========================================================================
## Step 2: read the series matrix to extract clinical data (group + FEV1)
## =========================================================================
cat("\n===== Step 2: Read clinical information =====\n")

matrix_file <- file.path(out_dir, "GSE57148_series_matrix.txt.gz")
matrix_url  <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE57nnn/GSE57148/matrix/GSE57148_series_matrix.txt.gz"
if (!file.exists(matrix_file)) {
  cat("No local series matrix, downloading...\n")
  system(paste0("curl -sk --max-time 300 -o \"", matrix_file, "\" \"", matrix_url, "\""))
}
gse <- getGEO(filename = matrix_file, getGPL = FALSE)
pd  <- pData(gse)
cat("pData column names:\n  ", paste(colnames(pd), collapse = "\n   "), "\n")

# Group assignment: COPD vs control
# Prefer searching for the "copd" keyword in title / source_name / characteristics
find_grp <- function(pd) {
  candidates <- c("title", "source_name_ch1")
  for (cn in intersect(candidates, colnames(pd))) {
    v <- as.character(pd[[cn]])
    if (any(grepl("copd", v, ignore.case = TRUE))) return(v)
  }
  for (cn in colnames(pd)) {
    v <- as.character(pd[[cn]])
    if (any(grepl("copd", v, ignore.case = TRUE))) return(v)
  }
  NULL
}
grp_raw <- find_grp(pd)
if (is.null(grp_raw)) stop("Cannot identify COPD groups, please check the pData column names above")
is_copd <- grepl("copd", grp_raw, ignore.case = TRUE)
cat("Groups: COPD =", sum(is_copd), " / control =", sum(!is_copd), "\n")

# FEV1 extraction (if present in the series matrix)
fev1 <- NULL
fev_col <- grep("fev1|fev_1|fev", colnames(pd), ignore.case = TRUE, value = TRUE)
fev_col <- fev_col[!grepl("fvc|ratio", fev_col, ignore.case = TRUE)]
if (length(fev_col) > 0) {
  fev1 <- suppressWarnings(as.numeric(sub("^[^0-9.]*([0-9.]+).*$", "\\1",
                                          as.character(pd[[fev_col[1]]]))))
  cat("FEV1 field:", fev_col[1], " (valid values:", sum(!is.na(fev1)), ")\n")
} else {
  cat("No numeric FEV1 field in the series matrix, skipping FEV1 correlation analysis\n")
}

# Sample name alignment: FPKM column names are sample aliases ("1017-NOR"/"578-COPD"), pData rownames are GSM IDs.
# Extract aliases from title ("Normal lung tissue 1017-NOR" / "COPD lung tissue 578-COPD") to align.
title_raw <- as.character(pd$title)
sample_alias <- vapply(strsplit(title_raw, " "), function(x) x[length(x)], character(1))
sample_alias <- gsub("^\"|\"$", "", sample_alias)

common <- intersect(colnames(expr), sample_alias)
cat("Sample alignment: shared samples between FPKM and title aliases:", length(common), "\n")
if (length(common) == 0) {
  cat("Alignment failed, printing the first 5 FPKM column names and the first 5 title aliases for checking:\n")
  cat("  FPKM:", paste(head(colnames(expr), 5), collapse = ", "), "\n")
  cat("  Aliases:", paste(head(sample_alias, 5), collapse = ", "), "\n")
  stop("Sample names cannot be aligned")
}
expr <- expr[, common, drop = FALSE]
m <- match(common, sample_alias)
is_copd <- is_copd[m]
if (!is.null(fev1)) fev1 <- fev1[m]
names(is_copd) <- common

## =========================================================================
## Step 3: signature scoring
## =========================================================================
cat("\n===== Step 3: Signature scoring =====\n")

score149 <- score_geneset(expr, shared149)
score_ifn <- score_geneset(expr, ifn_genes)
score_pdc <- score_geneset(expr, pdc_markers)
cat("149 signature genes matched:", sum(shared149 %in% rownames(expr)), "/", length(shared149), "\n")
cat("IFN signature genes matched:", sum(ifn_genes %in% rownames(expr)), "/", length(ifn_genes), "\n")
cat("pDC marker genes matched:", sum(pdc_markers %in% rownames(expr)), "/", length(pdc_markers), "\n")

## =========================================================================
## Step 4: COPD vs control differences (Wilcoxon)
## =========================================================================
cat("\n===== Step 4: COPD vs control signature differences =====\n")

res_lines <- c("=== GSE57148 lung-side bulk validation ===", "",
               paste0("Samples: COPD = ", sum(is_copd), " / control = ", sum(!is_copd)), "")

for (nm in c("149 shared myeloid","IFN","pDC markers")) {
  sc <- switch(nm, "149 shared myeloid" = score149, "IFN" = score_ifn, "pDC markers" = score_pdc)
  if (all(is.na(sc))) { next }
  w <- wilcox.test(sc[is_copd], sc[!is_copd])
  med_c <- median(sc[is_copd], na.rm = TRUE); med_n <- median(sc[!is_copd], na.rm = TRUE)
  line <- sprintf("%s: COPD median=%.3f vs control=%.3f, Wilcoxon p=%.3e",
                  nm, med_c, med_n, w$p.value)
  cat(" ", line, "\n"); res_lines <- c(res_lines, line)
}

## =========================================================================
## Step 5: signature vs FEV1 correlation (if available)
## =========================================================================
if (!is.null(fev1) && sum(!is.na(fev1)) > 30) {
  cat("\n===== Step 5: Signature vs FEV1 correlation =====\n")
  res_lines <- c(res_lines, "", "Signature vs FEV1% correlation (Spearman):")
  for (nm in c("149 shared myeloid","IFN","pDC markers")) {
    sc <- switch(nm, "149 shared myeloid" = score149, "IFN" = score_ifn, "pDC markers" = score_pdc)
    if (all(is.na(sc))) next
    ok <- !is.na(fev1) & !is.na(sc)
    if (sum(ok) < 30) next
    sp <- cor.test(sc[ok], fev1[ok], method = "spearman")
    line <- sprintf("%s: rho=%.3f, p=%.3e", nm, sp$estimate, sp$p.value)
    cat(" ", line, "\n"); res_lines <- c(res_lines, line)
  }
}

writeLines(res_lines, file.path(out_dir, "path23_gse57148_lung_stats.txt"))

## =========================================================================
## Step 6: visualization
## =========================================================================
cat("\n===== Step 6: Visualization =====\n")

df <- data.frame(
  group = ifelse(is_copd, "COPD", "Control"),
  shared149 = score149, IFN = score_ifn, pDC = score_pdc
)
df$group <- factor(df$group, levels = c("Control", "COPD"))

# Boxplots (three signatures, COPD vs control)
df_long <- reshape(df, direction = "long", varying = c("shared149","IFN","pDC"),
                   v.names = "score", timevar = "sig", times = c("149 shared myeloid","IFN","pDC"))
df_long$sig <- factor(df_long$sig, levels = c("149 shared myeloid","IFN","pDC"))
p1 <- ggplot(df_long, aes(x = group, y = score, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
  facet_wrap(~ sig, scales = "free_y") +
  scale_fill_manual(values = c("Control" = "#90A4AE", "COPD" = "#E53935")) +
  labs(x = NULL, y = "Signature score (z)", title = "COPD lung signatures (GSE57148)") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "path23_gse57148_lung_boxplot.pdf"), p1, width = 8, height = 4)
ggsave(file.path(out_dir, "path23_gse57148_lung_boxplot.png"), p1, width = 8, height = 4, dpi = 300)

# FEV1 correlation scatter plot (IFN signature, if available)
if (!is.null(fev1) && sum(!is.na(fev1)) > 30) {
  df2 <- data.frame(FEV1 = fev1, IFN = score_ifn, shared149 = score149)
  df2 <- df2[complete.cases(df2), ]
  p2 <- ggplot(df2, aes(x = FEV1, y = IFN)) +
    geom_point(alpha = 0.5, color = "#1976D2") + geom_smooth(method = "lm", se = TRUE, color = "black") +
    labs(x = "FEV1 % predicted", y = "IFN signature score",
         title = "IFN signature vs lung function (GSE57148)") +
    theme_minimal(base_size = 12)
  ggsave(file.path(out_dir, "path23_gse57148_fev1_cor.pdf"), p2, width = 6, height = 5)
  ggsave(file.path(out_dir, "path23_gse57148_fev1_cor.png"), p2, width = 6, height = 5, dpi = 300)
}

cat("\nGSE57148 lung-side validation finished\n")
cat("Output: path23_gse57148_lung_stats.txt / _boxplot.pdf/.png / (optional) _fev1_cor.pdf/.png\n")
cat("Interpretation: if the 149-shared / IFN / pDC signatures are significantly higher in COPD lung tissue than in controls,\n")
cat("      it supports the single-cell finding of "upregulated myeloid/IFN programs on the COPD lung side";\n")
cat("      a negative correlation between IFN and FEV1 suggests interferon signaling increases as lung function declines.\n")
