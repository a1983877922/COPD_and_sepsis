# =========================================================================
# 10_sepsis_bulk_validation.R
# Cross-validate the "blood-side DC exhaustion" conclusion using the sepsis bulk cohort
# GSE66099 (47 healthy + 199 sepsis, Sweeney 2015, GPL570 Affymetrix U133 Plus 2.0
# microarray, whole blood).
#
# Note: previously misused GSE65682 — that was the MARS consortium's 802 ICU patients
# (pneumonia grouping, no healthy controls), which cannot do sepsis vs healthy. Switched to
# GSE66099 (clearly healthy+sepsis+SIRS).
#
# Difference from 08 (GSE248493, RNA-seq):
#   - microarray: expression matrix (log2) is directly in the Series Matrix, no counts download needed;
#   - rownames are probe IDs, need to map to gene symbol via GPL570 fData annotation.
#
# Dependencies: GEOquery (download/parse)
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

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  stop("GEOquery is not installed. Please run first: BiocManager::install('GEOquery')")
}
suppressPackageStartupMessages(library(GEOquery))

cat("\n==============================================================\n")
cat("Sepsis bulk validation (GSE66099: DC signal in Sepsis vs healthy)\n")
cat("==============================================================\n")

## DC marker gene sets (consistent with single-cell annotation)
dc_pan  <- c("FCER1A","CLEC9A","CD1C","CLEC10A","IRF7","IRF8","ITGAX")
pdc     <- c("LILRA4","CLEC4C","IL3RA","TCF4","IRF7","GZMB")
c_dc1   <- c("CLEC9A","XCR1","BATF3","CADM1")
c_dc2   <- c("CD1C","FCER1A","CLEC10A","SIRPA")

## =========================================================================
## Step 1: download GSE66099 (microarray, Series Matrix contains expression + probe annotation)
## =========================================================================
cat("\n===== Step 1: download GSE66099 =====\n")

# Prefer a local series matrix (avoid 81MB download timeout), otherwise download online
matrix_file <- file.path(out_dir, "GSE66099_matrix.txt.gz")
if (file.exists(matrix_file)) {
  cat("Using local series matrix:", matrix_file, "\n")
  gse <- getGEO(filename = matrix_file, getGPL = TRUE)
} else {
  gse <- getGEO("GSE66099", GSEMatrix = TRUE, getGPL = TRUE)
  gse <- gse[[1]]
}
expr <- exprs(gse)
cat("Expression matrix dimensions:", nrow(expr), "probes x", ncol(expr), "samples\n")

pd <- pData(gse)
cat("pData column names:", paste(colnames(pd), collapse = ", "), "\n")

## =========================================================================
## Step 2: probe ID -> gene symbol mapping
## =========================================================================
cat("\n===== Step 2: probe -> symbol mapping =====\n")

fd <- fData(gse)
cat("fData column names:", paste(colnames(fd), collapse = ", "), "\n")

# Find the gene symbol column (Affymetrix common: "Gene Symbol"/"GeneSymbol"/"Symbol"/"gene_assignment")
sym_col <- NULL
for (cand in c("Gene Symbol", "GeneSymbol", "Symbol", "Gene.Symbol",
               "gene_symbol", "GENE_SYMBOL")) {
  if (cand %in% colnames(fd)) { sym_col <- cand; break }
}
if (is.null(sym_col)) {
  # Fallback: fuzzy match column names containing "symbol"
  m <- grep("symbol", colnames(fd), ignore.case = TRUE)
  if (length(m) > 0) sym_col <- colnames(fd)[m[1]]
}
if (is.null(sym_col)) stop("Cannot find gene symbol column in fData, please check column names and set sym_col manually")

cat("Using symbol column:", sym_col, "\n")
sym <- as.character(fd[[sym_col]])
names(sym) <- rownames(fd)   # probe ID
cat("Probe symbol examples:", paste(head(sym, 5), collapse = ", "), "\n")

# Keep only probes that have a symbol
keep <- sym != "" & !is.na(sym) & sym != "---"
expr <- expr[keep, , drop = FALSE]
sym  <- sym[keep]

# Multiple probes for the same symbol -> take mean (common microarray handling)
sym_uniq <- unique(sym)
expr_sym <- matrix(NA, nrow = length(sym_uniq), ncol = ncol(expr),
                   dimnames = list(sym_uniq, colnames(expr)))
for (g in sym_uniq) {
  idx <- which(sym == g)
  if (length(idx) == 1) expr_sym[g, ] <- expr[idx, ]
  else expr_sym[g, ] <- colMeans(expr[idx, , drop = FALSE])
}
expr <- expr_sym
cat("Matrix after mapping:", nrow(expr), "genes x", ncol(expr), "samples\n")

## =========================================================================
## Step 3: sample grouping (disease:ch1 -> Control / Sepsis / SepticShock / SIRS)
## =========================================================================
cat("\n===== Step 3: grouping =====\n")

# GSE66099 grouping field is disease:ch1, values: Control / Sepsis / SepticShock / SIRS
if (!"disease:ch1" %in% colnames(pd)) {
  stop("Cannot find disease:ch1 field, please check pData column names")
}
group <- as.character(pd[["disease:ch1"]])
names(group) <- rownames(pd)
cat("Raw grouping (disease:ch1):\n"); print(table(group, useNA = "ifany"))

# Merge SepticShock into Sepsis (sepsis spectrum)
group[group == "SepticShock"] <- "Sepsis"

# Alignment (series matrix column names = GSM accession, consistent with pd rownames)
common <- intersect(colnames(expr), names(group))
expr  <- expr[, common, drop = FALSE]
group <- group[common]

# Keep only Sepsis and Control (exclude SIRS/NA)
keep_samp <- group %in% c("Sepsis", "Control")
expr  <- expr[, keep_samp, drop = FALSE]
group <- group[keep_samp]
cat("After alignment + filtering (Sepsis vs Control):\n")
print(table(group))

## =========================================================================
## Step 4: DC marker gene set scoring (z-score average)
## =========================================================================
cat("\n===== Step 4: DC signal scoring =====\n")

score_geneset <- function(expr_mat, genes) {
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) < 2) return(rep(NA_real_, ncol(expr_mat)))
  sub <- expr_mat[genes, , drop = FALSE]
  z   <- t(scale(t(sub)))
  z[is.na(z)] <- 0
  colMeans(z)
}

cat("DC marker hit counts:",
    "pan=", sum(dc_pan %in% rownames(expr)),
    " pDC=", sum(pdc %in% rownames(expr)),
    " cDC1=", sum(c_dc1 %in% rownames(expr)),
    " cDC2=", sum(c_dc2 %in% rownames(expr)), "\n")

score_dc   <- score_geneset(expr, dc_pan)
score_pdc  <- score_geneset(expr, pdc)
score_cdc1 <- score_geneset(expr, c_dc1)
score_cdc2 <- score_geneset(expr, c_dc2)

res <- data.frame(sample = colnames(expr), group = group,
                  DC = score_dc, pDC = score_pdc,
                  cDC1 = score_cdc1, cDC2 = score_cdc2,
                  stringsAsFactors = FALSE)
write.csv(res, file.path(out_dir, "path10_sepsis_bulk_dc_scores.csv"), row.names = FALSE)

## =========================================================================
## Step 5: Sepsis vs healthy statistical test
## =========================================================================
cat("\n===== Step 5: statistical test (Sepsis vs Control) =====\n")

stats_lines <- character()
for (gs in c("DC", "pDC", "cDC1", "cDC2")) {
  v_sep <- res[[gs]][res$group == "Sepsis"]
  v_ctrl <- res[[gs]][res$group == "Control"]
  if (length(v_sep) < 2 || length(v_ctrl) < 2) { cat(" [", gs, "] insufficient samples, skipping\n"); next }
  tt <- wilcox.test(v_sep, v_ctrl)
  line <- sprintf("%s: Sepsis mean=%.3f vs Control mean=%.3f, wilcox p=%.3e",
                  gs, mean(v_sep, na.rm = TRUE), mean(v_ctrl, na.rm = TRUE), tt$p.value)
  cat(" ", line, "\n")
  stats_lines <- c(stats_lines, line)
}
writeLines(stats_lines, file.path(out_dir, "path10_sepsis_bulk_dc_stats.txt"))

## Visualization
library(ggplot2)
res_long <- reshape2::melt(res, id.vars = c("sample", "group"),
                           variable.name = "geneset", value.name = "score")
p <- ggplot(res_long, aes(x = group, y = score, fill = group)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.4) +
  facet_wrap(~ geneset, scales = "free_y") +
  scale_fill_manual(values = c("Control" = "#4c8bf5", "Sepsis" = "#d64545")) +
  labs(title = "Sepsis whole-blood DC signature (GSE66099)", y = "gene-set z-score") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "path10_sepsis_bulk_dc_boxplot.pdf"), p, width = 8, height = 6)

cat("\nSepsis bulk validation complete\n")
cat("Output: path10_sepsis_bulk_dc_scores.csv / path10_sepsis_bulk_dc_stats.txt / path10_sepsis_bulk_dc_boxplot.pdf\n")
cat("Interpretation: if the DC/pDC score in the Sepsis group is significantly lower than Control, this supports the 'blood-side DC exhaustion' conclusion\n")
