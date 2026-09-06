# =========================================================================
# 08_bulk_validation.R
# Cross-validate using the bulk PBMC cohort GSE248493 (26 COPD + 13 healthy, 3' QuantSeq RNA-seq):
#   score the single-cell-derived DC marker gene set in bulk expression to see whether
#   COPD patients' peripheral blood DC signal is lower than healthy controls —
#   cross-validating the "blood-side DC exhaustion" conclusion.
#
# Key: GSE248493 is RNA-seq; the Series Matrix only contains sample info, no expression values
#   (exprs() returns 0 rows). Expression data is in the NCBI-officially-generated raw counts matrix:
#     GSE248493_raw_counts_GRCh38.p13_NCBI.tsv.gz   (first column GeneID/ENTREZ ID)
#   Download requires the type=rnaseq_counts dedicated endpoint; GeneID->Symbol mapping uses org.Hs.eg.db.
#
# Dependencies: GEOquery (sample grouping), data.table (read matrix, loaded in 00_config),
#       org.Hs.eg.db (ENTREZ->Symbol mapping)
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
  stop("GEOquery is not installed. Please run first:\n",
       "  BiocManager::install('GEOquery')\n",
       "then re-source this script.")
}
suppressPackageStartupMessages(library(GEOquery))

cat("\n==============================================================\n")
cat("Bulk validation (GSE248493: DC signal in COPD vs healthy)\n")
cat("==============================================================\n")

## DC marker gene sets (from single-cell annotation; pan-DC + pDC/cDC1/cDC2 specific)
dc_pan  <- c("FCER1A","CLEC9A","CD1C","CLEC10A","IRF7","IRF8","ITGAX")
pdc     <- c("LILRA4","CLEC4C","IL3RA","TCF4","IRF7","GZMB")
c_dc1   <- c("CLEC9A","XCR1","BATF3","CADM1")
c_dc2   <- c("CD1C","FCER1A","CLEC10A","SIRPA")

## =========================================================================
## Step 1: download sample metadata (Series Matrix) to get grouping
## =========================================================================
cat("\n===== Step 1: download GSE248493 sample grouping =====\n")

gse <- getGEO("GSE248493", GSEMatrix = TRUE, getGPL = FALSE)
gse <- gse[[1]]
pd  <- pData(gse)
cat("Number of samples:", nrow(pd), "\n")
cat("pData column names:", paste(colnames(pd), collapse = ", "), "\n")

# Grouping: titles look like "Healthy control [DK20193B 09]" vs "COPD [DK21067B 06]"
# judging by title is most reliable (the disease field is the full name "Chronic Obstructive Pulmonary
# Disease", which does not contain the substring "COPD" and would be misjudged)
group <- ifelse(grepl("Healthy", pd$title, ignore.case = TRUE), "Control", "COPD")
names(group) <- rownames(pd)   # GSM accession
cat("Grouping (by title):\n"); print(table(group))

## =========================================================================
## Step 2: download NCBI raw counts matrix + map gene symbol
## =========================================================================
cat("\n===== Step 2: download raw counts matrix + map symbol =====\n")

# raw counts is an NCBI-generated matrix, requires the type=rnaseq_counts dedicated endpoint
# (the ordinary format=file&file=... would 404)
counts_url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts",
                     "&acc=GSE248493&format=file",
                     "&file=GSE248493_raw_counts_GRCh38.p13_NCBI.tsv.gz")

counts_file <- file.path(out_dir, "GSE248493_raw_counts.tsv.gz")
if (!file.exists(counts_file)) {
  cat("Downloading raw counts matrix ...\n")
  download.file(counts_url, counts_file, mode = "wb", quiet = TRUE)
}
cat("counts file:", counts_file, "\n")

cnt <- as.data.frame(data.table::fread(counts_file))
cat("raw counts dimensions:", nrow(cnt), "genes x", ncol(cnt) - 1, "samples\n")

# Determine first column: numbers = GeneID (ENTREZ), letters = symbol
gid        <- cnt[[1]]
gid_name   <- colnames(cnt)[1]
is_symbol  <- any(grepl("[A-Za-z]", as.character(gid)))
cat("First column name:", gid_name, "; is symbol:", is_symbol, "\n")

mat <- as.matrix(cnt[, -1, drop = FALSE])
sample_ids <- colnames(mat)

if (!is_symbol) {
  # ENTREZ GeneID -> Symbol mapping (use org.Hs.eg.db, no annotation table download needed)
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("org.Hs.eg.db is required for ENTREZ->Symbol mapping, please run first:\n",
         "  BiocManager::install('org.Hs.eg.db')")
  }
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  sym <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = as.character(gid),
                               column = "SYMBOL", keytype = "ENTREZID",
                               multiVals = "first")
  sym <- unname(sym)
  sym[is.na(sym)] <- ""
  cat("ENTREZ -> Symbol mapping done; empty symbols:", sum(sym == ""), "\n")
} else {
  sym <- as.character(gid)
}

# Deduplicate: keep non-empty symbols on first occurrence
keep  <- sym != "" & !duplicated(sym)
mat   <- mat[keep, , drop = FALSE]
rownames(mat) <- sym[keep]
cat("Matrix after mapping:", nrow(mat), "genes x", ncol(mat), "samples\n")

# CPM + log2 normalization (remove library size differences, then per-gene z-score)
libsize  <- colSums(mat)
cpm      <- sweep(mat, 2, libsize / 1e6, "/")
expr_log <- log2(cpm + 1)

# Sample alignment (counts column names GSM accession vs grouping names)
common <- intersect(sample_ids, names(group))
expr_log <- expr_log[, common, drop = FALSE]
group    <- group[common]
cat("Number of aligned samples:", length(common), "\n")
print(table(group))

## =========================================================================
## Step 3: DC marker gene set scoring (z-score average)
## =========================================================================
cat("\n===== Step 3: DC signal scoring =====\n")

score_geneset <- function(expr_mat, genes) {
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) < 2) return(rep(NA_real_, ncol(expr_mat)))
  sub <- expr_mat[genes, , drop = FALSE]
  z   <- t(scale(t(sub)))     # per-gene z-score across samples
  z[is.na(z)] <- 0
  colMeans(z)
}

cat("DC marker hit counts:",
    "pan=", sum(dc_pan %in% rownames(expr_log)),
    " pDC=", sum(pdc %in% rownames(expr_log)),
    " cDC1=", sum(c_dc1 %in% rownames(expr_log)),
    " cDC2=", sum(c_dc2 %in% rownames(expr_log)), "\n")

score_dc   <- score_geneset(expr_log, dc_pan)
score_pdc  <- score_geneset(expr_log, pdc)
score_cdc1 <- score_geneset(expr_log, c_dc1)
score_cdc2 <- score_geneset(expr_log, c_dc2)

res <- data.frame(sample = colnames(expr_log), group = group,
                  DC = score_dc, pDC = score_pdc,
                  cDC1 = score_cdc1, cDC2 = score_cdc2,
                  stringsAsFactors = FALSE)
write.csv(res, file.path(out_dir, "path8_bulk_dc_scores.csv"), row.names = FALSE)

## =========================================================================
## Step 4: COPD vs healthy statistical test
## =========================================================================
cat("\n===== Step 4: statistical test (COPD vs Control) =====\n")

stats_lines <- character()
for (gs in c("DC", "pDC", "cDC1", "cDC2")) {
  v_copd <- res[[gs]][res$group == "COPD"]
  v_ctrl <- res[[gs]][res$group == "Control"]
  if (length(v_copd) < 2 || length(v_ctrl) < 2) { cat(" [", gs, "] insufficient samples, skipping\n"); next }
  tt <- wilcox.test(v_copd, v_ctrl)
  line <- sprintf("%s: COPD mean=%.3f vs Control mean=%.3f, wilcox p=%.3e",
                  gs, mean(v_copd, na.rm = TRUE), mean(v_ctrl, na.rm = TRUE), tt$p.value)
  cat(" ", line, "\n")
  stats_lines <- c(stats_lines, line)
}
writeLines(stats_lines, file.path(out_dir, "path8_bulk_dc_stats.txt"))

## Visualization: boxplot
library(ggplot2)
res_long <- reshape2::melt(res, id.vars = c("sample", "group"),
                           variable.name = "geneset", value.name = "score")
p <- ggplot(res_long, aes(x = group, y = score, fill = group)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  facet_wrap(~ geneset, scales = "free_y") +
  scale_fill_manual(values = c("Control" = "#4c8bf5", "COPD" = "#d64545")) +
  labs(title = "Bulk PBMC DC signature (GSE248493)", y = "gene-set z-score") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "path8_bulk_dc_boxplot.pdf"), p, width = 8, height = 6)

cat("\nBulk validation complete\n")
cat("Output: path8_bulk_dc_scores.csv / path8_bulk_dc_stats.txt / path8_bulk_dc_boxplot.pdf\n")
cat("Interpretation: if the pDC/DC score in the COPD group is significantly lower than Control, this supports the 'blood-side DC exhaustion' conclusion\n")
