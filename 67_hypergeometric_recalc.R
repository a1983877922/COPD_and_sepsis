# =========================================================================
# 67_hypergeometric_background_recalc.R
# Response to reviewer comment C1: the hypergeometric test background for the
#   149 shared programs must not be hardcoded to 20000; recalculate P using
#   "the number of genes actually tested in common across the two tissues" as background.
# Method: read blood object + lung object, take rownames intersection as background N,
#   recalculate the hypergeometric P.
# Output: console print + path67_hypergeom_recalc.txt
# =========================================================================

.this_file <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
if (.this_file == "" || !file.exists(.this_file)) {
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1])
}
SCRIPT_DIR <- dirname(.this_file)
config_file <- file.path(SCRIPT_DIR, "00_config.R")
if (!file.exists(config_file)) stop("Cannot find 00_config.R: ", config_file)
source(config_file)

suppressPackageStartupMessages(library(Seurat))

cat("\n===== 67: Refining the hypergeometric background =====\n")

blood_file <- file.path(out_dir, "path1_sepsis_copd_integrated.rds")
lung_file  <- file.path(out_dir, "path2_copd_lung.rds")
if (!file.exists(blood_file)) stop("Cannot find ", blood_file)
if (!file.exists(lung_file))  stop("Cannot find ", lung_file)

blood <- readRDS(blood_file)
lung  <- readRDS(lung_file)
cat("Blood object gene count (rownames):", nrow(blood), "\n")
cat("Lung object gene count (rownames):", nrow(lung),  "\n")

## M2 cell count / donor count check
cat("\n--- M2 cell count / donor count ---\n")
cat("Blood object cell count (ncol):", ncol(blood), "\n")
cat("Lung object cell count (ncol):", ncol(lung),  "\n")
for (f in c("donor","donor_id","sample","orig.ident","patient","subject")) {
  if (f %in% colnames(blood@meta.data))
    cat("  Blood object field", f, "distinct:", length(unique(blood@meta.data[[f]])), "\n")
}
for (f in c("sample","donor","orig.ident")) {
  if (f %in% colnames(lung@meta.data))
    cat("  Lung object field", f, "distinct:", length(unique(lung@meta.data[[f]])), "\n")
}

bg <- length(intersect(rownames(blood), rownames(lung)))
cat("Number of genes tested in common across the two tissues N =", bg, "\n")

# DEG counts consistent with script 21
blood_up <- 1354
lung_up  <- 1270
shared   <- 149

res <- c(
  "===== 67 Refining the hypergeometric background =====",
  sprintf("Blood object gene count: %d", nrow(blood)),
  sprintf("Lung object gene count: %d", nrow(lung)),
  sprintf("Genes tested in common N: %d", bg),
  sprintf("blood_up=%d, lung_up=%d, shared_up=%d", blood_up, lung_up, shared)
)

for (N in bg) {
  exp_ov <- blood_up * lung_up / N
  p <- phyper(shared - 1, blood_up, N - blood_up, lung_up, lower.tail = FALSE)
  line <- sprintf("Background N=%d: expected overlap=%.1f, P(X>=%d)=%.4e (%s)",
                  N, exp_ov, shared, p,
                  ifelse(shared > exp_ov, "enriched", "depleted"))
  res <- c(res, line)
  cat(line, "\n")
}

# Robustness: also print N=15000/12000 for reference
for (N in c(20000, 15000, 12000)) {
  exp_ov <- blood_up * lung_up / N
  p <- phyper(shared - 1, blood_up, N - blood_up, lung_up, lower.tail = FALSE)
  cat(sprintf("  (reference) N=%d: expected=%.1f, P=%.3e", N, exp_ov, p), "\n")
}

writeLines(res, file.path(out_dir, "path67_hypergeom_recalc.txt"))
cat("\nOutput: path67_hypergeom_recalc.txt\n")
cat("===== Script 67 finished =====\n")
