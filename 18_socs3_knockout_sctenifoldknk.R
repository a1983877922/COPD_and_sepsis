# =========================================================================
# 18_socs3_knockout_sctenifoldknk.R (SOCS3 virtual knockout with scTenifoldKnk)
#
# Use scTenifoldKnk to perform an in-silico virtual knockout (vKO) of SOCS3:
#   build a single-cell gene regulatory network (scGRN) from all pDC single-cell
#   expression, set SOCS3's outdegree edges to 0 to simulate knockout, and
#   manifold-align WT vs KO networks to find significantly perturbed genes
#   (diffRegulation).
#   (Use all pDC rather than only Healthy: SOCS3 expression in Healthy is
#   extremely sparse ~2.5%, while the Sepsis group has high SOCS3 expression;
#   using all pDC ensures SOCS3 has sufficient signal in the GRN.)
#
# Biological purpose (reverse interpretation of 'overexpression'):
#   In sepsis pDC, SOCS3 is up-regulated (log2FC=+3.29, the strongest in the
#   whole dataset) = a natural 'overexpression' state. If, after SOCS3 KO,
#   ISG/STAT target genes are significantly 'de-repressed' (FC>1), this proves
#   SOCS3 is a repressive regulator of the IFN program
#   => SOCS3 overexpression = suppressed ISG = pDC immunoparalysis.
#
# Why use scTenifoldKnk instead of CellOracle:
#   CellOracle's GRN comes from TF-motif scanning (TF-centric) and can only
#   perturb TFs; SOCS3 is not a transcription factor and cannot be modeled.
#   scTenifoldKnk's GRN is built from PC regression + tensor decomposition,
#   yielding a 'gene-gene' network that does not distinguish TFs, allowing
#   direct knockout of any gene (including SOCS3), and requires no base GRN.
#
# Input: path11_pdc_subset.rds (pDC object saved by script 11, includes group column)
# Output: path18_* series of results
# =========================================================================

.this_file <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
if (.this_file == "" || !file.exists(.this_file)) {
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1])
}
SCRIPT_DIR  <- dirname(.this_file)
config_file <- file.path(SCRIPT_DIR, "00_config.R")
if (!file.exists(config_file)) stop("00_config.R not found: ", config_file)
source(config_file)

cat("\n==============================================================\n")
cat("SOCS3 virtual knockout (scTenifoldKnk) - reverse evidence that SOCS3 overexpression -> immunoparalysis\n")
cat("==============================================================\n")

## =========================================================================
## Step 0: Install dependencies (scTenifoldNet + scTenifoldKnk, both from GitHub)
## =========================================================================
cat("\n===== Step 0: Check/install dependencies =====\n")

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("scTenifoldNet", quietly = TRUE)) {
  cat("scTenifoldNet not detected, installing from GitHub...\n")
  remotes::install_github("cailab-tamu/scTenifoldNet", upgrade = "never")
}
if (!requireNamespace("scTenifoldNet", quietly = TRUE)) {
  stop("scTenifoldNet installation failed")
}
if (!requireNamespace("scTenifoldKnk", quietly = TRUE)) {
  cat("scTenifoldKnk not detected, installing from GitHub...\n")
  remotes::install_github("cailab-tamu/scTenifoldKnk", upgrade = "never")
}
if (!requireNamespace("scTenifoldKnk", quietly = TRUE)) {
  stop("scTenifoldKnk installation failed")
}

suppressPackageStartupMessages(library(scTenifoldKnk))
cat("scTenifoldKnk loaded\n")

## =========================================================================
## Step 1: Read pDC object, extract Healthy (WT) group counts
## =========================================================================
cat("\n===== Step 1: Extract Healthy-group pDC counts =====\n")

pdc_file <- file.path(out_dir, "path11_pdc_subset.rds")
if (!file.exists(pdc_file)) stop(pdc_file, " not found, please run script 11 first")
pdc <- readRDS(pdc_file)
cat("pDC object:", ncol(pdc), "cells\n")
print(table(pdc$group))

# Build the GRN from all pDC (rather than only Healthy).
# Reason: SOCS3 expression in Healthy pDC is extremely sparse (only ~2.5%);
# if only Healthy is used, SOCS3's regulatory edges in the GRN are too weak and
# the KO effect cannot be detected. In all pDC, the Sepsis group has high SOCS3
# expression (log2FC=+3.29), which ensures SOCS3 has sufficient signal in the GRN.
#       scTenifoldKnk's "WT" means "observational samples without experimental
#       knockout", which all pDC satisfy.
#
# Highly-variable-gene filtering (critical): scTenifoldNet's tensor
# decomposition builds a nGenes×nGenes×nNet dense array; if all ~25000 genes are
# used, the tensor is 25000×25000×10 ≈ 50GB dense, and CP decomposition's
# unfold/kronecker runs out of memory -> errors 'dims contain missing values'.
# So filter to top 3000 HVGs + forcibly keep SOCS3/ISG, reducing the tensor to
# ~3000×3000×10.
pdc <- FindVariableFeatures(pdc, nfeatures = 3000, verbose = FALSE)
hv_genes <- VariableFeatures(pdc)

# Key genes of the SOCS3-STAT-IFN axis, forcibly kept (even if not in top 3000 HVGs)
isg_genes <- c("SOCS3", "SOCS1", "STAT1", "STAT2", "IRF7", "IRF9",
               "ISG15", "ISG20", "MX1", "MX2", "OAS1", "OAS2", "OAS3",
               "IFIT1", "IFIT2", "IFIT3", "IFI6", "IFI44", "IFI44L",
               "RSAD2", "USP18", "GBP1", "PSMB8", "UBE2L6", "IFIH1", "DDX58")
keep_genes <- unique(c(hv_genes, isg_genes[isg_genes %in% rownames(pdc)]))
cat("HVG count:", length(hv_genes),
    "| ISG forcibly kept:", sum(isg_genes %in% rownames(pdc)),
    "| total kept genes:", length(keep_genes), "\n")

# Extract raw counts (keep only filtered genes)
counts <- GetAssayData(pdc, assay = "RNA", layer = "counts")
counts <- counts[keep_genes, , drop = FALSE]
counts <- as(counts, "CsparseMatrix")  # convert to dgCMatrix to save memory
rm(pdc); gc()
cat("counts matrix (all pDC):", nrow(counts), "genes ×", ncol(counts), "cells\n")

## =========================================================================
## Step 2: Check SOCS3 presence + manual QC (version-independent, not relying on
##         scTenifoldKnk's QC parameters)
## =========================================================================
cat("\n===== Step 2: SOCS3 expression check + manual QC =====\n")

gKO <- "SOCS3"
if (!gKO %in% rownames(counts)) {
  stop(gKO, " not in counts, cannot knockout")
}
socs3_pct <- mean(counts[gKO, ] > 0)
cat("SOCS3 expression fraction in all pDC:", round(socs3_pct * 100, 2), "%\n")

# Manual QC (does not depend on scTenifoldKnk version-specific qc_* parameter names):
#   1) Filter out cells with library size < 500
#   2) Filter out cells with mitochondrial ratio > 10% (gene names starting with MT-)
lib_size <- Matrix::colSums(counts)
mt_idx   <- grep("^MT-", rownames(counts), ignore.case = TRUE)
mt_ratio <- if (length(mt_idx) > 0) {
  Matrix::colSums(counts[mt_idx, , drop = FALSE]) / lib_size
} else {
  rep(0, ncol(counts))
}
keep_cell <- lib_size >= 500 & mt_ratio <= 0.1
counts <- counts[, keep_cell, drop = FALSE]
cat("Cell count after manual QC:", ncol(counts), "(filtered", sum(!keep_cell), "cells)\n")

# Re-confirm SOCS3 is still present after QC
if (!gKO %in% rownames(counts)) {
  stop(gKO, " lost after QC")
}

## =========================================================================
## Step 3: Run scTenifoldKnk - knockout SOCS3 (pass only core parameters for version compatibility)
## =========================================================================
cat("\n===== Step 3: Run virtual knockout (SOCS3) =====\n")

nCores <- min(16, parallel::detectCores())
cat("Number of cores used:", nCores, "\n")
cat("Parameters: nc_nNet=10, nc_nCells=500, nc_nComp=3, td_K=3\n")
cat("(runtime depends on genes × cells; tensor decomposition may take several to ~ten minutes, please be patient)\n")

# Only pass core parameters present in all versions; do not pass qc/qc_minLibSize/
# qc_minPCT etc., to avoid 'unused arguments' errors due to the scTenifoldKnk
# version installed on the server.
ko_result <- scTenifoldKnk(
  countMatrix     = counts,
  gKO             = gKO,
  nc_nNet         = 10,
  nc_nCells       = 500,
  nc_nComp        = 3,
  nc_scaleScores  = TRUE,
  td_K            = 3,
  ma_nDim         = 2,
  nCores          = nCores
)

cat("\nSOCS3 virtual knockout complete\n")

## =========================================================================
## Step 4: Parse diffRegulation results
## =========================================================================
cat("\n===== Step 4: Parse results =====\n")

dr <- ko_result$diffRegulation
cat("diffRegulation table dimension:", nrow(dr), "genes\n")
cat("Column names:", paste(colnames(dr), collapse = ", "), "\n")

# Save full results
write.csv(dr, file.path(out_dir, "path18_socs3ko_diffRegulation.csv"),
          row.names = FALSE)

# Significantly perturbed genes (p.adj < 0.05)
sig <- dr[!is.na(dr$p.adj) & dr$p.adj < 0.05, ]
cat("\nNumber of significantly perturbed genes (FDR<0.05):", nrow(sig), "\n")
cat("Top 30 significantly perturbed genes:\n")
if (nrow(sig) > 0) {
  print(head(sig[order(sig$p.adj), ], 30))
}

## =========================================================================
## Step 5: Focus on whether ISG / STAT target genes are significantly perturbed
##         (core biological interpretation)
## =========================================================================
cat("\n===== Step 5: ISG / STAT target gene perturbation analysis =====\n")

# Key genes of the SOCS3-STAT-IFN axis (ISG + STAT target genes)
isg_genes <- c("SOCS3", "SOCS1", "STAT1", "STAT2", "IRF7", "IRF9",
               "ISG15", "ISG20", "MX1", "MX2", "OAS1", "OAS2", "OAS3",
               "IFIT1", "IFIT2", "IFIT3", "IFI6", "IFI44", "IFI44L",
               "RSAD2", "USP18", "GBP1", "PSMB8", "UBE2L6", "IFIH1", "DDX58")

isg_in <- isg_genes[isg_genes %in% dr$gene]
cat("Number of ISG/STAT target genes present in diffRegulation:", length(isg_in), "\n\n")

isg_tab <- dr[match(isg_in, dr$gene), ]
isg_tab <- isg_tab[order(isg_tab$p.adj), ]
# FC>1 = after KO the gene's regulatory relationship is 'released' (up-regulated), supporting SOCS3 as a repressor
isg_tab$direction <- ifelse(isg_tab$FC > 1, "up-regulated (de-repressed)", "down-regulated")
print(isg_tab)

write.csv(isg_tab, file.path(out_dir, "path18_socs3ko_ISG_summary.csv"),
          row.names = FALSE)

# Count ISG genes that are significant (FDR<0.05) and FC>1 (de-repressed)
isg_sig_up <- sum(!is.na(isg_tab$p.adj) & isg_tab$p.adj < 0.05 & isg_tab$FC > 1)
isg_sig    <- sum(!is.na(isg_tab$p.adj) & isg_tab$p.adj < 0.05)
cat("\nISG/STAT target genes significantly perturbed:", isg_sig, "/", nrow(isg_tab), "\n")
cat("Of which FC>1 (de-repressed/up-regulated):", isg_sig_up, "\n")

# Write statistics summary
summ <- c(
  paste0("Knockout gene: ", gKO),
  paste0("Input cell count (all pDC, after QC): ", ncol(counts)),
  paste0("SOCS3 expression fraction in all pDC: ", round(socs3_pct*100,2), "%"),
  paste0("Total diffRegulation genes: ", nrow(dr)),
  paste0("Significantly perturbed genes (FDR<0.05): ", nrow(sig)),
  paste0("ISG/STAT target gene count: ", nrow(isg_tab)),
  paste0("ISG/STAT significantly perturbed: ", isg_sig),
  paste0("ISG/STAT FC>1 (de-repressed): ", isg_sig_up),
  "",
  "Interpretation: if ISG/STAT target genes are significant after SOCS3 KO with FC>1 (de-repressed),",
  "      then SOCS3 is a repressive regulator of the IFN program;",
  "      sepsis pDC's SOCS3 up-regulation (overexpression) => suppresses ISG => immunoparalysis."
)
writeLines(summ, file.path(out_dir, "path18_socs3ko_stats.txt"))
cat("\nStatistics summary saved: path18_socs3ko_stats.txt\n")

## =========================================================================
## Step 6: Visualization (plotKO draws the KO-gene-centered subnetwork)
## =========================================================================
cat("\n===== Step 6: Visualization =====\n")

# Ensure ggplot2 is available (for volcano plot; plotKO base plot handled by tryCatch fallback)
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(ggplot2))

tryCatch({
  p <- plotKO(ko_result, gKO = gKO)
  ggsave(file.path(out_dir, "path18_socs3ko_network.pdf"), p,
         width = 8, height = 7)
  cat("KO-centered subnetwork plot saved: path18_socs3ko_network.pdf\n")
}, error = function(e) {
  cat("plotKO plotting failed (skipped):", conditionMessage(e), "\n")
})

# Supplement: volcano-style scatter of diffRegulation (FC vs -log10 p.adj)
tryCatch({
  dr$log10_padj <- -log10(dr$p.adj + 1e-300)
  dr$sig <- ifelse(!is.na(dr$p.adj) & dr$p.adj < 0.05, "sig", "ns")
  dr$is_isg <- ifelse(dr$gene %in% isg_genes, "ISG", "other")
  p2 <- ggplot(dr, aes(x = log2(FC + 1e-6), y = log10_padj,
                       color = sig)) +
    geom_point(size = 0.6, alpha = 0.6) +
    scale_color_manual(values = c("sig" = "red", "ns" = "grey70")) +
    geom_point(data = dr[dr$is_isg == "ISG", ], aes(x = log2(FC + 1e-6),
                 y = log10_padj), color = "blue", size = 1.2) +
    geom_text(data = dr[dr$gene %in% c("SOCS3","SOCS1","STAT1","STAT2",
                                        "IRF7","ISG15","MX1","OAS2") &
                          dr$sig == "sig", ],
              aes(label = gene), size = 2.5, vjust = -0.5, color = "blue") +
    labs(x = "log2(FC) after SOCS3 KO",
         y = "-log10(adjusted p)",
         title = "SOCS3 virtual knockout: differentially regulated genes") +
    theme_minimal()
  ggsave(file.path(out_dir, "path18_socs3ko_volcano.pdf"), p2,
         width = 7, height = 6)
  cat("Volcano plot saved: path18_socs3ko_volcano.pdf\n")
}, error = function(e) {
  cat("Volcano plot drawing failed (skipped):", conditionMessage(e), "\n")
})

cat("\nSOCS3 virtual knockout analysis complete\n")
cat("Outputs: path18_socs3ko_diffRegulation.csv / path18_socs3ko_ISG_summary.csv\n")
cat("      path18_socs3ko_stats.txt / path18_socs3ko_network.pdf / path18_socs3ko_volcano.pdf\n")
