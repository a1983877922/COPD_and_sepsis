#!/usr/bin/env Rscript
# =============================================================================
# Script 34: hdWGCNA figure generation
# Sepsis + COPD comorbidity analysis — supplementary hdWGCNA figures (Fig 3 panels)
# =============================================================================
# Purpose: script 13 only wrote the module gene CSV/TXT files and produced no
#       figures. This script reads path13_sc_hdwgcna.rds
#       (or re-runs hdWGCNA) and produces the figures needed for Fig 3:
#         (1) dendrogram (DendrogramPlot, showing the module assignment)
#         (2) module network plot (ModuleNetworkPlot)
#         (3) module eigengene heatmap (MEs by group)
#
# Output:
#   path34_hdwgcna_dendrogram.pdf
#   path34_hdwgcna_network/          (module network plots)
#   path34_hdwgcna_MEs_heatmap.pdf
# =============================================================================

## ---- Dependencies ----
if (!requireNamespace("hdWGCNA", quietly = TRUE)) {
  install.packages("hdWGCNA", repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(Seurat); library(hdWGCNA); library(ggplot2); library(dplyr)
})

## ---- Paths ----
this_file <- commandArgs(trailingOnly = FALSE)
.f <- grep("--file=", this_file, value = TRUE)
if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1]) else .this_file <- "."
out_dir <- dirname(normalizePath(.this_file))
source(file.path(out_dir, "00_config.R"))

## ============================================================================
## Read the hdWGCNA object (re-run if it does not exist)
## ============================================================================
rds_file <- file.path(out_dir, "path13_sc_hdwgcna.rds")
if (!file.exists(rds_file)) {
  cat("Cannot find path13_sc_hdwgcna.rds, re-running hdWGCNA (reusing the core steps of script 13)...\n")
  blood_file <- file.path(out_dir, "path1_sepsis_copd_integrated.rds")
  if (!file.exists(blood_file)) stop("Cannot find path1_sepsis_copd_integrated.rds")
  seu <- readRDS(blood_file)
  seu <- subset(seu, subset = cell_type %in% c("Monocyte", "DC"))
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- SetupForWGCNA(seu, gene_select = "fraction", fraction = 0.05,
                       wgcna_name = "myeloid")
  seu <- ScaleData(seu, features = GetWGCNAGenes(seu), verbose = FALSE)
  seu <- MetacellsByGroups(seu, group.by = c("cell_type", "group"),
                           k = 25, max_shared = 10, ident.group = "group")
  if ("NormalizeMetacells" %in% getNamespaceExports("hdWGCNA")) {
    seu <- NormalizeMetacells(seu, wgcna_name = "myeloid")
  } else {
    m_obj <- GetMetacellObject(seu, wgcna_name = "myeloid")
    m_obj <- NormalizeData(m_obj, verbose = FALSE)
    seu <- SetMetacellObject(seu, m_obj, wgcna_name = "myeloid")
  }
  seu <- SetDatExpr(seu, group_name = "Monocyte", group.by = "cell_type")
  seu <- TestSoftPowers(seu, networkType = "signed")
  sp <- 8
  pt <- tryCatch(GetPowerTable(seu), error = function(e) NULL)
  if (!is.null(pt)) {
    ok <- which(pt$SFT.R.sq >= 0.85)
    if (length(ok) > 0) sp <- pt$Power[ok[1]]
  }
  cat("soft power =", sp, "\n")
  seu <- ConstructNetwork(seu, soft_power = sp, setDatExpr = FALSE,
                          overwrite_tom = TRUE)
  seu <- ModuleEigengenes(seu, group.by.vars = "group")
  saveRDS(seu, rds_file)
  cat("Re-run finished and saved:", rds_file, "\n")
} else {
  cat("Reading hdWGCNA object:", rds_file, "\n")
  seu <- readRDS(rds_file)
}

## ============================================================================
## Generate figures
## ============================================================================
# Add ModuleConnectivity (script 13 omitted this step; ModuleNetworkPlot needs kME)
seu <- tryCatch(ModuleConnectivity(seu), error = function(e) {
  cat("ModuleConnectivity failed:", conditionMessage(e), "\n"); seu
})

# (1) dendrogram (PlotDendrogram uses base R graphics, needs the pdf() device, not ggsave)
cat("\n===== Figure (1) dendrogram =====\n")
dend_ok <- tryCatch({
  pdf(file.path(out_dir, "path34_hdwgcna_dendrogram.pdf"), width = 10, height = 8)
  PlotDendrogram(seu, wgcna_name = "myeloid",
                 main = "Monocyte hdWGCNA co-expression modules")
  dev.off()
  TRUE
}, error = function(e) {
  cat("Dendrogram error:", conditionMessage(e), "\n")
  tryCatch(dev.off(), error = function(e2) NULL)
  FALSE
})
if (dend_ok) {
  cat("Saved: path34_hdwgcna_dendrogram.pdf\n")
} else {
  cat("Dendrogram plotting failed (skipped)\n")
}

# (2) module network plot
cat("\n===== Figure (2) module network plot =====\n")
tryCatch(ModuleNetworkPlot(seu, wgcna_name = "myeloid",
                           outdir = file.path(out_dir, "path34_hdwgcna_network")),
         error = function(e) cat("Network plot failed:", conditionMessage(e), "\n"))

# (3) module eigengene heatmap (MEs by group)
cat("\n===== Figure (3) module eigengene heatmap =====\n")
tryCatch({
  # harmonized=FALSE returns metacell-level MEs (row names "cell_type#group_idx");
  # harmonized=TRUE returns hMEs projected back onto single cells (row names = cell barcodes), so FALSE is used here
  MEs <- GetMEs(seu, harmonized = FALSE)
  modules <- colnames(MEs)
  cat("First 3 rownames(MEs):", paste(head(rownames(MEs), 3), collapse = ", "), "\n")
  # Parse group from the metacell name (MetacellsByGroups naming format "cell_type#group_idx")
  mc_parts <- strsplit(rownames(MEs), "#", fixed = TRUE)
  MEs$group <- vapply(mc_parts, function(x) {
    g <- if (length(x) >= 2) x[2] else x[1]
    sub("_[0-9]+$", "", g)
  }, character(1))
  cat("group values:", paste(unique(MEs$group), collapse = ", "), "\n")
  # Long table: mean of each module by group (base R, avoids the tidyr dependency)
  me_long <- do.call(rbind, lapply(modules, function(m) {
    data.frame(group = MEs$group, module = m, ME = MEs[[m]],
               stringsAsFactors = FALSE)
  }))
  me_long <- aggregate(ME ~ group + module, data = me_long, FUN = mean)
  p <- ggplot(me_long, aes(x = group, y = module, fill = ME)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
    labs(x = NULL, y = NULL, fill = "Module\neigengene") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(out_dir, "path34_hdwgcna_MEs_heatmap.pdf"), p,
         width = 7, height = 8)
  cat("Saved: path34_hdwgcna_MEs_heatmap.pdf\n")
}, error = function(e) cat("MEs heatmap failed:", conditionMessage(e), "\n"))

cat("\nScript 34 hdWGCNA figure generation done\n")
