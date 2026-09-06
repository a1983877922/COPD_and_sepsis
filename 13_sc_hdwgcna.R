# =========================================================================
# 13_sc_hdwgcna.R (single-cell hdWGCNA)
# Single-cell co-expression network analysis (hdWGCNA): aggregate monocytes
# (Monocyte, path1 object) into metacells, then build co-expression modules to
# address single-cell sparsity.
#
# Complements 12 (bulk whole-transcriptome WGCNA): 12 operates at the bulk
# level, this script at the single-cell level.
#
# Dependencies: hdWGCNA (required), Seurat, WGCNA; all must be installed
#   hdWGCNA: remotes::install_github('smorabit/hdWGCNA')
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

if (!requireNamespace("hdWGCNA", quietly = TRUE)) {
  stop("hdWGCNA not installed. Please run first: remotes::install_github('smorabit/hdWGCNA')")
}
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Seurat not installed")
}
suppressPackageStartupMessages(library(hdWGCNA))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

cat("\n==============================================================\n")
cat("Single-cell hdWGCNA (monocyte co-expression modules)\n")
cat("==============================================================\n")

## =========================================================================
## Step 1: Load path1 object + extract monocytes
## =========================================================================
cat("\n===== Step 1: Extract monocytes =====\n")

seu <- readRDS(file.path(out_dir, "path1_sepsis_copd_integrated.rds"))
cat("Single-cell object:", ncol(seu), "cells\n")

# Extract monocyte/myeloid (the core carrier of the 149 shared genes = blood Sepsis monocytes)
seu <- subset(seu, subset = cell_type %in% c("Monocyte", "DC"))
# Fix desynchronized @cells mapping after subset (Seurat v5 bug)
cc <- GetAssayData(seu, assay = "RNA", layer = "counts")
seu[["RNA"]] <- CreateAssay5Object(counts = cc); rm(cc); gc()
# Critical: add the data layer (hdWGCNA's SetDatExpr reads the data layer,
#       otherwise it errors "Layer 'data' is empty" + "undefined columns selected")
seu <- NormalizeData(seu, verbose = FALSE)
cat("Monocyte/myeloid cells:", ncol(seu), "\n")
print(table(seu$cell_type, seu$group))

## =========================================================================
## Step 2: SetupForWGCNA + metacell aggregation
## =========================================================================
cat("\n===== Step 2: metacell aggregation =====\n")

seu <- SetupForWGCNA(seu, gene_select = "fraction", fraction = 0.05,
                     wgcna_name = "myeloid")
cat("SetupForWGCNA done, number of genes selected:", length(GetWGCNAGenes(seu)), "\n")

# Critical: add the scale.data layer. ModuleEigengenes(group.by.vars) internally
# calls RunHarmony's ProjectDim, which needs scale.data for the
# "scale.data %*% harmony_embeddings" projection. After subset rebuilds the
# assay, scale.data is dropped -> errors "non-conformable arguments".
# Only ScaleData on WGCNA genes to save memory (whole-genome on hundreds of
# thousands of cells is too expensive).
seu <- ScaleData(seu, features = GetWGCNAGenes(seu), verbose = FALSE)
cat("ScaleData done (added scale.data layer)\n")

# metacell aggregation (k metacells per group, addressing single-cell sparsity)
seu <- MetacellsByGroups(seu, group.by = c("cell_type", "group"),
                         k = 25, max_shared = 10, ident.group = "group")
cat("metacell aggregation done\n")

# Critical: the new hdWGCNA MetacellsByGroups only aggregates counts, does not
# generate a data layer, while SetDatExpr reads the metacell object's data layer
# by default -> errors "Layer 'data' is empty". Must manually run
# NormalizeMetacells to add the metacell data layer.
if ("NormalizeMetacells" %in% getNamespaceExports("hdWGCNA")) {
  seu <- NormalizeMetacells(seu, wgcna_name = "myeloid")
  cat("NormalizeMetacells done (added metacell data layer)\n")
} else {
  # Old version has no NormalizeMetacells; manually run NormalizeData on the metacell object
  m_obj <- GetMetacellObject(seu, wgcna_name = "myeloid")
  m_obj <- NormalizeData(m_obj, verbose = FALSE)
  seu <- SetMetacellObject(seu, m_obj, wgcna_name = "myeloid")
  cat("Manual NormalizeData on metacell object done\n")
}

## =========================================================================
## Step 3: SetDatExpr (using monocyte expression)
## =========================================================================
cat("\n===== Step 3: SetDatExpr =====\n")

seu <- SetDatExpr(seu, group_name = "Monocyte", group.by = "cell_type")
cat("SetDatExpr done\n")

## =========================================================================
## Step 4: Soft threshold + network construction
## =========================================================================
cat("\n===== Step 4: soft threshold + network construction =====\n")

seu <- TestSoftPowers(seu, networkType = "signed")

# Select soft power from TestSoftPowers results: the smallest power where SFT.R.sq first >= 0.85
pt <- tryCatch(seu@misc[["myeloid"]]$power_tables$TestSoftPowers,
               error = function(e) NULL)
if (!is.null(pt) && "SFT.R.sq" %in% colnames(pt) && "Power" %in% colnames(pt)) {
  idx <- which(pt$SFT.R.sq >= 0.85)
  sp <- if (length(idx) > 0) pt$Power[min(idx)] else 8
} else {
  sp <- 8
}
cat("soft power =", sp, "\n")

# overwrite_tom = TRUE avoids the "TOM already exists" error on reruns
seu <- ConstructNetwork(seu, soft_power = sp, setDatExpr = FALSE,
                        overwrite_tom = TRUE)
cat("Network construction done\n")

## =========================================================================
## Step 5: Module eigengenes + module-disease association
## =========================================================================
cat("\n===== Step 5: Module eigengenes =====\n")

seu <- ModuleEigengenes(seu, group.by.vars = "group")

# Extract module genes
modules <- GetModules(seu)
write.csv(modules, file.path(out_dir, "path13_sc_gene_module.csv"), row.names = FALSE)
cat("Number of modules:", length(unique(modules$module)), "\n")
print(table(modules$module))

# Gene list for each module
for (mod in unique(modules$module)) {
  gs <- modules$gene_name[modules$module == mod]
  writeLines(gs, file.path(out_dir, paste0("path13_sc_module_", mod, ".txt")))
}

# Module eigengenes (metacell level)
MEs <- GetMEs(seu, harmonized = TRUE)
write.csv(MEs, file.path(out_dir, "path13_sc_MEs.csv"))

# Map the 149 shared genes to single-cell modules
genes_149 <- readLines(file.path(out_dir, "path2_shared_myeloid_genes.txt"))
genes_149 <- genes_149[genes_149 != ""]
mapped_149 <- modules[modules$gene_name %in% genes_149, ]
cat("149 shared genes mapped to single-cell modules:", nrow(mapped_149), "\n")
print(table(mapped_149$module))
write.csv(mapped_149, file.path(out_dir, "path13_sc_149gene_module.csv"),
          row.names = FALSE)

## =========================================================================
## Step 6: Visualization
## =========================================================================
cat("\n===== Step 6: Visualization =====\n")

p_dend <- tryCatch(DendrogramPlot(seu, main = "Monocyte WGCNA co-expression modules"),
                   error = function(e) NULL)
if (!is.null(p_dend)) {
  ggsave(file.path(out_dir, "path13_sc_dendrogram.pdf"), p_dend, width = 10, height = 8)
}

p_net <- tryCatch(ModuleNetworkPlot(seu, outdir = file.path(out_dir, "path13_sc_network")),
                  error = function(e) NULL)

## =========================================================================
## Step 7: Save
## =========================================================================
cat("\n===== Step 7: Save =====\n")

saveRDS(seu, file.path(out_dir, "path13_sc_hdwgcna.rds"))
cat("Saved hdWGCNA object: path13_sc_hdwgcna.rds\n")

cat("\nSingle-cell hdWGCNA complete\n")
cat("Outputs: path13_sc_gene_module.csv / path13_sc_module_*.txt / path13_sc_MEs.csv /\n")
cat("      path13_sc_149gene_module.csv / path13_sc_dendrogram.pdf / path13_sc_hdwgcna.rds\n")
cat("Interpretation: compare with the bulk modules from 12 to see whether the single-cell modules align with the bulk modules / four major pathways\n")
