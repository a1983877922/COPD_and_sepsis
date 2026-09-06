# =========================================================================
# 11_pdc_deep_analysis.R
# Separate downstream analysis of pDC (plasmacytoid dendritic cells):
#   1) Differential genes (DEG): pDC gene changes in sepsis/infection control/COPD vs healthy
#   2) Pathway enrichment: GO/KEGG of sepsis-specific pDC genes (optional, requires clusterProfiler)
#   3) Trajectory analysis: pDC activation pseudotime (Monocle3, auto-installs SeuratWrappers + monocle3, not skipped)
#
# Note: blood-side pDC cell counts Healthy 5143 / Sepsis 1268 / IC 763 / COPD 87 / SP 41,
#   so COPD and Sepsis_Pneumonia groups have too few pDCs; DEG is reliable only for the first four groups,
#   interpret COPD with caution.
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

cat("\n==============================================================\n")
cat("pDC deep analysis (differential genes + pathway enrichment + trajectory)\n")
cat("==============================================================\n")

dc_markers <- list(
  "cDC1" = c("CLEC9A", "XCR1", "BATF3", "CADM1"),
  "cDC2" = c("CD1C", "CLEC10A", "FCER1A", "ITGAX", "HLA-DQA1"),
  "pDC"  = c("LILRA4", "CLEC4C", "IL3RA", "TCF4", "IRF7")
)

## =========================================================================
## Step 1: read path1 -> subset DC -> subclustering -> dc_subtype
## =========================================================================
cat("\n===== Step 1: blood-side DC subclustering =====\n")

blood_seu <- readRDS(file.path(out_dir, "path1_sepsis_copd_integrated.rds"))
cat("Blood object:", ncol(blood_seu), "cells\n")

blood_dc <- subset(blood_seu, subset = cell_type == "DC")
rm(blood_seu); gc()
cat("Blood-side DC count:", ncol(blood_dc), "\n")

# Fix out-of-sync @cells mapping after subset (Seurat v5 bug)
dc_counts <- GetAssayData(blood_dc, assay = "RNA", layer = "counts")
blood_dc[["RNA"]] <- CreateAssay5Object(counts = dc_counts)
rm(dc_counts); gc()

old_res <- resolution
resolution <- 0.6
blood_dc <- integrate_and_annotate(blood_dc, marker_list = dc_markers,
                                   col_name = "dc_subtype")
resolution <- old_res
cat("\nBlood-side DC subtypes:\n")
print(table(blood_dc$dc_subtype, blood_dc$group))

## =========================================================================
## Step 2: extract pDC subset + save
## =========================================================================
cat("\n===== Step 2: extract pDC =====\n")

pdc <- subset(blood_dc, subset = dc_subtype == "pDC")
# Fix subset bug (rebuild assay), and re-run NormalizeData to generate the data layer
# (run_deg's FindMarkers needs the data layer; if only counts were kept it errors "Layer data is empty")
pdc_counts <- GetAssayData(pdc, assay = "RNA", layer = "counts")
pdc[["RNA"]] <- CreateAssay5Object(counts = pdc_counts)
rm(pdc_counts); gc()
pdc <- NormalizeData(pdc, verbose = FALSE)
# Add dimensionality reduction (trajectory analysis Monocle3 needs UMAP/PCA)
pdc <- FindVariableFeatures(pdc, nfeatures = 2000, verbose = FALSE)
pdc <- ScaleData(pdc, features = VariableFeatures(pdc), verbose = FALSE)
pdc <- RunPCA(pdc, features = VariableFeatures(pdc), npcs = 20, verbose = FALSE)
pdc <- RunUMAP(pdc, dims = 1:20, verbose = FALSE)

cat("pDC count:", ncol(pdc), "\n")
print(table(pdc$group))
saveRDS(pdc, file.path(out_dir, "path11_pdc_subset.rds"))
cat("Saved pDC object: path11_pdc_subset.rds\n")

## =========================================================================
## Step 3: pDC differential genes (DEG)
## =========================================================================
cat("\n===== Step 3: pDC DEG =====\n")

Idents(pdc) <- pdc$group

# pDC cell count per group (determines which comparisons are reliable)
n_by_grp <- table(pdc$group)

deg_results <- list()

# 3.1 Sepsis vs Healthy (core)
if (n_by_grp[["Sepsis"]] >= 10 && n_by_grp[["Healthy"]] >= 10) {
  cat("\n--- Sepsis vs Healthy ---\n")
  deg_results[["Sepsis_vs_Healthy"]] <- run_deg(pdc, "Sepsis", "Healthy")
  write.csv(deg_results[["Sepsis_vs_Healthy"]],
            file.path(out_dir, "path11_pdc_DEG_Sepsis_vs_Healthy.csv"))
}

# 3.2 Infection_Control vs Healthy
if (n_by_grp[["Infection_Control"]] >= 10 && n_by_grp[["Healthy"]] >= 10) {
  cat("\n--- Infection_Control vs Healthy ---\n")
  deg_results[["IC_vs_Healthy"]] <- run_deg(pdc, "Infection_Control", "Healthy")
  write.csv(deg_results[["IC_vs_Healthy"]],
            file.path(out_dir, "path11_pdc_DEG_InfectionControl_vs_Healthy.csv"))
}

# 3.3 Sepsis vs Infection_Control (sepsis-specific, relative to infection)
if (n_by_grp[["Sepsis"]] >= 10 && n_by_grp[["Infection_Control"]] >= 10) {
  cat("\n--- Sepsis vs Infection_Control (sepsis-specific) ---\n")
  deg_results[["Sepsis_vs_IC"]] <- run_deg(pdc, "Sepsis", "Infection_Control")
  write.csv(deg_results[["Sepsis_vs_IC"]],
            file.path(out_dir, "path11_pdc_DEG_Sepsis_vs_InfectionControl.csv"))
}

# 3.4 COPD vs Healthy (caution: only 87 COPD pDCs)
if (n_by_grp[["COPD"]] >= 10 && n_by_grp[["Healthy"]] >= 10) {
  cat("\n--- COPD vs Healthy (few COPD pDCs, caution) ---\n")
  deg_results[["COPD_vs_Healthy"]] <- run_deg(pdc, "COPD", "Healthy")
  write.csv(deg_results[["COPD_vs_Healthy"]],
            file.path(out_dir, "path11_pdc_DEG_COPD_vs_Healthy.csv"))
}

## =========================================================================
## Step 4: sepsis-specific + comorbidity-shared genes
## =========================================================================
cat("\n===== Step 4: specific/shared genes =====\n")

# 4.1 Sepsis-specific (Sepsis vs IC up-regulated)
if (!is.null(deg_results[["Sepsis_vs_IC"]])) {
  sep_spec <- rownames(deg_results[["Sepsis_vs_IC"]])[
    deg_results[["Sepsis_vs_IC"]]$avg_log2FC > 0]
  cat("Number of sepsis-specific pDC up-regulated genes:", length(sep_spec), "\n")
  cat("Top 50:", paste(head(sep_spec, 50), collapse = ", "), "\n")
  writeLines(sep_spec, file.path(out_dir, "path11_pdc_sepsis_specific_genes.txt"))
}

# 4.2 COPD and Sepsis shared (comorbidity-shared at the pDC level)
if (!is.null(deg_results[["Sepsis_vs_Healthy"]]) &&
    !is.null(deg_results[["COPD_vs_Healthy"]])) {
  up_sep  <- rownames(deg_results[["Sepsis_vs_Healthy"]])[
    deg_results[["Sepsis_vs_Healthy"]]$avg_log2FC > 0]
  up_copd <- rownames(deg_results[["COPD_vs_Healthy"]])[
    deg_results[["COPD_vs_Healthy"]]$avg_log2FC > 0]
  shared  <- intersect(up_sep, up_copd)
  cat("Number of COPD and Sepsis shared up-regulated genes at pDC level:", length(shared), "\n")
  cat("Shared up-regulated genes:", paste(head(shared, 50), collapse = ", "), "\n")
  writeLines(shared, file.path(out_dir, "path11_pdc_shared_COPD_Sepsis_genes.txt"))
}

## =========================================================================
## Step 5: pathway enrichment (sepsis-specific pDC genes, optional)
## =========================================================================
cat("\n===== Step 5: pathway enrichment =====\n")

if (!is.null(deg_results[["Sepsis_vs_IC"]]) &&
    requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  suppressPackageStartupMessages(library(clusterProfiler))
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  genes <- rownames(deg_results[["Sepsis_vs_IC"]])[
    deg_results[["Sepsis_vs_IC"]]$avg_log2FC > 0]
  if (length(genes) >= 5) {
    eg <- tryCatch(bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db), error = function(e) NULL)
    if (!is.null(eg) && nrow(eg) >= 5) {
      ego <- enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                      pvalueCutoff = 0.05, qvalueCutoff = 0.2)
      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        write.csv(as.data.frame(ego), file.path(out_dir, "path11_pdc_GO_BP.csv"))
        cat("pDC sepsis-specific GO_BP enrichment saved\n")
      } else cat("GO_BP no significant enrichment\n")
    } else cat("Insufficient ENTREZ mapping, skipping enrichment\n")
  } else cat("Fewer than 5 sepsis-specific genes, skipping enrichment\n")
} else {
  cat("clusterProfiler/org.Hs.eg.db not installed or no Sepsis_vs_IC result, skipping enrichment\n")
}

## =========================================================================
## Step 6: trajectory analysis (pDC activation pseudotime, Monocle3) -- auto-install deps, do not skip
## =========================================================================
cat("\n===== Step 6: trajectory analysis (Monocle3) =====\n")

# 6.1 Ensure SeuratWrappers is installed (CRAN)
if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
  cat("SeuratWrappers not detected, installing from CRAN...\n")
  install.packages("SeuratWrappers", repos = "https://cloud.r-project.org",
                   Ncpus = max(1, parallel::detectCores() %/% 2))
}
if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
  stop("SeuratWrappers installation failed, please install manually: install.packages('SeuratWrappers')")
}

# 6.2 Ensure monocle3 is installed (Bioconductor)
if (!requireNamespace("monocle3", quietly = TRUE)) {
  cat("monocle3 not detected, installing from Bioconductor...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install("monocle3", update = FALSE, ask = FALSE)
}
if (!requireNamespace("monocle3", quietly = TRUE)) {
  stop("monocle3 installation failed, please install manually: BiocManager::install('monocle3')")
}

suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(SeuratWrappers))

# 6.3 Trajectory analysis (do not skip)
cds <- as.cell_data_set(pdc)
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)

# In non-interactive mode order_cells must explicitly specify the root (otherwise errors
# "root_pr_nodes or root_cells must be provided").
# Use the Healthy-group pDCs as the pseudotime origin (healthy resting state -> disease activation state).
root_cells <- colnames(cds)[colData(cds)$group == "Healthy"]
if (length(root_cells) == 0) root_cells <- colnames(cds)[1]  # fallback
cds <- order_cells(cds, root_cells = root_cells, reduction_method = "UMAP")

# Pseudotime plot
p_traj <- plot_cells(cds, color_cells_by = "pseudotime",
                     label_cell_groups = FALSE, label_leaves = FALSE)
ggsave(file.path(out_dir, "path11_pdc_trajectory.pdf"), p_traj,
       width = 7, height = 6)

# Trajectory colored by group
p_grp <- plot_cells(cds, color_cells_by = "group")
ggsave(file.path(out_dir, "path11_pdc_trajectory_by_group.pdf"), p_grp,
       width = 7, height = 6)

# Pseudotime range
pseudo <- pseudotime(cds)
cat("Pseudotime range:", round(min(pseudo, na.rm = TRUE), 2), "->",
    round(max(pseudo, na.rm = TRUE), 2), "\n")
cat("Trajectory analysis complete\n")

cat("\npDC deep analysis complete\n")
cat("Output: path11_pdc_subset.rds / path11_pdc_DEG_*.csv / path11_pdc_*_genes.txt\n")
cat("      path11_pdc_GO_BP.csv / path11_pdc_trajectory.pdf / path11_pdc_trajectory_by_group.pdf\n")
