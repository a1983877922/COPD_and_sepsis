###############################################################################
# Sepsis + COPD comorbidity analysis — Path 2: lung-blood cross-tissue analysis
# File: 02_build_lung_cross_tissue.R (formerly 02_path2_cross_tissue.R)
#
# Purpose: address "how the systemic and pulmonary compartments interact"
#   COPD is a local lung disease (lung tissue), sepsis is systemic inflammation (peripheral blood).
#   Their cell types do not overlap, so they cannot be naively merged; this script adopts
#   "independent annotation + cross-tissue comparison":
#     Lung side: GSE136831 (COPD lung atlas, with built-in Manuscript_Identity annotation)
#     Blood side: re-read 5 PBMC core datasets from raw data, reintegrate + annotate
#           (SCP548 + GSE279452 (adult) + GSE151263/GSE167363/GSE175453,
#            5000 cells/sample, not relying on sepsis_integrated.rds)
#   Compare immune programs (monocyte/macrophage, T cells, etc.) between the two sides to find
#   molecular programs shared between systemic and lung compartments.
#
# Workflow:
#   Step 1: read GSE136831 lung tissue (filter Control + COPD, exclude IPF) + light QC
#   Step 2: read sepsis blood raw data + integrate + annotate
#   Step 3: cross-tissue cell composition overview
#   Step 4: monocyte/macrophage shared program (blood Sepsis DEG ∩ lung COPD DEG)
#   Step 5: T cell cross-tissue comparison (optional)
#   Step 6: save
#
# Prerequisite: 00_config.R has been run (source this script)
###############################################################################

# Auto-locate this script's directory and find 00_config.R in the same directory (no need to change this path locally or on the server)
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
cat("Path 2: lung-blood cross-tissue analysis (COPD lung tissue vs sepsis peripheral blood)\n")
cat("==============================================================\n")

## =========================================================================
## Step 1: read GSE136831 lung tissue (COPD + Control)
## =========================================================================

cat("\n===== Step 1: read COPD lung tissue =====\n")

lung_seu <- read_gse136831(file.path(copd_dir, "GSE136831"), gse_id = "GSE136831")
cat("Lung tissue:", ncol(lung_seu), "cells\n")
print(table(lung_seu$disease))

# Light QC (lung tissue already filtered by original authors; here only add mitochondrial/hemoglobin filtering)
lung_seu[["percent.mt"]] <- PercentageFeatureSet(lung_seu, pattern = "^MT-")
hb_genes <- grep("^HB[^PS]", rownames(lung_seu), value = TRUE)
lung_seu[["percent.hb"]] <- if (length(hb_genes) > 0)
  PercentageFeatureSet(lung_seu, features = hb_genes) else 0

lung_seu <- subset(lung_seu,
                   subset = nFeature_RNA > 200 & percent.mt < 20)
cat("Lung tissue after QC:", ncol(lung_seu), "cells\n")

# Fix Seurat v5 subset bug (#8407): after subset, the assay's @cells mapping may be out of sync,
# and subsequent NormalizeData errors with "Cannot add new cells with [[<-". Rebuild the RNA assay to fix.
lc_counts <- GetAssayData(lung_seu, assay = "RNA", layer = "counts")
lung_seu[["RNA"]] <- CreateAssay5Object(counts = lc_counts)
rm(lc_counts); gc()

# Generate the data layer (FindMarkers needs log-normalized data; read_gse136831 only built the counts layer)
lung_seu <- NormalizeData(lung_seu, verbose = FALSE)

# Built-in cell-type annotation overview
cat("\nLung tissue cell types (Manuscript_Identity):\n")
print(sort(table(lung_seu$cell_type), decreasing = TRUE))

## =========================================================================
## Step 2: read sepsis blood raw data + integrate + annotate
## =========================================================================

cat("\n===== Step 2: read sepsis blood raw data + integrate + annotate =====\n")

# 2.1 Read 5 PBMC core datasets from raw data (5000 cells/sample, adult subset)
blood_seu <- read_sepsis_pbmc()

# 2.2 Integrate + annotate (normalize → HVG → PCA → Harmony → clustering → UMAP → annotation)
blood_seu <- integrate_and_annotate(blood_seu)

cat("Blood cell types:", paste(names(table(blood_seu$cell_type)), collapse = ", "), "\n")

# Blood UMAP (by dataset / condition / cell_type)
p_b1 <- DimPlot(blood_seu, reduction = "umap", group.by = "dataset") + ggtitle("Blood UMAP by Dataset")
p_b2 <- DimPlot(blood_seu, reduction = "umap", group.by = "condition") + ggtitle("Blood UMAP by Condition")
p_b3 <- DimPlot(blood_seu, reduction = "umap", group.by = "cell_type", label = TRUE) + ggtitle("Blood UMAP by CellType")
ggsave(file.path(out_dir, "path2_blood_UMAP.pdf"), plot = p_b1 | p_b2 | p_b3,
       width = 24, height = 8)

# Save blood integrated object (for reuse)
saveRDS(blood_seu, file.path(out_dir, "path2_blood_integrated.rds"))

## =========================================================================
## Step 3: cross-tissue cell composition overview (no merging, only side-by-side comparison)
## =========================================================================

cat("\n===== Step 3: cross-tissue cell composition overview =====\n")

# Blood-side immune composition (by condition)
blood_comp <- table(blood_seu$cell_type, blood_seu$condition)
write.csv(blood_comp, file.path(out_dir, "path2_blood_celltype_by_condition.csv"))
cat("Blood-side cell composition (cell_type x condition):\n")
print(blood_comp)

# Lung-side composition (by disease)
lung_comp <- table(lung_seu$cell_type, lung_seu$disease)
write.csv(lung_comp, file.path(out_dir, "path2_lung_celltype_by_disease.csv"))
cat("\nLung-side cell composition (cell_type x disease):\n")
print(lung_comp)

## =========================================================================
## Step 4: monocyte/macrophage shared program (blood Sepsis vs lung COPD)
## =========================================================================

cat("\n===== Step 4: monocyte/macrophage shared program =====\n")

# 4.1 Lung side: myeloid cells (CellType_Category == "Myeloid")
lung_myeloid <- subset(lung_seu,
                       subset = cell_category == "Myeloid" &
                         disease %in% c("Control", "COPD"))
Idents(lung_myeloid) <- "disease"
cat("Lung myeloid cells:", ncol(lung_myeloid), "cells\n")
print(table(lung_myeloid$cell_type))

if (ncol(lung_myeloid) > 50 && length(unique(lung_myeloid$disease)) == 2) {
  lung_mye_deg <- FindMarkers(lung_myeloid, ident.1 = "COPD", ident.2 = "Control",
                              only.pos = FALSE, min.pct = 0.1,
                              logfc.threshold = 0.25, test.use = "wilcox")
  write.csv(lung_mye_deg, file.path(out_dir, "path2_lung_myeloid_DEG_COPD_vs_Control.csv"))
  cat("Lung myeloid COPD vs Control DEGs:", nrow(lung_mye_deg), "\n")
}

# 4.2 Blood side: monocytes (Monocyte), Sepsis vs Healthy
blood_mono <- subset(blood_seu, subset = cell_type == "Monocyte" &
                       condition %in% c("Sepsis", "Healthy"))
Idents(blood_mono) <- "condition"
cat("\nBlood monocytes:", ncol(blood_mono), "cells\n")

if (ncol(blood_mono) > 50 && length(unique(blood_mono$condition)) == 2) {
  blood_mono_deg <- FindMarkers(blood_mono, ident.1 = "Sepsis", ident.2 = "Healthy",
                                only.pos = FALSE, min.pct = 0.1,
                                logfc.threshold = 0.25, test.use = "wilcox")
  write.csv(blood_mono_deg, file.path(out_dir, "path2_blood_monocyte_DEG_Sepsis_vs_Healthy.csv"))
  cat("Blood monocyte Sepsis vs Healthy DEGs:", nrow(blood_mono_deg), "\n")
}

# 4.3 Shared program: lung COPD up-regulated ∩ blood Sepsis up-regulated
if (exists("lung_mye_deg") && exists("blood_mono_deg")) {
  lung_up  <- rownames(lung_mye_deg)[lung_mye_deg$avg_log2FC > 0]
  blood_up <- rownames(blood_mono_deg)[blood_mono_deg$avg_log2FC > 0]
  shared   <- intersect(lung_up, blood_up)
  cat("\nLung COPD myeloid ∩ blood Sepsis monocyte shared up-regulated genes:", length(shared), "\n")
  cat("Shared up-regulated genes:", paste(head(shared, 100), collapse = ", "), "\n")
  writeLines(shared, file.path(out_dir, "path2_shared_myeloid_genes.txt"))
}

## =========================================================================
## Step 5: T cell cross-tissue comparison (optional)
## =========================================================================

cat("\n===== Step 5: T cell cross-tissue comparison =====\n")

# Lung-side T cells (T-related among CellType_Category == "Lymphoid")
lung_t <- subset(lung_seu, subset = cell_category == "Lymphoid" &
                   grepl("T|CD4|CD8", cell_type) &
                   disease %in% c("Control", "COPD"))
cat("Lung T cells:", ncol(lung_t), "cells\n")

blood_t <- subset(blood_seu, subset = cell_type == "T_cell" &
                    condition %in% c("Sepsis", "Healthy"))
cat("Blood T cells:", ncol(blood_t), "cells\n")

if (ncol(lung_t) > 50 && length(unique(lung_t$disease)) == 2) {
  Idents(lung_t) <- "disease"
  lung_t_deg <- FindMarkers(lung_t, ident.1 = "COPD", ident.2 = "Control",
                            only.pos = FALSE, min.pct = 0.1,
                            logfc.threshold = 0.25, test.use = "wilcox")
  write.csv(lung_t_deg, file.path(out_dir, "path2_lung_T_DEG_COPD_vs_Control.csv"))
  cat("Lung T cell COPD vs Control DEGs:", nrow(lung_t_deg), "\n")
}

## =========================================================================
## Step 6: save
## =========================================================================

cat("\n===== Step 6: save =====\n")

saveRDS(lung_seu, file.path(out_dir, "path2_copd_lung.rds"))

cat("Saved:\n")
cat("  - path2_copd_lung.rds (COPD lung tissue object, with built-in annotation)\n")
cat("Path 2 analysis complete\n")
cat("\n[Note] How to use cross-tissue comparison conclusions:\n")
cat("  - path2_shared_myeloid_genes.txt = monocyte/macrophage program shared by systemic (blood) and local lung compartments\n")
cat("  - cross-validate with the shared genes from Path 1 to support the 'shared immune dysregulation in comorbidity' conclusion\n")
cat("Session Info:\n")
sessionInfo()
