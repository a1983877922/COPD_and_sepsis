###############################################################################
# Sepsis + COPD comorbidity analysis — Path 1: direct blood-blood merging
# File: 01_assemble_blood_integrated.R (formerly 01_path1_blood_blood_merge.R)
#
# Purpose: address "shared vs differentiated systemic circulating immune states"
#   Merge sepsis PBMC with COPD PBMC (GSE249584) into a single integrated space,
#   using a five-group design (Healthy / Infection_Control / COPD / Sepsis / Sepsis_Pneumonia)
#   to identify shared and differentiated immune dysregulation states and infer comorbidity mechanisms.
#
# Data:
#   Sepsis: re-read 5 PBMC 10x core datasets from raw data (not relying on merged_seu_qc.rds)
#           SCP548 + GSE279452 (adult subset) + GSE151263/GSE167363/GSE175453
#           (excluded: GSE252331/GSE216009/Kwok2023 whole blood, OMIX005600 GEXSCOPE,
#            GSE279451 overlaps with GSE279452 from the same study)
#   COPD: GSE249584 (peripheral blood PBMC, 10x 3' v3.1, 7 HC + 8 COPD)
#
# Workflow:
#   Step 1: read sepsis 5 datasets from raw data + read GSE249584 + QC
#   Step 2: merge + downsample
#   Step 3: normalize + highly variable genes + Scale + PCA
#   Step 4: Harmony integration (batch = dataset)
#   Step 5: clustering + UMAP
#   Step 6: cell annotation
#   Step 7: four-group comorbidity analysis (proportions / DEGs / shared states)
#   Step 8: CellChat comorbidity communication analysis
#   Step 9: save
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
cat("Path 1: direct blood-blood merging (sepsis PBMC + COPD PBMC)\n")
cat("==============================================================\n")

## =========================================================================
## Resumable checkpoint: continue from last successful step after an error, no need to rerun from scratch
##   cp1_merged.rds     = after Step 2 (read data + merge + downsample + grouping)
##   cp2_integrated.rds = after Step 5 (normalize + Harmony + clustering + UMAP, unannotated)
##   cp3_annotated.rds  = after Step 7 (annotation + FindAllMarkers + comorbidity analysis, with cell_type/group)
## =========================================================================
cp1_file <- file.path(out_dir, "cp1_merged.rds")
cp2_file <- file.path(out_dir, "cp2_integrated.rds")
cp3_file <- file.path(out_dir, "cp3_annotated.rds")

if (file.exists(cp3_file)) {
  cat("\n[checkpoint] Detected cp3_annotated.rds, continuing from Step 8 (CellChat)\n")
  merged <- readRDS(cp3_file)
  run_step1to2 <- FALSE; run_step3to5 <- FALSE; run_step6to7 <- FALSE
} else if (file.exists(cp2_file)) {
  cat("\n[checkpoint] Detected cp2_integrated.rds, continuing from Step 6 (annotation)\n")
  merged <- readRDS(cp2_file)
  run_step1to2 <- FALSE; run_step3to5 <- FALSE; run_step6to7 <- TRUE
} else if (file.exists(cp1_file)) {
  cat("\n[checkpoint] Detected cp1_merged.rds, continuing from Step 3 (normalization)\n")
  merged <- readRDS(cp1_file)
  run_step1to2 <- FALSE; run_step3to5 <- TRUE; run_step6to7 <- TRUE
} else {
  cat("\n[checkpoint] No checkpoint found, running from scratch\n")
  run_step1to2 <- TRUE; run_step3to5 <- TRUE; run_step6to7 <- TRUE
}

## =========================================================================
## Step 1-2: read data + merge (skip if cp1/cp2 exists)
## =========================================================================
if (run_step1to2) {

## =========================================================================
## Step 1: load sepsis object + read COPD PBMC
## =========================================================================

cat("\n===== Step 1: load sepsis + read COPD =====\n")

# 1.1 Sepsis: re-read 5 PBMC core datasets from raw data
#   keep only: SCP548 + GSE279452 (adult subset) + GSE151263/GSE167363/GSE175453
sepsis_seu <- read_sepsis_pbmc()

# 1.2 COPD: read GSE249584 (peripheral blood PBMC, 10x mtx)
copd_gse_dir <- file.path(copd_dir, "GSE249584")
copd_list <- read_copd_10x(copd_gse_dir, gse_id = "GSE249584",
                           max_cells = Inf)

# Verify condition mapping is fully filled in (a TODO means the user has not filled it yet)
copd_cond_vals <- unique(sapply(copd_list, function(x) unique(x$condition)))
if (any(grepl("TODO", copd_cond_vals))) {
  cat("\n[Warning] GSE249584 condition mapping incomplete, please edit copd_condition_map in 00_config.R\n")
}

# 1.3 COPD QC
copd_seu <- merge(copd_list[[1]], y = copd_list[-1],
                  add.cell.ids = names(copd_list), project = "COPD")
rm(copd_list); gc()
copd_seu <- mad_qc(copd_seu, nmads = 3)
cat("COPD cells after QC:", ncol(copd_seu), "cells\n")
print(table(copd_seu$condition))

## =========================================================================
## Step 2: merge sepsis + COPD
## =========================================================================

cat("\n===== Step 2: merge =====\n")

# Unify metadata columns (keep columns needed later)
common_meta <- intersect(colnames(sepsis_seu@meta.data),
                         colnames(copd_seu@meta.data))
cat("Shared metadata columns:", paste(common_meta, collapse = ", "), "\n")

merged <- merge(sepsis_seu, y = copd_seu,
                add.cell.ids = c("Sepsis", "COPD"),
                project = "Sepsis_COPD")
rm(sepsis_seu, copd_seu); gc()

# Key: after merge there are ~185 split layers (one per sample); subsequent NormalizeData/ScaleData
# would have huge memory overhead across many layers and cause OOM (std::bad_alloc). Immediately merge into a single counts layer.
merged <- JoinLayers(merged, assay = "RNA")
gc()

cat("Total cells after merge:", ncol(merged), "\n")
print(table(merged$dataset))

# 2.1 Downsampling (per-sample cap max_cells_per_sample)
# Note: read_sepsis_pbmc already downsampled sepsis to 5000 cells per sample; do it once more here
#       (mainly to constrain newly added COPD samples and prevent too many cells from a single sample after merge)
if (max_cells_per_sample > 0) {
  meta_df <- merged@meta.data
  cells_keep <- unlist(lapply(unique(meta_df$sample), function(s) {
    cells_s <- rownames(meta_df)[meta_df$sample == s]
    if (length(cells_s) > max_cells_per_sample)
      sample(cells_s, max_cells_per_sample) else cells_s
  }))
  rm(meta_df)
  if (length(cells_keep) < ncol(merged)) {
    set.seed(42)
    merged <- merged[, cells_keep]
    gc()
  }
  cat("After downsampling:", ncol(merged), "cells\n")
}

# 2.2 Define unified grouping (five-group comorbidity design)
#   Healthy           : healthy controls (sepsis healthy controls + COPD healthy controls)
#   Infection_Control : infection controls (infected but not meeting sepsis criteria; GSE279452 IC + SCP548 ICU-NoSEP/Leuk-UTI)
#   COPD              : COPD (GSE249584)
#   Sepsis            : adult sepsis (non-pulmonary source)
#   Sepsis_Pneumonia  : pulmonary sepsis (GSE151263, secondary to bacterial pneumonia, used as a proxy for COPD-related infection)
#   (excluded: Pediatric_* not included in core comparison)
merged$group <- NA_character_
merged$group[merged$dataset == "GSE249584" & merged$condition == "COPD"]    <- "COPD"
merged$group[merged$dataset == "GSE249584" & merged$condition == "Control"] <- "Healthy"
merged$group[merged$dataset == "GSE151263"]                                 <- "Sepsis_Pneumonia"
merged$group[merged$dataset != "GSE249584" & merged$dataset != "GSE151263" &
             merged$condition == "Healthy"]                                 <- "Healthy"
merged$group[merged$dataset != "GSE249584" & merged$dataset != "GSE151263" &
             merged$condition == "Sepsis"]                                  <- "Sepsis"
merged$group[merged$condition == "Infection_Control"]                       <- "Infection_Control"

cat("\nFive-group sample distribution (by dataset):\n")
print(table(merged$dataset, merged$group, useNA = "ifany"))

  # Save checkpoint 1 (after merge)
  saveRDS(merged, cp1_file)
  cat("[checkpoint] Saved cp1_merged.rds\n")
}

## =========================================================================
## Step 3-5: normalize + Harmony + clustering (skip if cp2 exists)
## =========================================================================
if (run_step3to5) {

## =========================================================================
## Step 3: normalization + highly variable genes + Scale + PCA
## =========================================================================

cat("\n===== Step 3: normalization / HVG / PCA =====\n")

merged <- NormalizeData(merged, verbose = FALSE)
merged <- FindVariableFeatures(merged, selection.method = "vst",
                               nfeatures = varGeneNum, verbose = FALSE)

write.table(VariableFeatures(merged),
            file = file.path(out_dir, "varGene_path1.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
merged <- ScaleData(merged, features = VariableFeatures(merged),
                    vars.to.regress = vars_to_regress, verbose = FALSE)
merged <- RunPCA(merged, features = VariableFeatures(merged),
                 npcs = computePCs, verbose = FALSE)

## =========================================================================
## Step 4: Harmony integration (batch = dataset)
## =========================================================================

cat("\n===== Step 4: Harmony batch correction =====\n")
cat(sprintf("RunHarmony (theta=%d, lambda=%d)...\n", harmony_theta, harmony_lambda))

merged <- RunHarmony(merged, group.by.vars = "dataset",
                     theta = harmony_theta, lambda = harmony_lambda,
                     verbose = FALSE)
cat("Harmony done\n")

## =========================================================================
## Step 5: clustering + UMAP
## =========================================================================

cat("\n===== Step 5: clustering + UMAP =====\n")

merged <- FindNeighbors(merged, reduction = "harmony", dims = pcadim)
merged <- FindClusters(merged, resolution = resolution, verbose = FALSE)
merged <- RunUMAP(merged, dims = pcadim, n.neighbors = umap_nneighbors,
                  reduction = "harmony", reduction.name = "umap", verbose = FALSE)

# Free memory: PCA/UMAP already computed; scale.data (3000 HVG x n_cells dense matrix ~8.6GB)
# is no longer needed, remove it to reduce memory and avoid OOM during subsequent subsetting
rna_layers <- merged@assays$RNA@layers
rna_layers[grep("scale\\.data", names(rna_layers))] <- NULL
merged@assays$RNA@layers <- rna_layers
rm(rna_layers)
gc()

cat("Number of clusters:", length(unique(merged$seurat_clusters)), "\n")

# Visualization: by dataset / group / cluster
p1 <- DimPlot(merged, reduction = "umap", group.by = "dataset", label = FALSE) +
  ggtitle("UMAP by Dataset")
p2 <- DimPlot(merged, reduction = "umap", group.by = "group", label = FALSE) +
  ggtitle("UMAP by Group (4 groups)")
p3 <- DimPlot(merged, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP by Cluster")
ggsave(file.path(out_dir, "path1_UMAP_dataset_group_cluster.pdf"),
       plot = p1 | p2 | p3, width = 24, height = 8)

# Batch correction effect: before vs after correction (PCA vs Harmony)
p_before <- DimPlot(merged, reduction = "pca", group.by = "dataset") +
  ggtitle("PCA (before correction)")
p_after  <- DimPlot(merged, reduction = "harmony", group.by = "dataset") +
  ggtitle("Harmony (after correction)")
ggsave(file.path(out_dir, "path1_batch_before_after.pdf"),
       plot = p_before | p_after, width = 16, height = 8)

  # Save checkpoint 2 (after integration + clustering, unannotated)
  saveRDS(merged, cp2_file)
  cat("[checkpoint] Saved cp2_integrated.rds\n")
}

## =========================================================================
## Step 6-7: annotation + comorbidity analysis (skip if cp3 exists)
## =========================================================================
if (run_step6to7) {

## =========================================================================
## Step 6: cell annotation (canonical marker genes)
## =========================================================================

cat("\n===== Step 6: cell annotation =====\n")

merged <- annotate_by_markers(merged, marker_list = blood_markers)

p_ann <- DimPlot(merged, reduction = "umap", group.by = "cell_type",
                 label = TRUE, label.size = 3) + ggtitle("Cell Type Annotation")
ggsave(file.path(out_dir, "path1_UMAP_celltype.pdf"), plot = p_ann,
       width = 10, height = 8)

# Annotation verification: DotPlot shows marker gene expression per cell type (only plot genes present in the object)
avail_markers <- unique(unlist(lapply(blood_markers,
                                      function(g) intersect(g, rownames(merged)))))
cat("Number of available marker genes:", length(avail_markers), "\n")
p_dot <- DotPlot(merged, features = avail_markers, group.by = "cell_type") +
  RotatedAxis() +
  ggtitle("Marker Gene Verification") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
ggsave(file.path(out_dir, "path1_dotplot_markers.pdf"), plot = p_dot,
       width = 16, height = 7)

# Annotation verification: stacked proportion of each cell type across groups (check for abnormal group enrichment)
p_comp <- ggplot(as.data.frame(table(merged$cell_type, merged$group)) %>%
                   rename(CellType = Var1, Group = Var2, Count = Freq),
                 aes(x = Group, y = Count, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "Fraction", title = "Cell Type Composition by Group") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "path1_celltype_by_group_verify.pdf"),
       plot = p_comp, width = 9, height = 6)

cat("Annotated cell types:", paste(names(table(merged$cell_type)), collapse = ", "), "\n")

## =========================================================================
## Step 6.5: all-cell-type markers + low-quality cluster detection
## =========================================================================

cat("\n===== Step 6.5: marker computation + low-quality cluster detection =====\n")

# 6.5.1 All-cell-type markers (FindAllMarkers, for low-quality detection + annotation validation)
all_markers <- FindAllMarkers(merged, only.pos = markerOnlyPos,
                              min.pct = markerMinPct,
                              logfc.threshold = markerThresh,
                              test.use = markerMethod, base = markerBase)
write.table(data.frame(gene = all_markers$gene, all_markers[, 1:6]),
            file = file.path(out_dir, "path1_AllMarkerGenes.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# 6.5.2 Low-quality cluster detection (high proportion of lncRNA/low-information genes among top markers → mark as LowQuality)
low_q_clusters <- identify_low_quality_clusters(all_markers)
if (length(low_q_clusters) > 0) {
  merged$cell_type[merged$seurat_clusters %in% low_q_clusters] <- "LowQuality"
  cat("Cells marked as LowQuality:", sum(merged$cell_type == "LowQuality"), "\n")
} else {
  cat("No low-quality clusters found\n")
}

## =========================================================================
## Step 7: five-group comorbidity analysis (proportions / DEGs / shared states)
## =========================================================================

cat("\n===== Step 7: five-group comorbidity analysis =====\n")

# Keep only the five core groups for comorbidity comparison
core_groups <- c("Healthy", "Infection_Control", "COPD", "Sepsis", "Sepsis_Pneumonia")
seu_core <- subset(merged, cells = colnames(merged)[merged$group %in% core_groups])
seu_core$group <- factor(seu_core$group, levels = core_groups)
# Exclude low-quality cells (marked LowQuality in Step 6.5)
seu_core <- subset(seu_core, subset = cell_type != "LowQuality")
cat("Five-group core cell count (after excluding low-quality):", ncol(seu_core), "\n")
print(table(seu_core$group))

# 7.1 Cell proportion analysis (Ro/e style: observed/expected ratio)
#     identify cell types jointly enriched/depleted in COPD and Sepsis
prop_tab <- table(seu_core$cell_type, seu_core$group)
prop_frac <- sweep(prop_tab, 2, colSums(prop_tab), "/")
# Expected = overall proportion of each cell type
expected <- rowSums(prop_tab) / sum(prop_tab)
roe <- sweep(prop_frac, 1, expected, "/")   # Ro/e matrix
write.csv(roe, file.path(out_dir, "path1_roe_celltype_by_group.csv"))

cat("\nCell type Ro/e (observed/expected, >1 enriched, <1 depleted):\n")
print(round(roe, 2))

# 7.2 Shared dysregulated cell types: cell types changing in the same direction in both COPD and Sepsis
#      (both Ro/e >1 or both <1 = immune composition changes shared by the comorbidity)
shared_up   <- rownames(roe)[roe[, "COPD"] > 1 & roe[, "Sepsis"] > 1]
shared_down <- rownames(roe)[roe[, "COPD"] < 1 & roe[, "Sepsis"] < 1]
cat("\nCell types jointly enriched in COPD and Sepsis:", paste(shared_up, collapse = ", "), "\n")
cat("Cell types jointly depleted in COPD and Sepsis:", paste(shared_down, collapse = ", "), "\n")

# Proportion stacked bar chart
p_prop <- ggplot(as.data.frame(prop_frac) %>% rename(CellType = Var1, Group = Var2, Frac = Freq),
                 aes(x = Group, y = Frac, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "Fraction", title = "Cell Type Proportion by Group") +
  theme_bw(base_size = 10)
ggsave(file.path(out_dir, "path1_celltype_proportion.pdf"),
       plot = p_prop, width = 8, height = 6)

# 7.3 Shared states of monocyte/macrophage subsets (focus of comorbidity study: inflammation/immune dysregulation)
#     compare DEGs among COPD / Sepsis / Healthy on monocytes to find shared vs specific programs
if ("Monocyte" %in% unique(seu_core$cell_type)) {
  cat("\n--- Monocyte differential analysis (shared vs specific) ---\n")
  seu_mono <- subset(seu_core, subset = cell_type == "Monocyte")
  Idents(seu_mono) <- "group"

  mono_copd <- run_deg(seu_mono, "COPD", "Healthy")
  mono_sep  <- run_deg(seu_mono, "Sepsis", "Healthy")
  write.csv(mono_copd, file.path(out_dir, "path1_monocyte_DEG_COPD_vs_Healthy.csv"))
  write.csv(mono_sep,  file.path(out_dir, "path1_monocyte_DEG_Sepsis_vs_Healthy.csv"))

  # Sepsis-specific (on top of infection) + infection-shared (infection vs healthy) layered comparison
  mono_sep_ic <- run_deg(seu_mono, "Sepsis", "Infection_Control")
  mono_ic_hc  <- run_deg(seu_mono, "Infection_Control", "Healthy")
  write.csv(mono_sep_ic, file.path(out_dir, "path1_monocyte_DEG_Sepsis_vs_InfectionControl.csv"))
  write.csv(mono_ic_hc,  file.path(out_dir, "path1_monocyte_DEG_InfectionControl_vs_Healthy.csv"))

  # Sepsis-specific genes (still up-regulated relative to infection controls) = core sepsis program after removing "general infection response"
  sepsis_specific <- rownames(mono_sep_ic)[mono_sep_ic$avg_log2FC > 0]
  cat("\nNumber of sepsis-specific up-regulated genes (Sepsis vs Infection_Control):", length(sepsis_specific), "\n")
  cat("Sepsis-specific up-regulated genes (top 50):", paste(head(sepsis_specific, 50), collapse = ", "), "\n")
  writeLines(sepsis_specific, file.path(out_dir, "path1_monocyte_sepsis_specific_genes.txt"))

  # Shared DEGs (same direction) = molecular program shared by the comorbidity
  genes_copd <- rownames(mono_copd)[mono_copd$avg_log2FC > 0]
  genes_sep  <- rownames(mono_sep)[mono_sep$avg_log2FC > 0]
  shared_deg <- intersect(genes_copd, genes_sep)
  cat("Number of shared up-regulated genes in COPD and Sepsis monocytes:", length(shared_deg), "\n")
  cat("Shared up-regulated genes (top 50):", paste(head(shared_deg, 50), collapse = ", "), "\n")
  writeLines(shared_deg, file.path(out_dir, "path1_monocyte_shared_genes.txt"))
}

  # Save checkpoint 3 (after annotation + comorbidity analysis)
  saveRDS(merged, cp3_file)
  cat("[checkpoint] Saved cp3_annotated.rds\n")
}

## =========================================================================
## Step 8: CellChat comorbidity communication analysis
## =========================================================================

# If resuming from cp3 (run_step6to7=FALSE), seu_core was not generated; rebuild it here
if (!run_step6to7) {
  core_groups <- c("Healthy", "Infection_Control", "COPD", "Sepsis", "Sepsis_Pneumonia")
  seu_core <- subset(merged, cells = colnames(merged)[merged$group %in% core_groups])
  seu_core$group <- factor(seu_core$group, levels = core_groups)
  seu_core <- subset(seu_core, subset = cell_type != "LowQuality")
  cat("Resumed from cp3, rebuilt seu_core:", ncol(seu_core), "cells\n")
}

cat("\n===== Step 8: CellChat comorbidity communication =====\n")

if (!requireNamespace("CellChat", quietly = TRUE)) {
  cat("CellChat not installed, skipping communication analysis\n")
  cat("Install: devtools::install_github('jinworks/CellChat')\n")
} else {
  library(CellChat)
  library(NMF)
  library(circlize)

  cellchat_dir <- file.path(out_dir, "path1_CellChat")
  dir.create(cellchat_dir, showWarnings = FALSE, recursive = TRUE)

  # 8.1 Define single-group CellChat construction function
  run_cellchat <- function(seu, group_name) {
    # Server with 512GB RAM: cap downsampling per group via cellchat_max_cells (default Inf = no downsampling)
    if (is.finite(cellchat_max_cells) && ncol(seu) > cellchat_max_cells) {
      set.seed(42)
      seu <- seu[, sample(colnames(seu), cellchat_max_cells)]
      cat("  [", group_name, "] downsampled to", cellchat_max_cells, "cells\n")
    }
    # After subset, @cells mapping is out of sync (Seurat v5 bug #8407); directly JoinLayers would error
    # "Cannot add new cells with [[<-". Instead rebuild the RNA assay (keep counts+data).
    cc_counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
    cc_data   <- GetAssayData(seu, assay = "RNA", layer = "data")
    seu[["RNA"]] <- CreateAssay5Object(counts = cc_counts, data = cc_data)
    rm(cc_counts, cc_data); gc()
    data_input <- GetAssayData(seu, assay = "RNA", layer = "data")
    meta_input <- data.frame(labels = seu$cell_type, row.names = colnames(seu))

    cc <- createCellChat(object = data_input, meta = meta_input, group.by = "labels")
    cc@DB <- CellChatDB.human
    cc <- subsetData(cc)
    set.seed(1234)
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1,
                            population.size = TRUE)
    # Note: newer CellChat removed the min.groups parameter, keeping only min.cells
    cc <- filterCommunication(cc, min.cells = 10)
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
    cat("  [", group_name, "] CellChat done\n")
    cc
  }

  # 8.2 Build the five groups separately (Healthy / Infection_Control / COPD / Sepsis / Sepsis_Pneumonia)
  cc_list <- list()
  for (grp in core_groups) {
    cells_grp <- colnames(seu_core)[seu_core$group == grp]
    if (length(cells_grp) < 500) { cat("  [", grp, "] too few cells, skipping\n"); next }
    seu_grp <- subset(seu_core, cells = cells_grp)
    cc_list[[grp]] <- run_cellchat(seu_grp, grp)
  }

  # 8.3 Disease/infection vs healthy differential communication
  for (dis in c("Infection_Control", "COPD", "Sepsis", "Sepsis_Pneumonia")) {
    if (!all(c(dis, "Healthy") %in% names(cc_list))) next
    cc_merged <- mergeCellChat(list(cc_list[["Healthy"]], cc_list[[dis]]),
                               add.names = c("Healthy", dis))

    pdf(file.path(cellchat_dir, paste0("diff_", dis, "_vs_Healthy.pdf")),
        width = 16, height = 8)
    par(mfrow = c(1, 2))
    netVisual_diffInteraction(cc_merged, weight.scale = TRUE,
                              title.name = c("Healthy", dis))
    dev.off()

    p_rank <- rankNet(cc_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
    ggsave(file.path(cellchat_dir, paste0("rankNet_", dis, "_vs_Healthy.pdf")),
           plot = p_rank, width = 10, height = 12)
    cat("  Saved ", dis, " vs Healthy differential communication\n")
  }

  # 8.4 Comorbidity-shared signals: pathways jointly enhanced in COPD and Sepsis
  #      (via rankNet information gain: pathways up-regulated relative to Healthy in both)
  if (all(c("COPD", "Sepsis", "Healthy") %in% names(cc_list))) {
    get_up_pathways <- function(cc_dis, cc_healthy) {
      # Approximate with interaction-strength difference: pathways enhanced in dis relative to healthy
      w_dis <- cc_dis@net$weight
      w_hc  <- cc_healthy@net$weight
      common <- intersect(rownames(w_dis), rownames(w_hc))
      diff <- rowSums(w_dis[common, , drop = FALSE]) -
              rowSums(w_hc[common, , drop = FALSE])
      names(diff)[diff > 0]
    }
    up_copd <- get_up_pathways(cc_list[["COPD"]], cc_list[["Healthy"]])
    up_sep  <- get_up_pathways(cc_list[["Sepsis"]], cc_list[["Healthy"]])
    shared_pathways <- intersect(up_copd, up_sep)
    cat("\nSignal pathways jointly enhanced in COPD and Sepsis:\n")
    print(shared_pathways)
    writeLines(shared_pathways, file.path(cellchat_dir, "shared_pathways_COPD_Sepsis.txt"))
  }

  saveRDS(cc_list, file.path(cellchat_dir, "cellchat_list.rds"))
}

## =========================================================================
## Step 9: save
## =========================================================================

cat("\n===== Step 9: save =====\n")

saveRDS(merged,   file.path(out_dir, "path1_sepsis_copd_integrated.rds"))
saveRDS(seu_core, file.path(out_dir, "path1_core_4groups.rds"))

cat("Saved:\n")
cat("  - path1_sepsis_copd_integrated.rds (full integrated object)\n")
cat("  - path1_core_4groups.rds (four-group core object)\n")
cat("Path 1 analysis complete\n")
cat("Session Info:\n")
sessionInfo()
