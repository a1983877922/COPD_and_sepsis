###############################################################################
# 05_reference_mapping.R — lung-to-blood reference mapping (label transfer)
#
# Purpose: use lung myeloid as reference and map blood monocytes onto it, to answer
#   "how many blood monocytes resemble lung-derived macrophages".
#   If the proportion of COPD/Sepsis blood monocytes predicted as "alveolar/macrophage" is
#   higher than in Healthy, this supports the causality of "imprint transmitted from lung to blood".
#
# Data:
#   Reference: path2_copd_lung.rds (lung myeloid, with built-in Manuscript_Identity annotation)
#   Query: path1_sepsis_copd_integrated.rds (blood monocytes, cell_type + group)
#
# Output:
#   path5_prediction_by_group.csv    (blood monocytes predicted lung myeloid label x group)
#   path5_lung_mac_prop.csv          (lung-derived macrophage proportion x group)
#   path5_prediction_barplot.pdf     (predicted label stacked plot)
#   path5_blood_mono_predicted.rds   (query object with prediction results)
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
cat("Lung-to-blood reference mapping (label transfer)\n")
cat("==============================================================\n")

## =========================================================================
## Step 1: build lung myeloid reference
## =========================================================================

cat("\n===== Step 1: build lung myeloid reference =====\n")

lung_seu <- readRDS(file.path(out_dir, "path2_copd_lung.rds"))
lung_myeloid <- subset(lung_seu, subset = cell_category == "Myeloid")
rm(lung_seu); gc()
cat("Lung myeloid:", ncol(lung_myeloid), "cells\n")

# Simplify labels: merge cMonocyte/ncMonocyte into Monocyte, cDC into DC,
# keep Macrophage_Alveolar / Macrophage / Mast original names
lung_myeloid$ref_label <- lung_myeloid$cell_type
lung_myeloid$ref_label[lung_myeloid$cell_type %in% c("cMonocyte", "ncMonocyte")] <- "Monocyte"
lung_myeloid$ref_label[lung_myeloid$cell_type %in% c("cDC1", "cDC2", "pDC", "DC_Mature", "DC_Langerhans")] <- "DC"
cat("Reference labels:\n")
print(table(lung_myeloid$ref_label))

# Rebuild assay (fix Seurat v5 bug with out-of-sync @cells mapping after subset),
# then re-run NormalizeData/ScaleData/PCA
mye_counts <- GetAssayData(lung_myeloid, assay = "RNA", layer = "counts")
lung_myeloid[["RNA"]] <- CreateAssay5Object(counts = mye_counts)
rm(mye_counts); gc()
lung_myeloid <- NormalizeData(lung_myeloid, verbose = FALSE)
lung_myeloid <- FindVariableFeatures(lung_myeloid, selection.method = "vst",
                                     nfeatures = 2000, verbose = FALSE)
lung_myeloid <- ScaleData(lung_myeloid, verbose = FALSE)
lung_myeloid <- RunPCA(lung_myeloid, npcs = 30, verbose = FALSE)

## =========================================================================
## Step 2: build blood monocyte query
## =========================================================================

cat("\n===== Step 2: build blood monocyte query =====\n")

blood_seu <- readRDS(file.path(out_dir, "path1_sepsis_copd_integrated.rds"))
blood_mono <- subset(blood_seu, subset = cell_type == "Monocyte" &
                       group %in% c("Healthy", "COPD", "Sepsis"))
rm(blood_seu); gc()
cat("Blood monocytes:", ncol(blood_mono), "cells\n")
print(table(blood_mono$group))

# Rebuild assay + process query using reference variable features
b_counts <- GetAssayData(blood_mono, assay = "RNA", layer = "counts")
blood_mono[["RNA"]] <- CreateAssay5Object(counts = b_counts)
rm(b_counts); gc()
blood_mono <- NormalizeData(blood_mono, verbose = FALSE)
# Use lung myeloid HVG, but first take the intersection (blood monocytes may lack some genes; direct indexing would go out of bounds)
feat_use <- intersect(VariableFeatures(lung_myeloid), rownames(blood_mono))
cat("Shared HVG:", length(feat_use), "/", length(VariableFeatures(lung_myeloid)), "\n")
blood_mono <- ScaleData(blood_mono, features = feat_use, verbose = FALSE)
blood_mono <- RunPCA(blood_mono, features = feat_use, npcs = 30, verbose = FALSE)

## =========================================================================
## Step 3: reference mapping (label transfer)
## =========================================================================

cat("\n===== Step 3: reference mapping (label transfer, may take 10-20 minutes) =====\n")

anchors <- FindTransferAnchors(reference = lung_myeloid, query = blood_mono,
                               dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = lung_myeloid$ref_label,
                            dims = 1:30)
blood_mono <- AddMetaData(blood_mono, metadata = predictions)
rm(anchors, predictions, lung_myeloid); gc()

cat("Prediction complete, prediction score distribution:\n")
print(summary(blood_mono$prediction.score.max))

## =========================================================================
## Step 4: statistics + visualization
## =========================================================================

cat("\n===== Step 4: statistics =====\n")

# 4.1 Predicted label x group contingency table
pred_tab <- table(blood_mono$group, blood_mono$predicted.id)
cat("\nPredicted label x group:\n")
print(pred_tab)
write.csv(pred_tab, file.path(out_dir, "path5_prediction_by_group.csv"))

# 4.2 Lung-derived macrophage proportion (Macrophage_Alveolar + Macrophage)
blood_mono$is_lung_mac <- blood_mono$predicted.id %in%
  c("Macrophage_Alveolar", "Macrophage")
lung_mac_tab <- table(blood_mono$group, blood_mono$is_lung_mac)
lung_mac_frac <- prop.table(lung_mac_tab, margin = 1)
cat("\nLung-derived macrophage proportion (predicted as alveolar/macrophage per group):\n")
print(round(lung_mac_frac, 4))
write.csv(lung_mac_frac, file.path(out_dir, "path5_lung_mac_prop.csv"))

# 4.3 Inter-group chi-square test
chi <- chisq.test(lung_mac_tab)
cat(sprintf("\nChi-square test of lung-derived macrophage proportion across three groups (Healthy/COPD/Sepsis): p = %.3e\n",
            chi$p.value))

# Pairwise Fisher test
pairs <- list(c("COPD", "Healthy"), c("Sepsis", "Healthy"), c("Sepsis", "COPD"))
for (pr in pairs) {
  sub <- lung_mac_tab[pr, , drop = FALSE]
  f <- fisher.test(sub)
  cat(sprintf("  %s vs %s: p = %.3e\n", pr[1], pr[2], f$p.value))
}

# 4.4 Predicted label stacked plot
pred_df <- as.data.frame(table(blood_mono$group, blood_mono$predicted.id))
colnames(pred_df) <- c("Group", "Predicted", "Count")
pred_df <- pred_df %>%
  group_by(Group) %>%
  mutate(Frac = Count / sum(Count)) %>%
  ungroup()
pred_df$Group <- factor(pred_df$Group, levels = c("Healthy", "COPD", "Sepsis"))

p_stack <- ggplot(pred_df, aes(x = Group, y = Frac, fill = Predicted)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "Fraction", x = NULL,
       title = "Blood monocytes mapped to lung myeloid types") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
ggsave(file.path(out_dir, "path5_prediction_barplot.pdf"),
       plot = p_stack, width = 8, height = 6)

# 4.5 Lung-derived macrophage proportion bar chart
mac_df <- as.data.frame(lung_mac_frac)
colnames(mac_df) <- c("Group", "is_lung_mac", "Frac")
mac_df <- mac_df[mac_df$is_lung_mac == "TRUE", ]
mac_df$Group <- factor(mac_df$Group, levels = c("Healthy", "COPD", "Sepsis"))
p_mac <- ggplot(mac_df, aes(x = Group, y = Frac, fill = Group)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("Healthy" = "#85B7EB", "COPD" = "#F0997B",
                               "Sepsis" = "#E24B4A")) +
  labs(y = "Lung-macrophage-like fraction", x = NULL,
       title = "Lung-macrophage-like cells in blood monocytes") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "path5_lung_mac_prop.pdf"),
       plot = p_mac, width = 6, height = 5)

# 4.6 Save query object (with prediction results)
saveRDS(blood_mono, file.path(out_dir, "path5_blood_mono_predicted.rds"))

cat("\nLung-to-blood reference mapping complete\n")
cat("Output:\n")
cat("  - path5_prediction_by_group.csv (predicted label x group)\n")
cat("  - path5_lung_mac_prop.csv (lung-derived macrophage proportion)\n")
cat("  - path5_prediction_barplot.pdf (stacked plot)\n")
cat("  - path5_lung_mac_prop.pdf (lung-derived macrophage proportion plot)\n")
cat("  - path5_blood_mono_predicted.rds (object with prediction results)\n")
