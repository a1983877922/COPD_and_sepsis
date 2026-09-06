# =============================================================================
# 44_Fig2_G_marker_dotplot.R
# Fig2 panel D (DotPlot version):
#   Note: the panel letters were re-ordered — cell identity verification (DotPlot) now comes
#   before the composition analysis, so G became D.
#   Reference, Fig2.G of the MASH paper:
#   Dot plot showing the expression of canonical marker genes in major cell types,
#   scaled by expression percentage (dot size) and average expression level (color).
#   x = major cell types, y = marker genes (grouped by lineage), size = percentage of
#   expressing cells, color = average expression.
# Criteria: consistent with the other Fig2 panels — read cp3_annotated.rds, keep only the five
#   core groups, drop LowQuality.
# Run: server (requires Seurat)
#   Rscript 44_Fig2_G_marker_dotplot.R [OUT_DIR]
# Outputs: Fig2D_marker_dotplot.pdf / .png
#          path44_marker_dotplot.csv
#          path44_marker_dotplot_plotdata.rds
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
script_dir <- tryCatch(dirname(normalizePath(sys.frames()[[1]]$ofile)),
                       error = function(e) getwd())
OUT_DIR <- ifelse(length(args) >= 1 && nzchar(args[1]), args[1], script_dir)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# shares the cell type order with Fig2F (single source of truth, avoids drift)
source(file.path(script_dir, "fig2_major_celltypes.R"))

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

## ---- Load the integrated object (same source as Fig2: cp3_annotated.rds) ----
obj_file <- file.path(OUT_DIR, "cp3_annotated.rds")
if (!file.exists(obj_file)) obj_file <- file.path(OUT_DIR, "path1_sepsis_copd_integrated.rds")
stopifnot(file.exists(obj_file))
seu <- readRDS(obj_file)
cat("Loaded:", basename(obj_file), "\n")

core_groups <- c("Healthy", "Infection_Control", "COPD", "Sepsis", "Sepsis_Pneumonia")
keep <- !is.na(seu@meta.data$cell_type) & seu@meta.data$cell_type != "LowQuality" &
        !is.na(seu@meta.data$group) & seu@meta.data$group %in% core_groups
seu <- subset(seu, cells = rownames(seu@meta.data)[keep])
cat("Cells retained:", ncol(seu), "\n")

## ---- Cell type list shared with Fig2F (single source of truth = the 8 major cell types in the Ro/e CSV) ----
roe_file <- file.path(OUT_DIR, "path1_roe_celltype_by_group.csv")
stopifnot(file.exists(roe_file))
roe_types <- rownames(read.csv(roe_file, row.names = 1, check.names = FALSE))

## ---- Canonical marker genes (1-3 per lineage, consistent with the blood cell markers) ----
marker_list <- list(
  T_cell     = c("CD3D", "CD8A", "IL7R"),
  B_cell     = c("MS4A1", "CD79A"),
  NK         = c("GNLY", "NKG7"),
  Monocyte   = c("LYZ", "S100A8", "CD14"),
  DC         = c("FCER1A", "ITGAX"),
  Neutrophil = c("FCGR3B"),
  Platelet   = c("PPBP"),
  RBC        = c("HBB"),
  Mast       = c("TPSAB1")
)

LINEAGE_COLS <- c(
  T_cell     = "#E64B35",
  B_cell     = "#4DBBD5",
  NK         = "#00A087",
  Monocyte   = "#3C5488",
  DC         = "#F39B7F",
  Neutrophil = "#8491B4",
  Platelet   = "#F0E442",
  RBC        = "#B09C85",
  Mast       = "#BB79C8"
)

## ---- Keep only cell types shared by F/G (the 8 types in the Ro/e CSV ∩ the cell types actually present in the object) ----
all_celltypes <- sort(unique(as.character(seu$cell_type)))
master_types  <- intersect(roe_types, all_celltypes)
marker_list   <- marker_list[names(marker_list) %in% master_types]
markers <- unlist(marker_list, use.names = FALSE)
lineage_vec <- rep(names(marker_list), times = sapply(marker_list, length))
names(lineage_vec) <- markers

avail_markers <- intersect(markers, rownames(seu))
lineage_vec <- lineage_vec[avail_markers]
cat("Available marker genes:", length(avail_markers), "/", length(markers), "\n")

## ---- Compute per-cell-type average expression (log-normalized data layer) and expression percentage ----
## Note (Seurat v5 pitfall):
##   row slicing of counts[g, ] returns an "unnamed" vector, so names(expr_count) has length 0
##   and cannot be used to obtain cell names (otherwise data.frame errors with
##   "differing number of rows: 0, N").
##   Always take cell names from colnames(seu).
suppressMessages(library(Matrix))

get_layer <- function(obj, layer) {
  tryCatch(GetAssayData(obj, assay = "RNA", layer = layer),
           error = function(e) GetAssayData(obj, assay = "RNA", slot = layer))
}

counts <- get_layer(seu, "counts")
data_l <- get_layer(seu, "data")

cells_sub <- colnames(seu)
ct_vec    <- as.character(seu@meta.data[cells_sub, "cell_type"])
keep_ct   <- !is.na(ct_vec) & ct_vec %in% master_types
counts <- counts[, cells_sub[keep_ct], drop = FALSE]
data_l <- data_l[, cells_sub[keep_ct], drop = FALSE]
ct_vec <- ct_vec[keep_ct]
cat("Cells used for G:", length(ct_vec), "\n")

## x-axis order: exactly the same as Fig2F
celltype_order <- FIG2_CELLTYPE_ORDER[FIG2_CELLTYPE_ORDER %in% master_types]
celltype_order <- celltype_order[celltype_order %in% ct_vec]

genes   <- avail_markers
avg_mat <- matrix(NA_real_, nrow = length(genes), ncol = length(celltype_order),
                  dimnames = list(genes, celltype_order))
pct_mat <- avg_mat

for (ct in celltype_order) {
  idx <- which(ct_vec == ct)
  if (!length(idx)) next
  d_sub <- data_l[genes, idx, drop = FALSE]
  avg_mat[, ct] <- as.numeric(Matrix::rowMeans(d_sub))
  c_sub <- counts[genes, idx, drop = FALSE]
  if (methods::is(c_sub, "sparseMatrix")) {
    c_sub@x <- as.numeric(c_sub@x > 0)          # expressed or not -> 0/1
  } else {
    c_sub <- (as.matrix(c_sub) > 0) * 1
  }
  pct_mat[, ct] <- as.numeric(Matrix::rowMeans(c_sub)) * 100
}

## Color scale convention (switchable):
##   DELOG = FALSE -> y = the direct log-normalized mean (recommended: smaller color-scale
##                    range, low expression is still visible)
##   DELOG = TRUE  -> y = expm1(mean(log1p(x))), the same convention as Seurat's own DotPlot
DELOG <- FALSE
if (DELOG) avg_mat <- expm1(avg_mat)

## Long table ----
res <- data.frame(
  Gene    = rep(rownames(avg_mat), times = ncol(avg_mat)),
  cell_type = rep(colnames(avg_mat), each  = nrow(avg_mat)),
  Lineage = lineage_vec[rep(rownames(avg_mat), times = ncol(avg_mat))],
  AvgExpr = as.numeric(avg_mat),
  PctExpr = as.numeric(pct_mat),
  stringsAsFactors = FALSE
)

## y order: by lineage order, then by gene name
res$cell_type <- factor(as.character(res$cell_type), levels = celltype_order)

lineage_level_order <- names(marker_list)
gene_order <- data.frame(Gene = avail_markers, Lineage = lineage_vec[avail_markers],
                         stringsAsFactors = FALSE) %>%
  mutate(Lineage = factor(Lineage, levels = lineage_level_order)) %>%
  arrange(Lineage, Gene) %>%
  pull(Gene)
res$Gene <- factor(res$Gene, levels = gene_order)

## Color vector (in y-axis gene order)
y_col_vec <- setNames(LINEAGE_COLS[lineage_vec[gene_order]], gene_order)

## ---- Save pre-plot data ----
plot_data <- list(
  res = res,
  core_groups = core_groups,
  LINEAGE_COLS = LINEAGE_COLS,
  marker_list = marker_list,
  gene_order = gene_order,
  celltype_order = levels(res$cell_type)
)
saveRDS(plot_data, file.path(OUT_DIR, "path44_marker_dotplot_plotdata.rds"))
write.csv(res, file.path(OUT_DIR, "path44_marker_dotplot.csv"), row.names = FALSE)
cat("Saved: path44_marker_dotplot.csv + path44_marker_dotplot_plotdata.rds\n")

## ---- Plot: DotPlot ----
dot_scale <- c(0, 25, 50, 75)
p_g <- ggplot(res, aes(x = cell_type, y = Gene, size = PctExpr, color = AvgExpr)) +
  geom_point(stroke = 0.2, shape = 16) +
  scale_color_gradientn(
    colors = c("#FFF5F0", "#FEE0D2", "#FB6A4A", "#CB181D", "#67000D"),
    name = "Mean expression\n(log-normalized)",
    limits = c(0, quantile(res$AvgExpr[res$AvgExpr > 0], 0.99))
  ) +
  scale_size_area(
    breaks = dot_scale,
    labels = paste0(dot_scale, "%"),
    name = "% cells",
    max_size = 6.5,
    limits = c(0, 100)
  ) +
  labs(
    x = NULL, y = NULL,
    title = "G  Canonical marker genes across major cell types",
    subtitle = "dot size = % expressing cells; color = average expression"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10, color = "grey20"),
    axis.text.y = element_text(size = 9, color = y_col_vec),
    axis.line = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(OUT_DIR, "Fig2D_marker_dotplot.pdf"), p_g,
       width = 7.5, height = 5.2)
cat("Saved: Fig2D_marker_dotplot.pdf\n")

tryCatch({
  grDevices::png(file.path(OUT_DIR, "Fig2D_marker_dotplot.png"),
                 width = 7.5, height = 5.2, units = "in", res = 300)
  print(p_g); grDevices::dev.off()
  cat("Saved: Fig2D_marker_dotplot.png\n")
}, error = function(e) cat("PNG writing skipped (regenerate from the PDF with PyMuPDF):", conditionMessage(e), "\n"))

cat("\n===== script 44 DotPlot version done =====\n")
