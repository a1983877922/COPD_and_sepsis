###############################################################################
# Sepsis + COPD comorbidity analysis — Script 38: Figure 2 publication-grade unified-style re-render (run on server)
# File (current name): 38_fig2_unified_rerender.R (original file name was in Chinese, meaning "unified-style re-render")
#
# Purpose: Re-render 5–6 panels of Figure 2 with a unified color palette.
#
# [Important · v2 revert] v1 changed UMAP to "subsampling + remove title/label + custom theme",
#   but the output was worse than the path1_* originals from script 01. Per user request, v2 **reverts the plotting**:
#     - UMAP uses all cells (no subsampling) + Seurat default theme
#     - Keep the original title and cell-type labels
#     - Proportion plot reverts to the original cell-pooled stacked bars
#     - QC plot reverts to the original two side-by-side DimPlots
#   **Only colors are unified** to this script's CELL_COLORS / GROUP_COLORS / DATASET_COLORS.
#   The four revert switches (SHOW_LABEL / SHOW_TITLE / STRAT_SAMPLE / DONOR_LEVEL)
#   are at the top; change the corresponding switch to get back the v1 effect.
#
# Data: read cp3_annotated.rds directly (the integrated+annotated object saved by script 01, containing
#       cell_type / group / dataset / seurat_clusters / PCA / Harmony / UMAP).
#       This object is in the server's 01script directory; no need to re-run integration.
#
# Usage: (server, UTF-8 locale, same as scripts 01~36)
#   cd /media/desk16/ysx5991/sepsis_copd/01script
#   Rscript 38_fig2_unified_rerender.R
#
# Output (to out_dir): each panel is saved as both PDF (vector) + PNG (300 dpi):
#   Fig2a_UMAP_celltype.{pdf,png}     Main plot: by cell type, legend shows n cells
#   Fig2b_UMAP_group.{pdf,png}        5-group design (including Sepsis_Pneumonia)
#   Fig2c_UMAP_dataset.{pdf,png}      6 datasets subsampled proportionally (prevent GSE279452 visual dominance)
#   (disabled) Fig2d_marker_dotplot —— duplicates Fig2D from script 44; keep the 44 version
#   Fig2e_celltype_proportion.{pdf,png}  Cell-pooled stacked bars (default original style)
#                                          DONOR_LEVEL=TRUE switches to donor-level stacked columns + scatter overlay
#   FigS1_batch_PCA_vs_Harmony.{pdf,png}  Integrated QC (proportional subsampling, technical QC -> supplementary figure)
#
#   Extra: intermediate products before plotting, for locally redrawing/tweaking without cp3_annotated.rds
#   path38_fig2_plotdata.rds    Pure plotting data: per-panel data.frame + palette + levels
#   path38_fig2_ggplots.rds    ggplot objects (a/b/c/e/f), can be ggsave'd locally directly
#
# Unified style conventions (hard-coded at the top of this script; change once to apply globally):
#   - Cell-type palette CELL_COLORS (12 colors, assigned by sorted type name, consistent across all plots)
#   - 5-group palette GROUP_COLORS (Healthy/Infection_Control/COPD/Sepsis/Sepsis_Pneumonia)
#   - Dataset palette DATASET_COLORS (6 colors)
#   - All UMAPs: no title, no axes, no grid; legend on the right, unified font size
#   - All ggsave: white background, no border, 300 dpi PNG
###############################################################################

## ---- Auto-locate script directory (works locally and on server) ----
.this_file <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
if (.this_file == "" || !file.exists(.this_file)) {
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1])
}
SCRIPT_DIR <- dirname(.this_file)
config_file <- file.path(SCRIPT_DIR, "00_config.R")
if (!file.exists(config_file)) stop("cannot find 00_config.R: ", config_file)
source(config_file)
out_dir <- SCRIPT_DIR
cat("Output directory:", out_dir, "\n")

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork); library(dplyr)
})
## tidyr: needed for Panel e's complete() (impute missing cell types within donor);
## not forcibly required; if missing, give a clear message instead of "could not find function".
if (!requireNamespace("tidyr", quietly = TRUE)) {
  stop("Missing tidyr package (needed for Panel e donor-level imputation). Please run on the server first:\n",
       "  install.packages('tidyr')\n",
       "or install it in the conda/R environment and re-run this script.")
}

## =========================== Unified style constants =================================
CELL_COLORS <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
                 "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
                 "#374E55", "#8A5082")
GROUP_LEVELS <- c("Healthy", "Infection_Control", "COPD", "Sepsis", "Sepsis_Pneumonia")
GROUP_COLORS <- c(Healthy = "#3C5488", Infection_Control = "#91D1C2",
                  COPD = "#F39B7F", Sepsis = "#E64B35", Sepsis_Pneumonia = "#8491B4")
DATASET_COLORS <- c("#4DBBD5", "#E64B35", "#00A087", "#3C5488", "#F39B7F", "#8491B4")
UMAP_THEME <- theme_void() +
  theme(legend.position = "right",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8, face = "plain"),
        legend.key.size = unit(0.35, "cm"))
SEED <- 123
DPI <- 300

## ---- Revert switches: default = original style from script 01, only colors unified to this script's palette ----
## Note: the first version of script 38 changed UMAP to subsampling + remove title/label + custom theme,
##       but the result was worse than the path1_* originals from script 01. Now reverting per user request, keeping only color unification.
SHOW_LABEL   <- TRUE    # Panel a shows cell-type labels (original DimPlot label = TRUE)
SHOW_TITLE   <- TRUE    # Keep title on each panel (original has it; set FALSE when assembling in Illustrator)
STRAT_SAMPLE <- FALSE   # FALSE = all cells (original); TRUE = proportional subsampling (prevent large-dataset visual dominance)
DONOR_LEVEL  <- FALSE   # FALSE = cell-pooled stacked bars (original); TRUE = donor-level stacked columns + scatter overlay
## ==========================================================================

cp_file <- file.path(out_dir, "cp3_annotated.rds")
if (!file.exists(cp_file)) stop("Cannot find ", cp_file, " (run script 01 first to generate checkpoint)")
cat("Loading cp3_annotated.rds ...\n")
seu <- readRDS(cp_file)
cat("Cells:", ncol(seu), " | Fields:", paste(intersect(c("cell_type", "group",
    "dataset", "seurat_clusters"), colnames(seu@meta.data)), collapse = ", "), "\n")

## ---- 1) Cleaning: remove LowQuality, fix 5-group factor ----
meta <- seu@meta.data
if ("LowQuality" %in% meta$cell_type) {
  cat("Removing LowQuality:", sum(meta$cell_type == "LowQuality"), "cells\n")
  seu <- subset(seu, subset = cell_type != "LowQuality")
}
if (!"group" %in% colnames(seu@meta.data)) stop("object has no group field")
seu$group <- factor(seu$group, levels = GROUP_LEVELS)
grp_tab <- table(seu$group)
cat("Cell counts per 5 groups:\n"); print(grp_tab)

ct_present <- sort(unique(as.character(seu$cell_type)))
if (length(ct_present) > length(CELL_COLORS)) {
  ct_cols <- setNames(hue_pal()(length(ct_present)), ct_present)
} else {
  ct_cols <- setNames(CELL_COLORS[seq_along(ct_present)], ct_present)
}
ds_present <- sort(unique(as.character(seu$dataset)))
ds_cols <- setNames(DATASET_COLORS[seq_along(ds_present)], ds_present)

## ---- 1b) Diagnostics: missing UMAP coordinates, the source of DimPlot warning "Removing N cells missing data" ----
## If many cells have NA UMAP, it means reductions were not synced after subset / some cells did not participate in dim reduction,
## and the plot will be silently cropped. Here we proactively report counts to judge whether to rebuild UMAP in script 01.
for (red in intersect(c("umap", "pca", "harmony"), names(seu@reductions))) {
  emb <- Embeddings(seu, red)
  k_na <- sum(!complete.cases(emb))
  cat(sprintf("[diag] %s coords contain NA: %d / %d cells\n", red, k_na, nrow(emb)))
  if (k_na > 0 && k_na < 200 && red == "umap") {
    cat("  Distribution by group:\n"); print(table(seu$group[!complete.cases(emb)]))
  }
}

## ---- Utility: proportional subsampling (dataset or group dimension), for display ----
strat_sample <- function(cells_vec, k = 1) {
  set.seed(SEED)
  unlist(lapply(unique(cells_vec), function(lv) {
    idx <- which(cells_vec == lv)
    if (length(idx) <= k) idx else sample(idx, k)
  }))
}

## =========================================================================
## Panel a: UMAP by cell type
##   [Revert to original style] all cells no subsampling + show cell-type label + keep title +
##   Seurat default theme; only colors unified to ct_cols.
## =========================================================================
if (STRAT_SAMPLE) {
  set.seed(SEED)
  disp_cells <- if (ncol(seu) > 3e5) sample(colnames(seu), 3e5) else colnames(seu)
  seu_d <- subset(seu, cells = disp_cells)
} else {
  seu_d <- seu                      # original: all cells
}
p_a <- DimPlot(seu_d, reduction = "umap", group.by = "cell_type",
               pt.size = 0.05, label = SHOW_LABEL, label.size = 3,
               cols = ct_cols) +
  labs(title = if (SHOW_TITLE) "Cell Type Annotation" else NULL)

## =========================================================================
## Panel b: UMAP by group (5 groups, including Sepsis_Pneumonia)
##   [Revert to original style] all cells + default theme + original title; colors unified to GROUP_COLORS
##   (original triad plot's group panel has no label; keep label = FALSE here)
## =========================================================================
if (STRAT_SAMPLE) {
  grp_sub <- strat_sample(seu$group, 2e4)
  seu_g <- subset(seu, cells = grp_sub)
} else {
  seu_g <- seu
}
p_b <- DimPlot(seu_g, reduction = "umap", group.by = "group",
               pt.size = 0.05, label = FALSE, cols = GROUP_COLORS) +
  labs(title = if (SHOW_TITLE) "UMAP by Group (5 groups)" else NULL)

## =========================================================================
## Panel c: UMAP by dataset
##   [Revert to original style] all cells + default theme + original title; colors unified to ds_cols
## =========================================================================
if (STRAT_SAMPLE) {
  ds_sub <- strat_sample(seu$dataset, 8e3)
  seu_c <- subset(seu, cells = ds_sub)
} else {
  seu_c <- seu
}
p_c <- DimPlot(seu_c, reduction = "umap", group.by = "dataset",
               pt.size = 0.05, label = FALSE, cols = ds_cols) +
  labs(title = if (SHOW_TITLE) "UMAP by Dataset" else NULL)

## =========================================================================
## Panel f: batch correction QC —— PCA (before) vs Harmony (after)
##   [Revert to original style] two side-by-side DimPlots + respective titles, colors unified to ds_cols.
##   df_qc below is still generated for path38_fig2_plotdata.rds to store data for local redraw.
## =========================================================================
emb_pca <- Embeddings(seu_c, "pca")[, 1:2]
emb_har <- Embeddings(seu_c, "harmony")[, 1:2]
df_qc <- data.frame(
  x = c(emb_pca[, 1], emb_har[, 1]),
  y = c(emb_pca[, 2], emb_har[, 2]),
  dataset = rep(seu_c$dataset, 2),
  step = rep(c("PCA (before Harmony)", "Harmony (after)"), each = nrow(emb_pca))
)
p_before <- DimPlot(seu_c, reduction = "pca", group.by = "dataset",
                    pt.size = 0.05, cols = ds_cols) +
  labs(title = if (SHOW_TITLE) "PCA (before correction)" else NULL)
p_after  <- DimPlot(seu_c, reduction = "harmony", group.by = "dataset",
                    pt.size = 0.05, cols = ds_cols) +
  labs(title = if (SHOW_TITLE) "Harmony (after correction)" else NULL)
p_f <- (p_before | p_after) + plot_layout(guides = "collect")

## =========================================================================
## Panel d: Dotplot marker genes —— [DISABLED]
##   Same type of plot as Fig2D_marker_dotplot from script 44, duplicated.
##   Keep the 44 version (shares FIG2_CELLTYPE_ORDER 8-type list and order with Fig2F/G,
##   and keeps the Mast row). This script no longer outputs this panel.
## =========================================================================
## ct_order <- ct_present
## if (exists("blood_markers", inherits = FALSE)) {
##   gm <- unlist(lapply(blood_markers, function(v) intersect(v, rownames(seu))))
## } else {
##   gm <- c("CD3D","CD3E","CD8A","IL7R","MS4A1","CD79A","GNLY","NKG7",
##           "CD14","LYZ","S100A8","FCER1A","CST3","CLEC9A","FCGR3B","CSF3R",
##           "PPBP","PF4","HBB","CPA3")
##   gm <- intersect(gm, rownames(seu))
## }
## gm <- gm[!duplicated(gm)]
## seu_d$cell_type <- factor(seu_d$cell_type, levels = ct_order)
## p_d <- DotPlot(seu_d, features = gm, group.by = "cell_type",
##                cols = c("lightgrey", "#185FA5")) +
##   RotatedAxis() +
##   labs(x = NULL, y = NULL) +
##   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
##         axis.text.y = element_text(size = 8),
##         legend.position = "right")

## =========================================================================
## Panel e: cell-type proportion —— donor-level stacked columns + donor scatter overlay if donor field exists, otherwise stacked bars
## =========================================================================
donor_field <- NULL
for (f in c("donor_id", "donor", "subject", "subject_id", "Subject_Identity", "sample")) {
  if (f %in% colnames(seu@meta.data) &&
      !all(is.na(seu@meta.data[[f]])) && !all(seu@meta.data[[f]] == "")) {
    donor_field <- f; break
  }
}
if (DONOR_LEVEL && !is.null(donor_field)) {
  cat("Detected donor field:", donor_field, " -> donor-level stacked columns + donor scatter overlay\n")
  ## Donor-level long table: per-donor per-cell-type frac + stacked segment midpoint y position (for scatter)
  df_dot <- seu@meta.data %>%
    filter(group %in% GROUP_LEVELS, cell_type %in% ct_present) %>%
    count(group, cell_type, .data[[donor_field]], name = "n") %>%
    group_by(group, .data[[donor_field]]) %>%
    tidyr::complete(cell_type = ct_present, fill = list(n = 0)) %>%  # impute missing cell types within donor = 0
    mutate(frac = n / sum(n)) %>%
    ungroup() %>%
    mutate(cell_type = factor(cell_type, levels = rev(ct_present))) %>%
    arrange(group, .data[[donor_field]], cell_type) %>%
    group_by(group, .data[[donor_field]]) %>%
    mutate(y_pos = cumsum(frac) - frac / 2) %>%                # segment midpoint, falls within the corresponding color block
    ungroup()
  n_don <- seu@meta.data %>% count(group, .data[[donor_field]]) %>% count(group, name = "n_donors")
  ## Within-group donor mean ± SE (stacked columns)
  df_don <- df_dot %>%
    group_by(group, cell_type) %>%
    summarise(mu = mean(frac),
              se = { s <- sd(frac); if (is.na(s)) 0 else s / sqrt(n()) },
              nd = n(), .groups = "drop") %>%
    left_join(n_don, by = "group")
  df_don$cell_type <- factor(as.character(df_don$cell_type), levels = rev(ct_present))
  lab_map <- setNames(paste0(n_don$group, "  (n=", n_don$n_donors, ")"), n_don$group)
  df_don$grp_lab <- unname(lab_map[as.character(df_don$group)])
  df_dot$grp_lab <- unname(lab_map[as.character(df_dot$group)])
  p_e <- ggplot() +
    geom_col(data = df_don, aes(grp_lab, mu, fill = cell_type), width = 0.7) +
    ## No error bars: (1) scatter already shows per-donor variation, error bars would be redundant;
    ## (2) if error bars were drawn at y=mu they would misalign with stacked blocks (cumulative height), visually "floating".
    ## se is still kept in df_don; add locally via path38_fig2_plotdata.rds when needed.
    geom_point(data = df_dot, aes(grp_lab, y_pos),
               shape = 21, fill = "black", color = "grey25",
               alpha = 0.30, size = 1.2,
               position = position_jitter(width = 0.12)) +
    scale_fill_manual(values = ct_cols, name = "Cell type") +
    labs(y = "Donor-level mean fraction", x = NULL,
         title = if (SHOW_TITLE) "Cell Type Composition by Group (donor level)" else NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank())
} else {
  ## [Revert to original style] cell-pooled stacked bars (consistent with path1_celltype_proportion from script 01),
  ## only colors unified to ct_cols; set DONOR_LEVEL to TRUE for donor-level stats to switch back.
  cat("-> cell-pooled stacked bars (original style; DONOR_LEVEL=TRUE switches to donor-level scatter version)\n")
  prop_df <- as.data.frame(table(seu$cell_type, seu$group)) %>%
    rename(CellType = Var1, Group = Var2, Count = Freq) %>%
    group_by(Group) %>% mutate(Frac = Count / sum(Count))
  prop_df$CellType <- factor(prop_df$CellType, levels = rev(ct_present))
  p_e <- ggplot(prop_df, aes(Group, Frac, fill = CellType)) +
    geom_col(position = "fill", width = 0.7) +
    scale_fill_manual(values = ct_cols, name = "Cell type") +
    labs(y = "Fraction", x = NULL,
         title = if (SHOW_TITLE) "Cell Type Composition by Group" else NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank())
}

## =========================================================================
## Save (unified 300 dpi)
## =========================================================================
save2 <- function(p, base, w, h) {
  ggsave(file.path(out_dir, paste0(base, ".pdf")), p, width = w, height = h)
  ggsave(file.path(out_dir, paste0(base, ".png")), p, width = w, height = h,
         dpi = DPI, device = grDevices::png)
  cat("Saved:", base, ".{pdf,png}\n")
}

## =========================================================================
## Pre-plot save: plotting data + ggplot objects
##   Purpose: locally redraw/tweak each panel without cp3_annotated.rds (the ~920k-cell large object),
##            no need to re-run on the server. The two switches can be turned off as needed.
## =========================================================================
SAVE_PLOTDATA <- TRUE   # light: pure data.frame + palette, can redraw from scratch with ggplot locally
SAVE_GGPLOTS  <- TRUE   # heavy: ggplot object includes data layer, can ggsave directly or +theme() tweak locally

umap_export <- function(p, nm) {
  d <- tryCatch(p$data, error = function(e) NULL)
  if (is.null(d) || !is.data.frame(d)) return(NULL)
  d$panel <- nm
  d
}
mb <- function(f) sprintf("%.1f MB", file.size(f) / 1024^2)

if (SAVE_PLOTDATA) {
  plotdata <- list(
    meta = list(
      script = "38_fig2_unified_rerender.R",
      created = Sys.time(),
      n_cells_total = ncol(seu),
      n_cells_display = list(
        a_celltype = if (is.null(p_a$data)) NA else nrow(p_a$data),
        b_group    = if (is.null(p_b$data)) NA else nrow(p_b$data),
        c_dataset  = if (is.null(p_c$data)) NA else nrow(p_c$data)
      ),
      group_cell_n = as.list(as.integer(table(seu$group))),
      group_cell_n_names = names(table(seu$group)),
      donor_field = if (exists("donor_field")) donor_field else NULL,
      celltype_levels = ct_present,
      dataset_levels = ds_present
    ),
    styles = list(
      cell_colors = ct_cols, group_colors = GROUP_COLORS,
      dataset_colors = ds_cols, group_levels = GROUP_LEVELS,
      umap_theme = UMAP_THEME, dpi = DPI
    ),
    umap_a = umap_export(p_a, "a_celltype"),
    umap_b = umap_export(p_b, "b_group"),
    umap_c = umap_export(p_c, "c_dataset"),
    qc_f = df_qc,
    prop_donor_mean = if (exists("df_don")) df_don else NULL,
    prop_donor_points = if (exists("df_dot")) df_dot else NULL,
    prop_pooled = if (exists("prop_df")) prop_df else NULL
  )
  f_pd <- file.path(out_dir, "path38_fig2_plotdata.rds")
  saveRDS(plotdata, f_pd)
  cat("Saved plotting data:", basename(f_pd), "(", mb(f_pd), ")\n")
  cat("  Local usage: pd <- readRDS('path38_fig2_plotdata.rds'); ",
      "str(pd, max.level = 2)\n")
}

if (SAVE_GGPLOTS) {
  gg_list <- list(a_celltype = p_a, b_group = p_b, c_dataset = p_c,
                  e_proportion = p_e, f_batch_QC = p_f)
  f_gg <- file.path(out_dir, "path38_fig2_ggplots.rds")
  saveRDS(gg_list, f_gg)
  cat("Saved ggplot objects:", basename(f_gg), "(", mb(f_gg), ")\n")
  cat("  Local usage: g <- readRDS('path38_fig2_ggplots.rds'); ",
      "ggsave('p.pdf', g$e_proportion + theme_bw(), w = 9, h = 6)\n")
}
## Sizes follow the original ratio from script 01 (a: 10x8, triad single ~8x7, QC double 16x8)
save2(p_a, "Fig2a_UMAP_celltype",   w = 10, h = 8)
save2(p_b, "Fig2b_UMAP_group",      w = 8,  h = 7)
save2(p_c, "Fig2c_UMAP_dataset",    w = 8,  h = 7)
## save2(p_d, "Fig2d_marker_dotplot",  w = 12, h = 6)   # disabled, see Panel d comment
save2(p_e, "Fig2e_celltype_proportion", w = 9, h = 6)
## Integrated QC is technical QC; the main 7-panel figure (A-G) has no room, so it becomes supplementary figure S1.
## To put it back in the main figure, rename the following line's file to "Fig2f_batch_PCA_vs_Harmony".
save2(p_f, "FigS1_batch_PCA_vs_Harmony", w = 16, h = 8)

cat("\n===== Script 38 finished =====\n")
cat("Note: if the cell_type legend n cells comes from the 300k subsampled subset rather than all cells, you can raise the 3e5 cap in Panel a or compute the full n from meta.\n")
