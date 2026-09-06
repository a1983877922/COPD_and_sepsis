#!/usr/bin/env Rscript
# =============================================================================
# 46_fig4c_socs3_celltype.R
# Fig4 Panel C: cell-type origin of SOCS3 and severity gradient (donor-level statistics)
# =============================================================================
# Panel intent (answers two questions following Fig4A/B):
#   C1  SOCS3 dot plot across 8 major cell types x 5 groups
#       size = fraction of expressing cells, color = mean expression -> "which cell types drive SOCS3"
#   C2  **Donor-level** five-group distribution of SOCS3 within the highest-expressing
#       cell type (boxplot + per-donor points)
#       Wilcoxon vs Healthy p-value annotations           -> "increases with severity, and the statistical unit is the donor"
#
# Convention consistency:
#   - Cell-type list/order reuses FIG2_CELLTYPE_ORDER from fig2_major_celltypes.R (consistent with Fig2F/G)
#   - Five-group colors reuse GROUP_COLS from Fig2E/F/G
#   - Reads the same integrated object cp3_annotated.rds, drops LowQuality, keeps only the five core groups
#
# Run: server (needs Seurat)
#   Rscript 46_fig4c_socs3_celltype.R [OUT_DIR]
#
# Output:
#   Fig4C_SOCS3_celltype_specificity.pdf / .png
#   Fig4C_SOCS3_celltype.csv        cell type x group aggregation (avg / pct / n_cells / n_donors)
#   Fig4C_SOCS3_donor.csv           donor-level means (for C2 and local rechecking)
#   Fig4C_SOCS3_donor_stats.txt     C2 Kruskal / Wilcoxon statistics
#   path46_fig4c_plotdata.rds       pre-plot data (re-plot locally without Seurat)
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
script_dir <- tryCatch(dirname(normalizePath(sys.frames()[[1]]$ofile)),
                       error = function(e) getwd())
OUT_DIR <- ifelse(length(args) >= 1 && nzchar(args[1]), args[1], script_dir)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## Cell-type order shared with Fig2/Fig4 (single source of truth)
source(file.path(script_dir, "fig2_major_celltypes.R"))

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr)
})
suppressMessages(library(Matrix))
has_pw <- requireNamespace("patchwork", quietly = TRUE)
if (!has_pw) {
  cat("!! patchwork not installed; C1/C2 will be output as two separate figures\n")
} else {
  suppressPackageStartupMessages(library(patchwork))
}

## =========================================================================
## Tunable parameters
## =========================================================================
GENES        <- "SOCS3"                       # can be extended to c("SOCS3","SOCS1","ISG15") etc.
GROUP_LEVELS <- c("Healthy", "Infection_Control", "COPD", "Sepsis", "Sepsis_Pneumonia")
GROUP_LABEL  <- c(Healthy = "Healthy", Infection_Control = "IC", COPD = "COPD",
                  Sepsis = "Sepsis", Sepsis_Pneumonia = "Sepsis+PNA")
GROUP_COLS   <- c(Healthy = "#3C5488", Infection_Control = "#91D1C2",
                  COPD = "#00A087", Sepsis = "#E64B35", Sepsis_Pneumonia = "#F39B7F")

## =========================================================================
## Step 1: load the object
## =========================================================================
obj_file <- file.path(OUT_DIR, "cp3_annotated.rds")
if (!file.exists(obj_file)) obj_file <- file.path(OUT_DIR, "path1_sepsis_copd_integrated.rds")
stopifnot(file.exists(obj_file))
seu <- readRDS(obj_file)
cat("loaded:", basename(obj_file), "\n")

keep <- !is.na(seu@meta.data$cell_type) & seu@meta.data$cell_type != "LowQuality" &
        !is.na(seu@meta.data$group)     & seu@meta.data$group %in% GROUP_LEVELS
seu <- subset(seu, cells = rownames(seu@meta.data)[keep])
cat("cells kept:", ncol(seu), "\n")

## Cell-type list = 8 types from the Ro/e CSV ∩ actually present ∩ FIG2_CELLTYPE_ORDER
roe_file <- file.path(OUT_DIR, "path1_roe_celltype_by_group.csv")
stopifnot(file.exists(roe_file))
roe_types <- rownames(read.csv(roe_file, row.names = 1, check.names = FALSE))
present   <- sort(unique(as.character(seu$cell_type)))
ct_order  <- FIG2_CELLTYPE_ORDER[
  FIG2_CELLTYPE_ORDER %in% intersect(roe_types, present)]
cat("cell types (", length(ct_order), "): ", paste(ct_order, collapse = ", "), "\n", sep = "")

## Donor field detection (same as script 38)
DONOR_CAND <- c("donor_id", "donor", "subject", "subject_id",
                "Subject_Identity", "sample")
donor_field <- NULL
for (f in DONOR_CAND) {
  if (f %in% colnames(seu@meta.data)) { donor_field <- f; break }
}
if (is.null(donor_field)) {
  cat("!! no donor field found; C2 degrades to a cell-level boxplot (no donor-level test)\n")
} else {
  cat("donor field:", donor_field, "\n")
}

## =========================================================================
## Step 2: extract expression (Seurat v5 pitfall: row slicing loses names; cell
##         names must be taken from colnames)
## =========================================================================
get_layer <- function(obj, layer) {
  tryCatch(GetAssayData(obj, assay = "RNA", layer = layer),
           error = function(e) GetAssayData(obj, assay = "RNA", slot = layer))
}
data_l <- get_layer(seu, "data")

genes_use <- intersect(GENES, rownames(seu))
if (!length(genes_use)) stop("none of these genes are in the object: ", paste(GENES, collapse = ", "))
cat("genes:", paste(genes_use, collapse = ", "), "\n")

cells <- colnames(seu)
meta  <- seu@meta.data[cells, , drop = FALSE]

## Long table: cells x genes
expr_list <- lapply(genes_use, function(g) {
  v <- as.numeric(data_l[g, , drop = TRUE])   # order matches colnames(seu)
  data.frame(cell = cells, gene = g, expr = v, stringsAsFactors = FALSE)
})
df_cell <- bind_rows(expr_list)

df_cell$cell_type <- as.character(meta[df_cell$cell, "cell_type"])
df_cell$group     <- as.character(meta[df_cell$cell, "group"])
df_cell$donor     <- if (is.null(donor_field)) NA_character_ else
                       as.character(meta[df_cell$cell, donor_field])

df_cell <- df_cell %>%
  filter(cell_type %in% ct_order, group %in% GROUP_LEVELS) %>%
  mutate(cell_type = factor(cell_type, levels = ct_order),
         group     = factor(group,     levels = GROUP_LEVELS))
cat("cell-gene records used for plotting:", nrow(df_cell), "\n")

## =========================================================================
## Step 3: aggregate — cell type x group (C1) and donor level (C2)
## =========================================================================
## C1: mean expression + expression percentage (log-normalized layer: expr>0 means counts>0)
df_ct <- df_cell %>%
  group_by(gene, cell_type, group) %>%
  summarise(avg = mean(expr),
            pct = mean(expr > 0) * 100,
            n_cells  = n(),
            n_donors = n_distinct(donor),
            .groups = "drop")

## C2: donor level (average within donor first, then summarize by group) — avoids pseudo-replication
if (!is.null(donor_field)) {
  df_donor <- df_cell %>%
    filter(!is.na(donor)) %>%
    group_by(gene, cell_type, group, donor) %>%
    summarise(donor_mean = mean(expr), .groups = "drop")
} else {
  df_donor <- df_cell %>%
    mutate(donor = cell) %>%
    group_by(gene, cell_type, group, donor) %>%
    summarise(donor_mean = expr, .groups = "drop")
}

## Pick the "highest-expressing cell type" as the C2 protagonist (by mean of the
## per-group average expression)
top_ct <- df_ct %>%
  filter(gene == genes_use[1]) %>%
  group_by(cell_type) %>%
  summarise(m = mean(avg), .groups = "drop") %>%
  arrange(desc(m)) %>%
  slice(1) %>% pull(cell_type) %>% as.character()
cat("highest-expressing cell type:", top_ct, " (C2 protagonist)\n")

## C2 statistics: Kruskal overall + each group vs Healthy (Wilcoxon)
stat_lines <- c(sprintf("=== Fig4C: %s donor-level expression statistics (cell type: %s) ===",
                        paste(genes_use, collapse = "/"), top_ct), "")
df_top <- df_donor %>% filter(cell_type == top_ct, gene == genes_use[1])
df_top$group <- factor(as.character(df_top$group), levels = GROUP_LEVELS)

grp_vals <- split(df_top$donor_mean, df_top$group)
## Note: do not use grp_vals[GROUP_LEVELS[len>0]] — the two vectors differ in
## length and would index misaligned. Correct approach: drop empty groups first,
## then order by GROUP_LEVELS.
grp_vals <- grp_vals[vapply(grp_vals, length, 1L) > 0]
grp_vals <- grp_vals[intersect(GROUP_LEVELS, names(grp_vals))]
if (length(grp_vals) >= 2) {
  kt <- kruskal.test(grp_vals)
  stat_lines <- c(stat_lines,
                  sprintf("Kruskal-Wallis (all five groups): chi2=%.3g, df=%d, p=%.3g",
                          kt$statistic, kt$parameter, kt$p.value), "")
}
if ("Healthy" %in% names(grp_vals)) {
  stat_lines <- c(stat_lines, "each group vs Healthy (Wilcoxon rank sum):")
  for (g in setdiff(names(grp_vals), "Healthy")) {
    wt <- suppressWarnings(wilcox.test(grp_vals[[g]], grp_vals[["Healthy"]]))
    stat_lines <- c(stat_lines,
      sprintf("  %-18s median=%.4f (n=%d) vs Healthy median=%.4f (n=%d), p=%.3g",
              g, median(grp_vals[[g]]), length(grp_vals[[g]]),
              median(grp_vals[["Healthy"]]), length(grp_vals[["Healthy"]]),
              wt$p.value))
  }
}
writeLines(stat_lines, file.path(OUT_DIR, "Fig4C_SOCS3_donor_stats.txt"))
cat("\n", paste(stat_lines, collapse = "\n"), "\n", sep = "")

## =========================================================================
## Step 4: save pre-plot data
## =========================================================================
plot_data <- list(
  df_cell   = df_cell,
  df_ct     = df_ct,
  df_donor  = df_donor,
  top_ct    = top_ct,
  ct_order  = ct_order,
  genes_use = genes_use,
  GROUP_LEVELS = GROUP_LEVELS, GROUP_LABEL = GROUP_LABEL, GROUP_COLS = GROUP_COLS,
  donor_field  = donor_field
)
saveRDS(plot_data, file.path(OUT_DIR, "path46_fig4c_plotdata.rds"))
write.csv(df_ct,    file.path(OUT_DIR, "Fig4C_SOCS3_celltype.csv"), row.names = FALSE)
write.csv(df_donor, file.path(OUT_DIR, "Fig4C_SOCS3_donor.csv"),    row.names = FALSE)
cat("\nsaved: Fig4C_SOCS3_celltype.csv / Fig4C_SOCS3_donor.csv / path46_fig4c_plotdata.rds\n")

## =========================================================================
## Step 5: C1 — cell type x group dot plot
## =========================================================================
df_c1 <- df_ct %>% filter(gene == genes_use[1])
df_c1$grp_lab <- factor(GROUP_LABEL[as.character(df_c1$group)],
                        levels = GROUP_LABEL[GROUP_LEVELS])

p_c1 <- ggplot(df_c1, aes(x = cell_type, y = grp_lab, size = pct, color = avg)) +
  geom_point(shape = 16, stroke = 0.2) +
  scale_color_gradientn(
    colors = c("#F7FBFF", "#DEEBF7", "#9ECAE1", "#3182BD", "#08519C"),
    name = "Mean expression\n(log-normalized)",
    limits = c(0, max(df_c1$avg) * 1.02)
  ) +
  scale_size_area(breaks = c(0, 25, 50, 75),
                  labels = c("0%", "25%", "50%", "75%"),
                  name = "% cells", max_size = 6.5, limits = c(0, 100)) +
  labs(x = NULL, y = NULL,
       title = "SOCS3 expression across major cell types and groups",
       subtitle = "dot size = % expressing cells; color = mean expression") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10, colour = "grey20"),
        axis.text.y = element_text(size = 10, colour = "grey20"),
        panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.position = "right", legend.box = "vertical",
        legend.title = element_text(size = 9), legend.text = element_text(size = 8),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 8.5, colour = "grey40"))

## =========================================================================
## Step 6: C2 — donor-level five-group distribution (boxplot + per-donor points)
## =========================================================================
df_c2 <- df_top %>%
  mutate(grp_lab = factor(GROUP_LABEL[as.character(group)],
                          levels = GROUP_LABEL[GROUP_LEVELS]))

p_c2 <- ggplot(df_c2, aes(x = grp_lab, y = donor_mean, colour = grp_lab)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, fill = NA,
               linewidth = 0.45, alpha = 0.85) +
  geom_jitter(position = position_jitter(width = 0.12, seed = 123),
              size = 1.6, alpha = 0.75, stroke = 0) +
  scale_colour_manual(values = GROUP_COLS, guide = "none") +
  labs(x = NULL,
       y = sprintf("%s mean expression%s", genes_use[1],
                   if (is.null(donor_field)) " (cell level)" else " (donor level)"),
       title = sprintf("%s-level %s in %s",
                       if (is.null(donor_field)) "Cell" else "Donor",
                       genes_use[1], top_ct),
       subtitle = sprintf("n = %d %s; Wilcoxon vs Healthy",
                          n_distinct(df_c2$donor),
                          if (is.null(donor_field)) "cells" else "donors")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10, colour = "grey20"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8.5, colour = "grey35"))

## =========================================================================
## Step 7: output — combined figure (when has_pw) + C1/C2 single panels
##         (always, for per-panel assembly of Fig4)
## =========================================================================
if (has_pw) {
  p_comb <- p_c1 + p_c2 + plot_layout(widths = c(1.05, 1), guides = "collect")
  ggsave(file.path(OUT_DIR, "Fig4C_SOCS3_celltype_specificity.pdf"), p_comb,
         width = 12.5, height = 5.4)
  cat("saved: Fig4C_SOCS3_celltype_specificity.pdf\n")
  tryCatch({
    grDevices::png(file.path(OUT_DIR, "Fig4C_SOCS3_celltype_specificity.png"),
                   width = 12.5, height = 5.4, units = "in", res = 300)
    print(p_comb); grDevices::dev.off()
    cat("saved: Fig4C_SOCS3_celltype_specificity.png\n")
  }, error = function(e) cat("PNG write skipped (can be recovered from PDF with PyMuPDF):",
                             conditionMessage(e), "\n"))
}

## ---- C1 / C2 single panels (for independent assembly of Fig4 panel c / d) ----
ggsave(file.path(OUT_DIR, "Fig4C1_SOCS3_celltype_dotplot.pdf"), p_c1,
       width = 7, height = 5.4)
ggsave(file.path(OUT_DIR, "Fig4C2_SOCS3_donor_boxplot.pdf"), p_c2,
       width = 6, height = 5.4)
tryCatch({
  grDevices::png(file.path(OUT_DIR, "Fig4C1_SOCS3_celltype_dotplot.png"),
                 width = 7, height = 5.4, units = "in", res = 300)
  print(p_c1); grDevices::dev.off()
  grDevices::png(file.path(OUT_DIR, "Fig4C2_SOCS3_donor_boxplot.png"),
                 width = 6, height = 5.4, units = "in", res = 300)
  print(p_c2); grDevices::dev.off()
  cat("saved: Fig4C1 / Fig4C2 single panels {pdf,png}\n")
}, error = function(e) cat("C1/C2 PNG write skipped (can be recovered from PDF with PyMuPDF):",
                           conditionMessage(e), "\n"))

cat("\n===== Script 46 done (Fig4C SOCS3 cell-type specificity) =====\n")
