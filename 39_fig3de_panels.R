###############################################################################
# Sepsis + COPD comorbidity analysis — Script 39: Figure 3 panels d/e (two new panels)
# File (current name): 39_fig3de_panels.R (originally named 39_Fig3d_e dual-tissue and leave-one-out panels.R)
#
# Background: Figure 3 uses the "6-panel full version" (confirmed by user):
#   a four-way Venn + hypergeometric / b 149 pathway enrichment / c 5-group heatmap (server) /
#   d blood-lung dual-tissue donor-level logFC scatter / e leave-one-dataset-out robustness / f GSE57148
# This script only outputs d and e (data are ready-made CSVs, can run locally or on server, no Seurat dependency).
#
#   d: x = blood monocyte donor-level logFC (Sepsis vs Healthy, edgeR)
#      y = lung myeloid donor-level logFC (COPD vs Control, edgeR)
#      all genes as grey base points, the original 149 shared genes highlighted red, key genes labeled;
#      extras: whole-genome Pearson/Spearman r, number of 149 genes in the upper-right quadrant (≈146/149).
#   e: number of 149 genes that remain directionally consistent (up-regulated on blood side) when leaving one dataset out, n=6 datasets;
#      dashed line 149 = removing no dataset; GSE279452 (largest sepsis cohort) marked red separately.
#
# Usage:
#   Server: Rscript 39_fig3de_panels.R
#           (auto-detect this directory as both data and output directory)
#   Windows local (non-UTF-8 locale): pass two directories via an ASCII junction:
#     Rscript 39_fig3de_panels.R <data_dir> <output_dir>
#
# Output (written to output directory):
#   Fig3d_cross_tissue_logFC_scatter.pdf / .png   (300 dpi)
#   Fig3e_leave_one_out_149.pdf / .png            (300 dpi)
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
mark_utf8 <- function(x) { if (!is.null(x)) Encoding(x) <- "UTF-8"; x }

## ---- Auto-locate script directory ----
.this_file <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
if (.this_file == "" || !file.exists(.this_file)) {
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1])
}
SCRIPT_DIR <- dirname(.this_file)

if (length(args) >= 1) DATA_DIR <- mark_utf8(args[1]) else DATA_DIR <- mark_utf8(SCRIPT_DIR)
if (length(args) >= 2) OUT_DIR  <- mark_utf8(args[2]) else OUT_DIR  <- mark_utf8(SCRIPT_DIR)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Data directory:", DATA_DIR, "\n")
cat("Output directory:", OUT_DIR, "\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

## =========================================================================
## d: blood-lung dual-tissue donor-level logFC scatter
## =========================================================================
cat("\n===== Panel d: blood-lung dual-tissue logFC scatter =====\n")

blood <- read.csv(file.path(DATA_DIR, "path28b_blood_mono_pseudobulk_deg.csv"),
                  stringsAsFactors = FALSE)
lung  <- read.csv(file.path(DATA_DIR, "path28b_lung_myeloid_pseudobulk_deg.csv"),
                  stringsAsFactors = FALSE)
stopifnot(all(c("gene", "logFC", "FDR") %in% colnames(blood)),
          all(c("gene", "logFC", "FDR") %in% colnames(lung)))

m <- merge(blood[, c("gene", "logFC", "FDR")],
           lung [, c("gene", "logFC", "FDR")],
           by = "gene", suffixes = c("_blood", "_lung"))
cat("Genes measured in both tissues:", nrow(m), "\n")

# original 149 shared genes (subject of the directional survival analysis)
f149 <- file.path(DATA_DIR, "path2_shared_myeloid_genes.txt")
if (!file.exists(f149)) f149 <- file.path(DATA_DIR, "path21_shared_up_genes.txt")
genes149 <- trimws(readLines(f149)); genes149 <- genes149[nzchar(genes149)]
cat("Original 149 genes:", length(genes149), "\n")

m$is149 <- m$gene %in% genes149
m$both_sig <- m$FDR_blood < 0.05 & m$FDR_lung < 0.05
r_pear  <- cor(m$logFC_blood, m$logFC_lung, method = "pearson")
r_spear <- cor(m$logFC_blood, m$logFC_lung, method = "spearman")
n149_q1 <- sum(m$is149 & m$logFC_blood > 0 & m$logFC_lung > 0)
r149    <- cor(m$logFC_blood[m$is149], m$logFC_lung[m$is149], method = "pearson")
cat(sprintf("Whole-genome Pearson r = %.3f | Spearman r = %.3f\n", r_pear, r_spear))
cat(sprintf("Within-149 Pearson r = %.3f | upper-right quadrant: %d / %d\n", r149, n149_q1, sum(m$is149)))

# Label genes: those named in the manuscript + the largest logFC product among both-side FDR<0.05 in the upper-right quadrant, capped at ~10
key_genes <- c("SERPINA1", "STAT3", "TREM1", "FCGR1A", "C1QA", "PDK4",
               "SOCS3", "IFIT3", "GBP1", "IFI6", "MT2A", "HP")
key_genes <- intersect(key_genes, m$gene)
cand <- m$is149 & m$both_sig & m$logFC_blood > 0 & m$logFC_lung > 0
extra <- m$gene[cand][order(m$logFC_blood[cand] * m$logFC_lung[cand],
                            decreasing = TRUE)]
extra <- setdiff(extra, key_genes)
labels <- c(key_genes, head(extra, max(0, 10 - length(key_genes))))
m$lab <- ifelse(m$gene %in% labels & m$is149, m$gene, "")

p_d <- ggplot(m, aes(logFC_blood, logFC_lung)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(data = m[!m$is149, ], colour = "grey75", size = 0.7, alpha = 0.5) +
  geom_point(data = m[m$is149, ],  colour = "#E64B35", size = 1.6, alpha = 0.9) +
  geom_text_repel(aes(label = lab), colour = "#791F1F", size = 2.6,
                  max.overlaps = 40, seed = 123, show.legend = FALSE,
                  min.segment.length = 0.2) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4, size = 3.3,
           colour = "grey20",
           label = sprintf("r over all genes = %.3f\nwithin 149 shared: r = %.3f, %d / %d co-upregulated",
                           r_pear, r149, n149_q1, sum(m$is149))) +
  scale_x_continuous(expand = expansion(mult = 0.04)) +
  scale_y_continuous(expand = expansion(mult = 0.04)) +
  labs(x = "Blood monocyte log2FC (sepsis vs healthy, donor-level edgeR)",
       y = "Lung myeloid log2FC (COPD vs control, donor-level edgeR)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "Fig3d_cross_tissue_logFC_scatter.pdf"), p_d,
       width = 6.6, height = 5.6)
ggsave(file.path(OUT_DIR, "Fig3d_cross_tissue_logFC_scatter.png"), p_d,
       width = 6.6, height = 5.6, dpi = 300, device = grDevices::png)
cat("Saved: Fig3d_cross_tissue_logFC_scatter.{pdf,png}\n")

## =========================================================================
## e: leave-one-dataset-out robustness bar chart
## =========================================================================
cat("\n===== Panel e: leave-one-dataset-out =====\n")

lo <- read.csv(file.path(DATA_DIR, "path31_lo_149_consistency.csv"),
               stringsAsFactors = FALSE)
stopifnot(all(c("dataset_removed", "n_149_ok") %in% colnames(lo)))
lo$pct <- round(100 * lo$n_149_ok / 149, 1)
lo$hl  <- lo$n_149_ok == min(lo$n_149_ok)      # the only notable outlier: GSE279452
lo$note <- ifelse(lo$hl, "largest sepsis cohort (93 donors)", "")
lo <- lo[order(lo$n_149_ok, decreasing = FALSE), ]
lo$dataset_removed <- factor(lo$dataset_removed, levels = lo$dataset_removed)

# Horizontal lollipop: one row per dataset; the value is "after removing that dataset, re-run donor-level edgeR,
# count of 149 genes still with log2FC>0"; red = the only outlier, with a reason note
p_e <- ggplot(lo, aes(n_149_ok, dataset_removed)) +
  geom_vline(xintercept = 147, colour = "grey55", linetype = "dashed",
             linewidth = 0.4) +
  geom_segment(aes(xend = 0, yend = dataset_removed, colour = hl),
               linewidth = 1.1) +
  geom_point(aes(colour = hl), size = 3.4) +
  scale_colour_manual(values = c("FALSE" = "#185FA5", "TRUE" = "#E64B35"),
                      guide = "none") +
  geom_text(aes(label = sprintf("%d/149  (%.1f%%)", n_149_ok, pct),
                x = n_149_ok + 6, colour = hl), hjust = 0, size = 3.1) +
  scale_x_continuous(limits = c(0, 196), breaks = c(0, 50, 100, 147, 150)) +
  labs(
    title = "Leave-one-dataset-out robustness",
    subtitle = "drop one blood dataset, re-run donor-level pseudobulk edgeR\nin blood monocytes (sepsis vs healthy); dashed = no removal (147/149)",
    x = "149-gene shared program retaining up-regulation (log2FC > 0)",
    y = "Dataset removed") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())

ggsave(file.path(OUT_DIR, "Fig3e_leave_one_out_149.pdf"), p_e,
       width = 7.2, height = 3.6)
ggsave(file.path(OUT_DIR, "Fig3e_leave_one_out_149.png"), p_e,
       width = 7.2, height = 3.6, dpi = 300, device = grDevices::png)
cat("Saved: Fig3e_leave_one_out_149.{pdf,png}\n")

cat("\n===== Script 39 finished =====\n")
