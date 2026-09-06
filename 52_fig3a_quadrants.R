# =========================================================================
# 52_Fig3A_four-quadrant overview.R
# Fig 3 panel a (v2): cross-tissue significance + directionality four-quadrant overview
#   Direction-agnostic "bilaterally significant" genes (lung myeloid COPD vs Ctrl & blood monocyte Sepsis vs Healthy
#   both FDR<0.05) decomposed into four quadrants: co-activated / co-inhibited / two reversal types.
#   Co-activated arm (137 double-significant) + up-regulated arm criterion (avg_log2FC>0) = 149-gene main program;
#   SOCS3 falls in the "lung-down / blood-up" reversal quadrant (highlighted, foreshadowing Fig4).
#
# Input: path2_lung_myeloid_DEG_COPD_vs_Control.csv
#       path2_blood_monocyte_DEG_Sepsis_vs_Healthy.csv
#       path21_shared_up_genes.txt (149, used for figure-legend corner note)
# Output: Fig3A_cross_tissue_quadrants.{pdf,png}
#       path52_quadrants.csv         (all genes x lung/blood logFC·FDR + quadrant label)
#       path52_fig3a_plotdata.rds    (pre-plot data: long table / stats / colors)
# Run: Rscript 52_Fig3A_four_quadrant_overview.R [OUT_DIR]
# =========================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})
has_repel <- requireNamespace("ggrepel", quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
if (length(args) >= 1 && nzchar(args[1])) script_dir <- args[1]
OUT_DIR <- if (length(args) >= 2 && nzchar(args[2])) args[2] else script_dir
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

f_lung  <- file.path(script_dir, "path2_lung_myeloid_DEG_COPD_vs_Control.csv")
f_blood <- file.path(script_dir, "path2_blood_monocyte_DEG_Sepsis_vs_Healthy.csv")
f_149   <- file.path(script_dir, "path21_shared_up_genes.txt")
stopifnot(file.exists(f_lung), file.exists(f_blood))

load_deg <- function(f) {
  d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  fc <- intersect(c("avg_log2FC", "log2FoldChange", "logFC"), colnames(d))[1]
  pa <- intersect(c("p_val_adj", "padj", "FDR"), colnames(d))[1]
  if (is.na(fc) || is.na(pa)) stop("Column names do not match: ", f)
  ## Gene-name column is empty string ("") or "gene"; always take the first column (position index)
  out <- data.frame(gene = as.character(d[[1]]),
                    logFC = as.numeric(d[[fc]]),
                    FDR   = as.numeric(d[[pa]]),
                    stringsAsFactors = FALSE)
  out[is.finite(out$logFC) & is.finite(out$FDR), ]
}
lung  <- load_deg(f_lung)
blood <- load_deg(f_blood)

## Union long table
df <- full_join(lung, blood, by = "gene", suffix = c("_lung", "_blood")) %>%
  filter(!is.na(logFC_lung) & !is.na(logFC_blood) &
         !is.na(FDR_lung)   & !is.na(FDR_blood))

## Classification
df <- df %>% mutate(
  sig_both = FDR_lung < 0.05 & FDR_blood < 0.05,
  quad = case_when(
    !sig_both                 ~ "ns",
    logFC_lung > 0 & logFC_blood > 0 ~ "co_up",
    logFC_lung < 0 & logFC_blood < 0 ~ "co_down",
    logFC_lung > 0 & logFC_blood < 0 ~ "lung_up_blood_down",
    TRUE                        ~ "lung_down_blood_up"
  ))

cnt <- df %>% count(quad)
print(cnt)
quad_cols <- c(
  co_up                 = "#A32D2D",
  co_down               = "#185FA5",
  lung_up_blood_down    = "#B5B0A6",
  lung_down_blood_up    = "#888780",
  ns                    = "#DCD8CF"
)

## 149 reference (used for corner note / annotation)
if (file.exists(f_149)) {
  g149 <- readLines(f_149, warn = FALSE)
  g149 <- trimws(g149[nzchar(g149)])
  n149 <- length(g149)
  in_co_up <- sum(g149 %in% df$gene[df$quad == "co_up"])
  cat(sprintf("Of the 149 genes, those in co_up (double-significant co-direction) quadrant: %d / %d\n", in_co_up, n149))
} else {
  n149 <- NA
}

## ===== Plot: cross-tissue four-quadrant scatter =====
df$quad <- factor(df$quad, levels = c("co_up", "co_down",
                                      "lung_up_blood_down", "lung_down_blood_up", "ns"))
bg <- df[df$quad == "ns", ]
fg <- df[df$quad != "ns", ]
hub <- df[df$gene == "SOCS3", ]

p <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_point(data = bg, aes(logFC_lung, logFC_blood),
             colour = quad_cols["ns"], size = 0.85, alpha = 0.32, stroke = 0) +
  geom_point(data = fg, aes(logFC_lung, logFC_blood, colour = quad),
             size = 1.8, alpha = 0.70, stroke = 0) +
  scale_colour_manual(
    values = quad_cols,
    breaks = c("co_up", "co_down", "lung_up_blood_down", "lung_down_blood_up"),
    labels = c(
      sprintf("Co-upregulated  (n = %d)", sum(df$quad == "co_up")),
      sprintf("Co-downregulated (n = %d)", sum(df$quad == "co_down")),
      sprintf("Lung-up / blood-down (n = %d)", sum(df$quad == "lung_up_blood_down")),
      sprintf("Lung-down / blood-up (n = %d)", sum(df$quad == "lung_down_blood_up"))),
    name = NULL)

## SOCS3 highlight
if (nrow(hub)) {
  p <- p +   geom_point(data = hub, aes(logFC_lung, logFC_blood),
                      shape = 21, size = 4.2, fill = "#7B2E00",
                      colour = "white", stroke = 0.7)
  if (has_repel) {
    p <- p + ggrepel::geom_text_repel(
      data = hub, aes(logFC_lung, logFC_blood, label = gene),
      colour = "#7B2E00", size = 4.0, fontface = "bold",
      box.padding = 0.5, min.segment.length = 0, nudge_y = -0.32,
      segment.colour = "#7B2E00", segment.size = 0.35)
  } else {
    p <- p + geom_text(data = hub, aes(logFC_lung, logFC_blood, label = gene),
                       colour = "#7B2E00", size = 4.0, fontface = "bold",
                       hjust = -0.2, vjust = -1.2)
  }
}

n_both <- sum(df$quad != "ns")
n_co_up <- sum(df$quad == "co_up")

## Diagnostics: points outside the ±2 display window (display clipped only, data not deleted)
lim <- 2
in_win <- with(df, abs(logFC_lung) <= lim & abs(logFC_blood) <= lim)
cat(sprintf("Display window ±%g: %d / %d points inside window (%.1f%%); %d points outside window\n",
            lim, sum(in_win), nrow(df), 100 * mean(in_win), sum(!in_win)))
if (any(!in_win)) {
  cat("  Genes outside window (significant genes): ",
      paste(head(df$gene[!in_win & df$quad != "ns"], 20), collapse = ", "), "\n", sep = "")
  cat("  Count outside window (significant genes): ", sum(!in_win & df$quad != "ns"), " / ",
      sum(!in_win), "\n", sep = "")
}
if (nrow(hub)) {
  cat(sprintf("  SOCS3 coordinates: lung = %.3f, blood = %.3f (inside window: %s)\n",
              hub$logFC_lung, hub$logFC_blood,
              ifelse(abs(hub$logFC_lung) <= lim & abs(hub$logFC_blood) <= lim, "yes", "no")))
}
sub_txt <- sprintf(
  "n = %d genes significant in both tissues; %d co-upregulated (the double-significant core of the 149-gene program); SOCS3 falls in the lung-down / blood-up quadrant",
  n_both, n_co_up)

p <- p +
  scale_x_continuous(expand = expansion(mult = 0.03)) +
  scale_y_continuous(expand = expansion(mult = 0.03)) +
  labs(title = "Cross-tissue significance and directionality",
       subtitle = sub_txt,
       x = "Lung myeloid log2FC (COPD vs control)",
       y = "Blood monocyte log2FC (sepsis vs healthy)") +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey93", linewidth = 0.25),
        legend.position = "right",
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"))

ggsave(file.path(OUT_DIR, "Fig3A_cross_tissue_quadrants.pdf"), p,
       width = 10, height = 6.2)
cat("Saved: Fig3A_cross_tissue_quadrants.pdf\n")
tryCatch({
  grDevices::png(file.path(OUT_DIR, "Fig3A_cross_tissue_quadrants.png"),
                 width = 10, height = 6.2, units = "in", res = 300)
  print(p); grDevices::dev.off()
  cat("Saved: Fig3A_cross_tissue_quadrants.png\n")
}, error = function(e) cat("PNG write skipped:", conditionMessage(e), "\n"))

## Save pre-plot data
plot_data <- list(
  df      = df,
  cnt     = cnt,
  quad_cols = quad_cols,
  hub_gene = if (nrow(hub)) hub$gene else NULL,
  n149    = n149,
  subtitle = sub_txt
)
saveRDS(plot_data, file.path(OUT_DIR, "path52_fig3a_plotdata.rds"))
write.csv(df, file.path(OUT_DIR, "path52_quadrants.csv"), row.names = FALSE)
cat("Saved: path52_quadrants.csv + path52_fig3a_plotdata.rds\n")

cat("\n===== Script 52 complete (Fig3A four-quadrant overview) =====\n")
