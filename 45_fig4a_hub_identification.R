#!/usr/bin/env Rscript
# =============================================================================
# 45_fig4a_hub_identification.R
# Fig4 Panel A: identify the hub gene SOCS3 from donor-level differential expression
# =============================================================================
# Background (important, relation to Fig3):
#   SOCS3 is **not** in the 149 shared myeloid programs (verified: neither
#   path2_shared_myeloid_genes.txt / path21_shared_up_genes.txt /
#   path28b_shared_donor_genes.txt contains SOCS3).
#   So this panel cannot use the narrative "SOCS3 locked from the 149" — that
#   would be wrong.
#   Correct narrative: SOCS3 emerges from the **whole-transcriptome donor-level
#   DEGs** as the single strong-signal gene passing the filter
#     "consistently significantly upregulated in all four groups (IC/COPD/Sepsis/SP)
#      + monotonically increasing with severity"
#   — an axis **independent of the 149 shared programs**.
#
# Verified key numbers (whole transcriptome, 15259 genes):
#   FDR<0.05 and logFC>0 in all four groups, with full monotonicity IC<COPD<Sepsis<SP
#   -> only 3 genes:
#     SOCS3   mean logFC 1.656  (0.911 / 1.440 / 1.587 / 2.686)   <- clear first by a wide margin
#     CLIC1   mean logFC 0.458
#     MFSD10  mean logFC 0.411
#
# Figures (two panels, patchwork):
#   A1 Global scatter: x = mean logFC across the four groups, y = -log10(strongest FDR),
#                grey points = whole transcriptome, hit candidates highlighted, SOCS3 red
#                with label; funnel counts of each filter step annotated in the corner.
#   A2 Trajectory lines: candidate genes x four groups (IC->COPD->Sepsis->SP) logFC lines,
#                SOCS3 in bold red, the rest as thin grey lines.
#
# Input (under the server OUT_DIR, produced by script 28):
#   path28_blood_mono_pseudobulk_{Infection_Control,COPD,Sepsis,Sepsis_Pneumonia}_vs_Healthy.csv
#   Columns: logFC, logCPM, F, PValue, FDR, gene
#
# Output:
#   Fig4A_hub_identification.pdf / .png
#   Fig4A_hub_candidates.csv          candidate gene details (for checking / writing the text)
#   path45_fig4a_plotdata.rds         pre-plot data (re-plot locally)
#
# Run: Rscript 45_fig4a_hub_identification.R [OUT_DIR]
# Depends on: ggplot2, dplyr, tidyr, patchwork (no Seurat dependency; runs in seconds)
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
script_dir <- tryCatch(dirname(normalizePath(sys.frames()[[1]]$ofile)),
                       error = function(e) getwd())
OUT_DIR <- ifelse(length(args) >= 1 && nzchar(args[1]), args[1], script_dir)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})
has_pw <- requireNamespace("patchwork", quietly = TRUE)
if (!has_pw) {
  cat("!! patchwork not installed; only A1/A2 single panels will be output (no combined figure)\n")
} else {
  suppressPackageStartupMessages(library(patchwork))
}
## Note: R's pkg::fn() errors out directly when the package is not installed
##       (it does not return NULL), so the check must happen before building the
##       ggplot; you cannot use the `+ { if (...) pkg::fn() }` pattern.
has_repel <- requireNamespace("ggrepel", quietly = TRUE)
if (!has_repel) cat("note: ggrepel not installed; gene labels use geom_text (may overlap slightly)\n")

## =========================================================================
## Tunable parameters
## =========================================================================
FDR_CUT   <- 0.05      # significance threshold
MONOTONE  <- "strict"  # "strict" = full monotonicity IC<COPD<Sepsis<SP
                       # "loose"  = Sepsis>IC and SP>Sepsis (COPD not forced in the middle)
TOPN_LINE <- 8         # A2 draws at most this many candidates (top N by mean logFC)
HUB_GENE  <- "SOCS3"   # highlighted gene

## =========================================================================
## Step 1: read the four donor-level DEG tables and merge
## =========================================================================
GROUPS <- c(Infection_Control = "Infection_Control_vs_Healthy",
            COPD              = "COPD_vs_Healthy",
            Sepsis            = "Sepsis_vs_Healthy",
            Sepsis_Pneumonia  = "Sepsis_Pneumonia_vs_Healthy")
GROUP_ORDER <- c("Infection_Control", "COPD", "Sepsis", "Sepsis_Pneumonia")
GROUP_LABEL <- c(Infection_Control = "IC", COPD = "COPD",
                 Sepsis = "Sepsis", Sepsis_Pneumonia = "Sepsis+PNA")

deg_list <- list()
for (g in names(GROUPS)) {
  f <- file.path(OUT_DIR, paste0("path28_blood_mono_pseudobulk_", GROUPS[g], ".csv"))
  if (!file.exists(f)) stop("missing input file: ", basename(f))
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(all(c("gene", "logFC", "FDR") %in% colnames(d)))
  deg_list[[g]] <- d[, c("gene", "logFC", "FDR")]
  cat(sprintf("read %-18s %d genes\n", g, nrow(d)))
}

## Wide merge: gene | logFC_xxx | FDR_xxx
wide <- deg_list[[names(GROUPS)[1]]]
names(wide) <- c("gene", paste0("logFC_", names(GROUPS)[1]), paste0("FDR_", names(GROUPS)[1]))
for (g in names(GROUPS)[-1]) {
  d <- deg_list[[g]]
  names(d) <- c("gene", paste0("logFC_", g), paste0("FDR_", g))
  wide <- inner_join(wide, d, by = "gene")
}
cat("genes after merge:", nrow(wide), "\n")

lf_cols <- paste0("logFC_", GROUP_ORDER)
fd_cols <- paste0("FDR_",   GROUP_ORDER)

wide$mean_logFC <- rowMeans(wide[, lf_cols], na.rm = TRUE)
wide$max_FDR    <- apply(wide[, fd_cols], 1, max, na.rm = TRUE)   # weakest group
wide$min_FDR    <- apply(wide[, fd_cols], 1, min, na.rm = TRUE)   # strongest evidence
wide$neglog10_FDR <- -log10(pmax(wide$min_FDR, 1e-300))

## =========================================================================
## Step 2: filter — concordantly significantly upregulated in all four groups + monotonic increase
## =========================================================================
lf <- as.matrix(wide[, lf_cols])
fd <- as.matrix(wide[, fd_cols])

sig_up <- rowSums(fd < FDR_CUT) == 4 & rowSums(lf > 0) == 4
cat("significantly upregulated in all four groups (FDR<", FDR_CUT, " and logFC>0): ", sum(sig_up), "\n", sep = "")

if (MONOTONE == "strict") {
  mono <- lf[, "logFC_Infection_Control"] < lf[, "logFC_COPD"] &
          lf[, "logFC_COPD"]              < lf[, "logFC_Sepsis"] &
          lf[, "logFC_Sepsis"]            < lf[, "logFC_Sepsis_Pneumonia"]
} else {
  mono <- lf[, "logFC_Infection_Control"] < lf[, "logFC_Sepsis"] &
          lf[, "logFC_Sepsis"]            < lf[, "logFC_Sepsis_Pneumonia"]
}
hit <- sig_up & mono
cat("of which monotonically increasing with severity (", MONOTONE, "): ", sum(hit), "\n", sep = "")

wide$hit <- hit
cand <- wide[hit, ] %>% arrange(desc(mean_logFC))
if (nrow(cand) == 0) stop("no hit candidate genes; relax MONOTONE or FDR_CUT")

cat("\ncandidate genes (descending by mean logFC):\n")
for (i in seq_len(min(nrow(cand), 15))) {
  cat(sprintf("  %2d. %-10s mean=%.3f  IC=%.3f COPD=%.3f Sepsis=%.3f SP=%.3f\n",
              i, cand$gene[i], cand$mean_logFC[i],
              cand$logFC_Infection_Control[i], cand$logFC_COPD[i],
              cand$logFC_Sepsis[i], cand$logFC_Sepsis_Pneumonia[i]))
}
hub_rank <- which(cand$gene == HUB_GENE)
if (length(hub_rank)) cat("\n>>> ", HUB_GENE, " rank: #", hub_rank[1], " / ", nrow(cand), "\n", sep = "")

## =========================================================================
## Step 3: save pre-plot data
## =========================================================================
line_genes <- cand$gene[seq_len(min(TOPN_LINE, nrow(cand)))]
cand_long <- cand %>%
  select(gene, all_of(lf_cols)) %>%
  pivot_longer(cols = all_of(lf_cols), names_to = "group", values_to = "logFC") %>%
  mutate(group = sub("^logFC_", "", group),
         group = factor(group, levels = GROUP_ORDER),
         gene  = factor(gene, levels = rev(cand$gene)))

plot_data <- list(
  wide_all  = wide[, c("gene", lf_cols, fd_cols, "mean_logFC", "max_FDR",
                       "min_FDR", "neglog10_FDR", "hit")],
  candidates = cand,
  cand_long  = cand_long,
  line_genes = line_genes,
  params = list(FDR_CUT = FDR_CUT, MONOTONE = MONOTONE, HUB_GENE = HUB_GENE),
  funnel = c(all = nrow(wide), sig_up = sum(sig_up), hit = sum(hit)),
  GROUP_ORDER = GROUP_ORDER, GROUP_LABEL = GROUP_LABEL
)
saveRDS(plot_data, file.path(OUT_DIR, "path45_fig4a_plotdata.rds"))
write.csv(cand, file.path(OUT_DIR, "Fig4A_hub_candidates.csv"), row.names = FALSE)
cat("\nsaved: Fig4A_hub_candidates.csv + path45_fig4a_plotdata.rds\n")

## =========================================================================
## Step 4: A1 — global scatter (whole-transcriptome background + hit highlights)
## =========================================================================
HUB_COL  <- "#B2182B"   # SOCS3 red
HIT_COL  <- "#2166AC"   # blue for other candidates
BG_COL   <- "grey75"

bg  <- wide[!wide$hit, ]
ht  <- wide[wide$hit & wide$gene != HUB_GENE, ]
hub <- wide[wide$gene == HUB_GENE, ]

funnel_lab <- sprintf("Screening: %s genes \U2192 %s concordant up \U2192 %s monotonic",
                      format(nrow(wide), big.mark = ","),
                      sum(sig_up), sum(hit))

p_a1 <- ggplot() +
  geom_point(data = bg, aes(mean_logFC, neglog10_FDR),
             colour = BG_COL, size = 0.45, alpha = 0.45, stroke = 0)

if (nrow(ht)) {
  p_a1 <- p_a1 + geom_point(data = ht, aes(mean_logFC, neglog10_FDR),
                            colour = HIT_COL, size = 2.8, alpha = 0.95, stroke = 0)
  ## Blue-point gene labels (nudged down to avoid SOCS3's top label)
  if (has_repel) {
    p_a1 <- p_a1 + ggrepel::geom_text_repel(
      data = ht, aes(mean_logFC, neglog10_FDR, label = gene),
      colour = HIT_COL, size = 3.6, fontface = "bold",
      box.padding = 0.4, min.segment.length = 0, nudge_y = -2,
      max.overlaps = 10, segment.colour = HIT_COL, segment.size = 0.3)
  } else {
    p_a1 <- p_a1 + geom_text(data = ht, aes(mean_logFC, neglog10_FDR, label = gene),
                             colour = HIT_COL, size = 3.6, fontface = "bold",
                             vjust = 1.4)
  }
}

if (nrow(hub)) {
  p_a1 <- p_a1 + geom_point(data = hub, aes(mean_logFC, neglog10_FDR),
                            colour = HUB_COL, size = 3.4, stroke = 0.6,
                            shape = 21, fill = HUB_COL)
  if (has_repel) {
    p_a1 <- p_a1 + ggrepel::geom_text_repel(
      data = hub, aes(mean_logFC, neglog10_FDR, label = gene),
      colour = HUB_COL, size = 4.0, fontface = "bold",
      box.padding = 0.5, min.segment.length = 0, nudge_y = 2.5,
      segment.colour = HUB_COL, segment.size = 0.4)
  } else {
    p_a1 <- p_a1 + geom_text(data = hub, aes(mean_logFC, neglog10_FDR, label = gene),
                             colour = HUB_COL, size = 4.0, fontface = "bold",
                             vjust = -0.9)
  }
}

p_a1 <- p_a1 +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.04))) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 14)) +
  labs(x = "Mean logFC across four comparisons (donor-level)",
       y = expression(-log[10](FDR) ~ "(best of four)")) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 8.5, colour = "grey35")) +
  labs(title = "Hub gene identification from donor-level DEGs",
       subtitle = funnel_lab) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.02, vjust = 1.6,
           label = sprintf("n = %d monotonic", sum(hit)),
           colour = HUB_COL, size = 3.6, fontface = "bold")

## =========================================================================
## Step 5: A2 — trajectory lines of candidate genes across the four groups
## =========================================================================
cand_long_line <- cand_long[cand_long$gene %in% line_genes, ]
cand_long_line$is_hub <- as.character(cand_long_line$gene) == HUB_GENE

## SOCS3 endpoint label data (pulled separately to avoid complex indexing)
hub_endpt <- cand_long_line[cand_long_line$gene == HUB_GENE &
                            cand_long_line$group == "Sepsis_Pneumonia", ]
hub_endpt <- hub_endpt[!is.na(hub_endpt$logFC), ]

p_a2 <- ggplot(cand_long_line,
               aes(x = group, y = logFC, group = gene,
                   colour = is_hub, size = is_hub, alpha = is_hub)) +
  geom_line(lineend = "round") +
  geom_point(shape = 21, fill = "white", stroke = 0.5) +
  scale_colour_manual(values = c("TRUE" = HUB_COL, "FALSE" = "grey55"), guide = "none") +
  scale_size_manual(values   = c("TRUE" = 1.0,     "FALSE" = 0.45),    guide = "none") +
  scale_alpha_manual(values  = c("TRUE" = 1.0,     "FALSE" = 0.75),    guide = "none") +
  scale_x_discrete(labels = GROUP_LABEL[GROUP_ORDER])

if (nrow(hub_endpt)) {
  hub_endpt$lab <- sprintf("%s  %.2f", hub_endpt$gene, hub_endpt$logFC)
  if (has_repel) {
    p_a2 <- p_a2 + ggrepel::geom_text_repel(
      data = hub_endpt, aes(label = lab),
      colour = HUB_COL, size = 3.8, fontface = "bold",
      nudge_x = 0.35, nudge_y = 0.12, segment.colour = HUB_COL, segment.size = 0.3)
  } else {
    p_a2 <- p_a2 + geom_text(data = hub_endpt, aes(label = lab),
                             colour = HUB_COL, size = 3.8, fontface = "bold",
                             hjust = -0.15, vjust = 1.4)
  }
}

p_a2 <- p_a2 +
  labs(x = NULL, y = "logFC vs Healthy (donor-level)",
       title = "Severity-graded trajectory of hub candidates",
       subtitle = sprintf("monotonically up across IC \U2192 COPD \U2192 Sepsis \U2192 Sepsis+PNA (n=%d)",
                          sum(hit))) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8.5, colour = "grey35"),
        axis.text.x = element_text(size = 10))

## =========================================================================
## Step 6: output — combined figure (when has_pw) + A1/A2 single panels
##         (always, for per-panel assembly of Fig4)
## =========================================================================
if (has_pw) {
  p_comb <- p_a1 + p_a2 + plot_layout(widths = c(1.15, 1), guides = "collect")
  ggsave(file.path(OUT_DIR, "Fig4A_hub_identification.pdf"), p_comb,
         width = 12.5, height = 5.4)
  cat("saved: Fig4A_hub_identification.pdf\n")
  tryCatch({
    grDevices::png(file.path(OUT_DIR, "Fig4A_hub_identification.png"),
                   width = 12.5, height = 5.4, units = "in", res = 300)
    print(p_comb); grDevices::dev.off()
    cat("saved: Fig4A_hub_identification.png\n")
  }, error = function(e) cat("PNG write skipped (can be recovered from PDF with PyMuPDF):",
                             conditionMessage(e), "\n"))
}

## ---- A1 / A2 single panels (for independent assembly of Fig4 panel a / b) ----
ggsave(file.path(OUT_DIR, "Fig4A1_hub_scatter.pdf"), p_a1, width = 7, height = 5.4)
ggsave(file.path(OUT_DIR, "Fig4A2_hub_trajectory.pdf"), p_a2, width = 6, height = 5.4)
tryCatch({
  grDevices::png(file.path(OUT_DIR, "Fig4A1_hub_scatter.png"),
                 width = 7, height = 5.4, units = "in", res = 300)
  print(p_a1); grDevices::dev.off()
  grDevices::png(file.path(OUT_DIR, "Fig4A2_hub_trajectory.png"),
                 width = 6, height = 5.4, units = "in", res = 300)
  print(p_a2); grDevices::dev.off()
  cat("saved: Fig4A1 / Fig4A2 single panels {pdf,png}\n")
}, error = function(e) cat("A1/A2 PNG write skipped (can be recovered from PDF with PyMuPDF):",
                           conditionMessage(e), "\n"))

cat("\n===== Script 45 done (Fig4A hub identification) =====\n")
