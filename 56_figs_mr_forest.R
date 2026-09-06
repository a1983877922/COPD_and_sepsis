# =============================================================================
# 56_FigS_MR_summary_forest_plot.R
# Purpose: Summarize path55_mr_results.csv (11 genes x 2 outcomes two-sample MR results) into
#       a supplementary forest plot (Fig S): one OR (95% CI) estimate per gene x outcome.
# Method: For each gene x outcome take the main method -- prefer IVW (>=2 SNP), otherwise Wald ratio (1 SNP).
# Expected result: all null (no significance after Bonferroni); plot used to show "expression MR does not suggest causality".
# Run: Rscript 56_FigS_MR_summary_forest_plot.R [OUT_DIR]   (pure CSV, seconds, no network needed)
# Output: FigS_MR_forest.pdf/.png + path56_mr_forest_data.csv + path56_mr_forest_plotdata.rds
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
if (length(args) >= 1 && nzchar(args[1])) script_dir <- args[1]
OUT_DIR <- if (length(args) >= 2 && nzchar(args[2])) args[2] else script_dir
out_dir <- OUT_DIR
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr)
})

csv_file <- file.path(out_dir, "path55_mr_results.csv")
stopifnot(file.exists(csv_file))
res <- read.csv(csv_file, stringsAsFactors = FALSE)
cat("Read MR results:", nrow(res), "rows\n")

## Gene order (consistent with 55 script GENES, SOCS3 main gene first)
gene_order <- c("SOCS3","STAT1","STAT3","JAK1","JAK2","TYK2",
                "IL6ST","ISG15","MX1","GBP1","IFI44L")
gene_order <- gene_order[gene_order %in% unique(res$gene)]

## For each gene x outcome take the main method
res$key <- paste(res$gene, res$outcome, sep = "__")
pick_main <- function(d) {
  if (any(d$method == "Inverse variance weighted")) {
    return(d[d$method == "Inverse variance weighted", ][1, ])
  }
  if (any(d$method == "Wald ratio")) {
    return(d[d$method == "Wald ratio", ][1, ])
  }
  d[which.max(abs(d$b)), ]
}
main <- do.call(rbind, lapply(split(res, res$key), pick_main))
rownames(main) <- NULL

main$OR <- exp(main$b)
main$lo <- exp(main$b - 1.96 * main$se)
main$hi <- exp(main$b + 1.96 * main$se)
main$p_lab <- formatC(main$pval, format = "e", digits = 2)
main$star <- ifelse(main$pval < 0.05, "*",
             ifelse(main$significant_bonf, "**", ""))
main$gene <- factor(main$gene, levels = rev(gene_order))
main$outcome <- factor(main$outcome,
                       levels = c("sepsis", "copd"),
                       labels = c("Sepsis (ieu-b-4980)", "COPD (finn-b-J10_COPD)"))

## Output data
write.csv(main[, c("gene","outcome","method","nsnp","b","se","pval",
                   "OR","lo","hi","significant_bonf")],
          file.path(out_dir, "path56_mr_forest_data.csv"), row.names = FALSE)

## Forest plot
p <- ggplot(main, aes(x = OR, y = gene)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey60", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.5,
                 colour = "grey35") +
  geom_point(aes(colour = method), size = 2.6) +
  geom_text(aes(label = p_lab), hjust = -0.2, vjust = -0.9, size = 2.6,
            colour = "grey30") +
  scale_colour_manual(values = c("Inverse variance weighted" = "#B2182B",
                                 "Wald ratio" = "#2166AC"),
                      name = "Estimator") +
  scale_x_log10() +
  facet_wrap(~ outcome, ncol = 2, scales = "free_y") +
  labs(title = "Two-sample MR: gene expression (eQTLGen) on disease risk",
       subtitle = "Points = OR per SD-increase in genetically-predicted expression; line = 95% CI; label = P.\nNo gene survived Bonferroni correction (threshold 0.0023); SOCS3 had a single cis-eQTL instrument.",
       x = "Odds ratio (95% CI), log scale", y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        plot.subtitle = element_text(size = 8, colour = "grey30"))

ggsave(file.path(out_dir, "FigS_MR_forest.pdf"), p, width = 9, height = 6.5)
cat("Saved: FigS_MR_forest.pdf\n")
tryCatch({
  ggsave(file.path(out_dir, "FigS_MR_forest.png"), p,
         width = 9, height = 6.5, dpi = 300, device = grDevices::png)
  cat("Saved: FigS_MR_forest.png\n")
}, error = function(e) cat("PNG write skipped:", conditionMessage(e), "\n"))

saveRDS(list(forest_data = main, plot = p),
        file.path(out_dir, "path56_mr_forest_plotdata.rds"))
cat("\n===== Script 56 complete =====\n")
