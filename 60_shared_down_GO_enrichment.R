# =============================================================================
# 60: GO enrichment of the shared down-regulated genes for Fig. S1
#      (v2: dual nominal + adjusted reporting, to settle the wording in Methods)
# Purpose: GO/KEGG of the blood-lung shared down-regulated genes (105). After BH
#       correction no pathway may remain significant (the Methods statement about
#       "basal transcription/RNA-processing" appears to come from the nominal
#       criterion), therefore:
#       * Run the full table with enrichGO (BP/CC/KEGG, pvalueCutoff=1, qvalueCutoff=1)
#       * Rank by nominal p and list the top hits, also giving the BH p.adjust, so
#         the authors can decide whether S1 shows "nominal enrichment concentrated
#         in housekeeping categories" or "no signal after correction".
# Output: path60_downGO_ranked_{BP,CC,KEGG}.csv (all terms, including pvalue/p.adjust)
#       FigS1_down105_GO_nominal.pdf/png  (top 10 nominal GO BP, annotated with post-correction status)
# Run:   Rscript 60_shared_down_GO_enrichment.R [OUT_DIR]  (needs clusterProfiler + org.Hs.eg.db)
# =============================================================================
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
if (length(args) >= 1 && nzchar(args[1])) script_dir <- args[1]
OUT_DIR <- if (length(args) >= 2 && nzchar(args[2])) args[2] else script_dir
out_dir <- OUT_DIR
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("===== 60 v2: shared down-regulated (105) GO/KEGG, nominal vs BH =====\n")
down <- trimws(readLines(file.path(out_dir, "path21_shared_down_genes.txt"),
                         warn = FALSE)); down <- down[down != ""]
cat("genes:", length(down), "\n")
map <- bitr(down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
eg <- unique(map$ENTREZID)
cat("ENTREZ:", length(eg), " (unmapped", length(down) - length(eg), ")\n")

rank_table <- function(res, ont) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
  df <- as.data.frame(res)
  df$ontology <- ont
  df <- df[, c("ontology", "ID", "Description", "GeneRatio", "BgRatio",
               "pvalue", "p.adjust", "qvalue", "geneID", "Count")]
  df[order(df$pvalue), ]
}
show_top <- function(df, label, k = 8) {
  if (is.null(df)) { cat(sprintf("[%s] no results\n", label)); return(NULL) }
  cat(sprintf("[%s] total terms %d; nominal p<0.05 %d; BH<0.05 %d\n", label,
              nrow(df), sum(df$pvalue < 0.05), sum(df$p.adjust < 0.05)))
  for (i in seq_len(min(k, nrow(df))))
    cat(sprintf("  top%02d p=%.2e padj=%.2e  %s\n", i, df$pvalue[i],
                df$p.adjust[i], substr(df$Description[i], 1, 62)))
  df
}

bp <- rank_table(tryCatch(enrichGO(eg, OrgDb = org.Hs.eg.db, ont = "BP",
                                   pAdjustMethod = "BH", pvalueCutoff = 1,
                                   qvalueCutoff = 1, readable = TRUE),
                          error = function(e) NULL), "BP")
cc <- rank_table(tryCatch(enrichGO(eg, OrgDb = org.Hs.eg.db, ont = "CC",
                                   pAdjustMethod = "BH", pvalueCutoff = 1,
                                   qvalueCutoff = 1, readable = TRUE),
                          error = function(e) NULL), "CC")
kg <- rank_table(tryCatch(enrichKEGG(eg, organism = "hsa", pvalueCutoff = 1,
                                     qvalueCutoff = 1),
                          error = function(e) NULL), "KEGG")
cat("\n")
show_top(bp, "GO BP")
cat("\n"); show_top(cc, "GO CC")
cat("\n"); show_top(kg, "KEGG")

if (!is.null(bp)) write.csv(bp, file.path(out_dir, "path60_downGO_ranked_BP.csv"),
                             row.names = FALSE)
if (!is.null(cc)) write.csv(cc, file.path(out_dir, "path60_downGO_ranked_CC.csv"),
                             row.names = FALSE)
if (!is.null(kg)) write.csv(kg, file.path(out_dir, "path60_downGO_ranked_KEGG.csv"),
                             row.names = FALSE)

## Top 10 nominal BP bars (colour indicates whether BH is passed), to see at a
## glance whether the signal is concentrated in housekeeping categories
if (!is.null(bp)) {
  top <- head(bp, 10)
  top$BH <- ifelse(top$p.adjust < 0.05, "BH<0.05", "n.s. after BH")
  top$Description <- factor(top$Description,
                            levels = rev(top$Description[order(top$pvalue)]))
  p <- ggplot(top, aes(-log10(pvalue), Description, fill = BH)) +
    geom_bar(stat = "identity", width = 0.6, alpha = 0.85) +
    scale_fill_manual(values = c("BH<0.05" = "#B2182B", "n.s. after BH" = "#4C8BF5")) +
    geom_text(aes(x = 0.05, label = Description), hjust = 0, size = 3.4) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.06))) +
    labs(x = "-Log10(nominal P)", y = "",
         title = "GO BP of shared down-regulated genes (n = 105)\n(nominal ranking; colour = status after BH correction)",
         fill = NULL) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_blank(), axis.ticks.length.y = unit(0, "cm"),
          axis.line.y = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(size = 11, hjust = 0.5))
  ggsave(file.path(out_dir, "FigS1_down105_GO_nominal.pdf"), p, width = 8, height = 6)
  ggsave(file.path(out_dir, "FigS1_down105_GO_nominal.png"), p, width = 8, height = 6, dpi = 300)
  cat("\nfigure written: FigS1_down105_GO_nominal.pdf/.png (top 10 nominal BP)\n")
}
cat("===== script 60 v2 finished =====\n")
