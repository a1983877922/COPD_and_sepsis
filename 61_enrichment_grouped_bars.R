# =============================================================================
# Script 61: grouped enrichment bar chart style   (following the grouped bar
# style of the "CD14/CD16 monocyte KEGG" plot)
# Key points of the plotting style (matching the example you sent): top_n
# pathways per group, grouped horizontal bars, two-colour fill, with the
#   Description text inset on the left of each bar (hjust=0) plus optional
#   coloured small geneID text, y-axis text left blank,
#   no legend, theme_classic + compact theme, saved with manually tuned width/height.
# Applied to:
#   A) Fig. S1 candidate: the 105 genes co-down-regulated in blood and lung
#        (path21_shared_down_genes.txt)
#       grouping = GO BP vs KEGG -> shows that the co-down-regulated pathways
#       concentrate in basal/housekeeping terms (echoing the
#       Methods "not pursued" rationale)
#   B) Comparison: the 149-gene up-regulated program (reads the ready-made
#        results in path6_GO_BP.csv + path6_KEGG.csv)
#       grouping = GO BP vs KEGG, same style -> lets the authors decide whether
#       to replace the old dotplot
# Output: FigS1_down105_GO_KEGG_bars.{pdf,png} / path61_*_bars_data.csv
# Run: Rscript this script [OUT_DIR]
#       (requires clusterProfiler, org.Hs.eg.db, ggplot2)
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

## ---------------------------------------------------------------- plotting function
enrich_bars <- function(d, cluster_col = "Cluster", top_n = 5,
                        metric = "p.adjust", colors = c("#6bb9d2", "#d55640"),
                        title = NULL, axis_label = NULL, show_geneid = TRUE,
                        geneid_col = "geneID", label_x = 0.05,
                        png = TRUE, pdf = TRUE, base = "path61_enrich_bars",
                        width = 6, height = NULL) {
  stopifnot(all(c(cluster_col, "Description", metric) %in% colnames(d)))
  d <- as.data.frame(d)
  d$m <- -log10(d[[metric]])
  d <- d[is.finite(d$m), ]
  cl <- unique(as.character(d[[cluster_col]]))
  if (length(colors) < length(cl)) colors <- rep(colors, length.out = length(cl))
  names(colors) <- cl
  ## Take the most significant top_n within each group; cluster blocks ordered by cl (block 1 on top)
  keep <- do.call(rbind, lapply(cl, function(cc) {
    x <- d[d[[cluster_col]] == cc, ]
    x <- x[order(x$m, decreasing = TRUE), ]
    x[seq_len(min(top_n, nrow(x))), , drop = FALSE]
  }))
  keep[[cluster_col]] <- factor(keep[[cluster_col]], levels = cl)
  ## y axis: within a group, from most to least significant top to bottom; Description de-duplicated, order preserved
  keep$Description <- factor(as.character(keep$Description),
                             levels = rev(unique(as.character(keep$Description))))
  nbar <- nrow(keep)
  if (is.null(height)) height <- 1.3 + 0.55 * nbar
  g <- ggplot(keep, aes(x = m, y = Description, fill = .data[[cluster_col]])) +
    geom_bar(stat = "identity", width = 0.5, alpha = 0.8) +
    scale_fill_manual(values = colors) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.06))) +
    geom_text(aes(x = label_x, label = Description), size = 3.8, hjust = 0)
  if (show_geneid && geneid_col %in% colnames(keep)) {
    keep$.gcol <- colors[as.character(keep[[cluster_col]])]
    g <- g + geom_text(data = keep,
                       aes(x = label_x, label = .data[[geneid_col]]),
                       size = 2.6, hjust = 0, vjust = 2.4,
                       colour = keep$.gcol)
  }
  if (is.null(axis_label))
    axis_label <- paste(rev(cl), collapse = "                 ")
  g <- g + labs(x = sprintf("-Log10(%s)", sub("^p\\.", "adjusted ", metric)),
                y = axis_label, title = title) +
    theme_classic(base_size = 11) +
    theme(axis.title = element_text(size = 13),
          axis.text = element_text(size = 11),
          axis.text.y = element_blank(),
          axis.ticks.length.y = unit(0, "cm"),
          axis.line.y = element_blank(),
          plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
          legend.position = "none",
          plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5))
  if (pdf) ggsave(file.path(out_dir, paste0(base, ".pdf")), g,
                  width = width, height = height)
  if (png) ggsave(file.path(out_dir, paste0(base, ".png")), g,
                  width = width, height = height, dpi = 300)
  write.csv(keep, file.path(out_dir, paste0(base, "_data.csv")),
            row.names = FALSE)
  cat("  [bar] rows:", nbar, "| groups:", paste(cl, collapse = " / "),
      "| output:", base, ".pdf/.png/.csv\n")
  g
}

## ---------------------------------------------------------------- A) down-regulated 105
cat("\n===== A) 105 genes co-down-regulated in blood and lung: GO BP vs KEGG (Fig S1 candidate) =====\n")
down <- trimws(readLines(file.path(out_dir, "path21_shared_down_genes.txt"),
                         warn = FALSE)); down <- down[down != ""]
cat("Number of genes:", length(down), "\n")
map <- bitr(down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
eg <- unique(map$ENTREZID)
cat("ENTREZ:", length(eg), "\n")
build_df <- function(res, cluster) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
  df <- as.data.frame(res); df$Cluster <- cluster
  df[, c("Cluster", "Description", "pvalue", "p.adjust", "geneID", "Count")]
}
bp_d  <- build_df(tryCatch(enrichGO(eg, OrgDb = org.Hs.eg.db, ont = "BP",
                                    pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.2, readable = TRUE),
                           error = function(e) NULL), "GO BP")
kg_d  <- build_df(tryCatch(enrichKEGG(eg, organism = "hsa", pvalueCutoff = 0.05,
                                      qvalueCutoff = 0.2),
                           error = function(e) NULL), "KEGG")
cat("BP significant:", if (is.null(bp_d)) 0 else nrow(bp_d),
    "| KEGG significant:", if (is.null(kg_d)) 0 else nrow(kg_d), "\n")
if (!is.null(bp_d)) { cat("  BP top5:\n"); print(head(bp_d$Description, 5)) }
if (!is.null(kg_d)) { cat("  KEGG top5:\n"); print(head(kg_d$Description, 5)) }
dA <- rbind(bp_d, kg_d)
if (!is.null(dA) && nrow(dA) > 0)
  enrich_bars(dA, top_n = 5, metric = "p.adjust",
              colors = c("GO BP" = "#6bb9d2", "KEGG" = "#d55640"),
              title = "Shared down-regulated genes (blood & lung, n = 105)",
              axis_label = "KEGG                  GO BP",
              base = "FigS1_down105_GO_KEGG_bars", width = 7, show_geneid = TRUE)

## ---------------------------------------------------------------- B) up-regulated 149 comparison
cat("\n===== B) 149-gene up-regulated program: GO BP vs KEGG (ready-made results from script 6, same style) =====\n")
f_bp <- file.path(out_dir, "path6_GO_BP.csv"); f_kg <- file.path(out_dir, "path6_KEGG.csv")
if (file.exists(f_bp) && file.exists(f_kg)) {
  bp_u <- read.csv(f_bp, stringsAsFactors = FALSE); bp_u$Cluster <- "GO BP"
  kg_u <- read.csv(f_kg, stringsAsFactors = FALSE); kg_u$Cluster <- "KEGG"
  bp_u <- bp_u[bp_u$p.adjust <= 0.05 & !is.na(bp_u$p.adjust), ]
  kg_u <- kg_u[kg_u$p.adjust <= 0.05 & !is.na(kg_u$p.adjust), ]
  dU <- rbind(bp_u[, c("Cluster", "Description", "pvalue", "p.adjust",
                       "geneID", "Count")],
              kg_u[, c("Cluster", "Description", "pvalue", "p.adjust",
                       "geneID", "Count")])
  cat("BP significant:", nrow(bp_u), "| KEGG significant:", nrow(kg_u), "\n")
  enrich_bars(dU, top_n = 5, metric = "p.adjust",
              colors = c("GO BP" = "#6bb9d2", "KEGG" = "#d55640"),
              title = "Shared up-regulated program (149 genes)",
              axis_label = "KEGG                  GO BP",
              base = "path61_up149_GO_KEGG_bars", width = 7, show_geneid = TRUE)
} else {
  cat("Missing path6_GO_BP.csv / path6_KEGG.csv, skipping B (A is unaffected)\n")
}
cat("\n===== Script 61 done =====\n")
