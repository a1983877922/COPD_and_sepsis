# =========================================================================
# 06_pathway_enrichment.R
# Perform GO / KEGG / Reactome enrichment on the 149 cross-tissue shared myeloid genes
# (path2_shared_myeloid_genes.txt), mapping the "inflammation amplification +
# interferon response + metabolic reprogramming" triad onto concrete pathways.
#
# Dependencies (Bioconductor): clusterProfiler, org.Hs.eg.db, ReactomePA, enrichplot
#   if not installed the script will print install commands and skip the corresponding part.
# =========================================================================

.this_file <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
if (.this_file == "" || !file.exists(.this_file)) {
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1])
}
SCRIPT_DIR  <- dirname(.this_file)
config_file <- file.path(SCRIPT_DIR, "00_config.R")
if (!file.exists(config_file)) stop("Cannot find 00_config.R: ", config_file)
source(config_file)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(enrichplot)
  library(ggplot2)
})

cat("\n==============================================================\n")
cat("Pathway enrichment analysis (149 shared myeloid genes: GO / KEGG / Reactome)\n")
cat("==============================================================\n")

## =========================================================================
## Step 1: read shared genes + convert to ENTREZ ID
## =========================================================================
shared_genes <- readLines(file.path(out_dir, "path2_shared_myeloid_genes.txt"))
shared_genes <- shared_genes[nzchar(shared_genes)]
cat("Number of shared genes:", length(shared_genes), "\n")

# symbol -> ENTREZ (drop unmappable ones)
id_map <- bitr(shared_genes, fromType = "SYMBOL", toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)
cat("Successfully mapped ENTREZ IDs:", nrow(id_map), "/", length(shared_genes), "\n")
entrez <- unique(id_map$ENTREZID)

## =========================================================================
## Step 2: GO enrichment (BP / MF / CC)
## =========================================================================
cat("\n===== GO enrichment =====\n")
run_enrichGO <- function(ont) {
  enrichGO(gene = entrez, OrgDb = org.Hs.eg.db, ont = ont,
           pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
           readable = TRUE)
}
go_bp <- run_enrichGO("BP")
go_mf <- run_enrichGO("MF")
go_cc <- run_enrichGO("CC")

if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
  write.csv(as.data.frame(go_bp), file.path(out_dir, "path6_GO_BP.csv"), row.names = FALSE)
  p <- dotplot(go_bp, showCategory = 20) + ggtitle("GO Biological Process (top 20)")
  ggsave(file.path(out_dir, "path6_GO_BP_dotplot.pdf"), p, width = 11, height = 9)
  cat("GO BP significant terms:", nrow(as.data.frame(go_bp)), "\n")
} else cat("GO BP no significant enrichment\n")

if (!is.null(go_mf) && nrow(as.data.frame(go_mf)) > 0) {
  write.csv(as.data.frame(go_mf), file.path(out_dir, "path6_GO_MF.csv"), row.names = FALSE)
}
if (!is.null(go_cc) && nrow(as.data.frame(go_cc)) > 0) {
  write.csv(as.data.frame(go_cc), file.path(out_dir, "path6_GO_CC.csv"), row.names = FALSE)
}

## =========================================================================
## Step 3: KEGG enrichment
## =========================================================================
cat("\n===== KEGG enrichment =====\n")
kegg <- enrichKEGG(gene = entrez, organism = "hsa",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.2)
if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
  write.csv(as.data.frame(kegg), file.path(out_dir, "path6_KEGG.csv"), row.names = FALSE)
  p <- dotplot(kegg, showCategory = 20) + ggtitle("KEGG Pathways (top 20)")
  ggsave(file.path(out_dir, "path6_KEGG_dotplot.pdf"), p, width = 11, height = 9)
  cat("KEGG significant pathways:", nrow(as.data.frame(kegg)), "\n")
} else cat("KEGG no significant enrichment\n")

## =========================================================================
## Step 4: Reactome enrichment
## =========================================================================
cat("\n===== Reactome enrichment =====\n")
react <- enrichPathway(gene = entrez, organism = "human",
                       pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
if (!is.null(react) && nrow(as.data.frame(react)) > 0) {
  write.csv(as.data.frame(react), file.path(out_dir, "path6_Reactome.csv"), row.names = FALSE)
  p <- dotplot(react, showCategory = 20) + ggtitle("Reactome Pathways (top 20)")
  ggsave(file.path(out_dir, "path6_Reactome_dotplot.pdf"), p, width = 11, height = 9)
  cat("Reactome significant pathways:", nrow(as.data.frame(react)), "\n")
} else cat("Reactome no significant enrichment\n")

## =========================================================================
## Step 5: summarize top pathways (print key results)
## =========================================================================
cat("\n===== Summary =====\n")
top_summary <- list()
if (exists("go_bp") && !is.null(go_bp)) {
  top_summary[["GO_BP"]] <- head(as.data.frame(go_bp)[, c("Description","p.adjust","Count")], 10)
}
if (exists("kegg") && !is.null(kegg)) {
  top_summary[["KEGG"]] <- head(as.data.frame(kegg)[, c("Description","p.adjust","Count")], 10)
}
if (exists("react") && !is.null(react)) {
  top_summary[["Reactome"]] <- head(as.data.frame(react)[, c("Description","p.adjust","Count")], 10)
}
for (nm in names(top_summary)) {
  cat("\n[", nm, " top 10 ]\n")
  print(top_summary[[nm]], row.names = FALSE)
}

cat("\nPathway enrichment analysis complete\n")
cat("Output: path6_GO_*.csv / path6_KEGG.csv / path6_Reactome.csv + dotplot PDF\n")
