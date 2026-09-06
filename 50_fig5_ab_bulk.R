# =========================================================================
# 50: bulk re-render of Figure 5 panels A/B
#
# Purpose: panels A / B of Figure 5 (bulk WGCNA side); reads the intermediate object
#   saved by 12_wgcna.R and re-renders them at a uniform size (base graphics, not ggplot):
#   A  Fig5A_bulk_dendrogram          bulk whole-transcriptome WGCNA gene dendrogram + module colour band
#   B  Fig5B_bulk_module_trait_heatmap  module eigengene x sepsis correlation heatmap (WGCNA convention)
#
# Requires: run the modified 12_wgcna.R first (it newly saves path12_wgcna_obj.rds), containing
#       geneTree / dynamicColors / MEs / moduleTraitCor / moduleTraitP
# Output:
#   Fig5A_bulk_dendrogram.pdf/.png
#   Fig5B_bulk_module_trait_heatmap.pdf/.png
# Server: Rscript 50_fig5_ab_bulk.R
# =========================================================================

.this_file <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
if (.this_file == "" || !file.exists(.this_file)) {
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1])
}
SCRIPT_DIR <- dirname(.this_file)
config_file <- file.path(SCRIPT_DIR, "00_config.R")
if (!file.exists(config_file)) stop("cannot find 00_config.R: ", config_file)
source(config_file)

suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library(ggplot2))

cat("\n==============================================================\n")
cat("Fig5 A/B: re-render of the bulk WGCNA panels\n")
cat("==============================================================\n")

obj_file <- file.path(out_dir, "path12_wgcna_obj.rds")
if (!file.exists(obj_file)) {
  stop("Missing path12_wgcna_obj.rds. Please run the modified 12_wgcna.R first to generate this checkpoint.")
}
w <- readRDS(obj_file)
geneTree <- w$geneTree; colors <- w$dynamicColors
MEs <- w$MEs; mTC <- w$moduleTraitCor; mTP <- w$moduleTraitP
cat("modules:", length(unique(colors)), "| genes:", ncol(w$datExpr), "\n")

## ---- A: bulk dendrogram ----
cat("\n===== A: dendrogram =====\n")
pdf(file.path(out_dir, "Fig5A_bulk_dendrogram.pdf"), width = 12, height = 6)
plotDendroAndColors(geneTree, colors, "Modules",
                    dendroLabels = FALSE, addGuide = TRUE, guideHang = 0.05,
                    main = "Bulk WGCNA modules (whole blood, GSE66099)")
dev.off()
png(file.path(out_dir, "Fig5A_bulk_dendrogram.png"),
    width = 12, height = 6, units = "in", res = 300)
plotDendroAndColors(geneTree, colors, "Modules",
                    dendroLabels = FALSE, addGuide = TRUE, guideHang = 0.05,
                    main = "Bulk WGCNA modules (whole blood, GSE66099)")
dev.off()
cat("saved: Fig5A_bulk_dendrogram.{pdf,png}\n")

## ---- B: module-trait heatmap ----
cat("\n===== B: module-trait heatmap =====\n")
mTC2 <- mTC[, 1, drop = FALSE]
colnames(mTC2) <- "Sepsis"
mTP2 <- mTP[, 1, drop = FALSE]
textM <- matrix(sprintf("%.2f", mTC2), nrow = nrow(mTC2),
                dimnames = list(rownames(mTC2), "Sepsis"))
# Significance markers: p<0.001 *** / <0.01 ** / <0.05 *
sig <- ifelse(mTP2 < 0.001, "***", ifelse(mTP2 < 0.01, "**",
              ifelse(mTP2 < 0.05, "*", "")))
textM <- matrix(paste0(sprintf("%.2f", mTC2), "\n", sig),
                nrow = nrow(mTC2), dimnames = dimnames(mTC2))

pdf(file.path(out_dir, "Fig5B_bulk_module_trait_heatmap.pdf"),
    width = 4.5, height = 7)
labeledHeatmap(Matrix = mTC2,
               xLabels = "Sepsis", yLabels = rownames(mTC2),
               colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textM, setStdMargins = FALSE,
               cex.text = 0.6, cex.lab = 0.9, cex.axis = 0.7,
               main = "Module eigengene correlation with sepsis\n(GSE66099 whole blood)")
dev.off()
png(file.path(out_dir, "Fig5B_bulk_module_trait_heatmap.png"),
    width = 4.5, height = 7, units = "in", res = 300)
labeledHeatmap(Matrix = mTC2,
               xLabels = "Sepsis", yLabels = rownames(mTC2),
               colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textM, setStdMargins = FALSE,
               cex.text = 0.6, cex.lab = 0.9, cex.axis = 0.7,
               main = "Module eigengene correlation with sepsis\n(GSE66099 whole blood)")
dev.off()
cat("saved: Fig5B_bulk_module_trait_heatmap.{pdf,png}\n")

## ---- data behind the figures ----
saveRDS(list(geneTree = geneTree, colors = colors,
             moduleTraitCor = mTC2, moduleTraitP = mTP2, softPower = w$softPower),
        file.path(out_dir, "path50_fig5ab_plotdata.rds"))
cat("saved path50_fig5ab_plotdata.rds\n")
cat("\n===== script 50 done =====\n")
