# =========================================================================
# 12_WGCNA.R  (bulk whole-transcriptome WGCNA + DEG mapping/intersection)
# Build WGCNA co-expression modules using all expressed genes (~5000 after low-variance
# filtering) in the sepsis bulk cohort (GSE66099, 246 samples), and map DEGs:
#   - module-sepsis (Sepsis vs Control) association
#   - map DEGs to modules -> "disease-related module ∩ significant DEG" hub key genes
#   - map the 149 shared myeloid genes to modules
#
# Dependencies: WGCNA (required), GEOquery (required), limma (DEG, optional, fallback wilcox)
# =========================================================================

.this_file <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
if (.this_file == "" || !file.exists(.this_file)) {
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1])
}
SCRIPT_DIR  <- dirname(.this_file)
config_file <- file.path(SCRIPT_DIR, "00_config.R")
if (!file.exists(config_file)) stop("cannot find 00_config.R: ", config_file)
source(config_file)

if (!requireNamespace("WGCNA", quietly = TRUE)) {
  stop("WGCNA is not installed. Please run first: BiocManager::install('WGCNA')")
}
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  stop("GEOquery is not installed. Please run first: BiocManager::install('GEOquery')")
}
suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library(GEOquery))
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

cat("\n==============================================================\n")
cat("bulk whole-transcriptome WGCNA + DEG mapping (GSE66099)\n")
cat("==============================================================\n")

## =========================================================================
## Step 1: read the GSE66099 full expression matrix + grouping
## =========================================================================
cat("\n===== Step 1: read data =====\n")

matrix_file <- file.path(out_dir, "GSE66099_matrix.txt.gz")
if (file.exists(matrix_file)) {
  gse <- getGEO(filename = matrix_file, getGPL = TRUE)
} else {
  gse <- getGEO("GSE66099", GSEMatrix = TRUE, getGPL = TRUE); gse <- gse[[1]]
}
expr <- exprs(gse); pd <- pData(gse)

fd <- fData(gse)
sym <- as.character(fd[["Gene Symbol"]]); names(sym) <- rownames(fd)
keep <- sym != "" & !is.na(sym) & sym != "---"
expr <- expr[keep, , drop = FALSE]; sym <- sym[keep]
sym_uniq <- unique(sym)
expr_sym <- matrix(NA, nrow = length(sym_uniq), ncol = ncol(expr),
                   dimnames = list(sym_uniq, colnames(expr)))
for (g in sym_uniq) {
  idx <- which(sym == g)
  expr_sym[g, ] <- if (length(idx) == 1) expr[idx, ] else colMeans(expr[idx, , drop = FALSE])
}
expr <- expr_sym
cat("Full expression matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")

group <- as.character(pd[["disease:ch1"]]); names(group) <- rownames(pd)
group[group == "SepticShock"] <- "Sepsis"
common <- intersect(colnames(expr), names(group))
expr <- expr[, common, drop = FALSE]; group <- group[common]
keep_s <- group %in% c("Sepsis", "Control")
expr <- expr[, keep_s, drop = FALSE]; group <- group[keep_s]
cat("Samples:", ncol(expr), "\n"); print(table(group))

## =========================================================================
## Step 2: filter low-variance genes (keep the top 5000 most variable, to limit TOM cost)
## =========================================================================
cat("\n===== Step 2: filter high-variance genes =====\n")

gene_var <- sort(apply(expr, 1, var), decreasing = TRUE)
top_genes <- names(gene_var)[1:min(5000, length(gene_var))]
datExpr <- t(expr[top_genes, , drop = FALSE])   # samples x genes
cat("Genes used for WGCNA:", ncol(datExpr), "\n")

## =========================================================================
## Step 3: WGCNA (soft threshold + TOM + module detection)
## =========================================================================
cat("\n===== Step 3: WGCNA =====\n")

powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
softPower <- 6
sft <- tryCatch(pickSoftThreshold(datExpr, powerVector = powers,
                                  networkType = "signed", verbose = 0),
                error = function(e) NULL)
if (!is.null(sft) && !is.na(sft$powerEstimate)) {
  softPower <- sft$powerEstimate
  cat("softPower =", softPower, "\n")
} else cat("using the default softPower = 6\n")

cat("Building the network + module detection (TOM, may take a few minutes)...\n")
adjacency <- adjacency(datExpr, power = softPower, type = "signed")
TOM <- TOMsimilarity(adjacency); rm(adjacency); gc()
dissTOM <- 1 - TOM; rm(TOM); gc()
geneTree <- hclust(as.dist(dissTOM), method = "average")

dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = 30)
dynamicColors <- labels2colors(dynamicMods)
cat("Number of modules:", length(unique(dynamicColors)), "\n")
print(table(dynamicColors))

MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- orderMEs(MEList$eigengenes)

## =========================================================================
## Step 4: module-trait association (Sepsis)
## =========================================================================
cat("\n===== Step 4: module-trait association =====\n")

trait <- ifelse(group == "Sepsis", 1, 0)
moduleTraitCor <- cor(MEs, trait, use = "p")
moduleTraitP <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
mt <- data.frame(module = sub("^ME", "", rownames(moduleTraitCor)),
                 cor_sepsis = moduleTraitCor[, 1],
                 pvalue = moduleTraitP[, 1])
mt <- mt[order(mt$pvalue), ]
print(mt)
write.csv(mt, file.path(out_dir, "path12_bulk_module_trait.csv"), row.names = FALSE)

sig_mods <- mt$module[mt$pvalue < 0.05]
cat("Modules significantly associated with sepsis (p<0.05):", paste(sig_mods, collapse = ", "), "\n")

gene_module <- data.frame(gene = colnames(datExpr), module = dynamicColors,
                          stringsAsFactors = FALSE)
write.csv(gene_module, file.path(out_dir, "path12_bulk_gene_module.csv"),
          row.names = FALSE)

## ---- [for Fig5 re-rendering] save the WGCNA intermediate object (dendrogram/heatmap can be re-rendered without re-running) ----
saveRDS(list(geneTree = geneTree, dynamicColors = dynamicColors,
             datExpr = datExpr, MEs = MEs,
             moduleTraitCor = moduleTraitCor, moduleTraitP = moduleTraitP,
             softPower = softPower),
        file.path(out_dir, "path12_wgcna_obj.rds"))
cat("Saved WGCNA object: path12_wgcna_obj.rds (for re-rendering Fig5 A/B)\n")

## =========================================================================
## Step 5: DEG mapping/intersection
## =========================================================================
cat("\n===== Step 5: DEG mapping =====\n")

if (requireNamespace("limma", quietly = TRUE)) {
  suppressPackageStartupMessages(library(limma))
  design <- model.matrix(~ group)
  fit <- lmFit(expr, design)
  fit <- eBayes(fit)
  deg_res <- topTable(fit, coef = "groupSepsis", number = Inf,
                      sort.by = "P", adjust.method = "BH")
  deg_res$gene <- rownames(deg_res)
  write.csv(deg_res, file.path(out_dir, "path12_bulk_DEG_Sepsis_vs_Control.csv"),
            row.names = FALSE)
  cat("limma DEG done, number of significant genes:", sum(deg_res$adj.P.Val < 0.05), "\n")
} else {
  cat("limma not installed, using wilcox as an approximate DEG test\n")
  pvals <- apply(expr, 1, function(g) {
    tryCatch(wilcox.test(g[group == "Sepsis"], g[group == "Control"])$p.value,
             error = function(e) NA)
  })
  deg_res <- data.frame(gene = names(pvals), pvalue = pvals)
  deg_res <- deg_res[!is.na(deg_res$pvalue), ]
  deg_res$adj.P.Val <- p.adjust(deg_res$pvalue, "BH")
  cat("wilcox DEG done, number of significant genes:", sum(deg_res$adj.P.Val < 0.05), "\n")
}

deg_sig <- deg_res$gene[deg_res$adj.P.Val < 0.05]
deg_module <- gene_module[gene_module$gene %in% deg_sig, ]
cat("Number of significant DEGs:", length(deg_sig), " of these inside WGCNA modules:", nrow(deg_module), "\n")

# disease-related module ∩ significant DEG = key hub genes
hub_genes <- gene_module[gene_module$module %in% sig_mods &
                         gene_module$gene %in% deg_sig, ]
cat("Number of key genes in disease-related module ∩ significant DEG:", nrow(hub_genes), "\n")
write.csv(hub_genes, file.path(out_dir, "path12_bulk_hub_genes.csv"), row.names = FALSE)
if (nrow(hub_genes) > 0) print(head(hub_genes, 40))

# map the 149 shared genes to modules
genes_149 <- readLines(file.path(out_dir, "path2_shared_myeloid_genes.txt"))
genes_149 <- genes_149[genes_149 != ""]
mapped_149 <- gene_module[gene_module$gene %in% genes_149, ]
cat("Of the 149 shared genes, mapped to modules:", nrow(mapped_149), "\n")
print(table(mapped_149$module))
write.csv(mapped_149, file.path(out_dir, "path12_bulk_149gene_module.csv"),
          row.names = FALSE)

## =========================================================================
## Step 6: visualization
## =========================================================================
cat("\n===== Step 6: visualization =====\n")

pdf(file.path(out_dir, "path12_bulk_dendrogram.pdf"), width = 14, height = 8)
plotDendroAndColors(geneTree, dynamicColors, "Modules",
                    dendroLabels = FALSE, addGuide = TRUE,
                    guideHang = 0.05, main = "bulk whole-transcriptome WGCNA modules")
dev.off()

pdf(file.path(out_dir, "path12_bulk_module_trait_heatmap.pdf"), width = 5, height = 9)
labeledHeatmap(Matrix = moduleTraitCor, xLabels = "Sepsis",
               yLabels = rownames(moduleTraitCor),
               colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = round(moduleTraitCor, 2),
               setStdMargins = FALSE, cex.text = 0.7,
               main = "Module-Sepsis correlation")
dev.off()

cat("\nbulk whole-transcriptome WGCNA done\n")
cat("Outputs: path12_bulk_module_trait.csv / path12_bulk_gene_module.csv /\n")
cat("         path12_bulk_DEG_*.csv / path12_bulk_hub_genes.csv /\n")
cat("         path12_bulk_149gene_module.csv / two PDFs\n")
cat("Interpretation: look at the hub genes in the disease-related module ∩ significant DEG set, and at which modules the 149 genes fall into\n")
