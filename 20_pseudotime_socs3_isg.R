# =========================================================================
# 20_pseudotime_socs3_isg.R (pseudotime SOCS3/ISG dynamics)
#
# Purpose: respond to the reviewer's criticism -- "you only place sepsis pDC at
# the activated-exhausted endpoint, without looking at SOCS3/ISG dynamics along
# the pseudotime axis, so you cannot prove SOCS3 is progressively up-regulated
# with the activation process".
#
# Approach: build a Monocle3 trajectory for pDC (root=Healthy), bin along
# pseudotime, characterize the expression dynamics of SOCS3 and ISG + Spearman
# correlation, presenting DYNAMIC evidence of the "activation-suppression
# paradox" (SOCS3 and ISG co-upregulate with activation).
#
# Input: path11_pdc_subset.rds (saved by script 11, includes counts+data+group)
# Output: path20_ series
# =========================================================================

.this_file <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
if (.this_file == "" || !file.exists(.this_file)) {
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1])
}
SCRIPT_DIR  <- dirname(.this_file)
config_file <- file.path(SCRIPT_DIR, "00_config.R")
if (!file.exists(config_file)) stop("00_config.R not found: ", config_file)
source(config_file)

suppressPackageStartupMessages(library(ggplot2))

cat("\n==============================================================\n")
cat("Pseudotime SOCS3/ISG dynamics analysis\n")
cat("==============================================================\n")

# Genes of interest: SOCS3/SOCS1 (brake) + ISG (IFN readiness) + TF (driver)
genes_of_interest <- c("SOCS3", "SOCS1",
                       "ISG15", "MX1", "MX2", "OAS2", "ISG20",
                       "IFI44L", "IFIT2", "IFIT3",
                       "IRF7", "STAT1", "IRF9")

## =========================================================================
## Step 1: Read pDC object
## =========================================================================
cat("\n===== Step 1: Read pDC object =====\n")

pdc_file <- file.path(out_dir, "path11_pdc_subset.rds")
if (!file.exists(pdc_file)) stop(pdc_file, " not found, please run script 11 first")
pdc <- readRDS(pdc_file)
cat("pDC object:", ncol(pdc), "cells\n")
print(table(pdc$group))

avail <- genes_of_interest[genes_of_interest %in% rownames(pdc)]
cat("Genes of interest available:", length(avail), "\n")
genes_of_interest <- avail

## =========================================================================
## Step 2: Monocle3 trajectory (root=Healthy)
## =========================================================================
cat("\n===== Step 2: Monocle3 trajectory =====\n")

# Ensure dependencies
if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
  install.packages("SeuratWrappers", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("monocle3", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("monocle3", update = FALSE, ask = FALSE)
}
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(SeuratWrappers))

# Reuse the trajectory already built in script 11: read the saved cds if present, otherwise rebuild
cds_file <- file.path(out_dir, "path11b_pdc_cds.rds")
if (file.exists(cds_file)) {
  cds <- readRDS(cds_file)
  cat("Reusing saved cds: path11b_pdc_cds.rds\n")
} else {
  cat("Rebuilding trajectory (as.cell_data_set -> cluster_cells -> learn_graph -> order_cells)\n")
  cds <- as.cell_data_set(pdc)
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cds <- learn_graph(cds)
  root_cells <- colnames(cds)[colData(cds)$group == "Healthy"]
  if (length(root_cells) == 0) root_cells <- colnames(cds)[1]
  cds <- order_cells(cds, root_cells = root_cells, reduction_method = "UMAP")
  saveRDS(cds, cds_file)
  cat("cds saved:", cds_file, "\n")
}

pseudo <- pseudotime(cds)
pseudo <- pseudo[!is.na(pseudo)]
cat("Pseudotime range:", round(min(pseudo), 2), "->", round(max(pseudo), 2),
    " (", length(pseudo), "cells )\n")

## =========================================================================
## Step 3: Bin along pseudotime, compute SOCS3/ISG expression dynamics
## =========================================================================
cat("\n===== Step 3: Binned expression dynamics =====\n")

cells <- names(pseudo)
# Extract the data layer (log-normalized)
expr <- GetAssayData(pdc, assay = "RNA", layer = "data")
expr <- as.matrix(expr[genes_of_interest, cells, drop = FALSE])

n_bins <- 10
bin_id <- cut(pseudo, breaks = n_bins, labels = FALSE, include.lowest = TRUE)
# Use sapply to iterate seq_len(n_bins) instead of tapply: tapply skips empty
# bins, returning length < n_bins, which mismatches bin_mean (always n_bins rows) and errors.
bin_center <- sapply(seq_len(n_bins), function(b) {
  idx <- which(bin_id == b)
  if (length(idx) == 0) return(NA_real_)
  median(pseudo[idx])
})

# Mean expression per bin
bin_mean <- sapply(seq_len(n_bins), function(b) {
  idx <- which(bin_id == b)
  if (length(idx) == 0) return(rep(NA, length(genes_of_interest)))
  rowMeans(expr[, idx, drop = FALSE])
})
bin_mean <- t(bin_mean)  # n_bins x n_genes
colnames(bin_mean) <- genes_of_interest
rownames(bin_mean) <- paste0("bin", seq_len(n_bins))

cat("Binned mean expression (10 bins x genes):\n")
print(round(bin_mean, 3))

write.csv(bin_mean, file.path(out_dir, "path20_pseudotime_bin_expression.csv"))

## =========================================================================
## Step 4: Spearman correlation (gene vs pseudotime; SOCS3 vs ISG)
## =========================================================================
cat("\n===== Step 4: Correlation =====\n")

# Spearman correlation between gene expression and pseudotime
gene_pt_cor <- sapply(genes_of_interest, function(g) {
  x <- expr[g, ]
  if (all(x == x[1])) return(NA_real_)
  cor(x, pseudo, method = "spearman")
})
cat("Gene vs pseudotime Spearman correlation:\n")
print(round(gene_pt_cor, 3))

# Spearman correlation between SOCS3 and each ISG (at single-cell level)
socs3_expr <- expr["SOCS3", ]
socs3_isg_cor <- sapply(setdiff(genes_of_interest, "SOCS3"), function(g) {
  cor(socs3_expr, expr[g, ], method = "spearman")
})
cat("\nSOCS3 vs each gene Spearman correlation:\n")
print(round(socs3_isg_cor, 3))

cor_df <- data.frame(
  gene = genes_of_interest,
  rho_with_pseudotime = gene_pt_cor
)
write.csv(cor_df, file.path(out_dir, "path20_gene_pseudotime_cor.csv"),
          row.names = FALSE)

## =========================================================================
## Step 5: Visualization
## =========================================================================
cat("\n===== Step 5: Visualization =====\n")

# 5.1 Dynamic curves of SOCS3 + representative ISGs along pseudotime (binned)
plot_genes <- intersect(c("SOCS3", "ISG15", "MX1", "OAS2", "IFI44L", "IRF7"),
                        genes_of_interest)
df_line <- data.frame(
  bin = rep(seq_len(n_bins), length(plot_genes)),
  bin_center = rep(bin_center, length(plot_genes)),
  gene = rep(plot_genes, each = n_bins),
  expr = as.vector(bin_mean[, plot_genes])
)
df_line$gene <- factor(df_line$gene, levels = plot_genes)
# Distinguish SOCS3 (brake) from ISG
df_line$type <- ifelse(df_line$gene == "SOCS3", "SOCS3 (brake)", "ISG / IFN")
p1 <- ggplot(df_line, aes(x = bin_center, y = expr, color = gene)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_color_manual(values = c("SOCS3" = "#D32F2F", "ISG15" = "#1976D2",
                                "MX1" = "#1976D2", "OAS2" = "#1976D2",
                                "IFI44L" = "#1976D2", "IRF7" = "#1976D2")) +
  labs(x = "Pseudotime (Healthy root -> activated/exhausted)",
       y = "Mean expression (log-normalized)",
       title = "SOCS3 and ISGs rise together along pDC activation",
       subtitle = "Activation-suppression paradox along pseudotime") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")
ggsave(file.path(out_dir, "path20_socs3_isg_pseudotime.pdf"), p1,
       width = 7.5, height = 5)

# 5.2 Heatmap of all genes along pseudotime (binned z-score)
heat_mat <- bin_mean
# z-score per column
heat_z <- scale(heat_mat)
heat_z[is.nan(heat_z)] <- 0
df_heat <- expand.grid(gene = genes_of_interest, bin = seq_len(n_bins))
df_heat$z <- as.vector(t(heat_z))
df_heat$gene <- factor(df_heat$gene, levels = genes_of_interest)
p2 <- ggplot(df_heat, aes(x = bin, y = gene, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "z-score") +
  labs(x = "Pseudotime bin (1 = Healthy root, 10 = exhausted)",
       y = NULL,
       title = "Gene expression dynamics along pDC pseudotime") +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir, "path20_pseudotime_heatmap.pdf"), p2,
       width = 6, height = 5)

cat("Visualizations saved: path20_socs3_isg_pseudotime.pdf / path20_pseudotime_heatmap.pdf\n")

## =========================================================================
## Step 6: Statistics summary
## =========================================================================
summ <- c(
  "=== Pseudotime SOCS3/ISG dynamics ===",
  paste0("Pseudotime range: ", round(min(pseudo),2), " -> ", round(max(pseudo),2)),
  paste0("SOCS3 vs pseudotime Spearman rho = ", round(gene_pt_cor["SOCS3"], 3)),
  paste0("ISG15 vs pseudotime Spearman rho = ", round(gene_pt_cor["ISG15"], 3)),
  paste0("MX1 vs pseudotime Spearman rho = ", round(gene_pt_cor["MX1"], 3)),
  paste0("IRF7 vs pseudotime Spearman rho = ", round(gene_pt_cor["IRF7"], 3)),
  "",
  "Interpretation: if both SOCS3 and ISG show significant positive correlation with pseudotime (co-upregulated),",
  "      this dynamically supports the 'activation-suppression paradox' -- SOCS3 is progressively up-regulated with pDC activation,",
  "      acting as the accompanying negative-feedback brake of the IFN program, rather than an independent event preceding activation."
)
writeLines(summ, file.path(out_dir, "path20_socs3_isg_stats.txt"))
cat("\nStatistics summary saved: path20_socs3_isg_stats.txt\n")

cat("\nPseudotime SOCS3/ISG dynamics analysis complete\n")
cat("Outputs: path20_pseudotime_bin_expression.csv / path20_gene_pseudotime_cor.csv\n")
cat("      path20_socs3_isg_stats.txt / path20_socs3_isg_pseudotime.pdf / path20_pseudotime_heatmap.pdf\n")
