# =========================================================================
# 11b_pdc_trajectory_analysis.R
# Extract "quantitative" conclusions of the pDC pseudotime trajectory from path11_pdc_subset.rds saved by 11:
#   1) pseudotime distribution per group + Kruskal-Wallis test (whether Sepsis is significantly skewed to the trajectory end)
#   2) graph_test: genes significantly changing along the pseudotime (find trajectory "driver genes", expected SOCS3/TLR9/ISG)
#   3) expression trend plots of key genes along pseudotime (SOCS3/TLR9/ISG20/MX1/IRF7/IRF4/CIITA...)
#   4) branch point detection (whether the trajectory bifurcates)
#
# Dependencies: monocle3 + SeuratWrappers + Seurat + ggplot2 (auto-install monocle3/SeuratWrappers)
# Usage: run after 11 finishes (requires path11_pdc_subset.rds to have been generated)
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

cat("\n==============================================================\n")
cat("pDC trajectory quantitative analysis (pseudotime distribution + graph_test + branches)\n")
cat("==============================================================\n")

## =========================================================================
## Step 0: dependencies (auto-install, do not skip)
## =========================================================================
if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
  cat("SeuratWrappers not detected, installing from CRAN...\n")
  install.packages("SeuratWrappers", repos = "https://cloud.r-project.org",
                   Ncpus = max(1, parallel::detectCores() %/% 2))
}
if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
  stop("SeuratWrappers installation failed")
}
if (!requireNamespace("monocle3", quietly = TRUE)) {
  cat("monocle3 not detected, installing from Bioconductor...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install("monocle3", update = FALSE, ask = FALSE)
}
if (!requireNamespace("monocle3", quietly = TRUE)) {
  stop("monocle3 installation failed")
}
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(SeuratWrappers))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))

## =========================================================================
## Step 1: read pDC object
## =========================================================================
cat("\n===== Step 1: read pDC object =====\n")
pdc_file <- file.path(out_dir, "path11_pdc_subset.rds")
if (!file.exists(pdc_file)) stop("Cannot find path11_pdc_subset.rds, please run 11 first")
pdc <- readRDS(pdc_file)
cat("pDC cell count:", ncol(pdc), "\n")
print(table(pdc$group))

## =========================================================================
## Step 2: convert to CDS + preprocessing
## =========================================================================
cat("\n===== Step 2: convert to CDS + preprocessing =====\n")
cds <- as.cell_data_set(pdc)
cds <- tryCatch(estimate_size_factors(cds), error = function(e) cds)
rm(pdc); gc()

## =========================================================================
## Step 3: reconstruct trajectory (consistent with 11, root = Healthy)
## =========================================================================
cat("\n===== Step 3: reconstruct trajectory =====\n")
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)
root_cells <- colnames(cds)[colData(cds)$group == "Healthy"]
if (length(root_cells) == 0) root_cells <- colnames(cds)[1]
cds <- order_cells(cds, root_cells = root_cells, reduction_method = "UMAP")
cat("Trajectory reconstruction complete\n")

pseudo <- pseudotime(cds)
cat("Pseudotime range:", round(min(pseudo, na.rm = TRUE), 3), "->",
    round(max(pseudo, na.rm = TRUE), 3), "\n")

## =========================================================================
## Step 4: pseudotime distribution per group + statistical test
## =========================================================================
cat("\n===== Step 4: pseudotime distribution per group =====\n")

grp <- as.character(colData(cds)$group)
df_pt <- data.frame(pseudotime = as.numeric(pseudo), group = grp,
                    stringsAsFactors = FALSE)
df_pt <- df_pt[!is.na(df_pt$pseudotime), ]

# Summary table
summ <- do.call(rbind, lapply(sort(unique(df_pt$group)), function(g) {
  x <- df_pt$pseudotime[df_pt$group == g]
  data.frame(group = g, n = length(x),
             median = median(x), mean = mean(x), sd = sd(x),
             q25 = quantile(x, 0.25), q75 = quantile(x, 0.75))
}))
print(summ, row.names = FALSE)
write.csv(summ, file.path(out_dir, "path11b_pdc_pseudotime_by_group.csv"),
          row.names = FALSE)

# Kruskal-Wallis test (non-parametric, pseudotime is non-normal)
sink(file.path(out_dir, "path11b_pdc_pseudotime_stats.txt"))
cat("pDC pseudotime statistics per group (root = Healthy)\n")
cat("=============================================\n\n")
kw <- kruskal.test(pseudotime ~ group, data = df_pt)
cat("Kruskal-Wallis test (all-group comparison):\n")
cat("  chi-squared =", kw$statistic, " df =", kw$parameter,
    " p =", kw$p.value, "\n\n")

# Pairwise wilcox (Sepsis vs Healthy etc.)
pairs <- list(c("Sepsis", "Healthy"), c("Sepsis", "Infection_Control"),
              c("COPD", "Healthy"), c("Infection_Control", "Healthy"))
cat("Pairwise Wilcoxon test:\n")
for (pr in pairs) {
  if (all(pr %in% df_pt$group)) {
    a <- df_pt$pseudotime[df_pt$group == pr[1]]
    b <- df_pt$pseudotime[df_pt$group == pr[2]]
    if (length(a) >= 3 && length(b) >= 3) {
      wt <- wilcox.test(a, b)
      cat(sprintf("  %s vs %s: p = %.3e\n", pr[1], pr[2], wt$p.value))
    }
  }
}
sink()
cat("Pseudotime statistics saved: path11b_pdc_pseudotime_by_group.csv / _stats.txt\n")

# Boxplot
p_box <- ggplot(df_pt, aes(x = group, y = pseudotime, fill = group)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
  theme_minimal(base_size = 12) +
  labs(x = NULL, y = "Pseudotime", title = "pDC pseudotime by group") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out_dir, "path11b_pdc_pseudotime_boxplot.pdf"), p_box,
       width = 7, height = 5)

## =========================================================================
## Step 5: graph_test -- genes significantly changing along pseudotime
## =========================================================================
cat("\n===== Step 5: graph_test (genes changing along pseudotime) =====\n")

pr_test <- tryCatch(
  graph_test(cds, neighbor_graph = "principal_graph",
             reduction_method = "UMAP", cores = 4, verbose = FALSE),
  error = function(e) NULL)

if (is.null(pr_test) || nrow(pr_test) == 0) {
  cat("graph_test failed or returned no results (possible version compatibility issue)\n")
} else {
  # Sort by q_value
  ord <- order(pr_test$q_value, na.last = TRUE)
  pr_test <- pr_test[ord, ]
  write.csv(pr_test, file.path(out_dir, "path11b_pdc_pseudotime_genes.csv"),
            row.names = FALSE)
  cat("Number of significant genes (q<0.05):", sum(pr_test$q_value < 0.05, na.rm = TRUE), "\n")
  cat("Top 20 genes changing along pseudotime:\n")
  top <- head(pr_test, 20)
  gene_names <- if ("gene_short_name" %in% colnames(top)) top$gene_short_name else rownames(top)
  print(data.frame(gene = gene_names,
                   morans_I = round(top$morans_I, 3),
                   q_value = format(top$q_value, scientific = TRUE, digits = 3)))
  cat("Saved: path11b_pdc_pseudotime_genes.csv\n")
}

## =========================================================================
## Step 6: expression trend plots of key genes along pseudotime
## =========================================================================
cat("\n===== Step 6: key genes along pseudotime =====\n")

key_genes <- c("SOCS3", "SOCS1", "TLR9", "MYD88", "ISG20", "MX1", "ISG15",
               "IFI44L", "IRF7", "IRF4", "CIITA", "HLA-DQA1", "CD74")
key_genes <- key_genes[key_genes %in% rownames(cds)]
cat("Key genes available for plotting (", length(key_genes), "):\n", paste(key_genes, collapse = ", "), "\n")

if (length(key_genes) >= 2) {
  p_genes <- tryCatch(
    plot_genes_in_pseudotime(cds[key_genes, ], color_cells_by = "group",
                             cell_size = 0.4),
    error = function(e) NULL)
  if (!is.null(p_genes)) {
    ggsave(file.path(out_dir, "path11b_pdc_key_genes_pseudotime.pdf"), p_genes,
           width = 11, height = ceiling(length(key_genes) / 2) * 2.8)
    cat("Key gene trend plot saved: path11b_pdc_key_genes_pseudotime.pdf\n")
  } else {
    cat("Key gene trend plot failed to render (skipped)\n")
  }
}

## =========================================================================
## Step 7: branch point detection
## =========================================================================
cat("\n===== Step 7: branch point detection =====\n")

branch_info <- tryCatch({
  if (requireNamespace("igraph", quietly = TRUE)) {
    g <- principal_graph(cds)[["UMAP"]]
    deg <- igraph::degree(g)
    branch_nodes <- names(deg)[deg >= 3]
    leaf_nodes   <- names(deg)[deg == 1]
    cat("Number of principal-graph nodes:", length(deg), "\n")
    cat("Number of branch points (degree>=3):", length(branch_nodes), "\n")
    cat("Number of leaf nodes (degree=1):", length(leaf_nodes), "\n")
    if (length(branch_nodes) > 0) {
      cat("Branch point nodes:", paste(branch_nodes, collapse = ", "), "\n")
    }
    data.frame(n_nodes = length(deg), n_branch = length(branch_nodes),
               n_leaf = length(leaf_nodes))
  } else {
    cat("igraph not installed, skipping branch point detection\n")
    NULL
  }
}, error = function(e) {
  cat("Branch point detection failed:", conditionMessage(e), "\n")
  NULL
})
if (!is.null(branch_info)) {
  write.csv(branch_info, file.path(out_dir, "path11b_pdc_branch_info.csv"),
            row.names = FALSE)
}

## =========================================================================
## Step 8: save CDS
## =========================================================================
cat("\n===== Step 8: save =====\n")
saveRDS(cds, file.path(out_dir, "path11b_pdc_cds.rds"))
cat("Saved CDS: path11b_pdc_cds.rds\n")

cat("\npDC trajectory quantitative analysis complete\n")
cat("Output: path11b_pdc_pseudotime_by_group.csv / _stats.txt / _boxplot.pdf\n")
cat("      path11b_pdc_pseudotime_genes.csv (graph_test)\n")
cat("      path11b_pdc_key_genes_pseudotime.pdf / _branch_info.csv\n")
