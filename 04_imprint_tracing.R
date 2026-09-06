###############################################################################
# 04_imprint_tracing.R — lung-to-blood immune imprint scoring
#
# Purpose: score cells with 97 shared imprint genes (lung COPD myeloid ∩ blood Sepsis monocyte)
#   to test the "imprint from lung to blood" gradient hypothesis:
#         Lung Control < Lung COPD   (lung-local chronic inflammation establishes the imprint)
#         Blood Healthy < Blood COPD  (imprint enters blood as monocytes migrate)
#         Blood Sepsis high           (imprint amplified in sepsis)
#
# Data:
#   Lung side: path2_copd_lung.rds (GSE136831, with built-in cell_category + disease)
#   Blood side: path1_sepsis_copd_integrated.rds (integrated object, cell_type + group)
#
# Output:
#   path4_imprint_score_summary.csv  (mean imprint score per group)
#   path4_imprint_boxplot.pdf        (boxplot of imprint score across five groups)
#   path4_imprint_stats.txt          (inter-group wilcox test)
###############################################################################

# Auto-locate this script's directory and find 00_config.R in the same directory (no need to change this path locally or on the server)
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
cat("Imprint tracing analysis (lung-to-blood immune imprint scoring)\n")
cat("==============================================================\n")

## =========================================================================
## Step 1: read imprint genes
## =========================================================================

cat("\n===== Step 1: read imprint genes =====\n")

imprint_file <- file.path(out_dir, "path2_shared_myeloid_genes.txt")
if (!file.exists(imprint_file)) {
  stop("Cannot find imprint gene file: ", imprint_file,
       "\nPlease run 02_build_lung_cross_tissue.R (formerly 02_path2_cross_tissue.R) to generate this file")
}
imprint_genes <- readLines(imprint_file)
imprint_genes <- imprint_genes[nzchar(imprint_genes)]
cat("Number of imprint genes:", length(imprint_genes), "\n")

## =========================================================================
## Step 2: blood-side monocyte imprint scoring (load in steps to avoid holding two large objects at once)
## =========================================================================

cat("\n===== Step 2: blood-side monocyte imprint scoring =====\n")

blood_seu <- readRDS(file.path(out_dir, "path1_sepsis_copd_integrated.rds"))
cat("Blood object:", ncol(blood_seu), "cells\n")

blood_mono <- subset(blood_seu, subset = cell_type == "Monocyte" &
                       group %in% c("Healthy", "COPD", "Sepsis"))
rm(blood_seu); gc()
cat("Blood monocytes:", ncol(blood_mono), "cells\n")
print(table(blood_mono$group))

# Compute imprint score (only take imprint genes present in the object)
genes_blood <- intersect(imprint_genes, rownames(blood_mono))
cat("Imprint genes available on blood side:", length(genes_blood), "/", length(imprint_genes), "\n")
if (length(genes_blood) < 10) {
  stop("Insufficient imprint genes on blood side, check whether gene naming is consistent")
}
blood_mono <- AddModuleScore(blood_mono, features = list(genes_blood),
                             name = "imprint")
blood_mono$imprint_score <- blood_mono$imprint1

blood_df <- data.frame(
  tissue = "Blood",
  group  = blood_mono$group,
  imprint = blood_mono$imprint_score,
  stringsAsFactors = FALSE
)
blood_df$label <- paste("Blood", blood_df$group)
rm(blood_mono); gc()

## =========================================================================
## Step 3: lung-side myeloid imprint scoring
## =========================================================================

cat("\n===== Step 3: lung-side myeloid imprint scoring =====\n")

lung_seu <- readRDS(file.path(out_dir, "path2_copd_lung.rds"))
cat("Lung object:", ncol(lung_seu), "cells\n")

lung_myeloid <- subset(lung_seu, subset = cell_category == "Myeloid" &
                         disease %in% c("Control", "COPD"))
rm(lung_seu); gc()
cat("Lung myeloid:", ncol(lung_myeloid), "cells\n")
print(table(lung_myeloid$disease))

genes_lung <- intersect(imprint_genes, rownames(lung_myeloid))
cat("Imprint genes available on lung side:", length(genes_lung), "/", length(imprint_genes), "\n")
if (length(genes_lung) < 10) {
  stop("Insufficient imprint genes on lung side, check whether gene naming is consistent")
}
lung_myeloid <- AddModuleScore(lung_myeloid, features = list(genes_lung),
                               name = "imprint")
lung_myeloid$imprint_score <- lung_myeloid$imprint1

lung_df <- data.frame(
  tissue = "Lung",
  group  = lung_myeloid$disease,
  imprint = lung_myeloid$imprint_score,
  stringsAsFactors = FALSE
)
lung_df$label <- paste("Lung", lung_df$group)
rm(lung_myeloid); gc()

## =========================================================================
## Step 4: merge + statistical test + visualization
## =========================================================================

cat("\n===== Step 4: merge + statistics =====\n")

plot_df <- rbind(blood_df, lung_df)
plot_df$label <- factor(plot_df$label,
                        levels = c("Lung Control", "Lung COPD",
                                   "Blood Healthy", "Blood COPD", "Blood Sepsis"))

# 4.1 Mean imprint score summary per group
summary_df <- aggregate(imprint ~ label, data = plot_df,
                        FUN = function(x) c(mean = mean(x), median = median(x),
                                            sd = sd(x), n = length(x)))
summary_df <- do.call(data.frame, summary_df)
colnames(summary_df) <- c("label", "mean", "median", "sd", "n")
summary_df <- summary_df[order(match(summary_df$label, levels(plot_df$label))), ]
write.csv(summary_df, file.path(out_dir, "path4_imprint_score_summary.csv"),
          row.names = FALSE)
cat("\nMean imprint score per group:\n")
print(summary_df)

# 4.2 Inter-group wilcox test (key comparisons)
stats_lines <- character(0)
run_test <- function(a, b, name) {
  p <- wilcox.test(a, b)$p.value
  line <- sprintf("%s: p = %.3e (mean %.3f vs %.3f)",
                  name, p, mean(a), mean(b))
  cat(line, "\n")
  stats_lines <<- c(stats_lines, line)
}
cat("\nInter-group wilcox test:\n")
run_test(blood_df$imprint[blood_df$group == "COPD"],
         blood_df$imprint[blood_df$group == "Healthy"],
         "Blood COPD vs Blood Healthy")
run_test(blood_df$imprint[blood_df$group == "Sepsis"],
         blood_df$imprint[blood_df$group == "Healthy"],
         "Blood Sepsis vs Blood Healthy")
run_test(lung_df$imprint[lung_df$group == "COPD"],
         lung_df$imprint[lung_df$group == "Control"],
         "Lung COPD vs Lung Control")
run_test(blood_df$imprint[blood_df$group == "COPD"],
         lung_df$imprint[lung_df$group == "COPD"],
         "Blood COPD vs Lung COPD")
writeLines(stats_lines, file.path(out_dir, "path4_imprint_stats.txt"))

# 4.3 Boxplot
p <- ggplot(plot_df, aes(x = label, y = imprint, fill = tissue)) +
  geom_boxplot(outlier.size = 0.2, alpha = 0.85) +
  scale_fill_manual(values = c("Lung" = "#F0997B", "Blood" = "#85B7EB")) +
  labs(x = NULL, y = "Imprint score",
       title = "Lung-to-blood imprint score gradient") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")
ggsave(file.path(out_dir, "path4_imprint_boxplot.pdf"),
       plot = p, width = 8, height = 6)

# 4.4 Violin plot (supplementary, showing distribution)
p2 <- ggplot(plot_df, aes(x = label, y = imprint, fill = tissue)) +
  geom_violin(alpha = 0.85, trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.size = 0.2, fill = "white") +
  scale_fill_manual(values = c("Lung" = "#F0997B", "Blood" = "#85B7EB")) +
  labs(x = NULL, y = "Imprint score",
       title = "Lung-to-blood imprint score distribution") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")
ggsave(file.path(out_dir, "path4_imprint_violin.pdf"),
       plot = p2, width = 8, height = 6)

cat("\nImprint tracing analysis complete\n")
cat("Output:\n")
cat("  - path4_imprint_score_summary.csv (mean imprint score per group)\n")
cat("  - path4_imprint_boxplot.pdf (boxplot)\n")
cat("  - path4_imprint_violin.pdf (violin plot)\n")
cat("  - path4_imprint_stats.txt (inter-group wilcox test)\n")
