#!/usr/bin/env Rscript
# =============================================================================
# Script 32: imprint gradient faceted by dataset
# Sepsis + COPD comorbidity analysis — visualising the imprint gradient with
# dataset colouring / faceting
# =============================================================================
# Purpose: make explicit the leave-one-dataset-out finding that the monotonic
#       imprint gradient partly depends on the large datasets:
#       (1) The high Sepsis imprint is contributed mainly by GSE279452
#       (2) The low Healthy imprint is contributed mainly by SCP548
#       Colouring + faceting by dataset shows, at a glance, the imprint
#       distribution of each dataset within a group.
#
# Imprint definition (consistent with scripts 28/31):
#   cell_type == "Monocyte", group in (Healthy, COPD, Sepsis)
#   imprint = donor-level mean of the 149-gene z-scores (per gene across cells)
#
# Output:
#   path32_imprint_donor_by_dataset.csv   long table (donor/imprint/group/dataset, for local plot editing)
#   path32_imprint_group_dataset.pdf      x=group, coloured by dataset
#   path32_imprint_facet.pdf              facet by group, x=dataset
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(ggplot2)
})

## ---- Paths ----
this_file <- commandArgs(trailingOnly = FALSE)
.f <- grep("--file=", this_file, value = TRUE)
if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1]) else .this_file <- "."
out_dir <- dirname(normalizePath(.this_file))
source(file.path(out_dir, "00_config.R"))

## ---- Donor field detection (patient preferred, sample as fallback) [same as 28b/31] ----
detect_donor <- function(seu) {
  meta <- seu@meta.data
  get <- function(f) if (f %in% colnames(meta)) as.character(meta[[f]]) else NULL
  patient <- get("patient"); sample <- get("sample")
  if (!is.null(patient)) {
    if (!is.null(sample)) {
      na_idx <- is.na(patient) | patient == ""
      patient[na_idx] <- sample[na_idx]
      return(patient)
    }
    return(patient)
  }
  for (cand in c("donor_id","donor","subject","subject_id","Subject_Identity","sample")) {
    v <- get(cand); if (!is.null(v)) return(v)
  }
  NULL
}

## ============================================================================
## Read blood object + original 149 genes
## ============================================================================
cat("\n########## Reading blood object + original 149 genes ##########\n")
blood_file <- file.path(out_dir, "path1_sepsis_copd_integrated.rds")
if (!file.exists(blood_file)) stop("Cannot find path1_sepsis_copd_integrated.rds")
blood <- readRDS(blood_file)
donor_vec <- detect_donor(blood)
cat("Donor field detection done\n")

ds_field <- if ("dataset" %in% colnames(blood@meta.data)) "dataset" else "orig.ident"
cat("Dataset field:", ds_field, "\n")

orig_file <- file.path(out_dir, "path2_shared_myeloid_genes.txt")
if (!file.exists(orig_file)) orig_file <- file.path(out_dir, "path21_shared_up_genes.txt")
orig_149 <- if (file.exists(orig_file)) readLines(orig_file) else character(0)
orig_149 <- orig_149[orig_149 != ""]
cat("Number of original 149 genes:", length(orig_149), "\n")

## ============================================================================
## Compute donor-level imprint (mean of 149 z-scores)
## ============================================================================
cat("\n########## Computing donor-level imprint ##########\n")
keep <- !is.na(blood$cell_type) & blood$cell_type == "Monocyte" &
        !is.na(blood$group) & blood$group %in% c("Healthy", "COPD", "Sepsis")

data_mat <- GetAssayData(blood, assay = "RNA", layer = "data")
genes <- orig_149[orig_149 %in% rownames(data_mat)]
data_mat <- data_mat[genes, , drop = FALSE]
data_z <- t(scale(t(as.matrix(data_mat))))   # per-gene z-score across cells
cells <- rownames(blood@meta.data)[keep]
data_z <- data_z[, cells, drop = FALSE]

d <- donor_vec[keep]
donors <- unique(d); donors <- donors[!is.na(donors) & donors != ""]
imprint <- vapply(donors, function(x)
  mean(colMeans(data_z[, d == x, drop = FALSE], na.rm = TRUE)), numeric(1))

# Assemble long table: donor, imprint, group, dataset
grp <- vapply(donors, function(x) {
  g <- unique(blood$group[donor_vec == x & keep]); if (length(g) > 0) g[1] else NA_character_
}, character(1))
ds <- vapply(donors, function(x) {
  v <- unique(blood@meta.data[[ds_field]][donor_vec == x & keep])
  if (length(v) > 0) v[1] else NA_character_
}, character(1))

plot_df <- data.frame(
  donor   = donors,
  imprint = as.numeric(imprint),
  group   = grp,
  dataset = ds,
  stringsAsFactors = FALSE)
plot_df <- plot_df[!is.na(plot_df$group) & !is.na(plot_df$dataset), ]
plot_df$group <- factor(plot_df$group, levels = c("Healthy", "COPD", "Sepsis"))

cat("Number of donor-level imprint samples:", nrow(plot_df), "\n")
print(table(plot_df$group, plot_df$dataset))

write.csv(plot_df, file.path(out_dir, "path32_imprint_donor_by_dataset.csv"), row.names = FALSE)

## ============================================================================
## Plotting
## ============================================================================
# Colour palette for 6 datasets (first 6 of RColorBrewer Dark2 + fallback)
ds_levels <- sort(unique(plot_df$dataset))
ds_cols <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02")
names(ds_cols) <- NULL
if (length(ds_levels) > 6) ds_cols <- colorRampPalette(ds_cols)(length(ds_levels))

# Figure 1: x=group, coloured by dataset
p1 <- ggplot(plot_df, aes(x = group, y = imprint, color = dataset)) +
  geom_boxplot(outlier.shape = NA, color = "grey40", width = 0.6, alpha = 0.3) +
  geom_jitter(width = 0.18, size = 1.6, alpha = 0.8) +
  scale_color_manual(values = ds_cols) +
  labs(x = NULL, y = "Inflammatory imprint (149-gene z-score)",
       title = "Inflammatory imprint by group, coloured by dataset") +
  theme_classic(base_size = 11) +
  theme(legend.position = "right",
        plot.title = element_text(size = 11, face = "bold"))
ggsave(file.path(out_dir, "path32_imprint_group_dataset.pdf"), p1, width = 7, height = 4.5)

# Figure 2: facet by group, x=dataset (most direct view of imprint differences between datasets within a group)
p2 <- ggplot(plot_df, aes(x = dataset, y = imprint, color = dataset)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  facet_wrap(~ group, scales = "free_x") +
  scale_color_manual(values = ds_cols) +
  labs(x = NULL, y = "Inflammatory imprint (149-gene z-score)",
       title = "Inflammatory imprint by dataset within each group") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(size = 11, face = "bold"))
ggsave(file.path(out_dir, "path32_imprint_facet.pdf"), p2, width = 8, height = 4)

cat("\nDone. Output:\n")
cat("  path32_imprint_donor_by_dataset.csv   (long table, for local plot editing)\n")
cat("  path32_imprint_group_dataset.pdf      (x=group, coloured by dataset)\n")
cat("  path32_imprint_facet.pdf              (facet by group, x=dataset)\n")
