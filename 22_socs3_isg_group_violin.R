# =========================================================================
# Script 22: violin plots of SOCS3 and ISGs across the five groups
#
# Purpose: draw violin/box plots of SOCS3 and representative ISGs in pDCs
#       across the five groups.
#       This is a more convincing mechanistic figure than the weak pseudotime
#       correlation (script 20, rho~0.1, negative): it directly shows that
#       SOCS3 is a between-group state feature of sepsis pDCs
#       (matching the DEG: SOCS3 log2FC=+3.29, p=1.3e-183, the strongest signal
#       overall), rather than a continuous gradient along the resting-to-
#       activated trajectory.
#
# Input: path11_pdc_subset.rds (saved by script 11, contains counts+data+group)
# Output: path22_* series
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

suppressPackageStartupMessages(library(ggplot2))

cat("\n==============================================================\n")
cat("SOCS3 / ISG violin plots across the five pDC groups\n")
cat("==============================================================\n")

# Genes of interest: SOCS3/SOCS1 (brakes) + IRF7/TLR9 (regulation/sensing) + ISGs (effectors)
genes_of_interest <- c("SOCS3", "SOCS1",
                       "IRF7", "TLR9",
                       "ISG15", "MX1", "OAS2", "IFI44L", "ISG20")

## =========================================================================
## Step 1: read pDC object + assemble long table
## =========================================================================
cat("\n===== Step 1: read pDC object =====\n")

pdc_file <- file.path(out_dir, "path11_pdc_subset.rds")
if (!file.exists(pdc_file)) stop("Cannot find ", pdc_file, " please run script 11 first")
pdc <- readRDS(pdc_file)
cat("pDC object:", ncol(pdc), "cells\n")
print(table(pdc$group))

# Order the five groups by disease severity
grp_order <- c("Healthy", "Infection_Control", "COPD",
               "Sepsis", "Sepsis_Pneumonia")
pdc$group <- factor(as.character(pdc$group), levels = grp_order)

genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(pdc)]
cat("Genes of interest available:", length(genes_of_interest), "\n")

# Extract the data layer (log-normalized)
expr <- GetAssayData(pdc, assay = "RNA", layer = "data")
expr <- as.matrix(expr[genes_of_interest, , drop = FALSE])

df_list <- lapply(genes_of_interest, function(g) {
  data.frame(gene = g,
             group = pdc$group,
             expr = as.numeric(expr[g, ]),
             stringsAsFactors = FALSE)
})
df <- do.call(rbind, df_list)
df$gene <- factor(df$gene, levels = genes_of_interest)

# SOCS3 highlight flag (red box outline)
df$is_socs3 <- ifelse(df$gene == "SOCS3", "SOCS3", "ISG / other")

## =========================================================================
## Step 2: statistical tests
## =========================================================================
cat("\n===== Step 2: statistical tests =====\n")

stat_lines <- c("=== SOCS3 / ISG expression statistics across the five pDC groups ===", "")

# Per-gene Kruskal-Wallis (overall difference across the five groups)
kw_tab <- sapply(genes_of_interest, function(g) {
  d <- df[df$gene == g, ]
  if (length(unique(d$group)) < 2) return(NA_real_)
  tryCatch(kruskal.test(expr ~ group, data = d)$p.value,
           error = function(e) NA_real_)
})
stat_lines <- c(stat_lines, "Kruskal-Wallis (overall difference across the five groups):")
cat("Kruskal-Wallis (overall difference across the five groups):\n")
for (g in genes_of_interest) {
  line <- sprintf("  %-8s p = %.3e", g, kw_tab[g])
  cat(line, "\n"); stat_lines <- c(stat_lines, line)
}

# SOCS3, each group vs Healthy (Wilcoxon) (core test)
stat_lines <- c(stat_lines, "", "SOCS3, each group vs Healthy (Wilcoxon):")
cat("\nSOCS3, each group vs Healthy (Wilcoxon):\n")
socs3_wt <- df$expr[df$gene == "SOCS3" & df$group == "Healthy"]
for (g in setdiff(levels(pdc$group), "Healthy")) {
  socs3_g <- df$expr[df$gene == "SOCS3" & df$group == g]
  if (length(socs3_g) < 3) {
    line <- sprintf("  %-20s too few cells (<3), skipped", g)
  } else {
    p <- tryCatch(wilcox.test(socs3_g, socs3_wt)$p.value,
                  error = function(e) NA_real_)
    line <- sprintf("  %-20s p = %.3e", g, p)
  }
  cat(line, "\n"); stat_lines <- c(stat_lines, line)
}

writeLines(stat_lines, file.path(out_dir, "path22_socs3_isg_stats.txt"))

## =========================================================================
## Step 3: visualisation
## =========================================================================
cat("\n===== Step 3: visualisation =====\n")

# Colour scheme for the five groups
grp_cols <- c("Healthy" = "#90A4AE",          # grey-blue
              "Infection_Control" = "#4DB6AC", # teal
              "COPD" = "#FFB300",              # amber
              "Sepsis" = "#E53935",            # red
              "Sepsis_Pneumonia" = "#8E24AA")  # purple

p <- ggplot(df, aes(x = group, y = expr, fill = group)) +
  geom_violin(scale = "width", alpha = 0.55, trim = TRUE, color = NA) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9,
               aes(color = is_socs3), linewidth = 0.4) +
  geom_jitter(width = 0.18, size = 0.25, alpha = 0.15) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = grp_cols, drop = FALSE) +
  scale_color_manual(values = c("SOCS3" = "#D32F2F", "ISG / other" = "grey30"),
                     guide = "none") +
  labs(x = NULL, y = "Expression (log-normalized)",
       title = "SOCS3 and interferon-stimulated genes across pDC states",
       subtitle = "SOCS3 (red box) is a group-level state feature, strongest in sepsis") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom",
        legend.title = element_blank())

ggsave(file.path(out_dir, "path22_socs3_isg_violin.pdf"), p,
       width = 9, height = 8)
ggsave(file.path(out_dir, "path22_socs3_isg_violin.png"), p,
       width = 9, height = 8, dpi = 300)

cat("Saved: path22_socs3_isg_violin.pdf / .png\n")
cat("Statistics saved: path22_socs3_isg_stats.txt\n")
cat("\nSOCS3 / ISG violin plots across the five groups done\n")
