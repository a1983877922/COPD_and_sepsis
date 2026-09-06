# =========================================================================
# 27_socs3_isg_monocyte_violin.R
#
# Purpose: Draw violin plots of SOCS3/ISG across the five "Monocyte" groups.
#   Replaces the original script 22 "pDC violin plot" — because the sensitivity
#   analysis (script 25) showed pDC annotation is biased by SCP548 contamination,
#   whereas SOCS3 in "blood monocytes" is a robust shared up-regulated gene
#   (COPD vs Healthy log2FC=2.03, Sepsis vs Healthy=1.67),
#   completely independent of pDC annotation.
#
# Key point: Additionally save the "long-table data before plotting" (CSV) for
#   easy local re-formatting of the figure.
#
# Input: path1_sepsis_copd_integrated.rds
# Output: path27_ series (including path27_plot_data.csv for local re-plotting)
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
cat("SOCS3 / ISG five-group monocyte (Monocyte) violin plot\n")
cat("==============================================================\n")

genes_of_interest <- c("SOCS3", "SOCS1",
                       "IRF7", "TLR9",
                       "ISG15", "MX1", "OAS2", "IFI44L", "ISG20")

grp_order <- c("Healthy", "Infection_Control", "COPD",
               "Sepsis", "Sepsis_Pneumonia")

## =========================================================================
## Step 1: Load blood object → subset monocytes
## =========================================================================
cat("\n===== Step 1: load blood object + subset Monocyte =====\n")

blood_file <- file.path(out_dir, "path1_sepsis_copd_integrated.rds")
if (!file.exists(blood_file)) stop("Cannot find ", blood_file, " please run script 01 first")
blood_seu <- readRDS(blood_file)
cat("blood object:", ncol(blood_seu), "cells\n")

mono <- subset(blood_seu, subset = cell_type == "Monocyte")
rm(blood_seu); gc()
cat("monocyte cell count:", ncol(mono), "\n")

mono$group <- factor(as.character(mono$group), levels = grp_order)
print(table(mono$group))

## =========================================================================
## Step 2: Extract data-layer expression → long table
## =========================================================================
cat("\n===== Step 2: extract expression =====\n")

genes <- genes_of_interest[genes_of_interest %in% rownames(mono)]
expr <- GetAssayData(mono, assay = "RNA", layer = "data")
expr <- as.matrix(expr[genes, , drop = FALSE])

df_full <- do.call(rbind, lapply(genes, function(g) {
  data.frame(cell = colnames(mono),
             gene = g,
             expr = as.numeric(expr[g, ]),
             group = as.character(mono$group),
             dataset = as.character(mono$dataset),
             stringsAsFactors = FALSE)
}))
df_full$gene <- factor(df_full$gene, levels = genes)

## =========================================================================
## Step 3: Statistical tests (full data)
## =========================================================================
cat("\n===== Step 3: statistical tests (all monocytes) =====\n")

stat_lines <- c("=== SOCS3 / ISG five-group monocyte expression statistics ===", "")

kw_tab <- sapply(genes, function(g) {
  d <- df_full[df_full$gene == g, ]
  if (length(unique(d$group)) < 2) return(NA_real_)
  tryCatch(kruskal.test(expr ~ group, data = d)$p.value, error = function(e) NA_real_)
})
stat_lines <- c(stat_lines, "Kruskal-Wallis (overall difference across five groups):")
cat("Kruskal-Wallis (overall difference across five groups):\n")
for (g in genes) {
  line <- sprintf("  %-8s p = %.3e", g, kw_tab[g])
  cat(line, "\n"); stat_lines <- c(stat_lines, line)
}

stat_lines <- c(stat_lines, "", "SOCS3 each group vs Healthy (Wilcoxon):")
cat("\nSOCS3 each group vs Healthy (Wilcoxon):\n")
socs3_healthy <- df_full$expr[df_full$gene == "SOCS3" & df_full$group == "Healthy"]
for (g in setdiff(grp_order, "Healthy")) {
  socs3_g <- df_full$expr[df_full$gene == "SOCS3" & df_full$group == g]
  if (length(socs3_g) < 3) {
    line <- sprintf("  %-20s too few cells (<3), skip", g)
  } else {
    p <- tryCatch(wilcox.test(socs3_g, socs3_healthy)$p.value,
                  error = function(e) NA_real_)
    med <- median(socs3_g, na.rm = TRUE); med_h <- median(socs3_healthy, na.rm = TRUE)
    line <- sprintf("  %-20s median=%.3f vs %.3f, p=%.3e", g, med, med_h, p)
  }
  cat(line, "\n"); stat_lines <- c(stat_lines, line)
}

writeLines(stat_lines, file.path(out_dir, "path27_socs3_isg_mono_stats.txt"))

## =========================================================================
## Step 4: Save the long-table data before plotting (sampled, for local re-plotting)
## =========================================================================
cat("\n===== Step 4: save plotting data (sampled) =====\n")

# At most 3000 cells per group, to keep CSV size manageable while representing every group (esp. COPD/SP with few cells)
set.seed(123)
MAX_PER_GROUP <- 3000
keep_idx <- unlist(lapply(grp_order, function(g) {
  idx <- which(df_full$group == g)
  if (length(idx) <= MAX_PER_GROUP) return(idx)
  sample(idx, MAX_PER_GROUP)
}))
df_plot <- df_full[keep_idx, , drop = FALSE]
cat("plotting data after sampling:", nrow(df_plot), "rows (original", nrow(df_full), "rows)\n")
cat("sampled cell count per group:\n")
print(table(df_plot$group[df_plot$gene == "SOCS3"]))

write.csv(df_plot, file.path(out_dir, "path27_plot_data.csv"), row.names = FALSE)
cat("Saved plotting data: path27_plot_data.csv (can re-plot directly with ggplot locally)\n")

## =========================================================================
## Step 5: Server-side plotting (basic version, can be tuned locally)
## =========================================================================
cat("\n===== Step 5: server-side plotting =====\n")

grp_cols <- c("Healthy" = "#90A4AE", "Infection_Control" = "#4DB6AC",
              "COPD" = "#FFB300", "Sepsis" = "#E53935", "Sepsis_Pneumonia" = "#8E24AA")

df_plot$is_socs3 <- ifelse(df_plot$gene == "SOCS3", "SOCS3", "ISG / other")

p <- ggplot(df_plot, aes(x = group, y = expr, fill = group)) +
  geom_violin(scale = "width", alpha = 0.55, trim = TRUE, color = NA) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9,
               aes(color = is_socs3), linewidth = 0.4) +
  geom_jitter(width = 0.18, size = 0.2, alpha = 0.12) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = grp_cols, drop = FALSE) +
  scale_color_manual(values = c("SOCS3" = "#D32F2F", "ISG / other" = "grey30"),
                     guide = "none") +
  labs(x = NULL, y = "Expression (log-normalized)",
       title = "SOCS3 and ISGs across monocyte states",
       subtitle = "SOCS3 (red box) is a robust monocyte-level group feature, strongest in COPD/sepsis") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom", legend.title = element_blank())

ggsave(file.path(out_dir, "path27_socs3_isg_mono_violin.pdf"), p, width = 9, height = 8)
ggsave(file.path(out_dir, "path27_socs3_isg_mono_violin.png"), p, width = 9, height = 8, dpi = 300)

cat("\nSOCS3/ISG monocyte violin plot complete\n")
cat("Output: path27_socs3_isg_mono_stats.txt / path27_plot_data.csv\n")
cat("      path27_socs3_isg_mono_violin.pdf/.png\n")
cat("Note: path27_plot_data.csv is the long table before plotting (columns: cell/gene/expr/group/dataset),\n")
cat("      after copying locally it can be re-plotted with ggplot with custom colors, fonts, facets, significance marks, etc.\n")
