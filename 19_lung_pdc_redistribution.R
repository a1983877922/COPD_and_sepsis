# =========================================================================
# 19_lung_pdc_redistribution.R (pulmonary vs non-pulmonary sepsis pDC redistribution)
#
# Purpose: answer the reviewer's core question -- is the pDC blood redistribution
# significant ONLY in pulmonary sepsis? Compare the proportions and Ro/e of DC
# subtypes (pDC/cDC1/cDC2) across 5 groups, focusing on:
#     Sepsis_Pneumonia (pulmonary sepsis) vs Sepsis (non-pulmonary/other sepsis)
#     vs Infection_Control (non-sepsis infection) vs Healthy
#
# Scientific question: along the COPD (chronic lung disease) -> pulmonary
# infection -> sepsis axis, is the pDC sequestration/depletion specifically
# driven by 'pulmonary' infection, or shared by all sepsis?
#
# Input: path1_sepsis_copd_integrated.rds (blood object, includes cell_type + group)
# Output: path19_ series
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
cat("Pulmonary vs non-pulmonary sepsis: pDC redistribution comparison\n")
cat("==============================================================\n")

dc_markers <- list(
  "cDC1" = c("CLEC9A", "XCR1", "BATF3", "CADM1"),
  "cDC2" = c("CD1C", "CLEC10A", "FCER1A", "ITGAX"),
  "pDC"  = c("LILRA4", "CLEC4C", "IL3RA", "TCF4", "IRF7")
)

## =========================================================================
## Step 1: Read blood object -> extract DC -> subclustering to obtain dc_subtype
## =========================================================================
cat("\n===== Step 1: Blood-side DC subclustering =====\n")

blood_file <- file.path(out_dir, "path1_sepsis_copd_integrated.rds")
if (!file.exists(blood_file)) stop(blood_file, " not found, please run script 01 first")
blood_seu <- readRDS(blood_file)
cat("Blood object:", ncol(blood_seu), "cells\n")

blood_dc <- subset(blood_seu, subset = cell_type == "DC")
rm(blood_seu); gc()
cat("Blood-side DC count:", ncol(blood_dc), "\n")

# Fix desynchronized @cells mapping after subset (Seurat v5 bug)
dc_counts <- GetAssayData(blood_dc, assay = "RNA", layer = "counts")
blood_dc[["RNA"]] <- CreateAssay5Object(counts = dc_counts)
rm(dc_counts); gc()

old_res <- resolution
resolution <- 0.6
blood_dc <- integrate_and_annotate(blood_dc, marker_list = dc_markers,
                                   col_name = "dc_subtype")
resolution <- old_res

cat("\nBlood-side DC subtype x five-group cell counts:\n")
cnt <- table(blood_dc$dc_subtype, blood_dc$group)
print(cnt)

## =========================================================================
## Step 2: Proportions + Ro/e (five groups)
## =========================================================================
cat("\n===== Step 2: Proportions and Ro/e =====\n")

subtypes <- c("pDC", "cDC1", "cDC2")
groups   <- c("Healthy", "Infection_Control", "COPD", "Sepsis", "Sepsis_Pneumonia")

cnt <- cnt[subtypes, groups, drop = FALSE]
cnt[is.na(cnt)] <- 0

# Within-group proportion (normalize each column)
prop <- sweep(cnt, 2, colSums(cnt), "/")
# Global proportion (each row)
global_prop <- rowSums(cnt) / sum(cnt)
# Ro/e = group proportion / global proportion
roe <- sweep(prop, 1, global_prop, "/")

cat("\nDC subtype proportions in each group:\n")
print(round(prop, 4))

cat("\nDC subtype Ro/e in each group:\n")
print(round(roe, 3))

write.csv(cnt,  file.path(out_dir, "path19_dc_subtype_count.csv"))
write.csv(prop, file.path(out_dir, "path19_dc_subtype_proportion.csv"))
write.csv(roe,  file.path(out_dir, "path19_dc_subtype_roe.csv"))

## =========================================================================
## Step 3: Statistical tests -- pDC proportion differences between key groups
##         (Fisher's exact test)
## =========================================================================
cat("\n===== Step 3: pDC proportion between-group statistics =====\n")

# 2x2 table: pDC vs other DC x group A vs group B
fisher_pdc <- function(g1, g2) {
  a <- cnt["pDC", g1]; b <- cnt["pDC", g2]
  c1 <- sum(cnt[, g1]) - a; d <- sum(cnt[, g2]) - b
  m <- matrix(c(a, b, c1, d), nrow = 2)
  f <- fisher.test(m)
  cat(sprintf("  %s vs %s: pDC proportion %.4f vs %.4f, OR=%.2f, p=%.3e\n",
              g1, g2, a/sum(cnt[,g1]), b/sum(cnt[,g2]),
              f$estimate, f$p.value))
  f$p.value
}

cat("pDC proportion Fisher's exact test:\n")
pvals <- list(
  SP_vs_Sepsis      = fisher_pdc("Sepsis_Pneumonia", "Sepsis"),
  SP_vs_IC          = fisher_pdc("Sepsis_Pneumonia", "Infection_Control"),
  Sepsis_vs_IC      = fisher_pdc("Sepsis", "Infection_Control"),
  Sepsis_vs_Healthy = fisher_pdc("Sepsis", "Healthy"),
  COPD_vs_Healthy   = fisher_pdc("COPD", "Healthy")
)

# Write statistics summary
summ <- c(
  "=== Pulmonary vs non-pulmonary sepsis pDC redistribution ===",
  paste0("pDC cell count: Healthy=", cnt["pDC","Healthy"],
         " / IC=", cnt["pDC","Infection_Control"],
         " / COPD=", cnt["pDC","COPD"],
         " / Sepsis=", cnt["pDC","Sepsis"],
         " / Sepsis_Pneumonia=", cnt["pDC","Sepsis_Pneumonia"]),
  paste0("pDC proportion: Healthy=", round(prop["pDC","Healthy"],4),
         " / IC=", round(prop["pDC","Infection_Control"],4),
         " / COPD=", round(prop["pDC","COPD"],4),
         " / Sepsis=", round(prop["pDC","Sepsis"],4),
         " / SP=", round(prop["pDC","Sepsis_Pneumonia"],4)),
  paste0("pDC Ro/e: Healthy=", round(roe["pDC","Healthy"],3),
         " / IC=", round(roe["pDC","Infection_Control"],3),
         " / COPD=", round(roe["pDC","COPD"],3),
         " / Sepsis=", round(roe["pDC","Sepsis"],3),
         " / SP=", round(roe["pDC","Sepsis_Pneumonia"],3)),
  paste0("Fisher SP_vs_Sepsis p=", format(pvals$SP_vs_Sepsis, scientific=TRUE, digits=3)),
  paste0("Fisher SP_vs_IC p=", format(pvals$SP_vs_IC, scientific=TRUE, digits=3)),
  "",
  "Interpretation: if pDC depletion shows no significant difference between Sepsis_Pneumonia and Sepsis,",
  "      the pDC blood redistribution is a sepsis-shared phenomenon, not pulmonary-specific;",
  "      if SP depletion is significantly worse, this supports pulmonary sepsis specifically driving pDC sequestration."
)
writeLines(summ, file.path(out_dir, "path19_pdc_redistribution_stats.txt"))
cat("\nStatistics summary saved: path19_pdc_redistribution_stats.txt\n")

## =========================================================================
## Step 4: Visualization
## =========================================================================
cat("\n===== Step 4: Visualization =====\n")

# 4.1 pDC proportion bar plot (five groups)
df_prop <- data.frame(
  group = factor(groups, levels = groups),
  pDC_prop = prop["pDC", ],
  pDC_n = cnt["pDC", ]
)
p1 <- ggplot(df_prop, aes(x = group, y = pDC_prop, fill = group)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0("n=", pDC_n)), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c("Healthy"="#4CAF50","Infection_Control"="#2196F3",
                               "COPD"="#FF9800","Sepsis"="#F44336",
                               "Sepsis_Pneumonia"="#9C27B0")) +
  labs(x = NULL, y = "pDC proportion (of DCs)",
       title = "pDC depletion in blood: pulmonary vs non-pulmonary sepsis") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")
ggsave(file.path(out_dir, "path19_pdc_proportion.pdf"), p1, width = 6.5, height = 4.5)

# 4.2 Ro/e dot plot (three subtypes x five groups)
df_roe <- expand.grid(subtype = subtypes, group = groups, stringsAsFactors = FALSE)
df_roe$roe <- as.vector(roe)
df_roe$subtype <- factor(df_roe$subtype, levels = subtypes)
df_roe$group <- factor(df_roe$group, levels = groups)
p2 <- ggplot(df_roe, aes(x = group, y = roe, color = subtype, group = subtype)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) + geom_line() +
  labs(x = NULL, y = "Ro/e (observed / expected)",
       title = "DC subset enrichment across groups") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out_dir, "path19_dc_subtype_roe.pdf"), p2, width = 7, height = 4.5)

cat("Visualizations saved: path19_pdc_proportion.pdf / path19_dc_subtype_roe.pdf\n")

cat("\nPulmonary vs non-pulmonary sepsis pDC redistribution analysis complete\n")
cat("Outputs: path19_dc_subtype_count.csv / _proportion.csv / _roe.csv\n")
cat("      path19_pdc_redistribution_stats.txt / path19_pdc_proportion.pdf / path19_dc_subtype_roe.pdf\n")
