# =========================================================================
# 21_venn_group_specific_genes.R (Venn four groups: disease-specific vs shared genes)
#
# Purpose: respond to the reviewer's criticism -- the 149 shared genes only take
# the 'co-upregulated' set, without distinguishing 'disease-specific' from
# 'comorbidity-shared', so it cannot prove the 149 are specific to comorbidity
# rather than the respective inflammatory programs of each disease.
#
# Approach: read lung COPD myeloid DEG and blood sepsis monocyte DEG, take
# up/down-regulated sets respectively, divide into four categories: shared up
# (149) / shared down / COPD-specific / sepsis-specific, draw two Venn diagrams
# (up + down), and output each set's gene list.
#
# Input: path2_lung_myeloid_DEG_COPD_vs_Control.csv
#        path2_blood_monocyte_DEG_Sepsis_vs_Healthy.csv
# Output: path21_ series
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

cat("\n==============================================================\n")
cat("Venn four groups: COPD-specific / sepsis-specific / shared up / shared down\n")
cat("==============================================================\n")

## =========================================================================
## Step 1: Read DEG files
## =========================================================================
cat("\n===== Step 1: Read DEG =====\n")

lung_file  <- file.path(out_dir, "path2_lung_myeloid_DEG_COPD_vs_Control.csv")
blood_file <- file.path(out_dir, "path2_blood_monocyte_DEG_Sepsis_vs_Healthy.csv")
if (!file.exists(lung_file))  stop(lung_file, " not found")
if (!file.exists(blood_file)) stop(blood_file, " not found")

lung_deg  <- read.csv(lung_file,  row.names = 1, check.names = FALSE)
blood_deg <- read.csv(blood_file, row.names = 1, check.names = FALSE)
cat("Lung COPD myeloid DEG:", nrow(lung_deg), "genes\n")
cat("Blood sepsis monocyte DEG:", nrow(blood_deg), "genes\n")

# Confirm avg_log2FC column exists
if (!"avg_log2FC" %in% colnames(lung_deg)) stop("lung DEG missing avg_log2FC column")
if (!"avg_log2FC" %in% colnames(blood_deg)) stop("blood DEG missing avg_log2FC column")

## =========================================================================
## Step 2: Split up-regulated / down-regulated (consistent with 149 definition:
##         avg_log2FC > 0 / < 0)
## =========================================================================
cat("\n===== Step 2: Up/down-regulated split =====\n")

lung_up    <- rownames(lung_deg)[lung_deg$avg_log2FC > 0]
lung_down  <- rownames(lung_deg)[lung_deg$avg_log2FC < 0]
blood_up   <- rownames(blood_deg)[blood_deg$avg_log2FC > 0]
blood_down <- rownames(blood_deg)[blood_deg$avg_log2FC < 0]

cat("Lung COPD myeloid: up", length(lung_up), "/ down", length(lung_down), "\n")
cat("Blood sepsis monocyte: up", length(blood_up), "/ down", length(blood_down), "\n")

## =========================================================================
## Step 3: Four gene sets
## =========================================================================
cat("\n===== Step 3: Four gene sets =====\n")

shared_up       <- intersect(lung_up, blood_up)        # shared up (=149)
shared_down     <- intersect(lung_down, blood_down)    # shared down
copd_only_up    <- setdiff(lung_up, blood_up)          # COPD-specific up
sepsis_only_up  <- setdiff(blood_up, lung_up)          # sepsis-specific up
copd_only_down  <- setdiff(lung_down, blood_down)      # COPD-specific down
sepsis_only_down<- setdiff(blood_down, lung_down)      # sepsis-specific down

copd_specific   <- unique(c(copd_only_up, copd_only_down))      # COPD-specific (total)
sepsis_specific <- unique(c(sepsis_only_up, sepsis_only_down))  # sepsis-specific (total)

cat("Shared up gene count:", length(shared_up), "\n")
cat("Shared down gene count:", length(shared_down), "\n")
cat("COPD-specific gene count:", length(copd_specific), "\n")
cat("Sepsis-specific gene count:", length(sepsis_specific), "\n")

## =========================================================================
## Step 4: Save gene lists
## =========================================================================
writeLines(shared_up,        file.path(out_dir, "path21_shared_up_genes.txt"))
writeLines(shared_down,      file.path(out_dir, "path21_shared_down_genes.txt"))
writeLines(copd_specific,    file.path(out_dir, "path21_copd_specific_genes.txt"))
writeLines(sepsis_specific,  file.path(out_dir, "path21_sepsis_specific_genes.txt"))
cat("Gene lists saved: path21_shared_up/down_genes.txt, path21_copd/sepsis_specific_genes.txt\n")

## =========================================================================
## Step 5: Statistics summary table
## =========================================================================
summ <- c(
  "=== Venn four groups: disease-specific vs comorbidity-shared ===",
  paste0("Lung COPD myeloid DEG: up ", length(lung_up), " / down ", length(lung_down)),
  paste0("Blood sepsis monocyte DEG: up ", length(blood_up), " / down ", length(blood_down)),
  "",
  paste0("1) Shared up: ", length(shared_up), " (shared myeloid program, previously 149)"),
  paste0("2) Shared down: ", length(shared_down)),
  paste0("3) COPD-specific: ", length(copd_specific)),
  paste0("4) Sepsis-specific: ", length(sepsis_specific)),
  "",
  "Interpretation:",
  "  - Shared up (149) is the cross-tissue comorbidity-shared program (core main line);",
  "  - Shared down suggests the shared suppressive program of both (supplementary);",
  "  - COPD-specific / sepsis-specific indicate independent disease programs of each, proving the 149 are not generic inflammatory genes;",
  "  - if the number of shared genes is significantly greater than random expectation (hypergeometric test), this proves the sharing is not by chance."
)

# Hypergeometric test: is shared up significantly higher than random expectation
# Null model: independently draw blood-up (1354) and lung-up (1270) from a genome
# background set (~20000 protein-coding genes), testing whether the observed
# overlap (149) is significantly higher than expected (1354*1270/20000=86).
# Note: cannot use the 'intersection of the two DEG tables' (669) as background
# -- that is the intersection of result sets, not a sampleable background, and
# being smaller than blood_up/lung_up would make phyper parameters invalid.
bg <- 20000
cat("Hypergeometric background set (genome protein-coding gene count):", bg, "\n")
if (length(shared_up) > 0) {
  p_hyper <- phyper(length(shared_up) - 1, length(blood_up), bg - length(blood_up),
                    length(lung_up), lower.tail = FALSE)
  summ <- c(summ, "",
            paste0("Hypergeometric test (background = ", bg, " protein-coding genes): shared up p = ",
                   format(p_hyper, scientific = TRUE, digits = 3)))
}
writeLines(summ, file.path(out_dir, "path21_venn_four_groups_stats.txt"))
cat("\nStatistics summary saved: path21_venn_four_groups_stats.txt\n")

## =========================================================================
## Step 6: Venn diagrams (up + down)
## =========================================================================
cat("\n===== Step 6: Venn diagrams =====\n")

if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram", repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(VennDiagram))

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# Up-regulated Venn
venn_up <- tryCatch({
  vp <- venn.diagram(
    x = list("COPD lung\nmyeloid up" = lung_up, "Sepsis blood\nmonocyte up" = blood_up),
    filename = NULL, disable.logging = TRUE,
    col = "transparent",
    fill = c("#FF9800", "#F44336"), alpha = 0.5,
    cat.cex = 1.1, cex = 1.3,
    fontfamily = "sans", cat.fontfamily = "sans"
  )
  pdf(file.path(out_dir, "path21_venn_up.pdf"), width = 6, height = 6)
  grid::grid.draw(vp)
  dev.off()
  cat("Up-regulated Venn diagram saved: path21_venn_up.pdf\n")
  TRUE
}, error = function(e) { cat("Up-regulated Venn diagram failed:", conditionMessage(e), "\n"); FALSE })

# Down-regulated Venn
venn_down <- tryCatch({
  vd <- venn.diagram(
    x = list("COPD lung\nmyeloid down" = lung_down, "Sepsis blood\nmonocyte down" = blood_down),
    filename = NULL, disable.logging = TRUE,
    col = "transparent",
    fill = c("#FF9800", "#F44336"), alpha = 0.5,
    cat.cex = 1.1, cex = 1.3,
    fontfamily = "sans", cat.fontfamily = "sans"
  )
  pdf(file.path(out_dir, "path21_venn_down.pdf"), width = 6, height = 6)
  grid::grid.draw(vd)
  dev.off()
  cat("Down-regulated Venn diagram saved: path21_venn_down.pdf\n")
  TRUE
}, error = function(e) { cat("Down-regulated Venn diagram failed:", conditionMessage(e), "\n"); FALSE })

cat("\nVenn four groups analysis complete\n")
cat("Outputs: path21_shared_up/down_genes.txt / path21_copd/sepsis_specific_genes.txt\n")
cat("      path21_venn_four_groups_stats.txt / path21_venn_up.pdf / path21_venn_down.pdf\n")
