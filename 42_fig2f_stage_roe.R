# =============================================================================
# 42_Fig2_F_stage_roe.R
# Intent: Fig2 panel F — Ro/e comparison between two disease-severity stages,
#          a faithful analogue of "eMASH vs aMASH" in the MASH paper.
#          This project has no early/advanced MASH; the closest match to the
#          "same-disease progression" logic is the intra-sepsis severity axis:
#          early  = Sepsis            (sepsis without pneumonia)
#          advanced = Sepsis_Pneumonia (sepsis with pneumonia, more severe)
# Method: take the fraction of each major cell type per group from the integrated object
#         -> compute Ro/e per group (vs the pooled whole-cohort baseline)
#         -> draw a lollipop plot of the log2(RoE_advanced / RoE_early) shift
#            (red = enriched with progression, blue = depleted).
# Criteria: consistent with the other Fig2 panels — use cp3_annotated.rds, keep only the five
#           core groups, drop LowQuality.
# Run: server (requires Seurat / reads the integrated object)
#   Rscript 42_Fig2_F_stage_roe.R [OUT_DIR]
# Outputs: Fig2F_stage_roe.pdf / .png
#          path42_stage_roe.csv        (pre-plot data: per-cell-type proportions / RoE / shift
#                                      for the two groups, for local debugging)
#          path42_stage_roe_plotdata.rds (= plotting input, re-renderable locally)
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
script_dir <- tryCatch(dirname(normalizePath(sys.frames()[[1]]$ofile)),
                       error = function(e) getwd())
OUT_DIR <- ifelse(length(args) >= 1 && nzchar(args[1]), args[1], script_dir)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

## ---- Load the integrated object (same source as Fig2: cp3_annotated.rds) ----
cp3_file <- file.path(OUT_DIR, "cp3_annotated.rds")
if (!file.exists(cp3_file)) cp3_file <- file.path(OUT_DIR, "path1_sepsis_copd_integrated.rds")
stopifnot(file.exists(cp3_file))
seu <- readRDS(cp3_file)
cat("Loaded:", basename(cp3_file), "\n")

meta <- seu@meta.data
core_groups <- c("Healthy", "Infection_Control", "COPD", "Sepsis", "Sepsis_Pneumonia")
meta <- meta[!is.na(meta$cell_type) & meta$cell_type != "LowQuality" &
             !is.na(meta$group) & meta$group %in% core_groups, ]

## ---- Fractions of the major cell types per group ----
tab <- table(meta$cell_type, meta$group)[, core_groups, drop = FALSE]
prop <- prop.table(tab, margin = 2)                 # normalize within column (fraction within each group)
expected <- rowSums(tab) / sum(tab)                 # pooled whole-cohort fraction (weighted by cells per group)
RoE <- sweep(prop, 1, expected, "/")                # Ro/e per group

early   <- "Sepsis"
adv     <- "Sepsis_Pneumonia"
RoE_e   <- RoE[, early]
RoE_a   <- RoE[, adv]
log2shift <- log2(RoE_a / RoE_e)                    # >0 = enriched with progression

res <- data.frame(
  CellType = rownames(RoE),
  Prop_Sepsis = as.numeric(prop[, early]),
  Prop_SepsisPneumonia = as.numeric(prop[, adv]),
  RoE_Sepsis = as.numeric(RoE_e),
  RoE_SepsisPneumonia = as.numeric(RoE_a),
  log2_shift = as.numeric(log2shift)
)
# keep only cell types present in both stages (avoids Inf/-Inf in log2 breaking the plot)
res <- res[is.finite(res$log2_shift), ]
# keep only the numeric columns, so that round() is not applied to the character CellType column
num_cols <- vapply(res, is.numeric, logical(1))
print(round(res[, num_cols], 3))

## =============================================================================
## Step: save pre-plot data (for local debugging and re-rendering; no longer needs Seurat or the integrated object)
## =============================================================================
GROUP_COLS <- c(Healthy = "#3C5488", Infection_Control = "#91D1C2",
                COPD = "#00A087", Sepsis = "#E64B35", Sepsis_Pneumonia = "#F39B7F")
plot_data <- list(
  res = res,
  GROUP_COLS = GROUP_COLS,
  early = early, adv = adv,
  ylim = max(abs(res$log2_shift), 0.1) * 1.25
)
saveRDS(plot_data, file.path(OUT_DIR, "path42_stage_roe_plotdata.rds"))
write.csv(res, file.path(OUT_DIR, "path42_stage_roe.csv"), row.names = FALSE)
cat("Saved: path42_stage_roe.csv + path42_stage_roe_plotdata.rds\n")

## =============================================================================
## Plot: log2 shift lollipop plot
## =============================================================================
res_o <- res[order(res$log2_shift, decreasing = TRUE), ]
res_o$CellType <- factor(res_o$CellType, levels = res_o$CellType)
res_o$sign <- ifelse(res_o$log2_shift >= 0, "enriched in advanced", "depleted in advanced")

p_f <- ggplot(res_o, aes(x = log2_shift, y = CellType, colour = sign)) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed", linewidth = 0.5) +
  geom_segment(aes(x = 0, xend = log2_shift, y = CellType, yend = CellType),
               colour = "grey70", linewidth = 0.7) +
  geom_point(size = 3.2) +
  scale_colour_manual(values = c("enriched in advanced" = "#B2182B",
                                 "depleted in advanced" = "#2166AC"),
                      name = NULL) +
  geom_text(aes(label = sprintf("%.2f", log2_shift),
                x = log2_shift + 0.06 * sign(log2_shift)),
            size = 3.2, colour = "grey20", hjust = 0.5,
            angle = 90, vjust = ifelse(res_o$log2_shift >= 0, -0.2, 1.2)) +
  coord_cartesian(xlim = c(-plot_data$ylim, plot_data$ylim)) +
  labs(x = expression(log[2] ~ (RoE[advanced] / RoE[early])),
       y = NULL,
       title = "F  Ro/e shift along sepsis severity",
       subtitle = "Sepsis (early)  ->  Sepsis_Pneumonia (advanced)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey40"))

ggsave(file.path(OUT_DIR, "Fig2F_stage_roe.pdf"), p_f, width = 6.0, height = 4.6)
ggsave(file.path(OUT_DIR, "Fig2F_stage_roe.png"), p_f, width = 6.0, height = 4.6, dpi = 300, device = grDevices::png)
cat("Saved: Fig2F_stage_roe.{pdf,png}\n")
cat("\n===== script 42 done =====\n")
