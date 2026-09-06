# =============================================================================
# 43_Fig2_F_slope.R
# Fig2 panel G (multi-group Ro/e slope version, aligned with and extending the MASH Fig2.F layout):
#   Note: panel letters were reshuffled (F -> G) to free the D slot for the marker-gene DotPlot.
#   MASH Fig2.F = x=cell type, y=Ro/e, 2 lines for the two stage groups.
#   This project has 5 disease groups, so this version draws 5 lines (one per group) across
#   all major cell types, reusing the Fig2 GROUP_COLS five-group palette; more informative
#   than plotting only 2 groups.
# Data: path1_roe_celltype_by_group.csv (8 cell types x 5 groups Ro/e, produced by script 01)
# Run: works locally or on the server (reads a CSV only, no Seurat dependency)
#   Rscript 43_Fig2_F_slope.R [DIR]
# Outputs: Fig2G_slope_roe.pdf / .png  (overwrites the old version)
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
script_dir <- tryCatch(dirname(normalizePath(sys.frames()[[1]]$ofile)),
                       error = function(e) getwd())
OUT_DIR <- ifelse(length(args) >= 1 && nzchar(args[1]), args[1], script_dir)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Shares the cell-type order with Fig2G (single source of truth, avoids drift)
source(file.path(script_dir, "fig2_major_celltypes.R"))

csv_file <- file.path(OUT_DIR, "path1_roe_celltype_by_group.csv")
stopifnot(file.exists(csv_file))
df <- read.csv(csv_file, row.names = 1, check.names = FALSE)

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

# wide -> long
df$CellType <- rownames(df)
long <- df %>% tidyr::pivot_longer(cols = -CellType, names_to = "Group", values_to = "RoE")

# Fixed order: based on the 8 major cell types in the CSV, ordered by FIG2_CELLTYPE_ORDER
celltype_order <- FIG2_CELLTYPE_ORDER[FIG2_CELLTYPE_ORDER %in% rownames(df)]
group_order    <- c("Healthy", "Infection_Control", "COPD", "Sepsis", "Sepsis_Pneumonia")
GROUP_COLS <- c(Healthy = "#3C5488", Infection_Control = "#91D1C2",
                COPD = "#00A087", Sepsis = "#E64B35", Sepsis_Pneumonia = "#F39B7F")

long$CellType <- factor(long$CellType, levels = celltype_order[celltype_order %in% long$CellType])
long$Group    <- factor(long$Group, levels = group_order)

# Save the data before plotting
saveRDS(list(res = long, celltype_order = celltype_order, group_order = group_order,
             GROUP_COLS = GROUP_COLS),
        file.path(OUT_DIR, "path43_fig2f_plotdata.rds"))

p_f <- ggplot(long, aes(x = CellType, y = RoE, colour = Group, group = Group)) +
  geom_hline(yintercept = 1, colour = "grey50", linetype = "dashed", linewidth = 0.5) +
  geom_line(linewidth = 1.1, alpha = 0.9) +
  geom_point(size = 2.4, stroke = 0.3, shape = 16) +
  scale_colour_manual(values = GROUP_COLS, name = NULL) +
  # Ro/e is a ratio; a log2 axis centred on 1 keeps enrichment/depletion symmetric and readable (for a linear axis: just delete this line)
  scale_y_continuous(trans = "log2",
                     breaks = c(0.25, 0.5, 1, 2, 4),
                     labels = c("0.25", "0.5", "1", "2", "4"),
                     limits = c(0.12, 4.5)) +
  scale_x_discrete(expand = c(0.04, 0.04)) +
  labs(x = NULL, y = "Ro/e (log2 scale; dashed = 1)",
       title = "F  Ro/e of major cell types across disease groups",
       subtitle = "one line per group; colour = disease severity") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey40"),
        axis.text.x = element_text(size = 10, angle = 25, hjust = 1))

ggsave(file.path(OUT_DIR, "Fig2G_slope_roe.pdf"), p_f, width = 7.2, height = 5.0)
cat("Saved: Fig2G_slope_roe.pdf\n")

tryCatch({
  grDevices::png(file.path(OUT_DIR, "Fig2G_slope_roe.png"),
                 width = 7.2, height = 5.0, units = "in", res = 300)
  print(p_f); grDevices::dev.off()
  cat("Saved: Fig2G_slope_roe.png\n")
}, error = function(e) cat("PNG writing skipped (can be backfilled from the PDF with PyMuPDF):", conditionMessage(e), "\n"))

cat("\n===== script 43 done (5-group Ro/e slope version) =====\n")
