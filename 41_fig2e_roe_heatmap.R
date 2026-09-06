# =============================================================================
# 41_Fig2_E_roe_heatmap.R
# Figure: Fig2 panel F — Ratio of observed to expected (Ro/e) for major cell types
#       Note: panel letters have been re-ordered (E -> F) to free slot D for the
#       marker-gene DotPlot.
#       across the five groups (Healthy / Infection_Control / COPD / Sepsis /
#       Sepsis_Pneumonia).
# Notes: the data are already computed in path1_roe_celltype_by_group.csv (script 01;
#       Ro/e = cell-type fraction within each group / pooled fraction over the whole
#       cohort). This script only renders them as a publication-grade heatmap.
#       Ro/e = 1 means consistent with the overall cohort; >1 enriched (red),
#       <1 depleted (blue). Follows Fig2E of the MASH paper.
# Run:   works locally or on the server (only depends on ggplot2, no Seurat needed)
#   Rscript 41_Fig2_E_roe_heatmap.R [DATA_DIR] [OUT_DIR]
#   DATA_DIR / OUT_DIR both default to the directory containing this script.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
script_dir <- tryCatch(dirname(normalizePath(sys.frames()[[1]]$ofile)),
                       error = function(e) getwd())
DATA_DIR <- ifelse(length(args) >= 1 && nzchar(args[1]), args[1], script_dir)
OUT_DIR  <- ifelse(length(args) >= 2 && nzchar(args[2]), args[2], script_dir)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

roe_file <- file.path(DATA_DIR, "path1_roe_celltype_by_group.csv")
stopifnot(file.exists(roe_file))
roe <- read.csv(roe_file, row.names = 1, check.names = FALSE)
grp_order <- c("Healthy", "Infection_Control", "COPD", "Sepsis", "Sepsis_Pneumonia")
ct_order  <- c("T_cell", "NK", "B_cell", "Monocyte", "DC", "Mast",
               "Platelet", "RBC")
roe <- roe[, intersect(grp_order, colnames(roe)), drop = FALSE]
ct_keep  <- intersect(ct_order, rownames(roe))
ct_keep  <- c(ct_keep, setdiff(rownames(roe), ct_keep))   # the rest appended at the end
roe <- roe[ct_keep, , drop = FALSE]

dat <- cbind(CellType = rownames(roe), as.data.frame(roe)) %>%
  pivot_longer(-CellType, names_to = "Group", values_to = "RoE") %>%
  mutate(CellType = factor(CellType, levels = ct_keep),
         Group    = factor(Group, levels = colnames(roe)))

# Unified five-group palette for the whole manuscript (consistent with scripts 38/40),
# used only to colour the column labels
GROUP_COLS <- c(Healthy = "#3C5488", Infection_Control = "#91D1C2",
                COPD = "#00A087", Sepsis = "#E64B35", Sepsis_Pneumonia = "#F39B7F")
y_max <- max(dat$RoE, na.rm = TRUE)
y_max <- ifelse(y_max < 4, 4, ceiling(y_max))

p_e <- ggplot(dat, aes(x = Group, y = CellType, fill = RoE)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", RoE)),
            colour = ifelse(dat$RoE > y_max * 0.55, "white", "grey20"),
            size = 3.6) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 1, limit = c(0, y_max),
                       name = "Ro/e") +
  scale_x_discrete(labels = function(l)
    sapply(strsplit(l, "_"), function(x) paste(toupper(substr(x,1,1)),
                                                tolower(substr(x,2,nchar(x))),
                                                sep = "", collapse = " "))) +
  labs(x = NULL, y = NULL,
       title = "E  Ro/e of major cell types by group",
       subtitle = "observed / expected (pooled baseline); red = enriched, blue = depleted") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9,
                                    colour = GROUP_COLS[levels(dat$Group)]),
         axis.text.y = element_text(size = 10),
         panel.grid = element_blank(),
         legend.title = element_text(size = 10),
         plot.title = element_text(size = 12, face = "bold"),
         plot.subtitle = element_text(size = 8, colour = "grey40"))

ggsave(file.path(OUT_DIR, "Fig2F_roe_by_group.pdf"), p_e,
       width = 6.4, height = 4.8)
ggsave(file.path(OUT_DIR, "Fig2F_roe_by_group.png"), p_e,
       width = 6.4, height = 4.8, dpi = 300, device = grDevices::png)
write.csv(dat, file.path(OUT_DIR, "path41_roe_tidy.csv"), row.names = FALSE)
cat("saved: Fig2F_roe_by_group.{pdf,png}  +  path41_roe_tidy.csv\n")
cat("Ro/e range:", sprintf("%.2f", min(dat$RoE)), "~", sprintf("%.2f", max(dat$RoE)), "\n")
