#!/usr/bin/env Rscript
# 62e: re-render the DC communication panel (Fig. S3)
#
# Task: re-render Fig. S3 (DC communication per DC cell) using the same 5-group
# palette as Fig. S4 (purple ME by group), fixing the inconsistency of the
# earlier Fig. S3, which used the ggplot2 default hue palette and therefore did
# not match the color scheme of the rest of the manuscript.
#
# How to run:
#   1) On the server (or any environment without encoding problems), cd into the
#      working directory and run:
#        Rscript server_version/62e_redraw_dc_communication.R
#   2) On Windows, when Rscript cannot parse source strings under a non-ASCII
#      path, see the local bridge script
#      C:/Users/Administrator/AppData/Local/Temp/run_s3_paint.py
#      or copy this script into a pure-ASCII temporary directory and run it there.
#
# Input:  working directory/server_version/path9_dc_norm_comm.csv
# Output: Figure/FigS3_dc_comm.{pdf,png}
#         server_version/FigS3_dc_comm.{pdf,png}
#         server_version/path9_dc_norm_comm_v2.pdf  (historical version, does not overwrite the original pdf)

args      <- commandArgs(trailingOnly = TRUE)
workdir   <- if (length(args) >= 1 && nzchar(args[1])) args[1] else getwd()
setwd(workdir)

csv_path   <- file.path(workdir, "path9_dc_norm_comm.csv")
out_pdf    <- file.path(workdir, "FigS3_dc_comm.pdf")
out_png    <- file.path(workdir, "FigS3_dc_comm.png")
hist_pdf   <- file.path(workdir, "path9_dc_norm_comm_v2.pdf")

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
})

## 5-group palette (identical to 48_Fig5_DEF_IFN_module.R / Fig. S4)
GROUP_LEVELS5 <- c("Healthy", "Infection_Control", "COPD", "Sepsis", "Sepsis_Pneumonia")
GROUP_COLS5   <- c(Healthy          = "#3C5488",
                   Infection_Control = "#91D1C2",
                   COPD             = "#00A087",
                   Sepsis           = "#E64B35",
                   Sepsis_Pneumonia = "#F39B7F")

## R reads "MHC-II_perDC" from the path9 csv as "MHC.II_perDC" by default
plot_mods <- c("MHC.II_perDC", "Costim_perDC")
panel_lbl <- c("MHC.II_perDC" = "MHC-II per DC",
               "Costim_perDC" = "Costim per DC")

norm_raw <- read.csv(csv_path, stringsAsFactors = FALSE)
cat("=== 62e: re-render the DC communication panel (Fig. S3) ===\n")
cat("  source rows:", nrow(norm_raw), "\n")

nl <- tidyr::pivot_longer(norm_raw[, c("group", plot_mods)],
                          cols      = plot_mods,
                          names_to  = "module",
                          values_to = "prob_perDC")
nl$group      <- factor(nl$group,  levels = GROUP_LEVELS5)
nl$module     <- factor(nl$module, levels = plot_mods)
nl$module_lbl <- factor(panel_lbl[as.character(nl$module)],
                        levels = unname(panel_lbl))

p2 <- ggplot(nl, aes(x = group, y = prob_perDC, fill = group)) +
  geom_col(alpha = 0.85) +
  facet_wrap(~ module_lbl, scales = "free_y") +
  scale_fill_manual(values = GROUP_COLS5, guide = "none") +
  labs(title = "Communication strength per DC cell (normalized)",
       y        = "prob / DC count") +
  theme_bw(base_size = 11) +
  theme(axis.text.x        = element_text(angle = 30, hjust = 1),
        legend.position    = "none",
        plot.title         = element_text(size = 11),
        strip.background   = element_rect(fill = "white"),
        panel.grid.minor   = element_blank())

## Same 9 x 5 in size as the original script 09
ggsave(hist_pdf, p2, width = 9, height = 5)
ggsave(out_pdf,  p2, width = 9, height = 5)
ggsave(out_png,  p2, width = 9, height = 5, dpi = 300, device = grDevices::png)

cat("\nwritten:\n")
cat("  ", hist_pdf, "\n")
cat("  ", out_pdf,  "\n")
cat("  ", out_png,  "\n")
cat("\ncolors (Fig. S4 style):\n")
for (g in GROUP_LEVELS5) cat(sprintf("  %-20s %s\n", g, GROUP_COLS5[g]))
