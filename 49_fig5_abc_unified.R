# =========================================================================
# Script 49: Fig5 A/B/C re-rendered in a unified style
#
# Purpose: re-render the three method panels A / B / C of Fig5 in a unified style
#   A  module-trait correlation bars: bulk WGCNA (GSE66099 whole blood),
#      13 modules x sepsis correlation
#      (data: path12_bulk_module_trait.csv) - highlighting blue / turquoise / magenta
#   B  distribution of the 149 shared genes across bulk modules (including unmapped)
#      (data: path12_bulk_149gene_module.csv) - faithfully shows the 82/149 mapping,
#      with the majority in blue
#   C  single-cell hdWGCNA co-expression module dendrogram (purple = ISG module)
#      (object: path13_sc_hdwgcna.rds, hdWGCNA::PlotDendrogram, base device)
#
# Output (Fig5A_* / Fig5B_* / Fig5C_*, PDF+PNG):
#   Fig5A_bulk_module_trait.pdf/.png
#   Fig5B_149_bulk_modules.pdf/.png
#   Fig5C_sc_hdwgcna_dendrogram.pdf/.png
#   path49_fig5abc_plotdata.rds       data behind the figures
#
# Server: Rscript this script
# =========================================================================

.this_file <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
if (.this_file == "" || !file.exists(.this_file)) {
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1])
}
SCRIPT_DIR <- dirname(.this_file)
config_file <- file.path(SCRIPT_DIR, "00_config.R")
if (!file.exists(config_file)) stop("Cannot find 00_config.R: ", config_file)
source(config_file)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr)
})
has_hdwgcna <- requireNamespace("hdWGCNA", quietly = TRUE)

## Standard WGCNA module colours (consistent with the dendrograms in 13/34)
MOD_COL <- c(turquoise = "#00CED1", blue = "#0000FF", brown = "#A52A2A",
             green = "#00CD00", yellow = "#FFD700", red = "#FF0000",
             grey = "#BEBEBE", black = "#000000", pink = "#FF69B4",
             magenta = "#FF00FF", tan = "#D2B48C", greenyellow = "#ADFF2F",
             salmon = "#FA8072", purple = "#A020F0")

sig_star <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
}

cat("\n==============================================================\n")
cat("Fig5 A/B/C re-rendered in a unified style\n")
cat("==============================================================\n")

## ===========================================================================
## A — bulk WGCNA module-trait correlation (13 modules x sepsis)
## ===========================================================================
cat("\n########## A: module-trait correlation ==========\n")
trait_file <- file.path(out_dir, "path12_bulk_module_trait.csv")
if (!file.exists(trait_file)) stop("Missing path12_bulk_module_trait.csv (run 12_WGCNA.R first)")
mt <- read.csv(trait_file, check.names = FALSE, stringsAsFactors = FALSE)
mt$sig <- sig_star(mt$pvalue)
mt$module <- factor(mt$module, levels = mt$module[order(mt$cor_sepsis)])
mt$fill <- unname(MOD_COL[as.character(mt$module)])
mt$fill[is.na(mt$fill)] <- "grey60"
mt$is_key <- mt$module %in% c("blue", "turquoise", "magenta")
mt$alpha <- ifelse(mt$is_key, 1, 0.55)

p_a <- ggplot(mt, aes(cor_sepsis, module)) +
  geom_col(fill = mt$fill, alpha = mt$alpha, width = 0.68) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f%s", cor_sepsis, sig)),
            hjust = ifelse(mt$cor_sepsis >= 0, -0.25, 1.25),
            size = 3, colour = "grey25") +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.28))) +
  labs(title = "Bulk WGCNA modules vs sepsis status",
       subtitle = "Whole blood (GSE66099); module eigengene correlation with sepsis. Blue = neutrophil/\nemergency-granulopoiesis (r = 0.68), turquoise = lymphocyte (r = -0.60), magenta = IFN module (r = 0.01).",
       x = expression(italic(r) ~ "(eigengene vs sepsis)"), y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_text(size = 9, colour = mt$fill[order(mt$cor_sepsis)]),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.subtitle = element_text(size = 7.5, colour = "grey35"))
ggsave(file.path(out_dir, "Fig5A_bulk_module_trait.pdf"), p_a, width = 7, height = 5)
ggsave(file.path(out_dir, "Fig5A_bulk_module_trait.png"), p_a, width = 7, height = 5,
       dpi = 300, device = grDevices::png)
cat("Saved: Fig5A_bulk_module_trait.{pdf,png}\n")

## ===========================================================================
## B — distribution of the 149 shared genes across bulk modules
## ===========================================================================
cat("\n########## B: 149 -> bulk module distribution ==========\n")
f149 <- file.path(out_dir, "path12_bulk_149gene_module.csv")
if (!file.exists(f149)) stop("Missing path12_bulk_149gene_module.csv")
b149 <- read.csv(f149, check.names = FALSE, stringsAsFactors = FALSE)
gcol <- intersect(c("gene", "gene_name", "gene_id"), colnames(b149))[1]
if (is.na(gcol) || !"module" %in% colnames(b149)) stop("Column names of the 149-bulk CSV do not match")
b149 <- b149[!grepl("///", b149[[gcol]]), ]
cnt <- as.data.frame(table(b149$module), stringsAsFactors = FALSE)
names(cnt) <- c("module", "n")
unmapped_n <- 149 - sum(cnt$n)
cnt <- rbind(cnt, data.frame(module = "unmapped", n = unmapped_n))
cnt$module <- factor(cnt$module,
                     levels = cnt$module[order(-cnt$n)])
cnt$fill <- unname(MOD_COL[as.character(cnt$module)])
cnt$fill[is.na(cnt$fill)] <- "grey60"
cnt$fill[cnt$module == "unmapped"] <- "grey75"

p_b <- ggplot(cnt, aes(n, module)) +
  geom_col(fill = cnt$fill, width = 0.68) +
  geom_text(aes(label = n), hjust = -0.3, size = 3.2, colour = "grey25") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(title = "Mapping of the 149 shared genes onto bulk modules",
       subtitle = sprintf("%d of 149 genes mapped to a GSE66099 co-expression module; the majority fall in the\nsepsis-correlated blue module (r = 0.68), the remainder in the IFN module (magenta, n = %d).",
                          149 - unmapped_n,
                          if ("magenta" %in% cnt$module) cnt$n[cnt$module == "magenta"] else 0),
       x = "Number of shared genes", y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.subtitle = element_text(size = 8, colour = "grey35"))
ggsave(file.path(out_dir, "Fig5B_149_bulk_modules.pdf"), p_b, width = 7, height = 4.6)
ggsave(file.path(out_dir, "Fig5B_149_bulk_modules.png"), p_b, width = 7, height = 4.6,
       dpi = 300, device = grDevices::png)
cat("Saved: Fig5B_149_bulk_modules.{pdf,png}\n")

## ===========================================================================
## C — single-cell hdWGCNA co-expression module dendrogram (base device)
## ===========================================================================
cat("\n########## C: single-cell hdWGCNA dendrogram ==========\n")
c_ok <- FALSE
rds_file <- file.path(out_dir, "path13_sc_hdwgcna.rds")
if (!has_hdwgcna) {
  cat("!! hdWGCNA is not installed, skipping panel C (path34_hdwgcna_dendrogram.pdf from script 34 can be used instead)\n")
} else if (!file.exists(rds_file)) {
  cat("!! Missing path13_sc_hdwgcna.rds (run the single-cell hdWGCNA script or 34 first), skipping panel C\n")
} else {
  seu <- tryCatch(readRDS(rds_file), error = function(e) {
    cat("!! Failed to read rds:", conditionMessage(e), "\n"); NULL
  })
  if (!is.null(seu)) {
    c_ok <- tryCatch({
      suppressPackageStartupMessages(library(hdWGCNA))
      pdf(file.path(out_dir, "Fig5C_sc_hdwgcna_dendrogram.pdf"),
          width = 10, height = 7.5)
      PlotDendrogram(seu, wgcna_name = "myeloid",
                     main = "Single-cell hdWGCNA of blood monocytes\n(co-expression modules)")
      dev.off()
      png(file.path(out_dir, "Fig5C_sc_hdwgcna_dendrogram.png"),
          width = 10, height = 7.5, units = "in", res = 300)
      PlotDendrogram(seu, wgcna_name = "myeloid",
                     main = "Single-cell hdWGCNA of blood monocytes\n(co-expression modules)")
      dev.off()
      TRUE
    }, error = function(e) {
      cat("!! Dendrogram failed:", conditionMessage(e), "\n")
      tryCatch(dev.off(), error = function(e2) NULL)
      FALSE
    })
    if (c_ok) cat("Saved: Fig5C_sc_hdwgcna_dendrogram.{pdf,png}\n")
  }
}

## ===========================================================================
## Data behind the figures
## ===========================================================================
plotdata <- list(
  a = list(trait = mt, gg = p_a),
  b = list(dist = cnt, gg = p_b),
  c = list(ok = c_ok, note = "hdWGCNA base dendrogram; see PDF/PNG directly"),
  mod_col = MOD_COL
)
saveRDS(plotdata, file.path(out_dir, "path49_fig5abc_plotdata.rds"))
sz <- round(file.info(file.path(out_dir, "path49_fig5abc_plotdata.rds"))$size / 1e6, 2)
cat("Saved path49_fig5abc_plotdata.rds (", sz, " MB )\n")
cat("\n===== Script 49 done =====\n")
