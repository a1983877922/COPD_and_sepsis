###############################################################################
# Sepsis + COPD comorbidity analysis — Script 40: Figure 3c — 149 shared genes five-group heatmap (server)
# File (current name): 40_fig3c_group_heatmap.R (originally: 40_Fig3c five-group heatmap.R)
#
# Purpose: add Figure 3 panel c — use a heatmap to visualize the overall gradient of the "149 shared myeloid
#       program across disease severity (Healthy → Infection_Control → COPD → Sepsis →
#       Sepsis_Pneumonia)".
#   - rows  = 149 shared genes (path2_shared_myeloid_genes.txt)
#   - cols  = blood monocyte 5 groups (donor-level)
#   - value = each gene first z-normalized across all donors, then donor mean taken per group (row z)
#   - left bar = single-cell hdWGCNA module assignment (path13_sc_149gene_module.csv, grey = unassigned)
#            — connects to the module story in Figure 5
#
# Data scope (exactly consistent with scripts 28/28b/31, to avoid "switching to yet another scope" objections):
#   - blood object: path1_sepsis_copd_integrated.rds (saved by script 01, contains cell_type/
#             group/dataset + patient/sample donor fields)
#   - cells:   cell_type == "Monocyte"; all five groups kept
#   - donor-level: detect_donor (patient first, sample fallback, same as 28b/31)
#   - expression:   per-donor pseudobulk counts → edgeR DGEList + TMM → logCPM
#   - 5-group palette: exactly the same as GROUP_COLORS in script 38 (Fig2), unified across the manuscript
#
# Usage: (server, same as scripts 01~39)
#   cd /media/desk16/ysx5991/sepsis_copd/01script
#   Rscript 40_fig3c_group_heatmap.R
#
# Output (to out_dir):
#   Fig3c_149_heatmap_five_groups.pdf / .png    heatmap (row z)
#   path40_149_group_z.csv          149 x 5 group z means (long table, with module)
#   path40_149_donor_logcpm.rds     donor-level logCPM (149 x donors) for verification
#   path40_donor_meta.csv           donor - group - cell count
#   path40_fig3c_plotdata.rds       plotting input saved before plotting (for local re-render/debug)
###############################################################################

## ---- Auto-locate script directory ----
.this_file <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
if (.this_file == "" || !file.exists(.this_file)) {
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1])
}
SCRIPT_DIR <- dirname(.this_file)
config_file <- file.path(SCRIPT_DIR, "00_config.R")
if (!file.exists(config_file)) stop("cannot find 00_config.R: ", config_file)
source(config_file)
out_dir <- SCRIPT_DIR
cat("Output directory:", out_dir, "\n")

suppressPackageStartupMessages({
  library(Seurat); library(edgeR); library(Matrix)
})

## ---- Detect donor field (patient first, sample fallback) [same as 28b/31] ----
detect_donor <- function(seu) {
  meta <- seu@meta.data
  get <- function(f) if (f %in% colnames(meta)) as.character(meta[[f]]) else NULL
  patient <- get("patient"); sample <- get("sample")
  if (!is.null(patient)) {
    if (!is.null(sample)) {
      na_idx <- is.na(patient) | patient == ""
      patient[na_idx] <- sample[na_idx]
      return(list(field = "patient(+sample)", vec = patient))
    }
    return(list(field = "patient", vec = patient))
  }
  for (cand in c("donor_id", "donor", "subject", "subject_id", "Subject_Identity", "sample")) {
    v <- get(cand)
    if (!is.null(v)) return(list(field = cand, vec = v))
  }
  NULL
}

GROUP_LEVELS <- c("Healthy", "Infection_Control", "COPD", "Sepsis", "Sepsis_Pneumonia")
GROUP_COLS   <- c(Healthy = "#3C5488", Infection_Control = "#91D1C2",
                  COPD = "#F39B7F", Sepsis = "#E64B35", Sepsis_Pneumonia = "#8491B4")

## =========================================================================
## Step 1: read blood object + monocytes + donor-level pseudobulk
## =========================================================================
blood_file <- file.path(out_dir, "path1_sepsis_copd_integrated.rds")
if (!file.exists(blood_file)) stop("Cannot find blood object: ", blood_file, " (run script 01 first)")
cat("Loading blood object:", blood_file, "\n")
blood <- readRDS(blood_file)
cat("Blood object meta columns:", paste(colnames(blood@meta.data), collapse = ", "), "\n")

donor_info <- detect_donor(blood)
if (is.null(donor_info)) stop("Blood object has no donor field; cannot do donor-level aggregation")
cat("Donor field:", donor_info$field, "\n")

meta <- blood@meta.data
keep_cell <- !is.na(meta$cell_type) & meta$cell_type == "Monocyte" &
             !is.na(meta$group) & meta$group %in% GROUP_LEVELS
cat("Monocytes (5 groups):", sum(keep_cell), "\n")
print(table(meta$group[keep_cell]))

donor_vec <- donor_info$vec[keep_cell]
grp_vec   <- as.character(meta$group[keep_cell])
counts    <- GetAssayData(blood, assay = "RNA", layer = "counts")
counts    <- counts[, rownames(meta)[keep_cell], drop = FALSE]
rm(blood, meta); gc()

donors <- unique(donor_vec); donors <- donors[!is.na(donors) & donors != ""]
cat("Number of donors:", length(donors), "\n")
cts <- do.call(cbind, lapply(donors, function(d) {
  Matrix::rowSums(counts[, donor_vec == d, drop = FALSE])
}))
colnames(cts) <- donors
donor_grp <- vapply(donors, function(d) unique(grp_vec[donor_vec == d])[1], character(1))
cat("Donors per group:\n"); print(table(donor_grp))
rm(counts); gc()

# Keep genes expressed in >=2 donors for TMM normalization (edgeR standard pipeline)
y <- DGEList(counts = cts, group = factor(donor_grp, levels = GROUP_LEVELS))
keep_g <- filterByExpr(y, min.count = 5, min.total.count = 10)
y <- y[keep_g, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
logcpm <- cpm(y, log = TRUE, prior.count = 1)     # gene x donor
cat("Donor-level logCPM:", nrow(logcpm), "genes x", ncol(logcpm), "donors\n")

## =========================================================================
## Step 2: 149 genes + z-normalization + group means
## =========================================================================
f149 <- file.path(out_dir, "path2_shared_myeloid_genes.txt")
if (!file.exists(f149)) f149 <- file.path(out_dir, "path21_shared_up_genes.txt")
genes149 <- trimws(readLines(f149)); genes149 <- genes149[nzchar(genes149)]
genes149 <- intersect(genes149, rownames(logcpm))
cat("Genes among 149 expressed at donor level:", length(genes149), "\n")

X <- logcpm[genes149, , drop = FALSE]             # 149 x donor
Xz <- t(scale(t(X)))                               # row (gene) z-normalized across donors
if (any(is.na(Xz))) Xz[is.na(Xz)] <- 0            # genes with many donors whose sd is constant 0 set to 0

Z <- do.call(cbind, lapply(GROUP_LEVELS, function(g) {
  idx <- which(donor_grp == g)
  if (length(idx) == 1) Xz[, idx] else rowMeans(Xz[, idx, drop = FALSE])
}))
colnames(Z) <- GROUP_LEVELS
rownames(Z) <- genes149

donor_n <- table(donor_grp)[GROUP_LEVELS]
donor_n[is.na(donor_n)] <- 0
donor_meta <- data.frame(donor = donors, group = donor_grp,
                         n_mono_cells = as.integer(table(donor_vec)[donors]),
                         stringsAsFactors = FALSE)

## =========================================================================
## Step 3: row annotation (single-cell hdWGCNA module, links to Fig5)
## =========================================================================
mod_file <- file.path(out_dir, "path13_sc_149gene_module.csv")
row_mod <- rep(NA_character_, length(genes149)); names(row_mod) <- genes149
if (file.exists(mod_file)) {
  mod <- read.csv(mod_file, stringsAsFactors = FALSE)
  hit <- match(genes149, mod$gene_name)
  row_mod[!is.na(hit)] <- mod$module[hit[!is.na(hit)]]
  row_mod[!is.na(row_mod) & row_mod == "grey"] <- NA   # grey = unassigned, not colored
}
cat("Genes with module annotation:", sum(!is.na(row_mod)), "/", length(genes149), "\n")

## =========================================================================
## Step 4: row ordering = first by module (annotated ones like ISG first, grouped by module) then cluster within group
## =========================================================================
ord_factor <- ifelse(is.na(row_mod), "zzz_unassigned", row_mod)
Z_sorted <- NULL
for (mm in sort(unique(ord_factor))) {
  idx <- which(ord_factor == mm)
  if (length(idx) >= 3) {
    sub_hc <- hclust(dist(Z[idx, , drop = FALSE]), method = "ward.D2")
    idx <- idx[sub_hc$order]
  }
  Z_sorted <- rbind(Z_sorted, Z[idx, , drop = FALSE])
}
if (is.null(Z_sorted)) stop("heatmap matrix is empty")
Z <- Z_sorted
row_mod_sorted <- row_mod[rownames(Z)]
genes_sorted   <- rownames(Z)

## =========================================================================
## Step 4b: save pre-plot data — for local re-render/debug (can draw without server object)
##
## Save an RDS that is "exactly the plotting input": sorted z matrix + row module annotation + palette +
## donors per group. Locally, just source the draw_heatmap logic from Step 5
## (or change to ggplot/ComplexHeatmap yourself); no longer depends on Seurat object or edgeR.
plotdata <- list(
  Z              = Z,                 # 149 x 5, sorted, rows = genes, cols = five groups
  row_module     = row_mod_sorted,    # hdWGCNA module corresponding to Z rows (NA = unassigned)
  gene_label     = intersect(c("IFI44L","IFI44","IFI6","IFIT3","GBP1","PSMB8",
                               "UBE2L6","LAP3","FCGR1A","C1QA","TREM1","CLEC4E",
                               "STAT3","SERPINA1","PDK4","S100A8","S100A9",
                               "MT2A","HP","SOCS3"), rownames(Z)),  # genes to label on plot
  GROUP_LEVELS   = GROUP_LEVELS,
  GROUP_COLS     = GROUP_COLS,
  donor_n        = donor_n,           # donors per group
  zlim           = 3,                 # color clipping range (consistent with plot)
  palette        = c("#2166AC", "#F7F7F7", "#B2182B")  # blue - white - red
)
saveRDS(plotdata, file.path(out_dir, "path40_fig3c_plotdata.rds"))
cat("Saved plotting data: path40_fig3c_plotdata.rds\n")
## Minimal local re-render example (commented out; uncomment when needed):
##   pd <- readRDS("path40_fig3c_plotdata.rds")
##   Z <- pd$Z; row_mod_sorted <- pd$row_module; donor_n <- pd$donor_n
##   GROUP_LEVELS <- pd$GROUP_LEVELS; GROUP_COLS <- pd$GROUP_COLS
##   # then copy the draw_heatmap() + pdf()/png() from Step 5 of this script

## =========================================================================
## Step 5: plotting (base graphics, no ComplexHeatmap/pheatmap dependency)
## =========================================================================
pdf_out <- file.path(out_dir, "Fig3c_149_heatmap_five_groups.pdf")
png_out <- file.path(out_dir, "Fig3c_149_heatmap_five_groups.png")

draw_heatmap <- function() {
  par(mar = c(5.2, 10, 2.2, 6.5), xpd = NA)
  zlim <- 3                                   # clip to [-3,3], white = 0
  Zc <- pmax(pmin(Z, zlim), -zlim)
  cols <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
  brk <- seq(-zlim, zlim, length.out = 102)
  nr <- nrow(Zc); nc <- ncol(Zc)
  vals <- as.vector(t(Zc))                    # row-major order, matches the i(row) j(col) double loop below
  ci  <- as.integer(cut(vals, breaks = brk, labels = FALSE))
  ci[is.na(ci)] <- 1L

  plot(NA, xlim = c(-6, nc + 1.6), ylim = c(-2.2, nr + 2), axes = FALSE,
       xlab = "", ylab = "")
  # main body
  for (i in seq_len(nr)) for (j in seq_len(nc))
    rect(j - 0.5, i - 0.5, j + 0.5, i + 0.5,
         col = cols[ci[(i - 1) * nc + j]], border = NA)
  # 5-group column header color blocks (top)
  for (j in seq_len(nc))
    rect(j - 0.48, nr + 0.28, j + 0.48, nr + 0.72,
         col = GROUP_COLS[GROUP_LEVELS[j]], border = NA)
  text(seq_len(nc), nr + 1.6, labels = paste0(GROUP_LEVELS, "  n=", donor_n[GROUP_LEVELS]),
       cex = 0.62, col = "grey10")
  # row gene names (only representative genes labeled) + module color bar (color bar flush to left edge of main body)
  sentinel <- c("IFI44L","IFI44","IFI6","IFIT3","GBP1","PSMB8","UBE2L6","LAP3",
                "FCGR1A","C1QA","TREM1","CLEC4E","STAT3","SERPINA1","PDK4",
                "S100A8","S100A9","MT2A","HP","SOCS3")
  lab <- intersect(sentinel, rownames(Z))
  pos <- match(lab, rownames(Z))
  text(-0.35, pos, labels = lab, adj = 1, cex = 0.42, col = "grey15")
  mod_col <- c(purple = "#800080", turquoise = "#00CCCC", red = "#FF0000",
               pink = "#FF69B4", magenta = "#FF00FF", tan = "#D2B48C",
               yellow = "#FFFF00", green = "#00CC00", blue = "#0000FF",
               brown = "#A52A2A", black = "#000000")
  for (i in seq_len(nr)) {
    m <- row_mod_sorted[i]
    rect(0.02, i - 0.5, 0.30, i + 0.5,
         col = if (!is.na(m) && m %in% names(mod_col)) mod_col[m] else "grey93",
         border = NA)
  }
  # color legend (right)
  lx0 <- nc + 1.35; lx1 <- lx0 + 0.4; ly0 <- nr * 0.25; ly1 <- nr * 0.75
  lseq <- seq(ly0, ly1, length.out = 101)
  for (k in seq_along(cols))
    rect(lx0, lseq[k], lx1, lseq[k] + (ly1 - ly0) / 101, col = cols[k], border = NA)
  text(lx1 + 0.18, ly1, labels = "3", cex = 0.55, adj = 0)
  text(lx1 + 0.18, ly0, labels = "-3", cex = 0.55, adj = 0)
  text(lx1 + 0.18, (ly0 + ly1) / 2, labels = "0", cex = 0.55, adj = 0)
  text((lx0 + lx1) / 2, nr + 1.6, labels = "row z", cex = 0.6, col = "grey30")
}

pdf(pdf_out, width = 7.2, height = 8.8)
draw_heatmap()
dev.off()
png(png_out, width = 7.2, height = 8.8, units = "in", res = 300)
draw_heatmap()
dev.off()
cat("Saved: Fig3c_149_heatmap_five_groups.{pdf,png}\n")

## =========================================================================
## Step 6: save data (for verification/assembly)
## =========================================================================
Z_df <- data.frame(gene = rownames(Z), module = row_mod_sorted, Z,
                   stringsAsFactors = FALSE, check.names = FALSE)
write.csv(Z_df, file.path(out_dir, "path40_149_group_z.csv"), row.names = FALSE)
saveRDS(list(logcpm = logcpm[genes149, , drop = FALSE], donor_grp = donor_grp),
        file.path(out_dir, "path40_149_donor_logcpm.rds"))
write.csv(donor_meta, file.path(out_dir, "path40_donor_meta.csv"), row.names = FALSE)
cat("Saved: path40_149_group_z.csv / path40_149_donor_logcpm.rds / path40_donor_meta.csv\n")

cat("\n===== Script 40 finished =====\n")
cat("To annotate row sources (e.g. manually annotated functional programs), just edit the row_mod section;\n")
cat("To reverse the 5-group order, change the GROUP_LEVELS order (colors update accordingly).\n")
