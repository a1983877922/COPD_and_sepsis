###############################################################################
# Sepsis + COPD comorbidity analysis — Script 37: lung-side (GSE136831) donor batch QC figure
# File: 37_lung_donor_batch_figure.R
#
# Purpose: answer a likely reviewer question — "is there any batch/donor effect
#   among the 46 lung-side donors?"
#   Path 1 (blood-blood merge) has path1_batch_before_after.pdf proving Harmony
#   works, but the lung side of Path 2 (GSE136831) never had any batch/donor
#   visualization — an explicit weak point of the paper (the lung myeloid DEGs
#   are half the source of the 149 shared programs, yet from a single dataset).
#   This script supplies that evidence:
#     A. Lung myeloid UMAP (colored by disease)
#     B. Lung myeloid UMAP (colored by donor, 46 donors)
#     C. Donor-level pseudobulk PCA (one point per donor) + PERMANOVA disease effect
#     D. Variance partition: pseudo-R2 of donor vs disease in cell-level PC space (same convention as script 25)
#
# Data: GSE136831 (Adams 2020 human lung atlas), only Control + COPD (IPF excluded),
#       then CellType_Category == "Myeloid" (115,557 cells, 46 donors).
#
# Memory strategy (important):
#   The raw mtx is 45947 x 312928 with 693 million nonzero entries; reading it
#   whole would need ~8-11 GB.
#   To run on a 32 GB machine, this script does NOT read it whole; instead it
#   streams through the mtx once, accomplishing two things at once:
#     (1) donor-level pseudobulk: all 115,557 myeloid cells summed per 46 donors
#     (2) cell-level matrix: proportion-stratified sampling of N_SUB cells per
#         donor (default 50,000) for PCA/UMAP
#   Measured scan speed ~3.75M rows/sec; the whole file takes ~3 minutes.
#
# Usage:
#   Server (Linux, locale = UTF-8):
#     Rscript 37_lung_donor_batch_figure.R
#   Windows local (non-UTF-8 locale; non-ASCII paths must be passed explicitly with encoding marked):
#     Rscript 37_lung_donor_batch_figure.R <GSE136831 dir> <output dir>
#
# Output:
#   path2_lung_UMAP_donor_disease.pdf / .png   four-panel figure
#   path37_lung_myeloid_umap.csv       cell-level UMAP coordinates (sampled cells)
#   path37_lung_donor_pca.csv          donor-level PCA coordinates + disease + cell counts
#   path37_lung_donor_cells.csv        myeloid cell count per donor
#   path37_lung_pc_variance.csv        donor/disease R2 per PC
#   path37_lung_variance_partition.txt variance partition + PERMANOVA results (can go straight into Methods/Results)
#   path37_lung_donor_logcpm.rds       donor-level pseudobulk log2CPM matrix (for recomputation)
#   path37_lung_myeloid_cache.rds      cache of streamed extraction results (skips the scan on reruns)
###############################################################################

args <- commandArgs(trailingOnly = TRUE)

## ---- Auto-locate script directory (supports command-line args; compatible with Windows non-UTF-8 locale) ----
.this_file <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
if (.this_file == "" || !file.exists(.this_file)) {
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1])
}
SCRIPT_DIR <- dirname(.this_file)

# When the Windows locale is not UTF-8, R cannot correctly parse Chinese string
# literals inside the script, and list.dirs() silently fails to discover non-ASCII
# directories (verified empirically). So when running locally, paths are passed on
# the command line and explicitly marked UTF-8 (R then accesses the filesystem
# via wide-char APIs).
mark_utf8 <- function(x) { if (!is.null(x)) Encoding(x) <- "UTF-8"; x }
if (length(args) >= 1) GSE_DIR_ARG <- mark_utf8(args[1]) else GSE_DIR_ARG <- NULL
if (length(args) >= 2) SCRIPT_DIR   <- mark_utf8(args[2])

config_file <- file.path(SCRIPT_DIR, "00_config.R")
if (!file.exists(config_file)) stop("cannot find 00_config.R: ", config_file)
source(config_file)

# out_dir is set by 00_config.R from SCRIPT_DIR; mark encoding again on Windows
out_dir <- mark_utf8(out_dir)
P <- function(...) mark_utf8(file.path(...))   # uniformly build UTF-8 output paths

gse_dir <- if (!is.null(GSE_DIR_ARG)) GSE_DIR_ARG else P(copd_dir, "GSE136831")
if (!dir.exists(gse_dir)) stop("GSE136831 dir not found: ", gse_dir)

cat("\n==============================================================\n")
cat("Script 37: lung-side (GSE136831) donor batch QC figure\n")
cat("==============================================================\n")
cat("GSE136831 dir:", gse_dir, "\n")
cat("Output dir    :", out_dir, "\n")

library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
suppressWarnings(suppressMessages(library(dplyr)))

## ========================== Tunable parameters ==========================
N_SUB   <- 50000    # cells sampled for cell-level PCA/UMAP (stratified by donor proportion)
NPC     <- 30       # number of PCs for cell-level PCA
NPC_VAR <- 20       # number of PCs used in variance partitioning
N_TOP   <- 2000     # number of highly variable genes for donor-level pseudobulk
NPERM   <- 999      # PERMANOVA permutations
SEED    <- 123
MIN_PER_DONOR <- 100  # minimum cells kept per donor when sampling
## ==========================================================================

cache_file <- P(out_dir, "path37_lung_myeloid_cache.rds")

## =========================================================================
## Step 1: metadata + donor mapping (lightweight, no matrix read)
## =========================================================================
cat("\n===== Step 1: read metadata =====\n")

meta_file <- P(gse_dir, "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz")
bc_file   <- P(gse_dir, "GSE136831_AllCells.cellBarcodes.txt.gz")
gene_file <- P(gse_dir, "GSE136831_AllCells.GeneIDs.txt.gz")
mtx_file  <- P(gse_dir, "GSE136831_RawCounts_Sparse.mtx.gz")

meta <- read.delim(gzfile(meta_file), stringsAsFactors = FALSE, check.names = FALSE)
barcodes <- readLines(gzfile(bc_file))
genes    <- read.delim(gzfile(gene_file), header = TRUE, stringsAsFactors = FALSE)

cat("metadata rows:", nrow(meta), " barcodes:", length(barcodes), " genes:", nrow(genes), "\n")
if (nrow(meta) != length(barcodes)) stop("metadata row count and barcode count disagree; cannot align")

# mtx column order = barcode file line order; metadata aligns on CellBarcode_Identity
if (!identical(barcodes, meta$CellBarcode_Identity)) {
  cat("metadata order disagrees with barcodes; reordering by barcode\n")
  meta <- meta[match(barcodes, meta$CellBarcode_Identity), ]
}
gene_sym <- make.unique(as.character(genes[[2]]))
G <- nrow(genes)

keep_lung <- meta$Disease_Identity %in% c("Control", "COPD") &
             meta$CellType_Category == "Myeloid"
donor_all <- as.character(meta$Subject_Identity)
donors    <- sort(unique(donor_all[keep_lung]))
nd        <- length(donors)

# Disease and cell counts per donor (all myeloid cells) — authoritative convention for donor-level analysis
donor_full <- as.data.frame(table(Subject  = donor_all[keep_lung],
                                  Disease  = meta$Disease_Identity[keep_lung]),
                            stringsAsFactors = FALSE)
donor_full <- donor_full[donor_full$Freq > 0, , drop = FALSE]
names(donor_full) <- c("donor", "disease", "n_cells")
donor_full$disease <- factor(donor_full$disease, levels = c("Control", "COPD"))

cat(sprintf("lung myeloid cells: %d (Control+COPD, Myeloid), donors: %d\n", sum(keep_lung), nd))

## =========================================================================
## Step 2: stream the mtx — produce donor pseudobulk and sampled cell-level matrix at once
## =========================================================================
if (file.exists(cache_file)) {
  cat("\n[cache] cache found, skipping streaming scan:", basename(cache_file), "\n")
  cache <- readRDS(cache_file)
  mat_sub <- cache$mat_sub; meta_sub <- cache$meta_sub
  donor_counts <- cache$donor_counts; gene_sub <- cache$gene_sub
} else {
  cat("\n===== Step 2: streaming scan of the sparse matrix (~3-8 minutes) =====\n")

  # Donor index per cell (0 = not included)
  donor_idx <- integer(length(barcodes))
  donor_idx[keep_lung] <- match(donor_all[keep_lung], donors)

  # Stratified sampling by donor proportion, at least MIN_PER_DONOR cells per donor
  set.seed(SEED)
  mye_pos <- which(keep_lung)
  alloc <- table(donor_all[mye_pos])
  n_target <- pmin(as.integer(alloc),
                   pmax(MIN_PER_DONOR, round(N_SUB * as.integer(alloc) / length(mye_pos))))
  names(n_target) <- names(alloc)          # fix: as.integer drops names, restore them
  sub_cells <- unlist(lapply(names(n_target), function(d) {
    p <- mye_pos[donor_all[mye_pos] == d]
    sample(p, min(n_target[[d]], length(p)))
  }))
  sub_idx <- integer(length(barcodes))
  sub_idx[sub_cells] <- seq_along(sub_cells)
  cat(sprintf("sampled cells: %d (from %d donors)\n", length(sub_cells), length(unique(donor_all[sub_cells]))))

  con <- gzcon(file(mtx_file, "rb"))
  hdr <- readLines(con, n = 1)
  repeat {
    l <- readLines(con, n = 1)
    if (!startsWith(l, "%") && nzchar(l)) break
  }
  dims <- as.numeric(strsplit(trimws(l), "\\s+")[[1]])
  cat("mtx dims:", dims[1], "x", dims[2], " nnz =", dims[3], "\n")
  stopifnot(dims[1] == G, dims[2] == length(barcodes))

  CHUNK <- 5e6
  cap   <- 170e6                     # preallocated (~110 million nonzero entries for the sampled cells)
  gi <- integer(cap); ci <- integer(cap); xv <- integer(cap); k <- 0L
  accum <- NULL
  nchunk <- 0L; t0 <- Sys.time()

  repeat {
    v <- scan(con, what = integer(), nmax = 3 * CHUNK, quiet = TRUE)
    if (length(v) < 3) break
    n <- length(v) %/% 3
    mm <- matrix(v[seq_len(3 * n)], ncol = 3, byrow = TRUE)
    g <- mm[, 1]; cc <- mm[, 2]; val <- mm[, 3]
    rm(mm, v)

    di <- donor_idx[cc]
    sel <- di > 0L
    if (any(sel)) {
      s <- sparseMatrix(i = g[sel], j = di[sel], x = as.numeric(val[sel]),
                        dims = c(G, nd), giveCsparse = TRUE)
      accum <- if (is.null(accum)) s else accum + s
    }

    sj <- sub_idx[cc]
    sel2 <- sj > 0L
    if (any(sel2)) {
      nn <- sum(sel2)
      if (k + nn > cap) {                       # grow capacity 1.5x if insufficient
        cap2 <- ceiling((k + nn) * 1.5)
        gi <- c(gi, integer(cap2 - length(gi)))
        ci <- c(ci, integer(cap2 - length(ci)))
        xv <- c(xv, integer(cap2 - length(xv)))
        cap <- cap2
        cat(sprintf("  [grow] cap -> %.0fM\n", cap / 1e6))
      }
      idx <- (k + 1L):(k + nn)
      gi[idx] <- g[sel2]; ci[idx] <- sj[sel2]; xv[idx] <- val[sel2]
      k <- k + nn
    }
    nchunk <- nchunk + 1L
    if (nchunk %% 20 == 0)
      cat(sprintf("  chunk %d: %.1fM rows accumulated, %.1fs\n", nchunk, nchunk * CHUNK / 1e6,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    if (n < CHUNK) break
  }
  close(con)
  cat(sprintf("scan done: %d chunks, %.1f minutes\n", nchunk,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  donor_counts <- as.matrix(accum)
  colnames(donor_counts) <- donors
  rownames(donor_counts) <- gene_sym
  cat("donor pseudobulk matrix:", dim(donor_counts)[1], "x", dim(donor_counts)[2], "\n")

  mat_sub <- sparseMatrix(i = gi[seq_len(k)], j = ci[seq_len(k)],
                          x = as.numeric(xv[seq_len(k)]),
                          dims = c(G, length(sub_cells)), giveCsparse = TRUE)
  rm(gi, ci, xv); gc()
  rownames(mat_sub) <- gene_sym
  colnames(mat_sub) <- barcodes[sub_cells]

  meta_sub <- data.frame(
    cell       = barcodes[sub_cells],
    donor      = donor_all[sub_cells],
    disease    = factor(meta$Disease_Identity[sub_cells], levels = c("Control", "COPD")),
    cell_type  = meta$Manuscript_Identity[sub_cells],
    stringsAsFactors = FALSE
  )
  gene_sub <- gene_sym

  saveRDS(list(mat_sub = mat_sub, meta_sub = meta_sub,
               donor_counts = donor_counts, gene_sub = gene_sub), cache_file)
  cat("[cache] saved:", basename(cache_file), "\n")
}

## =========================================================================
## Step 3: cell-level QC + normalization + PCA + UMAP
## =========================================================================
cat("\n===== Step 3: cell-level PCA + UMAP =====\n")

# Light QC (same as script 02, path2 cross-tissue: nFeature > 200 & percent.mt < 20)
mt_genes <- grep("^MT-", gene_sub)
tot_cnt  <- Matrix::colSums(mat_sub)
mt_pct   <- if (length(mt_genes) > 0)
  Matrix::colSums(mat_sub[mt_genes, , drop = FALSE]) / tot_cnt * 100 else rep(0, ncol(mat_sub))
nfeat    <- Matrix::colSums(mat_sub > 0)
qc_keep  <- nfeat > 200 & mt_pct < 20 & tot_cnt > 0
cat(sprintf("QC: kept %d / %d sampled cells (nFeature>200 & percent.mt<20)\n",
            sum(qc_keep), length(qc_keep)))
mat_sub  <- mat_sub[, qc_keep, drop = FALSE]
meta_sub <- meta_sub[qc_keep, , drop = FALSE]
meta_sub$percent.mt <- mt_pct[qc_keep]
rownames(meta_sub) <- meta_sub$cell

seu <- CreateSeuratObject(counts = mat_sub, meta.data = meta_sub,
                          project = "GSE136831_lung_myeloid",
                          min.cells = 3, min.features = 0)
rm(mat_sub); gc()
cat("Seurat object:", ncol(seu), "cells,", nrow(seu), "genes\n")

set.seed(SEED)
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = 3000, verbose = FALSE)
seu <- ScaleData(seu, vars.to.regress = "percent.mt", verbose = FALSE)
seu <- RunPCA(seu, npcs = NPC, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:min(30, NPC), n.neighbors = 30, verbose = FALSE)
cat("PCA + UMAP done\n")

umap_df <- data.frame(
  cell      = colnames(seu),
  donor     = seu$donor,
  disease   = seu$disease,
  cell_type = seu$cell_type,
  UMAP_1    = Embeddings(seu, "umap")[, 1],
  UMAP_2    = Embeddings(seu, "umap")[, 2],
  stringsAsFactors = FALSE
)
write.csv(umap_df, P(out_dir, "path37_lung_myeloid_umap.csv"), row.names = FALSE)
cat("saved: path37_lung_myeloid_umap.csv\n")

donor_cells <- donor_full          # per-donor counts use the "all myeloid cells" convention (not the sample)
write.csv(as.data.frame(donor_cells), P(out_dir, "path37_lung_donor_cells.csv"), row.names = FALSE)
cat("saved: path37_lung_donor_cells.csv\n")

## =========================================================================
## Step 4: donor-level pseudobulk PCA + PERMANOVA
## =========================================================================
cat("\n===== Step 4: donor-level pseudobulk PCA + PERMANOVA =====\n")

# CPM normalization (implemented here, no edgeR dependency; matches edgeR::cpm)
lib <- Matrix::colSums(donor_counts)
cpm <- sweep(donor_counts, 2, lib, "/") * 1e6
keep_g <- Matrix::rowSums(cpm > 1) >= 3
logcpm <- log2(cpm[keep_g, , drop = FALSE] + 1)
cat("donor pseudobulk: kept", nrow(logcpm), "genes (CPM>1 in >=3 donors)\n")

gv <- apply(logcpm, 1, var)
top_g <- names(sort(gv, decreasing = TRUE))[seq_len(min(N_TOP, length(gv)))]
X <- t(logcpm[top_g, , drop = FALSE])        # 46 x 2000
Xs <- scale(X)

grp <- factor(donor_full$disease[match(donors, donor_full$donor)],
              levels = c("Control", "COPD"))
names(grp) <- donors
cat("donor groups:"); print(table(grp))

set.seed(SEED)
pca_d <- prcomp(Xs, center = TRUE, scale. = FALSE)
pve <- summary(pca_d)$importance[2, seq_len(5)] * 100

# Hand-written PERMANOVA (no vegan dependency): Euclidean distance, 999 permutations
permanova <- function(Xm, g, nperm = NPERM, seed = SEED) {
  D <- as.matrix(dist(Xm)); n <- nrow(D); k <- nlevels(g)
  SST <- sum(D^2) / n
  ssw_of <- function(gg) {
    s <- 0
    for (lev in levels(gg)) {
      ii <- which(gg == lev)
      if (length(ii) > 1) s <- s + sum(D[ii, ii, drop = FALSE]^2) / length(ii)
    }
    s
  }
  ssw <- ssw_of(g); ssb <- SST - ssw
  Fobs <- (ssb / (k - 1)) / (ssw / (n - k))
  set.seed(seed)
  Fp <- replicate(nperm, {
    g2 <- sample(g); s2 <- ssw_of(g2)
    (SST - s2) / (k - 1) / (s2 / (n - k))
  })
  list(R2 = ssb / SST, F = Fobs, p = (1 + sum(Fp >= Fobs)) / (1 + nperm))
}
pm <- permanova(Xs, grp)
cat(sprintf("PERMANOVA (donor-level, %d donors): R2(disease) = %.4f, F = %.3f, p = %.4f\n",
            nrow(Xs), pm$R2, pm$F, pm$p))

pc1_p <- tryCatch(wilcox.test(pca_d$x[, 1] ~ grp)$p.value, error = function(e) NA)
pc2_p <- tryCatch(wilcox.test(pca_d$x[, 2] ~ grp)$p.value, error = function(e) NA)
cat(sprintf("PC1 disease difference Wilcoxon p = %.4g | PC2 p = %.4g\n", pc1_p, pc2_p))

donor_pca <- data.frame(
  donor    = donors,
  disease  = grp,
  n_cells  = as.integer(donor_cells$n_cells[match(donors, donor_cells$donor)]),
  PC1 = pca_d$x[, 1], PC2 = pca_d$x[, 2], PC3 = pca_d$x[, 3],
  stringsAsFactors = FALSE
)
write.csv(donor_pca, P(out_dir, "path37_lung_donor_pca.csv"), row.names = FALSE)
saveRDS(logcpm, P(out_dir, "path37_lung_donor_logcpm.rds"))
cat("saved: path37_lung_donor_pca.csv / path37_lung_donor_logcpm.rds\n")

## =========================================================================
## Step 5: variance partition (cell-level PC space, same convention as script 25)
## =========================================================================
cat("\n===== Step 5: variance partition (donor vs disease) =====\n")

emb  <- Embeddings(seu, "pca")[, seq_len(NPC_VAR), drop = FALSE]
sdev <- seu@reductions$pca@stdev
w_all <- sdev^2 / sum(sdev^2)
w <- w_all[seq_len(NPC_VAR)] / sum(w_all[seq_len(NPC_VAR)])

pc_tab <- do.call(rbind, lapply(seq_len(NPC_VAR), function(j) {
  d <- data.frame(PC = j, x = emb[, j], disease = seu$disease, donor = factor(seu$donor))
  data.frame(PC = j, var_prop = w_all[j],
             R2_disease = summary(lm(x ~ disease, d))$r.squared,
             R2_donor   = summary(lm(x ~ donor,   d))$r.squared)
}))
write.csv(pc_tab, P(out_dir, "path37_lung_pc_variance.csv"), row.names = FALSE)

R2_disease <- sum(w * pc_tab$R2_disease)
R2_donor   <- sum(w * pc_tab$R2_donor)
ratio      <- R2_disease / R2_donor
cat(sprintf("pseudo-R2 (weighted by PC variance, top %d PCs): donor = %.4f | disease = %.4f\n",
            NPC_VAR, R2_donor, R2_disease))
cat(sprintf("→ disease explains %.1f%% of between-donor variation\n", 100 * ratio))

vp_lines <- c(
  "===== Lung-side (GSE136831) myeloid variance partition =====",
  sprintf("Cells: %d sampled cells (stratified by donor proportion, %d myeloid cells in total), %d donors (Control %d / COPD %d)",
          ncol(seu), sum(keep_lung), nd, sum(grp == "Control"), sum(grp == "COPD")),
  sprintf("Cell-level PC-space pseudo-R2 (top %d PCs, variance-weighted):", NPC_VAR),
  sprintf("  donor   pseudo-R2 = %.4f", R2_donor),
  sprintf("  disease pseudo-R2 = %.4f", R2_disease),
  sprintf("  disease/donor explanatory ratio = %.4f  (disease accounts for %.1f%% of between-donor variation)", ratio, 100 * ratio),
  sprintf("  within-donor (cell-level) residual = %.4f", 1 - R2_donor),
  "",
  "Donor-level pseudobulk PERMANOVA (Euclidean, 999 permutations):",
  sprintf("  R2(disease) = %.4f, F = %.3f, p = %.4f", pm$R2, pm$F, pm$p),
  sprintf("  PC1 disease difference Wilcoxon p = %.4g | PC2 p = %.4g", pc1_p, pc2_p),
  "",
  "Interpretation:",
  "  1) donor pseudo-R2 captures how well \"donor identity\" explains cell transcriptional state,",
  "     including disease effect plus individual/technical differences;",
  "  2) disease pseudo-R2 is the part of it explained by COPD/Control; their ratio is the disease",
  "     share of between-donor variation;",
  "  3) if that share is low (e.g. <25%), lung myeloid variation is dominated by donor individual",
  "     differences — exactly why this paper switched to \"donor-level pseudobulk + edgeR\"",
  "     (scripts 28/28b) instead of cell-level DEGs."
)
writeLines(vp_lines, P(out_dir, "path37_lung_variance_partition.txt"))
cat("saved: path37_lung_variance_partition.txt\n")

## =========================================================================
## Step 6: figures
## =========================================================================
cat("\n===== Step 6: figures =====\n")

col_dis <- c(Control = "#4DBBD5", COPD = "#E64B35")
donors_use <- sort(unique(umap_df$donor))
pal_donor <- setNames(scales::hue_pal()(length(donors_use)), donors_use)

pt <- 0.35
pA <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = disease)) +
  geom_point(size = pt, alpha = 0.6, stroke = 0) +
  scale_colour_manual(values = col_dis, name = "Disease") +
  labs(title = "A  Lung myeloid UMAP — by disease",
       subtitle = sprintf("GSE136831, %d myeloid cells (Control+COPD), %d donors",
                          nrow(umap_df), length(donors_use)),
       x = "UMAP 1", y = "UMAP 2") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_bw(base_size = 11) + theme(panel.grid = element_blank())

pB <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = donor)) +
  geom_point(size = pt, alpha = 0.6, stroke = 0) +
  scale_colour_manual(values = pal_donor, name = "Donor") +
  labs(title = "B  Lung myeloid UMAP — by donor",
       subtitle = "each colour = one donor lung (46 donors)",
       x = "UMAP 1", y = "UMAP 2") +
  guides(colour = guide_legend(ncol = 3, override.aes = list(size = 2, alpha = 1),
                               keywidth = 0.4, keyheight = 0.4)) +
  theme_bw(base_size = 11) + theme(panel.grid = element_blank(),
                                   legend.key.size = unit(0.28, "cm"),
                                   legend.text = element_text(size = 6))

pC <- ggplot(donor_pca, aes(PC1, PC2, colour = disease, size = n_cells)) +
  geom_point(alpha = 0.85) +
  scale_colour_manual(values = col_dis, name = "Disease") +
  scale_size_continuous(range = c(2, 7), name = "Myeloid cells\n(all, n per donor)") +
  stat_ellipse(aes(group = disease), level = 0.68, linewidth = 0.5, inherit.aes = FALSE,
               data = donor_pca, mapping = aes(PC1, PC2, colour = disease)) +
  labs(title = "C  Donor-level pseudobulk PCA",
       subtitle = sprintf("one point = one donor (n=%d), log2CPM top %d variable genes",
                          nrow(donor_pca), length(top_g)),
       x = sprintf("PC1 (%.1f%%)", pve[1]), y = sprintf("PC2 (%.1f%%)", pve[2])) +
  { if (requireNamespace("ggrepel", quietly = TRUE))
      ggrepel::geom_text_repel(aes(label = donor), size = 2.4, colour = "grey25",
                               max.overlaps = 60, seed = SEED, show.legend = FALSE)
    else geom_text(aes(label = donor), size = 2.4, colour = "grey25",
                   vjust = -0.9, show.legend = FALSE) } +
  annotate("text", x = Inf, y = -Inf, hjust = 1.02, vjust = -0.6, size = 3.4, colour = "grey20",
           label = sprintf("PERMANOVA: R2(disease) = %.3f, p = %.3f\nPC1 Wilcoxon p = %.3g",
                           pm$R2, pm$p, pc1_p)) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_bw(base_size = 11) + theme(panel.grid = element_blank())

vp_df <- data.frame(
  comp = factor(c("Disease (COPD vs Control)", "Donor-specific (residual)", "Within-donor (cell-level)"),
                levels = c("Disease (COPD vs Control)", "Donor-specific (residual)", "Within-donor (cell-level)")),
  value = c(R2_disease, R2_donor - R2_disease, 1 - R2_donor)
)
pD <- ggplot(vp_df, aes(x = 1, y = value, fill = comp)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = sprintf("%.3f", value), colour = comp),
            position = position_stack(vjust = 0.5), size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Disease (COPD vs Control)" = "#E64B35",
                               "Donor-specific (residual)"  = "#F39B7F",
                               "Within-donor (cell-level)"  = "#B2B2B2"), name = NULL) +
  scale_colour_manual(values = c("Disease (COPD vs Control)" = "white",
                                 "Donor-specific (residual)"  = "white",
                                 "Within-donor (cell-level)"  = "grey25"), guide = "none") +
  coord_flip() + labs(
    title = "D  Variance partition (cell-level PC space)",
    subtitle = sprintf("disease explains %.1f%% of between-donor variance", 100 * ratio),
    x = NULL, y = "Proportion of variance") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), legend.position = "bottom")

lay <- "
AABB
CCDD
"
pfig <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "GSE136831 lung myeloid: donor / disease structure and variance partition",
    subtitle = sprintf("Path 2 (cross-tissue) control figure | %d donors (Control %d / COPD %d), %d myeloid cells subsampled from %d",
                       nd, sum(grp == "Control"), sum(grp == "COPD"), nrow(umap_df), sum(keep_lung)),
    theme = theme(plot.title = element_text(size = 15, face = "bold"),
                  plot.subtitle = element_text(size = 10, colour = "grey30")))

fig_pdf <- P(out_dir, "path2_lung_UMAP_donor_disease.pdf")
fig_png <- P(out_dir, "path2_lung_UMAP_donor_disease.png")
ggsave(fig_pdf, plot = pfig, width = 16, height = 13)
ggsave(fig_png, plot = pfig, width = 16, height = 13, dpi = 200)
cat("saved:", basename(fig_pdf), "/", basename(fig_png), "\n")

cat("\n===== Script 37 done =====\n")
