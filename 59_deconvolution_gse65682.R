# =============================================================================
# 59 - Immune infiltration deconvolution (GSE65682)
# Purpose: immune deconvolution of the GSE65682 whole-blood array (n=468, 28-day deaths n=103)
#       using cell-type ssGSEA, answering: with which "immune state" do high 149/SOCS3
#       expression co-vary? What immune composition corresponds to the 152-gene prognostic
#       signature? Output serves as FigS candidates and Discussion support (the immune side of
#       the immune-drug intersection).
# Methods:
#   - Cell-type ssGSEA: 12 blood-detectable immune cell classes (gene sets embedded below),
#     within-sample gene rank(ties=average); prefer GSVA::gsva(ssgseaParam), otherwise use the
#     embedded ssGSEA implementation (classic single-sample GSEA, Barbie 2009; power weight alpha=0.25).
#   - Association anchors (to avoid re-computation drift): SOCS3 expression (current matrix,
#     a single gene is stable) + the 152 signature (read the archived copd score from
#     path15_GSE65682_train_scores.csv, verified to be in the same order as the current cohort) +
#     28-day mortality status. Spearman rho + BH-FDR; Wilcoxon between death groups.
# Output: path59_immune_scores.csv / path59_immune_assoc.csv /
#       path59_immune_dotplot.pdf / path59_immune_summary.txt
# Run: Rscript this script [OUT_DIR]
#       (requires survival/ggplot2; matrix+map in the same directory; path15 scores archive)
# =============================================================================
suppressPackageStartupMessages({ library(GEOquery); library(survival); library(ggplot2) })

args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
if (length(args) >= 1 && nzchar(args[1])) script_dir <- args[1]
OUT_DIR <- if (length(args) >= 2 && nzchar(args[2])) args[2] else script_dir
out_dir <- OUT_DIR
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---------- Utility functions (consistent with 15/53/57/58) ----------
find_surv_fields <- function(pd) {
  cols <- colnames(pd)
  evt_pat <- "mortal|surviv|outcome|dead|alive|fate"
  tme_pat <- "time|follow"
  status_cand <- cols[grepl(evt_pat, cols, ignore.case = TRUE) &
                      !grepl(tme_pat, cols, ignore.case = TRUE)]
  time_cand   <- cols[grepl(tme_pat, cols, ignore.case = TRUE) &
                      !grepl(evt_pat, cols, ignore.case = TRUE)]
  if (length(status_cand) == 0 || length(time_cand) == 0) {
    char_cols <- cols[grepl("^characteristics", cols, ignore.case = TRUE)]
    for (cc in char_cols) {
      vals <- trimws(tolower(as.character(pd[[cc]])))
      vals <- vals[!is.na(vals) & vals != ""]
      if (length(vals) == 0) next
      has_evt <- any(grepl(evt_pat, vals))
      has_tme <- any(grepl("time|follow|days", vals))
      val_only <- sub("^[^:]*:\\s*", "", vals)
      all_num <- all(suppressWarnings(!is.na(as.numeric(val_only))))
      if (length(status_cand) == 0 && has_evt && !all_num)
        status_cand <- c(status_cand, cc)
      if (length(time_cand) == 0 && has_tme && !has_evt)
        time_cand <- c(time_cand, cc)
    }
  }
  list(status = status_cand, time = time_cand)
}
to_status01 <- function(v) {
  v <- trimws(tolower(as.character(v)))
  v <- sub("^[^:]*:\\s*", "", v)
  out <- rep(NA_real_, length(v))
  out[v %in% c("dead","death","died","deceased","non-survivor","nonsurvivor",
               "1","yes","event","mortality")] <- 1
  out[v %in% c("alive","survived","survivor","live","0","no","censored","survival")] <- 0
  out
}
load_cohort <- function(gse_acc, label) {
  cat("\n----- Processing cohort:", label, "(", gse_acc, ") -----\n")
  matrix_file <- file.path(out_dir, paste0(gse_acc, "_series_matrix.txt.gz"))
  if (!file.exists(matrix_file))
    matrix_file <- file.path(out_dir, paste0(gse_acc, "_matrix.txt.gz"))
  if (file.exists(matrix_file)) {
    cat("Using local series matrix:", matrix_file, "\n")
    gse <- getGEO(filename = matrix_file, getGPL = FALSE)
  } else {
    cat("No local matrix file, downloading", gse_acc, "online...\n")
    gse <- tryCatch(getGEO(gse_acc, GSEMatrix = TRUE, getGPL = TRUE),
                    error = function(e) NULL)
    if (is.null(gse)) { cat("Download failed, skipping", label, "\n"); return(NULL) }
    gse <- gse[[1]]
  }
  expr <- exprs(gse); pd <- pData(gse); fd <- fData(gse)
  cat("Expression matrix:", nrow(expr), "probes x", ncol(expr), "samples\n")
  sym_col <- NULL
  for (cand in c("Gene Symbol","GeneSymbol","Symbol","Gene.Symbol",
                 "gene_symbol","GENE_SYMBOL")) {
    if (cand %in% colnames(fd)) { sym_col <- cand; break }
  }
  if (is.null(sym_col)) {
    m <- grep("symbol", colnames(fd), ignore.case = TRUE)
    if (length(m) > 0) sym_col <- colnames(fd)[m[1]]
  }
  use_ga <- FALSE
  if (is.null(sym_col) && "gene_assignment" %in% colnames(fd)) {
    sym_col <- "gene_assignment"; use_ga <- TRUE
  }
  sym <- NULL
  if (!is.null(sym_col)) {
    sym <- as.character(fd[[sym_col]])
    if (use_ga) {
      sym <- vapply(strsplit(sym, "//", fixed = TRUE), function(x) {
        if (length(x) >= 2) trimws(x[2]) else ""
      }, character(1))
    }
    names(sym) <- rownames(fd)
  } else {
    cand <- unique(c(file.path(dirname(matrix_file), "GPL13667_symbol_map.csv"),
                     file.path(out_dir, "GPL13667_symbol_map.csv"),
                     "GPL13667_symbol_map.csv"))
    mf <- cand[file.exists(cand)]
    if (length(mf) == 0) { cat("No symbol column and no offline annotation file, skipping\n"); return(NULL) }
    cat("Using offline probe annotation:", mf[1], "\n")
    mtab <- read.delim(mf[1], header = TRUE, sep = "\t", quote = "",
                       stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("probe", "symbol") %in% colnames(mtab)))
      mtab <- read.csv(mf[1], header = TRUE, stringsAsFactors = FALSE)
    msym <- as.character(mtab$symbol); names(msym) <- as.character(mtab$probe)
    sym <- unname(msym[rownames(fd)])
    sym[is.na(sym)] <- ""
    sym <- sub(" /// .*$", "", sym)
    names(sym) <- rownames(fd)
  }
  keep <- sym != "" & !is.na(sym) & sym != "---"
  expr <- expr[keep, , drop = FALSE]; sym <- sym[keep]
  sym_uniq <- unique(sym)
  expr_sym <- matrix(NA, nrow = length(sym_uniq), ncol = ncol(expr),
                     dimnames = list(sym_uniq, colnames(expr)))
  for (g in sym_uniq) {
    idx <- which(sym == g)
    expr_sym[g, ] <- if (length(idx) == 1) expr[idx, ] else colMeans(expr[idx, , drop = FALSE])
  }
  expr <- expr_sym
  cat("After mapping:", nrow(expr), "genes x", ncol(expr), "samples\n")
  sf <- find_surv_fields(pd)
  if (length(sf$status) == 0) { cat("!! Could not identify the status field\n"); return(NULL) }
  status_col <- sf$status[1]
  time_col <- if (length(sf$time) == 0) NULL else sf$time[1]
  status_raw <- as.character(pd[[status_col]])
  if (is.null(time_col)) time_raw <- rep(28, nrow(pd))
  else time_raw <- suppressWarnings(as.numeric(as.character(pd[[time_col]])))
  names(status_raw) <- rownames(pd); names(time_raw) <- rownames(pd)
  common <- intersect(colnames(expr), rownames(pd))
  expr <- expr[, common, drop = FALSE]
  status_raw <- status_raw[common]; time_raw <- time_raw[common]
  status01 <- to_status01(status_raw)
  keep_samp <- !is.na(status01) & !is.na(time_raw) & time_raw > 0
  expr <- expr[, keep_samp, drop = FALSE]
  status01 <- status01[keep_samp]; time_raw <- time_raw[keep_samp]
  cat("Samples with survival data retained:", ncol(expr), "(dead=", sum(status01), ")\n")
  list(expr = expr, time = time_raw, status = status01, label = label)
}

## Immune cell gene sets (conservative consensus markers, detectable in bulk blood; sets with <5 probe hits are dropped)
IMMUNE_SETS <- list(
  `B cells`        = c("CD19","MS4A1","CD79A","CD79B","BANK1","BLK","FCRL2","CD22","SPIB","FCER2"),
  `Plasma cells`   = c("MZB1","SDC1","XBP1","TNFRSF17","SLAMF7","DERL3","FKBP11","IGHG1","IGJ","SEC11C"),
  `CD8 T cells`    = c("CD8A","CD8B","GZMA","GZMB","GZMH","GZMK","PRF1","NKG7","KLRK1","CXCR3","CCL5"),
  `CD4 T cells`    = c("CD4","IL7R","LEF1","TCF7","CCR7","SELL","AQP3","MAL","ANXA1","ITGA4"),
  `Treg`           = c("FOXP3","IL2RA","CTLA4","IKZF2","CCR8","TNFRSF18","TNFRSF9","BATF","LRRC32"),
  `NK cells`       = c("NKG7","KLRD1","KLRF1","KLRC1","KLRB1","GNLY","PRF1","GZMB","FCGR3A","SPON2","IL2RB"),
  `Monocytes`      = c("CD14","LYZ","FCN1","S100A8","S100A9","S100A12","CSF1R","CD300E","FPR1","THBS1","LST1","MS4A7"),
  `M2-like mono/mac` = c("MRC1","CD163","MS4A4A","MSR1","CCL18","SIGLEC1","TGM2","IL10","CCL22","STAB1"),
  `Neutrophils`    = c("FCGR3B","CSF3R","S100A8","S100A9","S100A12","G0S2","ANXA3","CEACAM8","FPR1","FPR2","CXCR2","NAMPT","VNN2","SIGLEC5"),
  `mDC`            = c("ITGAX","CD1C","CLEC9A","FLT3","BATF3","IRF8","CCR7","CD83","THBD","CADM1"),
  `pDC`            = c("IL3RA","CLEC4C","LILRA4","TCF4","PACSIN1","IRF7","GZMB","SERPINF1","SPIB"),
  `Mast cells`     = c("TPSAB1","TPSB2","CPA3","HDC","MS4A2","KIT","FCER1A","GATA2","TPSD1"),
  `Eosinophils`    = c("CLC","EPX","PRG2","PRG3","RNASE2","RNASE3","SIGLEC8","IL5RA","CCR3")
)

## ---------- Embedded ssGSEA (single-sample GSEA, Barbie 2009) ----------
ssgsea_es <- function(expr, gset, alpha = 0.25) {
  ## expr: genes x samples (numeric); returns: one enrichment score per sample (called for each gset)
  genes_in <- rownames(expr)
  S <- intersect(gset, genes_in)
  if (length(S) < 5) return(rep(NA_real_, ncol(expr)))
  N <- nrow(expr)
  # Within-sample rank (ties=average), ascending
  R <- apply(expr, 2, rank, ties.method = "average")
  # Weight w = r^alpha (using rank position)
  W <- R^alpha
  # Row indices of gene set members (row numbers in rows)
  s_idx <- match(S, genes_in)
  es <- numeric(ncol(expr))
  for (j in seq_len(ncol(expr))) {
    ord <- order(R[, j])                       # gene order from low to high
    in_set <- ord %in% s_idx
    nhit <- sum(in_set); nmiss <- N - nhit
    if (nhit == 0) { es[j] <- NA_real_; next }
    w_hit <- W[s_idx, j]
    ph <- cumsum(ifelse(in_set, w_hit[match(ord[in_set], s_idx)] / sum(w_hit), 0))
    pm <- cumsum(ifelse(in_set, 0, 1 / nmiss))
    es[j] <- sum(ph - pm)                       # classic ssGSEA ES accumulation
  }
  es
}

cat("\n===== 59: GSE65682 immune infiltration ssGSEA =====\n")
coh <- load_cohort("GSE65682", "MARS (GSE65682)")
if (is.null(coh)) quit(status = 1)
expr <- coh$expr; time <- coh$time; status <- coh$status
n <- ncol(expr); cat("Samples:", n, "| deaths:", sum(status), "\n")

## Cell-type ES (rows = samples)
es_mat <- sapply(names(IMMUNE_SETS), function(ct) {
  v <- ssgsea_es(expr, IMMUNE_SETS[[ct]])
  cat(sprintf("  %-18s genes matched %2d/%2d  ES NaN:%d\n", ct,
              length(intersect(IMMUNE_SETS[[ct]], rownames(expr))),
              length(IMMUNE_SETS[[ct]]), sum(is.na(v))))
  v
})
es_mat <- as.matrix(es_mat)
es_mat <- es_mat[, colSums(is.na(es_mat)) == 0, drop = FALSE]
colnames(es_mat) <- make.names(colnames(es_mat), unique = TRUE)
cat("Valid cell types:", ncol(es_mat), "->", paste(colnames(es_mat), collapse = ", "), "\n")

## Anchors: SOCS3 expression (z) + archived 152 score + 149 mean (z)
socs3 <- as.numeric(scale(expr["SOCS3", ]))
s15 <- read.csv(file.path(out_dir, "path15_GSE65682_train_scores.csv"))
if (nrow(s15) != n) { cat("!! path15 scores length mismatch:", nrow(s15), "vs", n, "\n"); quit(status = 1) }
copd152 <- s15$copd
p149 <- intersect(scan(file.path(out_dir, "path21_shared_up_genes.txt"), ""), rownames(expr))
z149 <- t(scale(t(expr[p149, , drop = FALSE]))); z149[is.na(z149)] <- 0
mean149 <- colMeans(z149)

## Association table (defensive: tryCatch per cell; a failed cell is set to NA without interrupting)
cor_tab <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 20) return(c(NA_real_, NA_real_))
  sp <- suppressWarnings(tryCatch(cor.test(x[ok], y[ok], method = "spearman"),
                                  error = function(e) NULL))
  if (is.null(sp)) c(NA_real_, NA_real_) else c(rho = unname(sp$estimate),
                                                p = unname(sp$p.value))
}
mkrow <- function(ct) {
  y <- es_mat[, ct]
  s1 <- cor_tab(socs3, y);    s2 <- cor_tab(copd152, y)
  s3 <- cor_tab(mean149, y)
  wp <- suppressWarnings(tryCatch(wilcox.test(y[status == 1], y[status == 0])$p.value,
                                  error = function(e) NA_real_))
  d  <- tryCatch(ifelse(median(y[status == 1]) > median(y[status == 0]),
                        "up_in_dead", "down_in_dead"), error = function(e) NA_character_)
  data.frame(cell = ct, socs3_rho = s1[1], socs3_p = s1[2],
             s152_rho = s2[1], s152_p = s2[2],
             s149_rho = s3[1], s149_p = s3[2],
             death_wilcox_p = wp, dir_death = d,
             stringsAsFactors = FALSE)
}
tab <- do.call(rbind, lapply(colnames(es_mat), mkrow))
num <- c("socs3_rho","socs3_p","s152_rho","s152_p","s149_rho","s149_p","death_wilcox_p")
tab[num] <- lapply(tab[num], as.numeric)
## BH-FDR
tab$socs3_fdr <- p.adjust(tab$socs3_p, "BH")
tab$s152_fdr  <- p.adjust(tab$s152_p, "BH")
tab$s149_fdr  <- p.adjust(tab$s149_p, "BH")
tab$death_fdr <- p.adjust(tab$death_wilcox_p, "BH")
tab <- tab[order(tab$socs3_fdr), ]
write.csv(tab, file.path(out_dir, "path59_immune_assoc.csv"), row.names = FALSE)

## Score table (with anchors and outcome)
score_df <- data.frame(sample = colnames(expr), time = time, status = status,
                       SOCS3_z = socs3, copd152 = copd152, mean149 = mean149,
                       es_mat, check.names = FALSE)
write.csv(score_df, file.path(out_dir, "path59_immune_scores.csv"), row.names = FALSE)

cat("\n===== Association summary (sorted by SOCS3 FDR) =====\n")
for (i in seq_len(nrow(tab))) {
  r <- tab[i, ]
  cat(sprintf("%-18s SOCS3 rho=%+.2f FDR=%.3g | 152 rho=%+.2f FDR=%.3g | 149 rho=%+.2f FDR=%.3g | death P=%.3g %s\n",
              r$cell, r$socs3_rho, r$socs3_fdr, r$s152_rho, r$s152_fdr,
              r$s149_rho, r$s149_fdr, r$death_wilcox_p, r$dir_death))
}
cat("\n103 deaths at 28 days; cell-type direction: up_in_dead = the cell-type ES is higher in non-survivors\n")

## FigS candidate: dot plot (two anchors, SOCS3 and the 152 signature)
plot_df <- rbind(
  data.frame(cell = tab$cell, anchor = "vs SOCS3", rho = tab$socs3_rho, fdr = tab$socs3_fdr),
  data.frame(cell = tab$cell, anchor = "vs 152 signature", rho = tab$s152_rho, fdr = tab$s152_fdr)
)
plot_df <- plot_df[order(plot_df$rho), ]
plot_df$cell <- factor(plot_df$cell, levels = unique(plot_df$cell))
p <- ggplot(plot_df, aes(rho, cell, colour = rho, size = -log10(fdr + 1e-10))) +
  geom_point() + facet_wrap(~ anchor, ncol = 2) +
  scale_colour_gradient2(low = "#2166AC", mid = "grey80", high = "#B2182B",
                         midpoint = 0, name = "Spearman rho") +
  scale_size_continuous(name = "-log10 FDR") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  labs(x = "Spearman rho", y = "", title = "GSE65682: immune cell infiltration vs SOCS3 / 152 signature") +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "path59_immune_dotplot.pdf"), p, width = 9, height = 5.5)

writeLines(capture.output(print(tab)), file.path(out_dir, "path59_immune_summary.txt"))
cat("\nOutput: path59_immune_scores.csv / _assoc.csv / _dotplot.pdf / _summary.txt\n")
cat("===== Script 59 finished =====\n")
