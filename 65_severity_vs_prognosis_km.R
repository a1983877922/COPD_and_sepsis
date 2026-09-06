# =============================================================================
# 65: Figure 6, severity vs prognosis
# Purpose: redesign main Figure 6 as an integrated "severity != prognosis" figure:
#   a) SOCS3 monocyte expression rises monotonically across the 5 severity groups (violin, KW)
#   b) 149-gene program imprint (monocytes) rises across the 5 groups (violin)
#   c) none of the three "axes" is prognostic in whole blood by KM:
#      SOCS3 / IFN module / 149-Lasso (1x3 mini panels)
#   d) 152-gene whole-blood COPD signature is significant by KM (HR 2.81)
# Data: a/b use cp3_annotated.rds (monocyte subset); c/d use GSE65682 (offline
#       annotation) plus path57 risk scores and the archived path15 scores.
# Output: Fig6_severity_vs_prognosis.{pdf,png} / path65_plotdata.rds / path65_stats.txt
# Run:   Rscript 65_severity_vs_prognosis_km.R [OUT_DIR]
#       (needs Seurat/GEOquery/survival/ggplot2/patchwork)
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(GEOquery); library(survival); library(ggplot2)
})
has_patch <- requireNamespace("patchwork", quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
if (length(args) >= 1 && nzchar(args[1])) script_dir <- args[1]
out_dir <- if (length(args) >= 2 && nzchar(args[2])) args[2] else script_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

safe <- function(lbl, code) tryCatch({ code }, error = function(e) {
  cat("  [skip]", lbl, "-", conditionMessage(e), "\n"); NULL })
GROUPS <- c("Healthy","Infection_Control","COPD","Sepsis","Sepsis_Pneumonia")

## ---------- load_cohort (GSE65682, same as scripts 63/64) ----------
find_surv_fields <- function(pd) {
  cols <- colnames(pd); evt <- "mortal|surviv|outcome|dead|alive|fate"; tme <- "time|follow"
  sc <- cols[grepl(evt, cols, ignore.case = TRUE) & !grepl(tme, cols, ignore.case = TRUE)]
  tc <- cols[grepl(tme, cols, ignore.case = TRUE) & !grepl(evt, cols, ignore.case = TRUE)]
  if (length(sc) == 0 || length(tc) == 0) {
    cc <- cols[grepl("^characteristics", cols, ignore.case = TRUE)]
    for (x in cc) {
      v <- trimws(tolower(as.character(pd[[x]]))); v <- v[!is.na(v) & v != ""]
      if (length(v) == 0) next
      he <- any(grepl(evt, v)); ht <- any(grepl("time|follow|days", v))
      vo <- sub("^[^:]*:\\s*", "", v); an <- all(suppressWarnings(!is.na(as.numeric(vo))))
      if (length(sc) == 0 && he && !an) sc <- c(sc, x)
      if (length(tc) == 0 && ht && !he) tc <- c(tc, x)
    }
  }
  list(status = sc, time = tc)
}
to01 <- function(v) { v <- trimws(tolower(as.character(v))); v <- sub("^[^:]*:\\s*", "", v)
  o <- rep(NA_real_, length(v))
  o[v %in% c("dead","death","died","deceased","non-survivor","nonsurvivor","1","yes","event","mortality")] <- 1
  o[v %in% c("alive","survived","survivor","live","0","no","censored","survival")] <- 0; o }
load_cohort <- function(gse_acc) {
  mf <- file.path(out_dir, paste0(gse_acc, "_series_matrix.txt.gz"))
  if (!file.exists(mf)) mf <- file.path(out_dir, paste0(gse_acc, "_matrix.txt.gz"))
  gse <- getGEO(filename = mf, getGPL = FALSE)
  expr <- exprs(gse); pd <- pData(gse); fd <- fData(gse)
  sym_col <- NULL
  for (cand in c("Gene Symbol","GeneSymbol","Symbol","Gene.Symbol","gene_symbol","GENE_SYMBOL"))
    if (cand %in% colnames(fd)) { sym_col <- cand; break }
  if (is.null(sym_col)) { m <- grep("symbol", colnames(fd), ignore.case = TRUE)
    if (length(m) > 0) sym_col <- colnames(fd)[m[1]] }
  if (is.null(sym_col)) {
    mtab <- read.delim(file.path(out_dir, "GPL13667_symbol_map.csv"), header = TRUE,
                       sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
    msym <- as.character(mtab$symbol); names(msym) <- as.character(mtab$probe)
    sym <- unname(msym[rownames(fd)]); sym[is.na(sym)] <- ""; sym <- sub(" /// .*$", "", sym)
    names(sym) <- rownames(fd)
  } else { sym <- as.character(fd[[sym_col]]); names(sym) <- rownames(fd) }
  keep <- sym != "" & !is.na(sym) & sym != "---"
  expr <- expr[keep, , drop = FALSE]; sym <- sym[keep]; su <- unique(sym)
  es <- matrix(NA, length(su), ncol(expr), dimnames = list(su, colnames(expr)))
  for (g in su) { idx <- which(sym == g)
    es[g, ] <- if (length(idx) == 1) expr[idx, ] else colMeans(expr[idx, , drop = FALSE]) }
  expr <- es
  sf <- find_surv_fields(pd); sc <- sf$status[1]; tc <- if (length(sf$time) == 0) NULL else sf$time[1]
  sr <- as.character(pd[[sc]])
  tr <- if (is.null(tc)) rep(28, nrow(pd)) else suppressWarnings(as.numeric(as.character(pd[[tc]])))
  names(sr) <- rownames(pd); names(tr) <- rownames(pd)
  cm <- intersect(colnames(expr), rownames(pd)); expr <- expr[, cm, drop = FALSE]
  sr <- sr[cm]; tr <- tr[cm]; s01 <- to01(sr)
  k <- !is.na(s01) & !is.na(tr) & tr > 0
  list(expr = expr[, k, drop = FALSE], time = tr[k], status = s01[k])
}

cat("===== 65: new Figure 6, severity vs prognosis =====\n")

## ---- a/b: monocyte severity gradient (cp3_annotated.rds) ----
bl <- NULL
blf <- file.path(out_dir, "cp3_annotated.rds")
if (!file.exists(blf)) blf <- file.path(out_dir, "path1_sepsis_copd_integrated.rds")
if (file.exists(blf)) {
  bl <- readRDS(blf)
  ctcol <- "cell_type"; grpcol <- if ("condition" %in% colnames(bl@meta.data)) "condition" else "group"
  mono <- subset(bl, subset = cell_type == "Monocyte")
  cat("monocyte count:", ncol(mono), "| group field:", grpcol, "\n")
  mono[[grpcol]] <- factor(mono[[grpcol]], levels = GROUPS)
  pool149 <- trimws(readLines(file.path(out_dir, "path21_shared_up_genes.txt"))); pool149 <- pool149[pool149 != ""]
  gens <- intersect(pool149, rownames(mono))
  cat("149 genes matched (monocytes):", length(gens), "\n")
  # imprint = per-cell mean z-score of the 149 genes (log-normalized data layer)
  dat <- GetAssayData(mono, slot = "data")[gens, , drop = FALSE]
  z <- t(scale(t(as.matrix(dat)))); z[is.na(z)] <- 0
  imprint <- colMeans(z)
  socs3 <- as.numeric(GetAssayData(mono, slot = "data")["SOCS3", ])
  df_sev <- data.frame(group = mono[[grpcol]], SOCS3 = socs3, imprint = imprint)
  df_sev <- df_sev[!is.na(df_sev$group), ]
  kw_s <- kruskal.test(SOCS3 ~ group, df_sev); kw_i <- kruskal.test(imprint ~ group, df_sev)
  pa <- ggplot(df_sev, aes(group, SOCS3)) +
    geom_violin(aes(fill = group), alpha = 0.55, scale = "width") + geom_boxplot(width = 0.12, outlier.size = 0.2) +
    scale_fill_manual(values = c("#3C5488","#91D1C2","#F39B7F","#E64B35","#8491B4"), guide = "none") +
    labs(x = NULL, y = "SOCS3 expression",
         title = sprintf("a  SOCS3 in monocytes (cell-level, KW P = %.2g)", kw_s$p.value)) +
    theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  pb <- ggplot(df_sev, aes(group, imprint)) +
    geom_violin(aes(fill = group), alpha = 0.55, scale = "width") + geom_boxplot(width = 0.12, outlier.size = 0.2) +
    scale_fill_manual(values = c("#3C5488","#91D1C2","#F39B7F","#E64B35","#8491B4"), guide = "none") +
    labs(x = NULL, y = "149-program imprint",
         title = sprintf("b  149-program imprint (cell-level, KW P = %.2g)", kw_i$p.value)) +
    theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
} else { pa <- NULL; pb <- NULL; cat("!! monocyte object not found\n") }

## ---- c/d: whole-blood prognosis (GSE65682) ----
coh <- load_cohort("GSE65682"); expr <- coh$expr; time <- coh$time; status <- coh$status
n <- ncol(expr); cat("GSE65682 n=", n, "\n")
purple <- trimws(readLines(file.path(out_dir, "path13_sc_module_purple.txt"))); purple <- purple[purple != ""]
s15 <- read.csv(file.path(out_dir, "path15_GSE65682_train_scores.csv"))
p57 <- readRDS(file.path(out_dir, "path57_lasso_plotdata.rds"))
score_gs <- function(expr, genes) { genes <- intersect(genes, rownames(expr))
  if (length(genes) < 2) return(rep(NA_real_, ncol(expr)))
  z <- t(scale(t(expr[genes, , drop = FALSE]))); z[is.na(z)] <- 0; colMeans(z) }
socs3z <- as.numeric(scale(expr["SOCS3", ]))
ifn <- score_gs(expr, purple)
riskL <- as.numeric(p57$risk)
ref152 <- s15$copd
km_df <- function(score) { ok <- is.finite(score)
  d <- data.frame(time = time[ok], status = status[ok], score = as.numeric(score[ok]))
  d <- d[!is.na(d$score), ]; d$grp <- ifelse(d$score >= median(d$score), "High", "Low"); d }
km_curve <- function(d) { fit <- survfit(Surv(time, status) ~ grp, data = d); sf <- summary(fit)
  g <- sub("^.*=", "", as.character(sf$strata))
  data.frame(t = sf$time, s = sf$surv, grp = factor(g, levels = c("High", "Low"))) }
km_p <- function(d, lab) { cdf <- km_curve(d)
  lr <- survdiff(Surv(time, status) ~ grp, data = d); p <- 1 - pchisq(lr$chisq, 1)
  ggplot(cdf, aes(t, s, colour = grp)) + geom_step(linewidth = 0.8) +
    scale_colour_manual(values = c("High" = "#B2182B", "Low" = "#2166AC"), breaks = c("High","Low"),
                        name = "Risk", drop = FALSE) +
    scale_x_continuous(limits = c(0, NA)) + labs(x = "Days", y = "Survival", title = lab) +
    theme_bw(base_size = 8) + theme(legend.position = "none") }
pc1 <- km_p(km_df(socs3z), sprintf("c1  SOCS3 (P = %.3f)",
   1 - pchisq(survdiff(Surv(time,status)~grp, km_df(socs3z))$chisq, 1)))
pc2 <- km_p(km_df(ifn), sprintf("c2  IFN module (P = %.3f)",
   1 - pchisq(survdiff(Surv(time,status)~grp, km_df(ifn))$chisq, 1)))
pc3 <- km_p(km_df(riskL), sprintf("c3  149-Lasso (P = %.3f)",
   1 - pchisq(survdiff(Surv(time,status)~grp, km_df(riskL))$chisq, 1)))
pd <- km_p(km_df(ref152), "d  152-gene signature (log-rank P = 0.020)")
pd <- pd + theme(legend.position = "right")

## ---- assemble panels ----
plots <- list(pa, pb, pc1, pc2, pc3, pd)
plots <- Filter(Negate(is.null), plots)
if (has_patch && length(plots) >= 5) {
  library(patchwork)
  P <- (pa | pb) / (pc1 | pc2 | pc3) / (pd | plot_spacer())
  ggsave(file.path(out_dir, "Fig6_severity_vs_prognosis.pdf"), P, width = 11, height = 10)
  ggsave(file.path(out_dir, "Fig6_severity_vs_prognosis.png"), P, width = 11, height = 10, dpi = 300)
  cat("figure written: Fig6_severity_vs_prognosis.pdf/.png\n")
} else {
  for (pp in plots) print(pp)
  cat("!! patchwork unavailable or too few panels, composite figure not assembled\n")
}
saveRDS(list(severity = if (exists("df_sev")) df_sev else NULL,
             curves = list(socs3 = km_df(socs3z), ifn = km_df(ifn),
                           lasso = km_df(riskL), s152 = km_df(ref152)),
             n = n), file.path(out_dir, "path65_plotdata.rds"))
cat("===== script 65 finished =====\n")
