# =============================================================================
# 63_axis_prognosis_socs3_ifn.R
# Purpose: closed-loop test - directly test in GSE65682 whole blood (n=468) whether the
#       "single-cell mechanism axes" carry a 28-day mortality prognostic signal:
#         A. SOCS3 (single-gene z)
#         B. IFN/purple module score (path13_sc_module_purple.txt, 65 genes)
#         C. 149 shared program score (path21_shared_up_genes.txt) [control, corroborates FigS6]
#         D. 152 whole-blood signature (reads archived path15 scores, positive reference, should reproduce HR~2.81/0.020)
#       Expected: A/B/C not significant, D significant -> confirms that the IFN/SOCS3 axis
#       marks severity rather than outcome; prognosis is carried by an independent
#       whole-blood inflammatory fingerprint. Results[46] can honestly be rewritten accordingly.
# Outputs: path63_axis_prognosis_stats.txt / path63_IFN_SOCS3_KM.{pdf,png} /
#       path63_plotdata.rds
# Run: Rscript this script [OUT_DIR]  (requires GEOquery/survival/ggplot2)
# =============================================================================
suppressPackageStartupMessages({
  library(GEOquery); library(survival); library(ggplot2)
})
has_pROC <- requireNamespace("pROC", quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
if (length(args) >= 1 && nzchar(args[1])) script_dir <- args[1]
OUT_DIR <- if (length(args) >= 2 && nzchar(args[2])) args[2] else script_dir
out_dir <- OUT_DIR
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---------- Reuse load_cohort (same criteria as 15/53/57/58/59) ----------
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
      vals <- trimws(tolower(as.character(pd[[cc]]))); vals <- vals[!is.na(vals) & vals != ""]
      if (length(vals) == 0) next
      has_evt <- any(grepl(evt_pat, vals)); has_tme <- any(grepl("time|follow|days", vals))
      val_only <- sub("^[^:]*:\\s*", "", vals)
      all_num <- all(suppressWarnings(!is.na(as.numeric(val_only))))
      if (length(status_cand) == 0 && has_evt && !all_num) status_cand <- c(status_cand, cc)
      if (length(time_cand) == 0 && has_tme && !has_evt) time_cand <- c(time_cand, cc)
    }
  }
  list(status = status_cand, time = time_cand)
}
to_status01 <- function(v) {
  v <- trimws(tolower(as.character(v))); v <- sub("^[^:]*:\\s*", "", v)
  out <- rep(NA_real_, length(v))
  out[v %in% c("dead","death","died","deceased","non-survivor","nonsurvivor","1","yes","event","mortality")] <- 1
  out[v %in% c("alive","survived","survivor","live","0","no","censored","survival")] <- 0
  out
}
score_geneset <- function(expr_mat, genes) {
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) < 2) return(rep(NA_real_, ncol(expr_mat)))
  z <- t(scale(t(expr_mat[genes, , drop = FALSE]))); z[is.na(z)] <- 0
  colMeans(z)
}
load_cohort <- function(gse_acc, label) {
  cat("\n----- processing cohort:", label, "(", gse_acc, ") -----\n")
  matrix_file <- file.path(out_dir, paste0(gse_acc, "_series_matrix.txt.gz"))
  if (!file.exists(matrix_file)) matrix_file <- file.path(out_dir, paste0(gse_acc, "_matrix.txt.gz"))
  if (file.exists(matrix_file)) {
    cat("Using the local series matrix:", matrix_file, "\n")
    gse <- getGEO(filename = matrix_file, getGPL = FALSE)
  } else {
    cat("No local matrix, downloading online", gse_acc, "...\n")
    gse <- tryCatch(getGEO(gse_acc, GSEMatrix = TRUE, getGPL = TRUE), error = function(e) NULL)
    if (is.null(gse)) { cat("Download failed\n"); return(NULL) }; gse <- gse[[1]]
  }
  expr <- exprs(gse); pd <- pData(gse); fd <- fData(gse)
  cat("Expression matrix:", nrow(expr), "x", ncol(expr), "\n")
  sym_col <- NULL
  for (cand in c("Gene Symbol","GeneSymbol","Symbol","Gene.Symbol","gene_symbol","GENE_SYMBOL"))
    if (cand %in% colnames(fd)) { sym_col <- cand; break }
  if (is.null(sym_col)) { m <- grep("symbol", colnames(fd), ignore.case = TRUE)
    if (length(m) > 0) sym_col <- colnames(fd)[m[1]] }
  use_ga <- FALSE
  if (is.null(sym_col) && "gene_assignment" %in% colnames(fd)) { sym_col <- "gene_assignment"; use_ga <- TRUE }
  sym <- NULL
  if (!is.null(sym_col)) {
    sym <- as.character(fd[[sym_col]])
    if (use_ga) sym <- vapply(strsplit(sym, "//", fixed = TRUE), function(x)
      if (length(x) >= 2) trimws(x[2]) else "", character(1))
    names(sym) <- rownames(fd)
  } else {
    cand <- unique(c(file.path(dirname(matrix_file), "GPL13667_symbol_map.csv"),
                     file.path(out_dir, "GPL13667_symbol_map.csv"), "GPL13667_symbol_map.csv"))
    mf <- cand[file.exists(cand)]
    if (length(mf) == 0) { cat("No symbol column and no offline annotation\n"); return(NULL) }
    cat("Using offline probe annotation:", mf[1], "\n")
    mtab <- read.delim(mf[1], header = TRUE, sep = "\t", quote = "",
                       stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("probe","symbol") %in% colnames(mtab))) mtab <- read.csv(mf[1], stringsAsFactors = FALSE)
    msym <- as.character(mtab$symbol); names(msym) <- as.character(mtab$probe)
    sym <- unname(msym[rownames(fd)]); sym[is.na(sym)] <- ""; sym <- sub(" /// .*$", "", sym)
    names(sym) <- rownames(fd)
  }
  keep <- sym != "" & !is.na(sym) & sym != "---"
  expr <- expr[keep, , drop = FALSE]; sym <- sym[keep]
  sym_uniq <- unique(sym)
  expr_sym <- matrix(NA, length(sym_uniq), ncol(expr), dimnames = list(sym_uniq, colnames(expr)))
  for (g in sym_uniq) { idx <- which(sym == g)
    expr_sym[g, ] <- if (length(idx) == 1) expr[idx, ] else colMeans(expr[idx, , drop = FALSE]) }
  expr <- expr_sym
  cat("After mapping:", nrow(expr), "genes\n")
  sf <- find_surv_fields(pd)
  if (length(sf$status) == 0) { cat("!! no status field\n"); return(NULL) }
  status_col <- sf$status[1]; time_col <- if (length(sf$time) == 0) NULL else sf$time[1]
  status_raw <- as.character(pd[[status_col]])
  if (is.null(time_col)) time_raw <- rep(28, nrow(pd)) else
    time_raw <- suppressWarnings(as.numeric(as.character(pd[[time_col]])))
  names(status_raw) <- rownames(pd); names(time_raw) <- rownames(pd)
  common <- intersect(colnames(expr), rownames(pd))
  expr <- expr[, common, drop = FALSE]; status_raw <- status_raw[common]; time_raw <- time_raw[common]
  status01 <- to_status01(status_raw)
  keep_samp <- !is.na(status01) & !is.na(time_raw) & time_raw > 0
  expr <- expr[, keep_samp, drop = FALSE]; status01 <- status01[keep_samp]; time_raw <- time_raw[keep_samp]
  cat("Samples retained with survival data:", ncol(expr), "(dead=", sum(status01), ")\n")
  list(expr = expr, time = time_raw, status = status01, label = label)
}

cat("===== 63: IFN/SOCS3 whole-blood prognostic closed-loop test =====\n")
coh <- load_cohort("GSE65682", "MARS (GSE65682)")
if (is.null(coh)) quit(status = 1)
expr <- coh$expr; time <- coh$time; status <- coh$status; n <- ncol(expr)

purple <- trimws(readLines(file.path(out_dir, "path13_sc_module_purple.txt")))
purple <- purple[purple != ""]
pool149 <- trimws(readLines(file.path(out_dir, "path21_shared_up_genes.txt")))
pool149 <- pool149[pool149 != ""]
cat("purple genes matched:", sum(purple %in% rownames(expr)), "/", length(purple),
    "| 149 genes matched:", sum(pool149 %in% rownames(expr)), "/", length(pool149), "\n")

## Axis scores
socs3_z <- as.numeric(scale(expr["SOCS3", ]))
ifn     <- score_geneset(expr, purple)
prog149 <- score_geneset(expr, pool149)
## 152-gene reference (archived, sample order already verified)
s15 <- read.csv(file.path(out_dir, "path15_GSE65682_train_scores.csv"))
if (nrow(s15) != n) { cat("!! path15 scores length mismatch\n"); quit(status = 1) }
ref152 <- s15$copd

axes <- list(
  `SOCS3 (z)`   = socs3_z,
  `IFN module (65)` = ifn,
  `149 program`  = prog149,
  `152 signature (ref)` = ref152
)

test_axis <- function(score, lab) {
  ok <- is.finite(score)
  df <- data.frame(time = time[ok], status = status[ok], score = as.numeric(score[ok]))
  df <- df[!is.na(df$score), ]
  df$grp <- ifelse(df$score >= median(df$score), "High", "Low")
  fit <- coxph(Surv(time, status) ~ score, data = df)
  cs <- summary(fit)$conf.int
  lr <- survdiff(Surv(time, status) ~ grp, data = df)
  auc <- if (has_pROC) as.numeric(pROC::auc(pROC::roc(status ~ score, data = df,
                                                     levels = c(0, 1), direction = "<"))) else NA
  list(HR = unname(cs[1]), lo = unname(cs[3]), hi = unname(cs[4]),
       p_cox = summary(fit)$coefficients[5], p_lr = 1 - pchisq(lr$chisq, 1),
       auc = auc, df = df, lab = lab)
}

res <- lapply(names(axes), function(a) test_axis(axes[[a]], a))
names(res) <- names(axes)

cat("\n===== axis prognosis summary (GSE65682, n=", n, ", dead=", sum(status), ") =====\n", sep = "")
lines <- c()
for (a in names(res)) {
  r <- res[[a]]
  cat(sprintf("%-20s HR=%.2f (%.2f-%.2f) CoxP=%.3g | log-rank P=%.3g | AUC=%.3f\n",
              a, r$HR, r$lo, r$hi, r$p_cox, r$p_lr, r$auc))
  lines <- c(lines, sprintf("%s\tHR=%.3f\tCI=%.3f-%.3f\tCoxP=%.3e\tlogrankP=%.3e\tAUC=%.3f",
                            a, r$HR, r$lo, r$hi, r$p_cox, r$p_lr, r$auc))
}
writeLines(lines, file.path(out_dir, "path63_axis_prognosis_stats.txt"))

## KM two panels: SOCS3 + IFN module
km_curve <- function(df) {
  fit <- survfit(Surv(time, status) ~ grp, data = df); sf <- summary(fit)
  g <- sub("^.*=", "", as.character(sf$strata))
  data.frame(t = sf$time, s = sf$surv, grp = factor(g, levels = c("High", "Low")))
}
c_socs <- km_curve(res[["SOCS3 (z)"]]$df); c_socs$sig <- sprintf("SOCS3 (log-rank P = %.3f)", res[["SOCS3 (z)"]]$p_lr)
c_ifn  <- km_curve(res[["IFN module (65)"]]$df); c_ifn$sig <- sprintf("IFN module 65 genes (P = %.3f)", res[["IFN module (65)"]]$p_lr)
c_all <- rbind(c_socs, c_ifn)
p <- ggplot(c_all, aes(t, s, colour = grp)) +
  geom_step(linewidth = 0.9) +
  facet_wrap(~ sig, ncol = 2) +
  scale_colour_manual(values = c("High" = "#B2182B", "Low" = "#2166AC"),
                      breaks = c("High", "Low"), name = "Risk group", drop = FALSE) +
  scale_x_continuous(limits = c(0, NA)) +
  labs(x = "Days", y = "Survival probability") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "path63_IFN_SOCS3_KM.pdf"), p, width = 10, height = 4.8)
ggsave(file.path(out_dir, "path63_IFN_SOCS3_KM.png"), p, width = 10, height = 4.8, dpi = 300)

saveRDS(list(axes = axes, res = res, n = n), file.path(out_dir, "path63_plotdata.rds"))
cat("\nOutputs: path63_axis_prognosis_stats.txt / path63_IFN_SOCS3_KM.pdf/.png / path63_plotdata.rds\n")
cat("===== script 63 done =====\n")
