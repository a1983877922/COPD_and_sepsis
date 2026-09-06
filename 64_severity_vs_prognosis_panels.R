# =============================================================================
# 64: Figure S6, four panels of severity vs prognosis
# Purpose: expand Fig. S6 to 4 panels, one figure conveying "severity biology !=
#   prognostic fingerprint":
#   a) 152-gene whole-blood COPD signature (positive reference)  HR~2.81, log-rank~0.020
#   b) 149-gene program, Lasso 12-gene signature                 HR~1.07, log-rank~0.206
#   c) SOCS3 (single-gene z score)                               HR~1.02, log-rank~0.855
#   d) IFN/purple module (65 genes) score                        HR~0.80, log-rank~0.217
#   Same cohort GSE65682 n=468, same median cut-off, same style; reconstructed from
#   the closed-loop results of script 63 plus archived scores.
# Data sources: a reads path15_GSE65682_train_scores.csv; b reads path57_lasso_plotdata.rds;
#          c/d are recomputed here from the current offline-annotated expr matrix
#          (SOCS3 single gene / purple module).
# Output: FigS6_axis_vs_fingerprint_KM.{pdf,png} / path64_stats.txt / path64_plotdata.rds
# Run:    Rscript 64_severity_vs_prognosis_panels.R [OUT_DIR]
# =============================================================================
suppressPackageStartupMessages({
  library(GEOquery); library(survival); library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
if (length(args) >= 1 && nzchar(args[1])) script_dir <- args[1]
out_dir <- if (length(args) >= 2 && nzchar(args[2])) args[2] else script_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---------- load_cohort (same as script 63) ----------
find_surv_fields <- function(pd) {
  cols <- colnames(pd)
  evt_pat <- "mortal|surviv|outcome|dead|alive|fate"; tme_pat <- "time|follow"
  status_cand <- cols[grepl(evt_pat, cols, ignore.case = TRUE) & !grepl(tme_pat, cols, ignore.case = TRUE)]
  time_cand   <- cols[grepl(tme_pat, cols, ignore.case = TRUE) & !grepl(evt_pat, cols, ignore.case = TRUE)]
  if (length(status_cand) == 0 || length(time_cand) == 0) {
    char_cols <- cols[grepl("^characteristics", cols, ignore.case = TRUE)]
    for (cc in char_cols) {
      vals <- trimws(tolower(as.character(pd[[cc]]))); vals <- vals[!is.na(vals) & vals != ""]
      if (length(vals) == 0) next
      has_evt <- any(grepl(evt_pat, vals)); has_tme <- any(grepl("time|follow|days", vals))
      val_only <- sub("^[^:]*:\\s*", "", vals); all_num <- all(suppressWarnings(!is.na(as.numeric(val_only))))
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
load_cohort <- function(gse_acc) {
  matrix_file <- file.path(out_dir, paste0(gse_acc, "_series_matrix.txt.gz"))
  if (!file.exists(matrix_file)) matrix_file <- file.path(out_dir, paste0(gse_acc, "_matrix.txt.gz"))
  gse <- getGEO(filename = matrix_file, getGPL = FALSE)
  expr <- exprs(gse); pd <- pData(gse); fd <- fData(gse)
  sym_col <- NULL
  for (cand in c("Gene Symbol","GeneSymbol","Symbol","Gene.Symbol","gene_symbol","GENE_SYMBOL"))
    if (cand %in% colnames(fd)) { sym_col <- cand; break }
  if (is.null(sym_col)) { m <- grep("symbol", colnames(fd), ignore.case = TRUE)
    if (length(m) > 0) sym_col <- colnames(fd)[m[1]] }
  if (is.null(sym_col)) {
    mf <- file.path(out_dir, "GPL13667_symbol_map.csv")
    mtab <- read.delim(mf, header = TRUE, sep = "\t", quote = "",
                       stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("probe","symbol") %in% colnames(mtab))) mtab <- read.csv(mf, stringsAsFactors = FALSE)
    msym <- as.character(mtab$symbol); names(msym) <- as.character(mtab$probe)
    sym <- unname(msym[rownames(fd)]); sym[is.na(sym)] <- ""; sym <- sub(" /// .*$", "", sym)
    names(sym) <- rownames(fd)
  } else {
    sym <- as.character(fd[[sym_col]]); names(sym) <- rownames(fd)
  }
  keep <- sym != "" & !is.na(sym) & sym != "---"
  expr <- expr[keep, , drop = FALSE]; sym <- sym[keep]
  sym_uniq <- unique(sym)
  expr_sym <- matrix(NA, length(sym_uniq), ncol(expr), dimnames = list(sym_uniq, colnames(expr)))
  for (g in sym_uniq) { idx <- which(sym == g)
    expr_sym[g, ] <- if (length(idx) == 1) expr[idx, ] else colMeans(expr[idx, , drop = FALSE]) }
  expr <- expr_sym
  sf <- find_surv_fields(pd)
  status_col <- sf$status[1]; time_col <- if (length(sf$time) == 0) NULL else sf$time[1]
  status_raw <- as.character(pd[[status_col]])
  time_raw <- if (is.null(time_col)) rep(28, nrow(pd)) else suppressWarnings(as.numeric(as.character(pd[[time_col]])))
  names(status_raw) <- rownames(pd); names(time_raw) <- rownames(pd)
  common <- intersect(colnames(expr), rownames(pd))
  expr <- expr[, common, drop = FALSE]; status_raw <- status_raw[common]; time_raw <- time_raw[common]
  status01 <- to_status01(status_raw)
  keep_samp <- !is.na(status01) & !is.na(time_raw) & time_raw > 0
  expr <- expr[, keep_samp, drop = FALSE]; status01 <- status01[keep_samp]; time_raw <- time_raw[keep_samp]
  list(expr = expr, time = time_raw, status = status01)
}

cat("===== 64: FigS6 four panels (severity != prognosis) =====\n")
coh <- load_cohort("GSE65682"); expr <- coh$expr; time <- coh$time; status <- coh$status
n <- ncol(expr); cat("n=", n, " dead=", sum(status), "\n")

purple <- trimws(readLines(file.path(out_dir, "path13_sc_module_purple.txt"))); purple <- purple[purple != ""]
s15 <- read.csv(file.path(out_dir, "path15_GSE65682_train_scores.csv"))
p57 <- readRDS(file.path(out_dir, "path57_lasso_plotdata.rds"))

socs3  <- as.numeric(scale(expr["SOCS3", ]))
ifn    <- score_geneset(expr, purple)
ref152 <- s15$copd
riskL  <- as.numeric(p57$risk)

km_stat <- function(score) {
  ok <- is.finite(score)
  df <- data.frame(time = time[ok], status = status[ok], score = as.numeric(score[ok]))
  df <- df[!is.na(df$score), ]
  df$grp <- ifelse(df$score >= median(df$score), "High", "Low")
  fit <- coxph(Surv(time, status) ~ score, data = df); cs <- summary(fit)$conf.int
  lr <- survdiff(Surv(time, status) ~ grp, data = df)
  list(HR = unname(cs[1]), p_lr = 1 - pchisq(lr$chisq, 1), df = df)
}
rA <- km_stat(ref152); rB <- km_stat(riskL); rC <- km_stat(socs3); rD <- km_stat(ifn)

curve_of <- function(r) {
  fit <- survfit(Surv(time, status) ~ grp, data = r$df); sf <- summary(fit)
  g <- sub("^.*=", "", as.character(sf$strata))
  data.frame(t = sf$time, s = sf$surv, grp = factor(g, levels = c("High", "Low")))
}
labs <- c(
  sprintf("a  152-gene COPD signature\nlog-rank P = %.3f | HR = %.2f", rA$p_lr, rA$HR),
  sprintf("b  149-gene program, Lasso (12 genes)\nlog-rank P = %.3f | HR = %.2f", rB$p_lr, rB$HR),
  sprintf("c  SOCS3\nlog-rank P = %.3f | HR = %.2f", rC$p_lr, rC$HR),
  sprintf("d  IFN module (65 genes)\nlog-rank P = %.3f | HR = %.2f", rD$p_lr, rD$HR)
)
curves <- rbind(cbind(curve_of(rA), sig = factor(labs[1], levels = labs)),
                cbind(curve_of(rB), sig = factor(labs[2], levels = labs)),
                cbind(curve_of(rC), sig = factor(labs[3], levels = labs)),
                cbind(curve_of(rD), sig = factor(labs[4], levels = labs)))
p <- ggplot(curves, aes(t, s, colour = grp)) +
  geom_step(linewidth = 0.8) +
  facet_wrap(~ sig, ncol = 2) +
  scale_colour_manual(values = c("High" = "#B2182B", "Low" = "#2166AC"),
                      breaks = c("High", "Low"), name = "Risk group", drop = FALSE) +
  scale_x_continuous(limits = c(0, NA)) +
  labs(x = "Days", y = "Survival probability") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(size = 8.5, lineheight = 0.95),
        legend.position = "bottom")
ggsave(file.path(out_dir, "FigS6_axis_vs_fingerprint_KM.pdf"), p, width = 8.5, height = 7)
ggsave(file.path(out_dir, "FigS6_axis_vs_fingerprint_KM.png"), p, width = 8.5, height = 7, dpi = 300)

stats <- sprintf(
  "FigS6 four panels (GSE65682 n=%d): a 152 HR=%.2f p_lr=%.3f | b 149-Lasso HR=%.2f p_lr=%.3f | c SOCS3 HR=%.2f p_lr=%.3f | d IFN HR=%.2f p_lr=%.3f\n",
  n, rA$HR, rA$p_lr, rB$HR, rB$p_lr, rC$HR, rC$p_lr, rD$HR, rD$p_lr)
writeLines(stats, file.path(out_dir, "path64_stats.txt"))
saveRDS(list(curves = curves, rA = rA, rB = rB, rC = rC, rD = rD, n = n),
        file.path(out_dir, "path64_plotdata.rds"))
cat(stats)
cat("===== script 64 finished =====\n")
