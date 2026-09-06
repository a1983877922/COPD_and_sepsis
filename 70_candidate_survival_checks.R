# =========================================================================
# 53 four-gene severity / prognostic test
# Question: using only the 4 genes that were consistently up-regulated across the four
#           blood-side groups (SOCS3 / CLIC1 / MFSD10 / ITGB2), can we reproduce a prognosis
#           similar to the 152-gene signature in GSE65682 (MARS, n=479, 28-day follow-up)?
# Method: 4-gene mean z-score -> univariate Cox + age-adjusted Cox + KM (log-rank)
#         + 28-day mortality AUC; High/Low = median split (same criteria as script 15)
#
# Reuses the load_cohort/score_geneset logic of script 15 (same cohort, same criteria).
# Outputs: path53_GSE65682_4gene_cox_stats.txt
#          path53_GSE65682_4gene_scores.csv
# Run: Rscript 70_candidate_survival_checks.R [OUT_DIR]
# =========================================================================
suppressPackageStartupMessages({
  library(GEOquery); library(survival)
})
has_pROC <- requireNamespace("pROC", quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
if (length(args) >= 1 && nzchar(args[1])) script_dir <- args[1]
OUT_DIR <- if (length(args) >= 2 && nzchar(args[2])) args[2] else script_dir
out_dir <- OUT_DIR
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---------- Utility functions (same as script 15, so the criteria are comparable) ----------
score_geneset <- function(expr_mat, genes) {
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) < 2) return(rep(NA_real_, ncol(expr_mat)))
  sub <- expr_mat[genes, , drop = FALSE]
  z   <- t(scale(t(sub)))
  z[is.na(z)] <- 0
  colMeans(z)
}

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
  cat("\n----- processing cohort:", label, "(", gse_acc, ") -----\n")
  matrix_file <- file.path(out_dir, paste0(gse_acc, "_series_matrix.txt.gz"))
  if (!file.exists(matrix_file))
    matrix_file <- file.path(out_dir, paste0(gse_acc, "_matrix.txt.gz"))
  if (file.exists(matrix_file)) {
    cat("Using the local series matrix:", matrix_file, "\n")
    gse <- getGEO(filename = matrix_file, getGPL = TRUE)
  } else {
    cat("No local matrix file, downloading online", gse_acc, "...\n")
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
  if (is.null(sym_col)) { cat("symbol column not found, skipping\n"); return(NULL) }
  sym <- as.character(fd[[sym_col]])
  if (use_ga) {
    sym <- vapply(strsplit(sym, "//", fixed = TRUE), function(x) {
      if (length(x) >= 2) trimws(x[2]) else ""
    }, character(1))
  }
  names(sym) <- rownames(fd)
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
  if (length(sf$status) == 0) { cat("!! status field not identified\n"); return(NULL) }
  status_col <- sf$status[1]
  time_col <- if (length(sf$time) == 0) NULL else sf$time[1]
  cat("Using status field:", status_col, " / time field:",
      if (is.null(time_col)) "fixed 28 days" else time_col, "\n")
  status_raw <- as.character(pd[[status_col]])
  if (is.null(time_col)) {
    time_raw <- rep(28, nrow(pd))
  } else {
    time_raw <- suppressWarnings(as.numeric(as.character(pd[[time_col]])))
  }
  names(status_raw) <- rownames(pd); names(time_raw) <- rownames(pd)
  common <- intersect(colnames(expr), rownames(pd))
  expr <- expr[, common, drop = FALSE]
  status_raw <- status_raw[common]; time_raw <- time_raw[common]
  age_col <- grep("age", colnames(pd), ignore.case = TRUE, value = TRUE)
  age_col <- age_col[!grepl("stage|agent|percent", age_col, ignore.case = TRUE)]
  age <- if (length(age_col) > 0) suppressWarnings(as.numeric(as.character(pd[[age_col[1]]][common]))) else NULL
  status01 <- to_status01(status_raw)
  keep_samp <- !is.na(status01) & !is.na(time_raw) & time_raw > 0
  expr <- expr[, keep_samp, drop = FALSE]
  status01 <- status01[keep_samp]; time_raw <- time_raw[keep_samp]
  if (!is.null(age)) age <- age[keep_samp]
  cat("Samples retained with survival data:", ncol(expr), "(dead=", sum(status01), ")\n")
  list(expr = expr, time = time_raw, status = status01, label = label, age = age)
}

## ---------- Main analysis (v2: three-gene composite + four-gene control + single genes one by one, for finalizing the R2 wording) ----------
HUBS <- c("SOCS3", "CLIC1", "MFSD10")
ALL4 <- c(HUBS, "ITGB2")

cat("===== candidate severity gene prognostic test (GSE65682, reusing the criteria of 53) =====\n")
cohort <- load_cohort("GSE65682", "GSE65682_hub")
if (is.null(cohort)) quit(status = 1)

res_lines <- character()
run_score <- function(tag, genes) {
  hits <- intersect(genes, rownames(cohort$expr))
  cat("Genes matched:", tag, " ->", paste(hits, collapse = ","), " (", length(hits), "/", length(genes), ")\n", sep = "")
  if (length(hits) < 2) { cat(" !! ", tag, " fewer than 2 matched, skipping\n"); return(invisible(NULL)) }
  score <- score_geneset(cohort$expr, genes)
  df <- data.frame(time = cohort$time, status = cohort$status, score = score)
  df$grp <- ifelse(df$score > median(df$score, na.rm = TRUE), "High", "Low")
  fit <- coxph(Surv(time, status) ~ score, data = df); s <- summary(fit)
  hr <- s$conf.int[1, "exp(coef)"]; lo <- s$conf.int[1, "lower .95"]; hi <- s$conf.int[1, "upper .95"]
  pv <- s$coefficients[1, "Pr(>|z|)"]
  lr <- survdiff(Surv(time, status) ~ grp, data = df); lrp <- 1 - pchisq(lr$chisq, df = 1)
  auc <- tryCatch(as.numeric(pROC::roc(df$status, df$score, quiet = TRUE)$auc), error = function(e) NA)
  line <- sprintf("[%s] composite score, univariate Cox: HR=%.3f (%.3f-%.3f), p=%.3e | log-rank p=%.3e | AUC=%.3f",
                  tag, hr, lo, hi, pv, lrp, auc)
  cat(" ", line, "\n"); res_lines <<- c(res_lines, line)
}
run_gene <- function(gene) {
  if (!gene %in% rownames(cohort$expr)) { cat("!! ", gene, " is not in the expression matrix\n"); return(invisible(NULL)) }
  z <- scale(as.numeric(cohort$expr[gene, ]))[, 1]
  df <- data.frame(time = cohort$time, status = cohort$status, score = z)
  df$grp <- ifelse(df$score > median(df$score, na.rm = TRUE), "High", "Low")
  fit <- coxph(Surv(time, status) ~ score, data = df); s <- summary(fit)
  hr <- s$conf.int[1, "exp(coef)"]; lo <- s$conf.int[1, "lower .95"]; hi <- s$conf.int[1, "upper .95"]
  pv <- s$coefficients[1, "Pr(>|z|)"]
  lr <- survdiff(Surv(time, status) ~ grp, data = df); lrp <- 1 - pchisq(lr$chisq, df = 1)
  auc <- tryCatch(as.numeric(pROC::roc(df$status, df$score, quiet = TRUE)$auc), error = function(e) NA)
  line <- sprintf("[%s single-gene z] univariate Cox: HR=%.3f (%.3f-%.3f), p=%.3e | log-rank p=%.3e | AUC=%.3f",
                  gene, hr, lo, hi, pv, lrp, auc)
  cat(" ", line, "\n"); res_lines <<- c(res_lines, line)
}

run_score("3-hub-gene composite (SOCS3/CLIC1/MFSD10)", HUBS)
run_score("4-candidate-gene composite (including ITGB2, control of 53)", ALL4)
for (g in ALL4) run_gene(g)

writeLines(res_lines, file.path(out_dir, "path70_hub_survival_stats.txt"))
cat("\nSaved: path70_hub_survival_stats.txt\n")
cat("===== 70 done =====\n")
