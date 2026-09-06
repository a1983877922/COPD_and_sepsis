# =============================================================================
# 57_lasso_cox_signature.R
# Purpose: following the "machine-learning refined prognostic genes" idea from
#       online bioinformatics tutorials — select prognostic genes from a
#       candidate pool in a data-driven manner, rather than a fixed gene set
#       (152). Question: using the 149 shared programs as the candidate pool,
#       can Lasso-Cox select a prognostically valuable refined signature in
#       GSE65682?
# Candidate pool: 149 shared genes ∩ measurable genes in GSE65682 (GPL570) (hit count printed)
# Baseline control: the 152 whole-blood signature in that cohort had log-rank P=0.020, Cox HR=2.81 (95%CI 0.99-7.97)
# Method: cv.glmnet(family="cox", alpha=1), 10-fold CV; signature = nonzero-coefficient
#       genes at lambda.1se (falls back to lambda.min if fully shrunk);
#       risk score = linear combination of standardized expression
# Output: path57_lasso_{hits,genes,coef,cox_stats}.csv|txt + figures (path/CV/KM)
# Run: Rscript 57_lasso_cox_signature.R [OUT_DIR]   (needs glmnet; local GSE65682 matrix)
# =============================================================================
suppressPackageStartupMessages({
  library(GEOquery); library(survival); library(glmnet)
})
has_pROC <- requireNamespace("pROC", quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
if (length(args) >= 1 && nzchar(args[1])) script_dir <- args[1]
OUT_DIR <- if (length(args) >= 2 && nzchar(args[2])) args[2] else script_dir
out_dir <- OUT_DIR
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---------- Utility functions (same as scripts 15/53, conventions comparable) ----------
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
    # getGPL=FALSE: avoid triggering a GPL annotation download from NCBI (the server
    # often times out accessing NCBI).
    # The embedded table of an Affy matrix (GPL13667/HG-U219) usually has no Gene
    # Symbol column; probe->gene mapping is done by the offline
    # GPL13667_symbol_map.csv fallback below.
    gse <- getGEO(filename = matrix_file, getGPL = FALSE)
  } else {
    cat("No local matrix file; downloading", gse_acc, "online...\n")
    gse <- tryCatch(getGEO(gse_acc, GSEMatrix = TRUE, getGPL = TRUE),
                    error = function(e) NULL)
    if (is.null(gse)) { cat("download failed, skipping", label, "\n"); return(NULL) }
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
    ## --- Offline probe->gene annotation fallback -----------------------------
    ## The embedded table of an Affy series matrix (GPL13667/HG-U219) has no Gene
    ## Symbol column, and server access to NCBI often times out; instead read the
    ## offline mapping table deployed alongside the script:
    ##   GPL13667_symbol_map.csv  (header: probe\tsymbol, one row per probe)
    cand <- unique(c(file.path(dirname(matrix_file), "GPL13667_symbol_map.csv"),
                     file.path(out_dir, "GPL13667_symbol_map.csv"),
                     "GPL13667_symbol_map.csv"))
    mf <- cand[file.exists(cand)]
    if (length(mf) == 0) {
      cat("no symbol column and no offline annotation file; skipping\n"); return(NULL)
    }
    cat("using offline probe annotation:", mf[1], "\n")
    mtab <- read.delim(mf[1], header = TRUE, sep = "\t", quote = "",
                       stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("probe", "symbol") %in% colnames(mtab)))
      mtab <- read.csv(mf[1], header = TRUE, stringsAsFactors = FALSE)
    msym <- as.character(mtab$symbol); names(msym) <- as.character(mtab$probe)
    sym <- unname(msym[rownames(fd)])
    sym[is.na(sym)] <- ""
    sym <- sub(" /// .*$", "", sym)      # multi-gene probes: keep only the first symbol to avoid dirty symbols
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
  cat("after mapping:", nrow(expr), "genes x", ncol(expr), "samples\n")

  sf <- find_surv_fields(pd)
  if (length(sf$status) == 0) { cat("!! could not identify a status field\n"); return(NULL) }
  status_col <- sf$status[1]
  time_col <- if (length(sf$time) == 0) NULL else sf$time[1]
  cat("using status field:", status_col, " / time field:",
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
  status01 <- to_status01(status_raw)
  keep_samp <- !is.na(status01) & !is.na(time_raw) & time_raw > 0
  expr <- expr[, keep_samp, drop = FALSE]
  status01 <- status01[keep_samp]; time_raw <- time_raw[keep_samp]
  cat("samples retained with survival data:", ncol(expr), "(dead=", sum(status01), ")\n")
  list(expr = expr, time = time_raw, status = status01, label = label)
}

## =============================================================================
## Main workflow
## =============================================================================
cat("\n===== Lasso-Cox refined prognosis (candidate pool = 149 shared programs) =====\n")

## Candidate pool: 149 shared genes
gene_file <- file.path(out_dir, "path21_shared_up_genes.txt")
if (!file.exists(gene_file)) gene_file <- file.path(out_dir, "path2_shared_myeloid_genes.txt")
pool149 <- if (file.exists(gene_file)) {
  readLines(gene_file, warn = FALSE)
} else { character(0) }
pool149 <- trimws(pool149[pool149 != ""])
cat("149 candidate genes read:", length(pool149), "\n")

## Load GSE65682
coh <- load_cohort("GSE65682", "MARS (GSE65682)")
if (is.null(coh)) quit(status = 1)
expr <- coh$expr; time <- coh$time; status <- coh$status
cat("samples:", ncol(expr), "| deaths:", sum(status), "| 28-day follow-up:", all(time == 28), "\n")

## Hit check
hits <- intersect(pool149, rownames(expr))
cat("\n149 genes hit in GSE65682:", length(hits), "/", length(pool149), "\n")
if (length(hits) > 0) cat("hit genes:", paste(sort(hits), collapse = ", "), "\n")
writeLines(hits, file.path(out_dir, "path57_lasso_149_hits.txt"))

if (length(hits) < 10) {
  cat("!! too few hits (<10); Lasso is meaningless. Conclusion: the 149 single-cell-specific programs are barely\n",
      "   measurable in whole-blood arrays, supporting the need for the whole-blood-defined broad fingerprint (152)\n",
      "   for prognosis. Script ends early.\n")
  quit(status = 0)
}

## Standardized expression + Lasso-Cox
## glmnet convention: x = samples x genes (nrow must equal n = sample count); z-score each column (gene).
X <- expr[hits, , drop = FALSE]     # genes x samples
X <- t(X)                            # samples x genes
X <- scale(X)                        # standardize per column (gene)
X[is.na(X)] <- 0
set.seed(2026)
cvfit <- tryCatch(
  cv.glmnet(X, Surv(time, status), family = "cox", alpha = 1, nfolds = 10),
  error = function(e) { cat("!! cv.glmnet failed:", conditionMessage(e), "\n"); NULL })
if (is.null(cvfit)) quit(status = 1)

l1 <- coef(cvfit, s = "lambda.1se")
lmin <- coef(cvfit, s = "lambda.min")
nz1 <- hits[as.numeric(l1) != 0]
nzmin <- hits[as.numeric(lmin) != 0]
cat("genes selected at lambda.1se:", length(nz1), "| selected at lambda.min:", length(nzmin), "\n")

## If 1se fully shrinks, fall back to lambda.min; if both empty, end early
sig_genes <- if (length(nz1) > 0) nz1 else nzmin
s_used <- if (length(nz1) > 0) "lambda.1se" else "lambda.min"
if (length(sig_genes) == 0) {
  cat("!! Lasso selected no genes (the 149 pool carries no prognostic signal in GSE65682).\n",
      "Conclusion: compared with the 152 whole-blood signature (HR 2.81, log-rank P=0.020), the 149\n",
      "     single-cell programs carry no measurable whole-blood prognostic signal — further support\n",
      "     for the necessity of the 152. Script ends (figures skipped).\n")
  saveRDS(list(hits = hits, cvfit = cvfit, sig_genes = character(0)),
          file.path(out_dir, "path57_lasso_plotdata.rds"))
  quit(status = 0)
}

## risk score = linear combination of standardized expression (one row per sample)
b_hat <- as.numeric(coef(cvfit, s = s_used))[match(sig_genes, hits)]
names(b_hat) <- sig_genes
risk <- rowSums(X[, sig_genes, drop = FALSE] * b_hat)

## Cox + KM + AUC
cox <- summary(coxph(Surv(time, status) ~ risk))
hr <- cox$conf.int[1]; hr_lo <- cox$conf.int[3]; hr_hi <- cox$conf.int[4]
p_cox <- cox$coefficients[5]
med <- median(risk)
grp <- ifelse(risk > med, "High", "Low")
lr <- survdiff(Surv(time, status) ~ grp)
p_lr <- 1 - pchisq(lr$chisq, 1)
auc <- if (has_pROC) {
  pROC::auc(pROC::roc(status ~ risk, levels = c(0, 1), direction = "<"))
} else { NA }

cat("\n===== Lasso-Cox refined signature results (GSE65682) =====\n")
cat("signature genes (n =", length(sig_genes), ",", s_used, "):",
    paste(sig_genes, collapse = ", "), "\n")
cat("Cox: HR =", round(hr, 3), "(", round(hr_lo, 3), "-", round(hr_hi, 3),
    ") P =", format(p_cox, digits = 3), "\n")
cat("KM High vs Low log-rank P =", format(p_lr, digits = 3), "\n")
cat("28-day death AUC =", round(as.numeric(auc), 3), "\n")
cat("reference, 152 signature: log-rank P = 0.020, HR = 2.81\n")

## Output
coef_df <- data.frame(gene = names(b_hat), coef = as.numeric(b_hat),
                      stringsAsFactors = FALSE)
write.csv(coef_df, file.path(out_dir, "path57_lasso_selected_genes.csv"), row.names = FALSE)
stats_txt <- sprintf(
  "GSE65682 Lasso-Cox refined signature (pool=149, hits=%d, selected=%d [%s])\nCox HR=%.3f (%.3f-%.3f) P=%.3e\nlog-rank P=%.3e (High deaths/total, see KM)\nAUC28d=%.3f\nBaseline 152: log-rank P=0.020 HR=2.81\n",
  length(hits), length(sig_genes), s_used, hr, hr_lo, hr_hi, p_cox, p_lr, auc)
writeLines(stats_txt, file.path(out_dir, "path57_lasso_cox_stats.txt"))

## Figures: lasso CV + coefficient path + KM
suppressPackageStartupMessages(library(ggplot2))
tryCatch({
  pdf(file.path(out_dir, "path57_lasso_cv.pdf"), width = 7, height = 5)
  plot(cvfit); title("Lasso-Cox 10-fold CV (149-gene candidate pool)")
  dev.off()
}, error = function(e) cat("CV figure failed:", conditionMessage(e), "\n"))

## KM: extract via survfit summary instead of hard-coding strata naming
## (compatible with "grp=High"/"High" etc.), avoiding empty strata subsets →
## colour scale "No shared levels" warnings + blank curves.
df_km  <- data.frame(time = time, status = status, grp = grp)
fit_km <- survfit(Surv(time, status) ~ grp, data = df_km)
sf     <- summary(fit_km)
km_grp <- sub("^.*=", "", as.character(sf$strata))   # "grp=High" -> "High"
km_df  <- data.frame(t = sf$time, s = sf$surv,
                     grp = factor(km_grp, levels = c("High", "Low")))
cat("KM strata:", paste(unique(as.character(sf$strata)), collapse = " | "),
    "| curve points:", nrow(km_df), "| High points:", sum(km_grp == "High"),
    "| Low points:", sum(km_grp == "Low"), "\n")
p_km <- ggplot(km_df, aes(t, s, colour = grp)) +
  geom_step(linewidth = 0.9) +
  scale_colour_manual(
    values = c("High" = "#B2182B", "Low" = "#2166AC"),
    breaks = c("High", "Low"), name = "Risk group", drop = FALSE) +
  labs(title = sprintf("Lasso signature (%d genes): log-rank P = %.3f",
                       length(sig_genes), p_lr),
       x = "Days", y = "Survival probability") +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "path57_lasso_KM.pdf"), p_km, width = 6, height = 4.5)

saveRDS(list(hits = hits, cvfit = cvfit, sig_genes = sig_genes, coef = b_hat,
             risk = risk, cox_p = p_cox, logrank_p = p_lr, auc = auc),
        file.path(out_dir, "path57_lasso_plotdata.rds"))
cat("\n===== Script 57 done =====\n")
