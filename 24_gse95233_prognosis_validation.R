# =========================================================================
# Script 24: GSE95233 prognosis validation
#
# Purpose: external prognostic validation in a third sepsis cohort (GSE95233),
#          complementing GSE65682 (training) + E-MTAB-4451 (external leukocyte).
#
# Data: GSE95233 (whole-blood array, GPL570 = Affymetrix HG-U133 Plus 2.0)
#   51 sepsis patients (longitudinal sampling, 2 time points per patient) + 22 healthy controls
#   with 28-day survival (17 deaths / 34 survivors)
#
# Reuses load_cohort / score_geneset / find_surv_fields / to_status01 from script 14,
# adds a COPD inflammation signature (reads path15_copd_signature_genes.txt), and runs
# Cox + KM + ROC for the three signatures.
#
# Note: on a whole-blood array the pDC abundance signature is expected to be diluted
#       (as in GSE65682), so a non-significant result is normal; the focus is on the
#       direction of the COPD inflammation and immunoparalysis signatures.
# Limitation: GSE95233 is longitudinally sampled (two time points per patient share the
#       outcome) while Cox treats samples as independent; n = 51 is small, so conclusions
#       are stated as "consistency of direction" rather than "strong significance".
#
# Requires: survival, GEOquery, ggplot2 (mandatory); survminer, pROC (optional)
# =========================================================================

.this_file <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
if (.this_file == "" || !file.exists(.this_file)) {
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1])
}
SCRIPT_DIR  <- dirname(.this_file)
config_file <- file.path(SCRIPT_DIR, "00_config.R")
if (!file.exists(config_file)) stop("Cannot find 00_config.R: ", config_file)
source(config_file)

for (pkg in c("GEOquery", "survival")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Not installed: ", pkg)
}
suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(ggplot2))

cat("\n==============================================================\n")
cat("GSE95233 external prognostic validation (28-day mortality)\n")
cat("==============================================================\n")

## Signatures
pdc_markers <- c("LILRA4","CLEC4C","IL3RA","TCF4","IRF7","IRF4","IRF8",
                 "GZMB","SPIB","BCL11A","RUNX2")         # pDC abundance
pdc_suppress <- c("SOCS3","SOCS1","TLR9","ISG20","MX1","ISG15","OAS2",
                  "IFI44L","MYD88","IRAK1")              # pDC immunoparalysis

copd_file <- file.path(out_dir, "path15_copd_signature_genes.txt")
copd_sig <- if (file.exists(copd_file)) trimws(readLines(copd_file)) else character(0)
cat("COPD inflammation signature gene count:", length(copd_sig), "\n")

## =========================================================================
## Utility functions (reused from script 14)
## =========================================================================

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
      if (length(status_cand) == 0 && has_evt && !all_num) {
        status_cand <- c(status_cand, cc)
        cat("  [detected] characteristics column", cc, "assigned to status (examples:",
            paste(head(vals, 3), collapse = "; "), ")\n")
      }
      if (length(time_cand) == 0 && has_tme && !has_evt) {
        time_cand <- c(time_cand, cc)
        cat("  [detected] characteristics column", cc, "assigned to time (examples:",
            paste(head(vals, 3), collapse = "; "), ")\n")
      }
    }
  }
  list(status = status_cand, time = time_cand)
}

to_status01 <- function(v) {
  v <- trimws(tolower(as.character(v)))
  v <- sub("^[^:]*:\\s*", "", v)
  v <- gsub("\\s+", "", v)   # strip spaces: "non survivor" -> "nonsurvivor"
  out <- rep(NA_real_, length(v))
  out[v %in% c("dead","death","died","deceased","non-survivor","nonsurvivor",
               "1","yes","event","mortality")] <- 1
  out[v %in% c("alive","survived","survivor","live","0","no","censored","survival")] <- 0
  out
}

load_cohort <- function(gse_acc, label, gpl_acc = "GPL570") {
  cat("\n----- Processing cohort:", label, "(", gse_acc, ") -----\n")
  matrix_file <- file.path(out_dir, paste0(gse_acc, "_series_matrix.txt.gz"))
  if (!file.exists(matrix_file)) {
    cat("No local series matrix; downloading...\n")
    url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/",
                  substr(gse_acc, 1, nchar(gse_acc) - 3), "nnn/", gse_acc,
                  "/matrix/", gse_acc, "_series_matrix.txt.gz")
    system(paste0("curl -sk --max-time 600 -o \"", matrix_file, "\" \"", url, "\""))
  }
  if (!file.exists(matrix_file)) stop("No series matrix: ", matrix_file)
  cat("Using series matrix:", matrix_file, "\n")

  # Manual parsing (works around a GEOquery parsing bug).
  # Key point: the GSE95233 series matrix is a multi-member gzip file (495 concatenated
  # gzip members); R's gzfile()/fread reads only the first member (~first 2383 lines),
  # so the begin marker is found but the end marker is lost.
  # Decompress fully with the system gunzip -c (which concatenates all members), then read
  # the plain text with readLines.
  txt_file <- sub("\\.gz$", "", matrix_file)
  if (!file.exists(txt_file) || file.info(txt_file)$size < 1000) {
    cat("Decompressing series matrix (multi-member gzip)...\n")
    system(paste0("gunzip -c '", matrix_file, "' > '", txt_file, "'"))
  }
  lines <- readLines(txt_file, warn = FALSE)

  # ---- Table (probe expression matrix) ----
  # Use fixed=TRUE byte-level matching to avoid regex grep failures caused by non-ASCII
  # bytes in the file (UTF-8 author names / special symbols) under a non-UTF-8 locale.
  begin <- grep("!series_matrix_table_begin", lines, fixed = TRUE)
  end   <- grep("!series_matrix_table_end", lines, fixed = TRUE)
  cat("Diagnostics: total lines", length(lines), "| begin marker:", paste(begin, collapse = ","),
      "| end marker:", paste(end, collapse = ","), "\n")
  if (length(begin) != 1 || length(end) != 1) stop("Abnormal series matrix table markers")
  tab <- data.table::fread(text = paste(lines[(begin + 1):(end - 1)], collapse = "\n"),
                           sep = "\t", header = TRUE, quote = "",
                           na.strings = c("NA", "null", "NULL", "Null"))
  probe_ids <- gsub("^\"|\"$", "", as.character(tab[[1]]))  # strip quotes (fread quote="" keeps them)
  expr <- as.matrix(tab[, -1, drop = FALSE])
  storage.mode(expr) <- "double"
  rownames(expr) <- probe_ids
  cat("Expression matrix:", nrow(expr), "probes x", ncol(expr), "samples\n")

  # ---- Sample metadata ----
  sample_lines <- lines[grepl("!Sample_", lines, fixed = TRUE)]
  sp <- strsplit(sample_lines, "\t", fixed = TRUE)
  fields <- vapply(sp, function(x) sub("!Sample_", "", x[1], fixed = TRUE), character(1))
  # characteristics_ch1 appears 4 times (gender/age/time point/survival); make.unique de-duplicates
  fields_uniq <- make.unique(fields)
  val_mat <- do.call(rbind, lapply(sp, function(x) x[-1]))
  gsm <- val_mat[which(fields == "geo_accession")[1], ]
  pd <- as.data.frame(val_mat, stringsAsFactors = FALSE)
  colnames(pd) <- gsm
  rownames(pd) <- fields_uniq
  pd <- as.data.frame(t(pd), stringsAsFactors = FALSE)
  pd[] <- lapply(pd, function(x) gsub("^\"|\"$", "", x))
  cat("Sample metadata:", nrow(pd), "samples x", ncol(pd), "fields\n")

  # ---- Probe -> symbol (GPL annotation) ----
  gpl <- tryCatch(getGEO(gpl_acc), error = function(e) NULL)
  if (is.null(gpl)) stop("Cannot retrieve GPL annotation ", gpl_acc)
  # The Table generic is not exported by Biobase, so access the dataTable@table slot of the
  # GPL object directly, falling back to the Table() generic (available once GEOquery is loaded).
  gtab <- tryCatch(gpl@dataTable@table, error = function(e) NULL)
  if (is.null(gtab)) gtab <- tryCatch(Table(gpl), error = function(e) NULL)
  if (is.null(gtab)) stop("Cannot access the GPL annotation table")
  sym_col <- intersect(c("Gene Symbol","GeneSymbol","Symbol","GENE_SYMBOL"), colnames(gtab))
  if (length(sym_col) == 0) stop("GPL annotation has no symbol column")
  sym_map <- setNames(as.character(gtab[[sym_col[1]]]), as.character(gtab$ID))
  cat("Diagnostics: GPL annotation columns:", paste(head(colnames(gtab), 6), collapse = ", "), "...\n")
  cat("Diagnostics: example probe IDs:", paste(head(probe_ids, 3), collapse = ", "),
      " | example GPL IDs:", paste(head(names(sym_map), 3), collapse = ", "),
      " | overlap:", sum(probe_ids %in% names(sym_map)), "/", length(probe_ids), "\n")
  sym <- sym_map[probe_ids]
  keep <- !is.na(sym) & sym != "" & sym != "---"
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

  # ---- Extract survival (4 characteristics fields: gender/age/time point/survival) ----
  surv_col <- colnames(pd)[sapply(pd, function(col) any(grepl("survival", col, fixed = TRUE)))][1]
  if (length(surv_col) == 0 || is.na(surv_col)) stop("No survival field found")
  status_raw <- pd[[surv_col]]
  status01 <- to_status01(status_raw)
  time_raw <- rep(28, nrow(pd))   # fixed 28-day observation window

  # ---- Extract age (for multivariable Cox adjustment) ----
  age_col <- colnames(pd)[sapply(pd, function(col) any(grepl("age:", col, fixed = TRUE)))][1]
  age <- if (length(age_col) == 1 && !is.na(age_col)) {
    suppressWarnings(as.numeric(sub("age:\\s*", "", pd[[age_col]])))
  } else NULL

  # ---- Take the Day1 baseline (time point = D01; two longitudinal time points per patient
  #      share the outcome, avoiding duplicated records in Cox) ----
  tp_col <- colnames(pd)[sapply(pd, function(col) any(grepl("time point", col, fixed = TRUE)))][1]
  if (length(tp_col) == 1 && !is.na(tp_col)) {
    tp <- pd[[tp_col]]
    is_day1 <- grepl("D01", tp, fixed = TRUE)
  } else {
    src <- pd[[grep("source_name", colnames(pd), value = TRUE)[1]]]
    is_day1 <- grepl("Day1", src, fixed = TRUE)
  }
  cat("Day1 baseline samples:", sum(is_day1), " (longitudinal Day2/3:", sum(!is_day1), ")\n")

  # Align expression column names (GSM) with pd row names
  common <- intersect(colnames(expr), rownames(pd))
  expr <- expr[, common, drop = FALSE]
  status01 <- status01[match(common, rownames(pd))]
  time_raw <- time_raw[match(common, rownames(pd))]
  is_day1  <- is_day1[match(common, rownames(pd))]
  if (!is.null(age)) age <- age[match(common, rownames(pd))]

  # Keep only Day1 samples with survival data
  keep_samp <- is_day1 & !is.na(status01) & time_raw > 0
  expr <- expr[, keep_samp, drop = FALSE]
  status01 <- status01[keep_samp]; time_raw <- time_raw[keep_samp]
  if (!is.null(age)) age <- age[keep_samp]

  cat("Samples retained with survival data:", ncol(expr), "(dead=", sum(status01), ")\n")
  if (!is.null(age) && sum(!is.na(age)) > 10) {
    cat("  Valid age values:", sum(!is.na(age)), "(usable for multivariable Cox adjustment)\n")
  }
  list(expr = expr, time = time_raw, status = status01, label = label, age = age)
}

## =========================================================================
## Three signatures: Cox + KM + ROC
## =========================================================================
run_three <- function(cohort, out_tag) {
  if (is.null(cohort)) return(NULL)
  expr <- cohort$expr; time <- cohort$time; status <- cohort$status

  score_pdc <- score_geneset(expr, pdc_markers)
  score_sup <- score_geneset(expr, pdc_suppress)
  score_copd <- score_geneset(expr, copd_sig)
  cat("\n[", cohort$label, "] genes matched: pDC", sum(pdc_markers %in% rownames(expr)),
      "| immunoparalysis", sum(pdc_suppress %in% rownames(expr)),
      "| COPD", sum(copd_sig %in% rownames(expr)), "\n")

  df <- data.frame(time = time, status = status,
                   pdc = score_pdc, sup = score_sup, copd = score_copd)
  write.csv(df, file.path(out_dir, paste0("path24_", out_tag, "_scores.csv")),
            row.names = FALSE)

  res_lines <- character()

  # Univariable Cox (three signatures)
  for (nm in c("pdc", "sup", "copd")) {
    if (all(is.na(df[[nm]]))) next
    fit <- tryCatch(coxph(Surv(time, status) ~ df[[nm]], data = df),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      s <- summary(fit)
      hr <- s$conf.int[1, "exp(coef)"]
      lo <- s$conf.int[1, "lower .95"]; hi <- s$conf.int[1, "upper .95"]
      pv <- s$coefficients[1, "Pr(>|z|)"]
      line <- sprintf("[%s] %s univariable Cox: HR=%.3f (%.3f-%.3f), p=%.3e",
                      cohort$label, nm, hr, lo, hi, pv)
      cat(" ", line, "\n"); res_lines <- c(res_lines, line)
    }
  }

  # Multivariable Cox (adjusted for age)
  if (!is.null(cohort$age) && sum(!is.na(cohort$age)) > 50) {
    df$age <- cohort$age
    for (nm in c("pdc", "sup", "copd")) {
      if (all(is.na(df[[nm]]))) next
      df$score <- df[[nm]]
      fit_m <- tryCatch(coxph(Surv(time, status) ~ score + age, data = df),
                        error = function(e) NULL)
      if (!is.null(fit_m)) {
        s <- summary(fit_m)
        hr <- s$conf.int[1, "exp(coef)"]
        pv <- s$coefficients[1, "Pr(>|z|)"]
        line <- sprintf("[%s] %s multivariable Cox (age-adjusted): HR=%.3f, p=%.3e",
                        cohort$label, nm, hr, pv)
        cat(" ", line, "\n"); res_lines <- c(res_lines, line)
      }
    }
  }

  # KM (COPD signature High vs Low, log-rank)
  df$copd_grp <- ifelse(df$copd >= median(df$copd, na.rm = TRUE), "High", "Low")
  fit_km <- survfit(Surv(time, status) ~ copd_grp, data = df)
  lr <- survdiff(Surv(time, status) ~ copd_grp, data = df)
  lr_p <- 1 - pchisq(lr$chisq, df = 1)
  line <- sprintf("[%s] COPD signature High vs Low log-rank p=%.3e", cohort$label, lr_p)
  cat(" ", line, "\n"); res_lines <- c(res_lines, line)

  if (requireNamespace("survminer", quietly = TRUE)) {
    suppressPackageStartupMessages(library(survminer))
    p <- tryCatch(ggsurvplot(fit_km, data = df, pval = TRUE, risk.table = FALSE,
                             legend.title = "COPD signature", palette = c("#d64545","#4c8bf5")),
                  error = function(e) NULL)
    if (!is.null(p)) {
      pdf(file.path(out_dir, paste0("path24_", out_tag, "_KM.pdf")), width = 6, height = 5)
      print(p); dev.off()
    }
  }

  # ROC (COPD signature predicting 28-day death)
  if (requireNamespace("pROC", quietly = TRUE)) {
    suppressPackageStartupMessages(library(pROC))
    roc_obj <- tryCatch(roc(df$status, df$copd, quiet = TRUE), error = function(e) NULL)
    if (!is.null(roc_obj)) {
      line <- sprintf("[%s] COPD signature predicting 28-day death AUC=%.3f", cohort$label,
                      as.numeric(roc_obj$auc))
      cat(" ", line, "\n"); res_lines <- c(res_lines, line)
    }
  }

  writeLines(res_lines, file.path(out_dir, paste0("path24_", out_tag, "_cox_stats.txt")))
  df
}

## =========================================================================
## Main workflow
## =========================================================================
cat("\n===== GSE95233 external prognostic validation =====\n")
gse95233 <- load_cohort("GSE95233", "GSE95233_valid")
run_three(gse95233, "GSE95233_valid")

cat("\nGSE95233 prognosis validation complete\n")
cat("Outputs: path24_*_scores.csv / path24_*_cox_stats.txt / path24_*_KM.pdf\n")
cat("Interpretation: a COPD signature HR>1 (higher score = higher mortality) with a direction consistent with\n")
cat("      GSE65682/E-MTAB-4451 supports the comorbidity hypothesis; the pDC abundance signature is expected to\n")
cat("      be diluted on a whole-blood array, so a non-significant result is normal (as in GSE65682).\n")
