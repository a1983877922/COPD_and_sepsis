# =========================================================================
# 14_pdc_prognosis_cox.R (pDC prognosis Cox)
# 28-day mortality prognosis analysis of the pDC depletion/immunoparalysis
# signature (true survival Cox)
#
# Training cohort: GSE65682 (MARS, Netherlands, GPL13667/HG-U219)
#   760 sepsis + 42 healthy controls, of which 479 sepsis cases have 28-day follow-up
#   (365 alive + 114 dead, mortality 23.8%; follow-up 23.19±9.20 days)
# Validation cohort: GSE134347 (Netherlands, 156 cases, 28-day, 77 dead, GPL17586)
#   (optional) E-MTAB-4451 (UK, 106 cases, 28-day, GPL10558, ArrayExpress)
#
# Hypothesis (from single-cell discovery): blood-side pDC depletion is a marker
#   of sepsis immunoparalysis
#   -> low abundance of pDC marker genes in bulk blood -> higher 28-day mortality
#   -> high pDC immunoparalysis genes (SOCS3/ISG) in bulk blood -> higher 28-day mortality
#
# Dependencies: survival, GEOquery, ggplot2 (required); survminer, pROC (optional, tryCatch)
# =========================================================================

.this_file <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
if (.this_file == "" || !file.exists(.this_file)) {
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1])
}
SCRIPT_DIR  <- dirname(.this_file)
config_file <- file.path(SCRIPT_DIR, "00_config.R")
if (!file.exists(config_file)) stop("00_config.R not found: ", config_file)
source(config_file)

for (pkg in c("GEOquery", "survival")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop(pkg, " not installed. Please run: BiocManager::install('", pkg, "')")
}
suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(ggplot2))

cat("\n==============================================================\n")
cat("pDC prognosis Cox (28-day mortality)\n")
cat("==============================================================\n")

## pDC gene signature (from single-cell discovery)
pdc_markers <- c("LILRA4","CLEC4C","IL3RA","TCF4","IRF7","IRF4","IRF8",
                 "GZMB","SPIB","BCL11A","RUNX2")         # pDC abundance/identity
pdc_suppress <- c("SOCS3","SOCS1","TLR9","ISG20","MX1","ISG15","OAS2",
                  "IFI44L","MYD88","IRAK1")              # pDC immunoparalysis

## =========================================================================
## Utility functions
## =========================================================================

# Geneset z-score scoring (reusing logic from 10)
score_geneset <- function(expr_mat, genes) {
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) < 2) return(rep(NA_real_, ncol(expr_mat)))
  sub <- expr_mat[genes, , drop = FALSE]
  z   <- t(scale(t(sub)))
  z[is.na(z)] <- 0
  colMeans(z)
}

# Automatically identify survival fields (status and time), print candidates for verification
# Critical: outcome-field keywords (mortality/survival/dead...) and time-field
# keywords (time/follow) are mutually exclusive. Note that GSE65682's actual
# fields are "mortality_event_28days" (status) and "time_to_event_28days" (time)
# -- both contain "event", so "event" cannot be used to distinguish them; the
# status field also contains "days" and must not be treated as a time field, so
# we exclude it and keep only time|follow.
find_surv_fields <- function(pd) {
  cols <- colnames(pd)
  evt_pat <- "mortal|surviv|outcome|dead|alive|fate"
  tme_pat <- "time|follow"
  status_cand <- cols[grepl(evt_pat, cols, ignore.case = TRUE) &
                      !grepl(tme_pat, cols, ignore.case = TRUE)]
  time_cand   <- cols[grepl(tme_pat, cols, ignore.case = TRUE) &
                      !grepl(evt_pat, cols, ignore.case = TRUE)]

  # When column names cannot be matched, scan the value content of the characteristics columns
  # (GSE134347's survival info is hidden in the values of characteristics_ch1.x, whose column names contain no keywords)
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
        cat("  [detected] characteristics column", cc, "identified as status (example:",
            paste(head(vals, 3), collapse = "; "), ")\n")
      }
      if (length(time_cand) == 0 && has_tme && !has_evt) {
        time_cand <- c(time_cand, cc)
        cat("  [detected] characteristics column", cc, "identified as time (example:",
            paste(head(vals, 3), collapse = "; "), ")\n")
      }
    }
  }
  list(status = status_cand, time = time_cand)
}

# Convert status values to 0/1 (Alive/survived=0, Dead/died=1)
to_status01 <- function(v) {
  v <- trimws(tolower(as.character(v)))
  # Handle "key: value" format (e.g. "survival: dead" -> "dead")
  v <- sub("^[^:]*:\\s*", "", v)
  out <- rep(NA_real_, length(v))
  out[v %in% c("dead","death","died","deceased","non-survivor","nonsurvivor",
               "1","yes","event","mortality")] <- 1
  out[v %in% c("alive","survived","survivor","live","0","no","censored","survival")] <- 0
  out
}

# Process a cohort into a unified format (expression + survival data)
load_cohort <- function(gse_acc, label) {
  cat("\n----- Processing cohort:", label, "(", gse_acc, ") -----\n")
  # Prefer local series matrix (official filename, also accepts old names)
  matrix_file <- file.path(out_dir, paste0(gse_acc, "_series_matrix.txt.gz"))
  if (!file.exists(matrix_file)) {
    matrix_file <- file.path(out_dir, paste0(gse_acc, "_matrix.txt.gz"))
  }
  if (file.exists(matrix_file)) {
    cat("Using local series matrix:", matrix_file, "\n")
    gse <- getGEO(filename = matrix_file, getGPL = TRUE)
  } else {
    cat("No local matrix file, downloading", gse_acc, "...\n")
    gse <- tryCatch(getGEO(gse_acc, GSEMatrix = TRUE, getGPL = TRUE),
                    error = function(e) NULL)
    if (is.null(gse)) { cat("Download failed, skipping", label, "\n"); return(NULL) }
    gse <- gse[[1]]
  }
  expr <- exprs(gse)
  pd   <- pData(gse)
  fd   <- fData(gse)
  cat("Expression matrix:", nrow(expr), "probes x", ncol(expr), "samples\n")

  # Probe -> symbol
  sym_col <- NULL
  for (cand in c("Gene Symbol","GeneSymbol","Symbol","Gene.Symbol",
                 "gene_symbol","GENE_SYMBOL")) {
    if (cand %in% colnames(fd)) { sym_col <- cand; break }
  }
  if (is.null(sym_col)) {
    m <- grep("symbol", colnames(fd), ignore.case = TRUE)
    if (length(m) > 0) sym_col <- colnames(fd)[m[1]]
  }
  # Still no symbol column -> try the gene_assignment column (Affymetrix HTA/ST
  # platforms, e.g. GPL17586, format "---" or "// symbol // description // ...",
  # take the 2nd segment split by //)
  use_ga <- FALSE
  if (is.null(sym_col) && "gene_assignment" %in% colnames(fd)) {
    sym_col <- "gene_assignment"; use_ga <- TRUE
  }
  if (is.null(sym_col)) { cat("No symbol column found, skipping\n"); return(NULL) }
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

  # Identify survival fields
  sf <- find_surv_fields(pd)
  cat("pData column names:\n  ", paste(colnames(pd), collapse = "\n   "), "\n")
  cat("Candidate status fields:", paste(sf$status, collapse = ", "), "\n")
  cat("Candidate time fields:", paste(sf$time, collapse = ", "), "\n")

  if (length(sf$status) == 0) {
    cat("!! Could not identify status field, please verify the pData column names above and specify manually\n")
    return(NULL)
  }
  status_col <- sf$status[1]
  # time field is optional: if not found, use a fixed 28-day observation window
  # (standard practice for sepsis 28-day mortality)
  if (length(sf$time) == 0) {
    cat("No time field found, using fixed 28-day observation window\n")
    time_col <- NULL
  } else {
    time_col <- sf$time[1]
  }
  cat("Using status field:", status_col, " / time field:",
      if (is.null(time_col)) "fixed 28 days" else time_col, "\n")

  status_raw <- as.character(pd[[status_col]])
  if (is.null(time_col)) {
    time_raw <- rep(28, nrow(pd))
  } else {
    time_raw <- suppressWarnings(as.numeric(as.character(pd[[time_col]])))
  }
  names(status_raw) <- rownames(pd); names(time_raw) <- rownames(pd)

  # Align
  common <- intersect(colnames(expr), rownames(pd))
  expr <- expr[, common, drop = FALSE]
  status_raw <- status_raw[common]; time_raw <- time_raw[common]

  # Identify covariates (age/gender, if fields exist, for multivariate Cox)
  age_col <- grep("age", colnames(pd), ignore.case = TRUE, value = TRUE)
  age_col <- age_col[!grepl("stage|agent|percent", age_col, ignore.case = TRUE)]
  age <- if (length(age_col) > 0) suppressWarnings(as.numeric(as.character(pd[[age_col[1]]][common]))) else NULL
  gender_col <- grep("gender|sex", colnames(pd), ignore.case = TRUE, value = TRUE)
  gender <- if (length(gender_col) > 0) as.character(pd[[gender_col[1]]][common]) else NULL

  status01 <- to_status01(status_raw)
  keep_samp <- !is.na(status01) & !is.na(time_raw) & time_raw > 0
  expr <- expr[, keep_samp, drop = FALSE]
  status01 <- status01[keep_samp]; time_raw <- time_raw[keep_samp]
  if (!is.null(age)) age <- age[keep_samp]
  if (!is.null(gender)) gender <- gender[keep_samp]

  cat("Number of samples with survival data retained:", ncol(expr), "(dead=", sum(status01), ")\n")
  if (!is.null(age)) cat("  age field:", age_col[1],
                         " | gender field:", if (length(gender_col) > 0) gender_col[1] else "none", "\n")

  list(expr = expr, time = time_raw, status = status01, label = label,
       age = age, gender = gender)
}

# Single-cohort Cox + KM + ROC
run_cohort <- function(cohort, out_tag) {
  if (is.null(cohort)) return(NULL)
  expr <- cohort$expr; time <- cohort$time; status <- cohort$status

  score_pdc <- score_geneset(expr, pdc_markers)
  score_sup <- score_geneset(expr, pdc_suppress)
  cat("\n[", cohort$label, "] pDC signature genes hit:",
      sum(pdc_markers %in% rownames(expr)), "/", length(pdc_markers),
      " | immunoparalysis hits:", sum(pdc_suppress %in% rownames(expr)), "\n")

  df <- data.frame(time = time, status = status,
                   pdc = score_pdc, sup = score_sup)

  # Grouping (median)
  df$pdc_grp <- ifelse(df$pdc >= median(df$pdc, na.rm = TRUE), "High", "Low")
  df$sup_grp <- ifelse(df$sup >= median(df$sup, na.rm = TRUE), "High", "Low")

  write.csv(df, file.path(out_dir, paste0("path14_", out_tag, "_scores.csv")),
            row.names = FALSE)

  res_lines <- character()

  # Univariate Cox
  for (nm in c("pdc", "sup")) {
    if (all(is.na(df[[nm]]))) next
    fit <- tryCatch(coxph(Surv(time, status) ~ df[[nm]], data = df),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      s <- summary(fit)
      hr <- s$conf.int[1, "exp(coef)"]
      lo <- s$conf.int[1, "lower .95"]; hi <- s$conf.int[1, "upper .95"]
      pv <- s$coefficients[1, "Pr(>|z|)"]
      line <- sprintf("[%s] %s univariate Cox: HR=%.3f (%.3f-%.3f), p=%.3e",
                      cohort$label, nm, hr, lo, hi, pv)
      cat(" ", line, "\n"); res_lines <- c(res_lines, line)
    }
  }

  # Multivariate Cox (adjust for age if available; also adjust gender if binary)
  if (!is.null(cohort$age) && sum(!is.na(cohort$age)) > 50) {
    df$age <- cohort$age
    for (nm in c("pdc", "sup")) {
      if (all(is.na(df[[nm]]))) next
      df$score <- df[[nm]]
      fit_m <- tryCatch(coxph(Surv(time, status) ~ score + age, data = df),
                        error = function(e) NULL)
      if (!is.null(fit_m)) {
        s <- summary(fit_m)
        hr <- s$conf.int[1, "exp(coef)"]
        pv <- s$coefficients[1, "Pr(>|z|)"]
        line <- sprintf("[%s] %s multivariate Cox (adjusted for age): HR=%.3f, p=%.3e",
                        cohort$label, nm, hr, pv)
        cat(" ", line, "\n"); res_lines <- c(res_lines, line)
      }
    }
  }

  # KM (pDC abundance High vs Low, log-rank)
  fit_km <- survfit(Surv(time, status) ~ pdc_grp, data = df)
  lr <- survdiff(Surv(time, status) ~ pdc_grp, data = df)
  lr_p <- 1 - pchisq(lr$chisq, df = 1)
  line <- sprintf("[%s] pDC High vs Low log-rank p=%.3e", cohort$label, lr_p)
  cat(" ", line, "\n"); res_lines <- c(res_lines, line)

  # KM plot (survminer optional)
  if (requireNamespace("survminer", quietly = TRUE)) {
    suppressPackageStartupMessages(library(survminer))
    p <- tryCatch(ggsurvplot(fit_km, data = df, pval = TRUE, risk.table = FALSE,
                             legend.title = "pDC signature", palette = c("#d64545","#4c8bf5")),
                  error = function(e) NULL)
    if (!is.null(p)) {
      pdf(file.path(out_dir, paste0("path14_", out_tag, "_KM.pdf")), width = 6, height = 5)
      print(p); dev.off()
    }
  }

  # ROC (pDC signature predicts 28-day death)
  if (requireNamespace("pROC", quietly = TRUE)) {
    suppressPackageStartupMessages(library(pROC))
    roc_obj <- tryCatch(roc(df$status, df$pdc, quiet = TRUE), error = function(e) NULL)
    if (!is.null(roc_obj)) {
      line <- sprintf("[%s] pDC signature predicts 28-day death AUC=%.3f", cohort$label,
                      as.numeric(roc_obj$auc))
      cat(" ", line, "\n"); res_lines <- c(res_lines, line)
    }
  }

  writeLines(res_lines, file.path(out_dir, paste0("path14_", out_tag, "_cox_stats.txt")))
  df
}

## =========================================================================
## Training cohort: GSE65682
## =========================================================================
cat("\n===== Training cohort GSE65682 =====\n")
train <- load_cohort("GSE65682", "GSE65682_train")
run_cohort(train, "GSE65682_train")

## =========================================================================
## Internal validation: GSE65682 random split (train 70% + test 30%)
## =========================================================================
# GSE134347's survival data is not in the GEO series matrix (characteristics
# columns are all age), so we use GSE65682 internal random split validation
# instead (standard practice in sepsis prognosis literature; MARS itself has
# discovery/validation splits)
if (!is.null(train) && ncol(train$expr) >= 50) {
  cat("\n===== Internal validation: GSE65682 random split =====\n")
  set.seed(123)
  n <- ncol(train$expr)
  idx <- sample(n, round(0.7 * n))
  make_split <- function(idx_vec, lbl) {
    list(expr = train$expr[, idx_vec, drop = FALSE],
         time = train$time[idx_vec], status = train$status[idx_vec],
         label = lbl,
         age = if (!is.null(train$age)) train$age[idx_vec] else NULL)
  }
  run_cohort(make_split(idx, "GSE65682_discovery"), "GSE65682_discovery")
  run_cohort(make_split(setdiff(seq_len(n), idx), "GSE65682_validation"), "GSE65682_validation")
}

cat("\npDC prognosis Cox analysis complete\n")
cat("Outputs: path14_*_scores.csv / path14_*_cox_stats.txt / path14_*_KM.pdf\n")
cat("Interpretation: if the pDC abundance signature has HR<1 (high score = better survival) or the immunoparalysis signature has HR>1 (high score = higher death),\n")
cat("      this supports the single-cell finding that 'pDC depletion/immunoparalysis -> higher 28-day mortality'\n")
cat("Note: GSE134347 survival is not on GEO; switched to GSE65682 internal split validation;\n")
cat("    for true external validation, manually curate GSE134347's original paper supplement table, or use ArrayExpress E-MTAB-4451.\n")
