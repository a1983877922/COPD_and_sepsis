# =========================================================================
# 16_emtab4451_external_validation.R (E-MTAB-4451 external validation)
# Use E-MTAB-4451 (Davenport 2016, Lancet Resp Med; UK ICU sepsis, Illumina
# HumanHT-12 v4 = GPL10558) for TRUE external validation of the pDC signature +
# COPD inflammatory pre-activation signature.
#
# Input (user has downloaded, place in the E-MTAB-4451/ subfolder in this
# script's directory):
#   E-MTAB-4451.sdrf.txt                       clinical info (includes 28 day survival)
#   Davenport_sepsis_Feb2016_normalised_106.txt  normalized expression (106 samples, ILMN probes)
#   E-MTAB-4451.idf.txt                        experiment description (not used by this script)
#
# Key points:
#   - SDRF column "Characteristics[28 day survival]": survivor / non survivor
#   - Expression matrix first column Probe_ID (ILMN_xxx), needs GPL10558
#     annotation mapping to gene symbol
#   - No follow-up time -> use fixed 28-day observation window (time=28)
#
# Dependencies: GEOquery (download GPL10558 annotation), survival, ggplot2;
#   survminer/pROC optional
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

EMTAB_DIR <- file.path(SCRIPT_DIR, "E-MTAB-4451")
if (!dir.exists(EMTAB_DIR)) stop("E-MTAB-4451 directory not found: ", EMTAB_DIR)

cat("\n==============================================================\n")
cat("E-MTAB-4451 external validation (pDC + COPD signature -> 28-day mortality)\n")
cat("==============================================================\n")

## pDC signature (same as 14)
pdc_markers <- c("LILRA4","CLEC4C","IL3RA","TCF4","IRF7","IRF4","IRF8",
                 "GZMB","SPIB","BCL11A","RUNX2")
pdc_suppress <- c("SOCS3","SOCS1","TLR9","ISG20","MX1","ISG15","OAS2",
                  "IFI44L","MYD88","IRAK1")
# COPD signature (if 15 has run, read its output)
copd_file <- file.path(out_dir, "path15_copd_signature_genes.txt")
copd_sig <- if (file.exists(copd_file)) readLines(copd_file) else NULL
copd_sig <- copd_sig[copd_sig != ""]
cat("COPD signature gene count:", length(copd_sig), "\n")

score_geneset <- function(expr_mat, genes) {
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) < 2) return(rep(NA_real_, ncol(expr_mat)))
  sub <- expr_mat[genes, , drop = FALSE]
  z   <- t(scale(t(sub)))
  z[is.na(z)] <- 0
  colMeans(z)
}

## =========================================================================
## Step 1: Read SDRF clinical info
## =========================================================================
cat("\n===== Step 1: Read SDRF clinical info =====\n")
sdrf <- read.table(file.path(EMTAB_DIR, "E-MTAB-4451.sdrf.txt"),
                   header = TRUE, sep = "\t", quote = "", comment.char = "",
                   stringsAsFactors = FALSE, check.names = FALSE)
cat("SDRF sample count:", nrow(sdrf), "\n")

surv_col <- grep("28 day survival", colnames(sdrf), ignore.case = TRUE)
if (length(surv_col) == 0) stop("'28 day survival' column not found in SDRF")
survival_raw <- as.character(sdrf[[surv_col[1]]])
status <- ifelse(grepl("non.?survivor|dead|died|deceased", survival_raw, ignore.case = TRUE), 1, 0)
cat("28 day survival distribution:\n"); print(table(survival_raw, useNA = "ifany"))

age_col <- grep("age", colnames(sdrf), ignore.case = TRUE)
age <- if (length(age_col) > 0) suppressWarnings(as.numeric(as.character(sdrf[[age_col[1]]]))) else NULL
sex_col <- grep("sex", colnames(sdrf), ignore.case = TRUE)
sex <- if (length(sex_col) > 0) as.character(sdrf[[sex_col[1]]]) else NULL

sample_id <- as.character(sdrf[["Source Name"]])
names(status) <- sample_id
time <- rep(28, length(status)); names(time) <- sample_id
if (!is.null(age)) names(age) <- sample_id

## =========================================================================
## Step 2: Read expression matrix
## =========================================================================
cat("\n===== Step 2: Read expression matrix =====\n")
expr_file <- file.path(EMTAB_DIR, "Davenport_sepsis_Feb2016_normalised_106.txt")
expr_raw <- read.table(expr_file, header = TRUE, sep = "\t", row.names = 1,
                       check.names = FALSE, stringsAsFactors = FALSE)
cat("Expression matrix:", nrow(expr_raw), "probes x", ncol(expr_raw), "samples\n")

## =========================================================================
## Step 3: Download GPL10558 annotation, probe -> symbol
## =========================================================================
cat("\n===== Step 3: GPL10558 annotation, probe -> symbol =====\n")
gpl <- getGEO("GPL10558")
gpl_table <- Table(gpl)
cat("GPL10558 annotation column names:", paste(colnames(gpl_table), collapse = ", "), "\n")

sym_col <- NULL
for (cand in c("Symbol","ILMN_Gene","Gene symbol","GeneSymbol","gene_symbol",
               "GENE_SYMBOL","Gene.Symbol")) {
  if (cand %in% colnames(gpl_table)) { sym_col <- cand; break }
}
if (is.null(sym_col)) {
  m <- grep("symbol|gene", colnames(gpl_table), ignore.case = TRUE)
  if (length(m) > 0) sym_col <- colnames(gpl_table)[m[1]]
}
if (is.null(sym_col)) stop("symbol column not found in GPL10558 annotation, please specify manually")
cat("Using symbol column:", sym_col, "\n")

probe2sym <- as.character(gpl_table[[sym_col]])
names(probe2sym) <- as.character(gpl_table$ID)
# An Illumina probe may map to multiple genes (separated by | or ///), take the first
probe2sym <- sub("\\|.*$", "", probe2sym)
probe2sym <- sub("///.*$", "", probe2sym)
probe2sym <- trimws(probe2sym)

# Map the expression matrix probe rows to gene symbol (take mean for multiple probes)
probes <- rownames(expr_raw)
sym <- probe2sym[probes]
sym[is.na(sym)] <- ""; sym[sym == ""] <- "---"
keep <- sym != "---"
mat <- as.matrix(expr_raw[keep, , drop = FALSE]); sym <- sym[keep]
sym_uniq <- unique(sym)
expr <- matrix(NA, nrow = length(sym_uniq), ncol = ncol(mat),
               dimnames = list(sym_uniq, colnames(mat)))
for (g in sym_uniq) {
  idx <- which(sym == g)
  expr[g, ] <- if (length(idx) == 1) mat[idx, ] else colMeans(mat[idx, , drop = FALSE])
}
cat("After mapping:", nrow(expr), "genes x", ncol(expr), "samples\n")

## =========================================================================
## Step 4: Sample alignment + scoring
## =========================================================================
cat("\n===== Step 4: Sample alignment + scoring =====\n")
common <- intersect(colnames(expr), sample_id)
expr <- expr[, common, drop = FALSE]
status <- status[common]; time <- time[common]
if (!is.null(age)) age <- age[common]
cat("Number of samples after alignment:", length(common), "(dead=", sum(status), ")\n")

scores <- data.frame(sample = common, time = time, status = status,
                     pdc = score_geneset(expr, pdc_markers),
                     sup = score_geneset(expr, pdc_suppress),
                     stringsAsFactors = FALSE)
if (length(copd_sig) >= 2) {
  scores$copd <- score_geneset(expr, copd_sig)
} else {
  scores$copd <- NA_real_
}
if (!is.null(age)) scores$age <- age
write.csv(scores, file.path(out_dir, "path16_EMTAB4451_scores.csv"), row.names = FALSE)

cat("Signature genes hit: pDC=", sum(pdc_markers %in% rownames(expr)),
    " immunoparalysis=", sum(pdc_suppress %in% rownames(expr)),
    " COPD=", if (length(copd_sig) >= 2) sum(copd_sig %in% rownames(expr)) else 0, "\n")

## =========================================================================
## Step 5: Cox + KM + ROC
## =========================================================================
cat("\n===== Step 5: Cox + KM + ROC =====\n")
res_lines <- character()

for (nm in c("pdc", "sup", "copd")) {
  if (all(is.na(scores[[nm]]))) next
  df <- data.frame(time = scores$time, status = scores$status, score = scores[[nm]])
  fit <- tryCatch(coxph(Surv(time, status) ~ score, data = df), error = function(e) NULL)
  if (!is.null(fit)) {
    s <- summary(fit)
    hr <- s$conf.int[1, "exp(coef)"]
    lo <- s$conf.int[1, "lower .95"]; hi <- s$conf.int[1, "upper .95"]
    pv <- s$coefficients[1, "Pr(>|z|)"]
    line <- sprintf("[E-MTAB-4451] %s univariate Cox: HR=%.3f (%.3f-%.3f), p=%.3e",
                    nm, hr, lo, hi, pv)
    cat(" ", line, "\n"); res_lines <- c(res_lines, line)
  }
}

# Multivariate Cox (adjust for age)
if (!is.null(scores$age) && sum(!is.na(scores$age)) > 50) {
  for (nm in c("pdc", "sup", "copd")) {
    if (all(is.na(scores[[nm]]))) next
    df <- data.frame(time = scores$time, status = scores$status,
                     score = scores[[nm]], age = scores$age)
    fit_m <- tryCatch(coxph(Surv(time, status) ~ score + age, data = df),
                      error = function(e) NULL)
    if (!is.null(fit_m)) {
      s <- summary(fit_m)
      hr <- s$conf.int[1, "exp(coef)"]
      pv <- s$coefficients[1, "Pr(>|z|)"]
      line <- sprintf("[E-MTAB-4451] %s multivariate Cox (adjusted for age): HR=%.3f, p=%.3e",
                      nm, hr, pv)
      cat(" ", line, "\n"); res_lines <- c(res_lines, line)
    }
  }
}

# KM (pDC abundance, High vs Low)
df_km <- data.frame(time = scores$time, status = scores$status, pdc = scores$pdc)
df_km$grp <- ifelse(df_km$pdc >= median(df_km$pdc, na.rm = TRUE), "High", "Low")
fit_km <- survfit(Surv(time, status) ~ grp, data = df_km)
lr <- survdiff(Surv(time, status) ~ grp, data = df_km)
lr_p <- 1 - pchisq(lr$chisq, df = 1)
line <- sprintf("[E-MTAB-4451] pDC High vs Low log-rank p=%.3e", lr_p)
cat(" ", line, "\n"); res_lines <- c(res_lines, line)

if (requireNamespace("survminer", quietly = TRUE)) {
  suppressPackageStartupMessages(library(survminer))
  p <- tryCatch(ggsurvplot(fit_km, data = df_km, pval = TRUE, risk.table = FALSE,
                           legend.title = "pDC signature", palette = c("#d64545","#4c8bf5")),
                error = function(e) NULL)
  if (!is.null(p)) {
    pdf(file.path(out_dir, "path16_EMTAB4451_KM.pdf"), width = 6, height = 5)
    print(p); dev.off()
  }
}

if (requireNamespace("pROC", quietly = TRUE)) {
  suppressPackageStartupMessages(library(pROC))
  for (nm in c("pdc", "sup", "copd")) {
    if (all(is.na(scores[[nm]]))) next
    roc_obj <- tryCatch(roc(scores$status, scores[[nm]], quiet = TRUE), error = function(e) NULL)
    if (!is.null(roc_obj)) {
      line <- sprintf("[E-MTAB-4451] %s predicts 28-day death AUC=%.3f", nm, as.numeric(roc_obj$auc))
      cat(" ", line, "\n"); res_lines <- c(res_lines, line)
    }
  }
}

writeLines(res_lines, file.path(out_dir, "path16_EMTAB4451_cox_stats.txt"))

cat("\nE-MTAB-4451 external validation complete\n")
cat("Outputs: path16_EMTAB4451_scores.csv / path16_EMTAB4451_cox_stats.txt / path16_EMTAB4451_KM.pdf\n")
cat("Interpretation: if the pDC abundance signature has HR<1 (high score = better survival) or the COPD/immunoparalysis signature has HR>1 (high score = higher death),\n")
cat("      an independent cohort validates that 'pDC depletion / COPD inflammatory pre-activation -> higher 28-day mortality'.\n")
