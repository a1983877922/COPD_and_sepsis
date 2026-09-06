#!/usr/bin/env Rscript
# =============================================================================
# 36 three-cohort ROC rebuild
# Sepsis + COPD comorbidity analysis — Fig6c three-cohort ROC (pure R download + manual parsing)
# =============================================================================
# Background: the Fig6c part of script 35 downloaded data online with getGEO, but getGEO no
#             longer works on the server. This script parses manually with curl + gunzip and
#             does not depend on getGEO.
#
# Three cohorts:
#   GSE95233    : read the existing path24_GSE95233_valid_scores.csv (no download needed)
#   E-MTAB-4451 : read the E-MTAB-4451/ subdirectory + download GPL10558.annot.gz (Illumina HT-12)
#   GSE65682    : download the series matrix + download GPL13667 family soft (Affymetrix U219)
#                 note: the GSE65682 platform is GPL13667 (U219, probes 117xxxxx_at), not GPL570
#
# Outputs: Fig6c_ROC_scores.csv / Fig6c_ROC_curve.csv / Fig6c_ROC.pdf/.png
# =============================================================================

## ---- Paths ----
this_file <- commandArgs(trailingOnly = FALSE)
.f <- grep("--file=", this_file, value = TRUE)
if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1]) else .this_file <- "."
out_dir <- dirname(normalizePath(.this_file))
setwd(out_dir)
cat("Working directory:", out_dir, "\n")

## ---- Dependencies ----
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
suppressPackageStartupMessages(library(ggplot2))
has_pROC <- requireNamespace("pROC", quietly = TRUE)
if (has_pROC) suppressPackageStartupMessages(library(pROC))

## ---- COPD signature ----
copd_file <- file.path(out_dir, "path15_copd_signature_genes.txt")
copd_sig  <- if (file.exists(copd_file)) trimws(readLines(copd_file)) else character(0)
copd_sig  <- copd_sig[copd_sig != ""]
cat("COPD signature gene count:", length(copd_sig), "\n")

## ---- Utility functions ----
score_geneset <- function(expr_mat, genes) {
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) < 2) return(rep(NA_real_, ncol(expr_mat)))
  sub <- expr_mat[genes, , drop = FALSE]
  z <- t(scale(t(sub))); z[is.na(z)] <- 0
  colMeans(z)
}

# download .gz + gunzip + read text (curl with -L to follow redirects; corrupted files are re-downloaded automatically)
fetch_gz_lines <- function(url, gz, tag) {
  txt <- sub("\\.gz$", "", gz)
  try_gunzip <- function() {
    if (file.exists(txt) && file.info(txt)$size > 1000) return(TRUE)
    if (!file.exists(gz) || file.info(gz)$size < 1000) return(FALSE)
    code <- system(paste0("gunzip -c \"", gz, "\" > \"", txt, "\""))
    code == 0 && file.exists(txt) && file.info(txt)$size > 1000
  }
  if (try_gunzip()) {
    cat("  [", tag, "] reusing local file:", basename(gz), "\n", sep = "")
  } else {
    # gz missing or corrupted (e.g. unexpected end of file), delete and download again
    if (file.exists(gz)) { cat("  [", tag, "] local file missing/corrupted, downloading again\n", sep = ""); file.remove(gz) }
    if (file.exists(txt)) file.remove(txt)
    cat("  [", tag, "] downloading:", url, "\n", sep = "")
    code <- system(paste0("curl -skL --retry 3 --retry-delay 3 --max-time 900 -o \"", gz, "\" \"", url, "\""))
    if (code != 0 || !file.exists(gz) || file.info(gz)$size < 1000) {
      cat("  [", tag, "] download failed\n", sep = ""); return(NULL)
    }
    if (!try_gunzip()) { cat("  [", tag, "] decompression failed\n", sep = ""); return(NULL) }
  }
  lines <- tryCatch(readLines(txt, warn = FALSE), error = function(e) NULL)
  if (is.null(lines) || length(lines) == 0) {
    cat("  [", tag, "] failed to read the text\n", sep = ""); return(NULL)
  }
  lines
}

# extract the probe -> symbol mapping from a SOFT-format annot file
parse_annot <- function(lines, tag) {
  begin <- grep("!platform_table_begin", lines, fixed = TRUE)
  end   <- grep("!platform_table_end", lines, fixed = TRUE)
  annot <- if (length(begin) == 1 && length(end) == 1) {
    tryCatch(read.delim(text = paste(lines[(begin + 1):(end - 1)], collapse = "\n"),
                        header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE,
                        check.names = FALSE), error = function(e) NULL)
  } else {
    tryCatch(read.delim(text = paste(lines, collapse = "\n"),
                        header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE,
                        check.names = FALSE), error = function(e) NULL)
  }
  if (is.null(annot)) { cat("  [", tag, "] failed to read the annot file\n", sep = ""); return(NULL) }
  id_col  <- intersect(c("ID", "Probe_ID", "probe_id", "ID_REF"), colnames(annot))
  sym_col <- intersect(c("Gene symbol", "GeneSymbol", "Symbol", "Gene Symbol",
                         "GENE_SYMBOL", "ILMN_Gene", "gene_symbol"), colnames(annot))
  if (length(id_col) == 0 || length(sym_col) == 0) {
    cat("  [", tag, "] annot is missing the ID/symbol column, column names:",
        paste(head(colnames(annot), 8), collapse = ", "), "\n", sep = "")
    return(NULL)
  }
  probe2sym <- as.character(annot[[sym_col[1]]])
  names(probe2sym) <- as.character(annot[[id_col[1]]])
  probe2sym <- sub("\\|.*$", "", probe2sym)
  probe2sym <- sub("///.*$", "", probe2sym)
  probe2sym <- trimws(probe2sym)
  probe2sym[is.na(probe2sym) | probe2sym == ""] <- "---"
  cat("  [", tag, "] mapped probes:", length(probe2sym), "\n", sep = "")
  probe2sym
}

collapse_to_symbol <- function(expr, probe2sym) {
  sym <- probe2sym[rownames(expr)]
  sym[is.na(sym)] <- "---"
  keep <- sym != "---"
  mat <- expr[keep, , drop = FALSE]; sym <- sym[keep]
  sym_uniq <- unique(sym)
  out <- matrix(NA, nrow = length(sym_uniq), ncol = ncol(mat),
                dimnames = list(sym_uniq, colnames(mat)))
  for (g in sym_uniq) {
    idx <- which(sym == g)
    out[g, ] <- if (length(idx) == 1) mat[idx, ] else colMeans(mat[idx, , drop = FALSE])
  }
  out
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
  v <- gsub("\\s+", "", v)
  out <- rep(NA_real_, length(v))
  out[v %in% c("dead","death","died","deceased","non-survivor","nonsurvivor",
               "1","yes","event","mortality")] <- 1
  out[v %in% c("alive","survived","survivor","live","0","no","censored","survival")] <- 0
  out
}

## ============================================================================
## Cohort 1: GSE95233 (read the existing scores.csv)
## ============================================================================
get_gse95233 <- function() {
  f <- file.path(out_dir, "path24_GSE95233_valid_scores.csv")
  if (!file.exists(f)) { cat("  [GSE95233] skipped: no scores.csv\n"); return(NULL) }
  d <- read.csv(f, stringsAsFactors = FALSE)
  if (!all(c("status", "copd") %in% colnames(d))) return(NULL)
  cat("  [GSE95233] succeeded:", nrow(d), "samples\n")
  data.frame(cohort = "GSE95233", status = d$status, score = d$copd, stringsAsFactors = FALSE)
}

## ============================================================================
## Cohort 2: E-MTAB-4451 (subdirectory + GPL10558 annot)
## ============================================================================
get_emtab4451 <- function() {
  emtab_dir <- file.path(out_dir, "E-MTAB-4451")
  expr_file <- file.path(emtab_dir, "Davenport_sepsis_Feb2016_normalised_106.txt")
  sdrf_file <- file.path(emtab_dir, "E-MTAB-4451.sdrf.txt")
  if (!file.exists(expr_file) || !file.exists(sdrf_file)) {
    zip_file <- file.path(out_dir, "E-MTAB-4451.zip")
    if (!file.exists(zip_file)) { cat("  [E-MTAB-4451] skipped: no subdirectory and no zip\n"); return(NULL) }
    if (!dir.exists(emtab_dir)) dir.create(emtab_dir)
    cat("  [E-MTAB-4451] unzipping ...\n")
    unzip(zip_file, exdir = emtab_dir)
  }
  if (!file.exists(expr_file) || !file.exists(sdrf_file)) return(NULL)

  sdrf <- read.table(sdrf_file, header = TRUE, sep = "\t", quote = "",
                     comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
  surv_col <- grep("28 day survival", colnames(sdrf), ignore.case = TRUE)
  if (length(surv_col) == 0) { cat("  [E-MTAB-4451] skipped: no survival field in the SDRF\n"); return(NULL) }
  survival_raw <- as.character(sdrf[[surv_col[1]]])
  status <- ifelse(grepl("non.?survivor|dead|died|deceased", survival_raw, ignore.case = TRUE), 1, 0)
  sample_id <- as.character(sdrf[["Source Name"]])
  names(status) <- sample_id

  expr_raw <- read.table(expr_file, header = TRUE, sep = "\t", row.names = 1,
                         check.names = FALSE, stringsAsFactors = FALSE)
  url <- "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz"
  gz  <- file.path(out_dir, "GPL10558.annot.gz")
  lines <- fetch_gz_lines(url, gz, "E-MTAB-4451")
  if (is.null(lines)) { cat("  [E-MTAB-4451] skipped: GPL10558 download failed\n"); return(NULL) }
  probe2sym <- parse_annot(lines, "E-MTAB-4451")
  if (is.null(probe2sym)) { cat("  [E-MTAB-4451] skipped: GPL10558 parsing failed\n"); return(NULL) }
  expr <- collapse_to_symbol(as.matrix(expr_raw), probe2sym)

  common <- intersect(colnames(expr), sample_id)
  expr <- expr[, common, drop = FALSE]; status <- status[common]
  score <- score_geneset(expr, copd_sig)
  cat("  [E-MTAB-4451] succeeded:", length(common), "samples; genes matched",
      sum(copd_sig %in% rownames(expr)), "\n")
  data.frame(cohort = "E-MTAB-4451", status = status, score = score, stringsAsFactors = FALSE)
}

## ============================================================================
## Cohort 3: GSE65682 (series matrix + GPL570 annot)
## ============================================================================
get_gse65682 <- function() {
  url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65682/matrix/GSE65682_series_matrix.txt.gz"
  gz  <- file.path(out_dir, "GSE65682_series_matrix.txt.gz")
  lines <- fetch_gz_lines(url, gz, "GSE65682")
  if (is.null(lines)) { cat("  [GSE65682] skipped: series matrix download failed\n"); return(NULL) }

  cat("  [GSE65682] parsing the series matrix ...\n")
  begin <- grep("!series_matrix_table_begin", lines, fixed = TRUE)
  end   <- grep("!series_matrix_table_end", lines, fixed = TRUE)
  if (length(begin) != 1 || length(end) != 1) {
    cat("  [GSE65682] skipped: abnormal table markers (begin=", length(begin), " end=", length(end), ")\n", sep = "")
    return(NULL)
  }
  tab <- data.table::fread(text = paste(lines[(begin + 1):(end - 1)], collapse = "\n"),
                           sep = "\t", header = TRUE, quote = "",
                           na.strings = c("NA", "null", "NULL", "Null"))
  probe_ids <- gsub("^\"|\"$", "", as.character(tab[[1]]))
  expr <- as.matrix(tab[, -1, drop = FALSE]); storage.mode(expr) <- "double"
  rownames(expr) <- probe_ids
  cat("  [GSE65682] expression matrix:", nrow(expr), "probes x", ncol(expr), "samples\n")

  # sample metadata
  sample_lines <- lines[grepl("!Sample_", lines, fixed = TRUE)]
  sp <- strsplit(sample_lines, "\t", fixed = TRUE)
  fields <- vapply(sp, function(x) sub("!Sample_", "", x[1], fixed = TRUE), character(1))
  fields_uniq <- make.unique(fields)
  val_mat <- do.call(rbind, lapply(sp, function(x) x[-1]))
  gsm <- val_mat[which(fields == "geo_accession")[1], ]
  pd <- as.data.frame(val_mat, stringsAsFactors = FALSE)
  colnames(pd) <- gsm; rownames(pd) <- fields_uniq
  pd <- as.data.frame(t(pd), stringsAsFactors = FALSE)
  pd[] <- lapply(pd, function(x) gsub("^\"|\"$", "", x))

  # GPL13667 (Affymetrix Human Genome U219) annotation — the platform used by GSE65682
  # there is no annot file (404), only the family soft file; U219 probes are 117xxxxx_at;
  # the column name is "Gene Symbol" (capital S)
  url2 <- "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13667/soft/GPL13667_family.soft.gz"
  gz2  <- file.path(out_dir, "GPL13667_family.soft.gz")
  lines2 <- fetch_gz_lines(url2, gz2, "GSE65682")
  if (is.null(lines2)) { cat("  [GSE65682] skipped: GPL13667 download failed\n"); return(NULL) }
  probe2sym <- parse_annot(lines2, "GSE65682")
  if (is.null(probe2sym)) { cat("  [GSE65682] skipped: GPL13667 parsing failed\n"); return(NULL) }
  expr <- collapse_to_symbol(expr, probe2sym)

  # survival
  sf <- find_surv_fields(pd)
  if (length(sf$status) == 0) {
    cat("  [GSE65682] skipped: status not identified; pData column names:",
        paste(head(colnames(pd), 10), collapse = ", "), "\n")
    return(NULL)
  }
  status01 <- to_status01(as.character(pd[[sf$status[1]]]))
  names(status01) <- rownames(pd)

  common <- intersect(colnames(expr), names(status01))
  expr <- expr[, common, drop = FALSE]; status01 <- status01[common]
  keep_s <- !is.na(status01)
  expr <- expr[, keep_s, drop = FALSE]; status01 <- status01[keep_s]

  score <- score_geneset(expr, copd_sig)
  cat("  [GSE65682] succeeded:", sum(keep_s), "samples (dead=", sum(status01),
      "); genes matched", sum(copd_sig %in% rownames(expr)), "\n", sep = "")
  data.frame(cohort = "GSE65682", status = status01, score = score, stringsAsFactors = FALSE)
}

## ============================================================================
## Collect + ROC
## ============================================================================
try_get <- function(tag, fun) {
  res <- tryCatch(fun(), error = function(e) {
    cat("  [", tag, "] error:", conditionMessage(e), "\n", sep = ""); NULL })
  if (is.null(res)) cat("  [", tag, "] no scores obtained\n", sep = "")
  res
}
roc_list <- list()
roc_list$GSE95233  <- try_get("GSE95233", get_gse95233)
roc_list$EMTAB4451 <- try_get("E-MTAB-4451", get_emtab4451)
roc_list$GSE65682  <- try_get("GSE65682", get_gse65682)
roc_list <- roc_list[!vapply(roc_list, is.null, logical(1))]
cat("Cohorts successfully collected:", length(roc_list), "\n")

if (length(roc_list) == 0) stop("scores are unavailable for all three cohorts")

roc_scores <- do.call(rbind, roc_list)
roc_scores <- roc_scores[!is.na(roc_scores$score), ]
write.csv(roc_scores, file.path(out_dir, "Fig6c_ROC_scores.csv"), row.names = FALSE)
cat("Saved: Fig6c_ROC_scores.csv (", nrow(roc_scores), " samples)\n", sep = "")

calc_roc <- function(status, score) {
  if (length(unique(status)) < 2 || any(table(status) < 2)) return(NULL)
  if (has_pROC) {
    r <- pROC::roc(status, score, quiet = TRUE, direction = "auto")
    df <- data.frame(fpr = 1 - r$specificities, tpr = r$sensitivities)
    df <- df[order(df$fpr, df$tpr), ]
    list(df = df, auc = as.numeric(r$auc))
  } else {
    n_pos <- sum(status == 1); n_neg <- sum(status == 0)
    auc <- (sum(rank(score)[status == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    ord <- order(score, decreasing = TRUE); s <- status[ord]
    df <- data.frame(fpr = c(0, cumsum(s == 0) / n_neg, 1),
                     tpr = c(0, cumsum(s == 1) / n_pos, 1))
    list(df = df, auc = auc)
  }
}

curve_list <- lapply(unique(roc_scores$cohort), function(cq) {
  sub <- roc_scores[roc_scores$cohort == cq, ]
  r <- calc_roc(sub$status, sub$score)
  if (is.null(r)) return(NULL)
  data.frame(cohort = cq, fpr = r$df$fpr, tpr = r$df$tpr, auc = r$auc,
             stringsAsFactors = FALSE)
})
curve_list <- curve_list[!vapply(curve_list, is.null, logical(1))]
roc_curve <- do.call(rbind, curve_list)
write.csv(roc_curve, file.path(out_dir, "Fig6c_ROC_curve.csv"), row.names = FALSE)
cat("Saved: Fig6c_ROC_curve.csv\n")

auc_lab <- unique(roc_curve[, c("cohort", "auc")])
auc_lab$label <- sprintf("%s  AUC=%.3f", auc_lab$cohort, auc_lab$auc)
cohort_col <- c(GSE65682 = "#d62728", `E-MTAB-4451` = "#1f77b4", GSE95233 = "#2ca02c")

p <- ggplot(roc_curve, aes(x = fpr, y = tpr, color = cohort)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = cohort_col) +
  coord_equal() +
  labs(x = "1 - Specificity (False positive rate)",
       y = "Sensitivity (True positive rate)",
       title = "COPD signature predicts 28-day mortality", color = "Cohort") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position.inside = c(0.75, 0.25))
for (i in seq_len(nrow(auc_lab))) {
  p <- p + annotate("text", x = 0.6, y = 0.30 - 0.06 * (i - 1),
                    label = auc_lab$label[i], size = 3.4, hjust = 0,
                    color = cohort_col[auc_lab$cohort[i]])
}
ggsave(file.path(out_dir, "Fig6c_ROC.pdf"), p, width = 6, height = 6)
ggsave(file.path(out_dir, "Fig6c_ROC.png"), p, width = 6, height = 6, dpi = 300)
cat("Saved: Fig6c_ROC.pdf/.png\n")

cat("\n36 three-cohort ROC done\n")
