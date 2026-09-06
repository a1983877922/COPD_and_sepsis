#!/usr/bin/env Rscript
# =============================================================================
# 35_auxiliary_panels.R
# Sepsis + COPD comorbidity analysis — three missing panels (Fig4a / Fig5b / Fig6c)
# =============================================================================
# Figure structure (6 main figures):
#   Fig1 workflow | Fig2 single-cell annotation | Fig3 149 shared programs | Fig4 SOCS3 |
#   Fig5 interferon WGCNA | Fig6 prognosis
#
# The 3 panels this script fills in:
#   Fig4a  SOCS3 donor-level logFC barplot (monotonically increasing: IC < COPD < Sepsis < SP)
#   Fig5b  149 shared programs -> hdWGCNA module mapping barplot (enriched in the purple interferon module)
#   Fig6c  Three-cohort COPD signature ROC curves (predicting 28-day mortality)
#
# Each figure saves a clean CSV data table BEFORE plotting, for easy local edits
# (data-driven re-plotting):
#   Fig4a_SOCS3_logFC.csv         SOCS3 logFC/FDR/PValue per group
#   Fig5b_149_module_mapping.csv  149 genes -> module stats + gene lists
#   Fig6c_ROC_scores.csv          three cohorts (cohort, status, score)
#   Fig6c_ROC_curve.csv           ROC curve point coordinates (FPR/TPR), editable for re-plotting
#
# Depends on: ggplot2 (required); pROC / GEOquery / survival (only needed to rebuild
#             Fig6c cohorts; degrades gracefully if missing)
# Note: this script only reads CSVs / local raw data and does not source 00_config.R
#       (to avoid the heavy Seurat dependency).
# =============================================================================

## ---- Path location (same as script 28, no 00_config dependency) ----
this_file <- commandArgs(trailingOnly = FALSE)
.f <- grep("--file=", this_file, value = TRUE)
if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1]) else .this_file <- "."
out_dir <- dirname(normalizePath(.this_file))
setwd(out_dir)
cat("Working directory:", out_dir, "\n")

## ---- Dependencies ----
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(ggplot2))

has_pROC      <- requireNamespace("pROC", quietly = TRUE)
has_GEOquery  <- requireNamespace("GEOquery", quietly = TRUE)
has_survival  <- requireNamespace("survival", quietly = TRUE)
if (has_pROC)      suppressPackageStartupMessages(library(pROC))
if (has_GEOquery)  suppressPackageStartupMessages(library(GEOquery))
if (has_survival)  suppressPackageStartupMessages(library(survival))

## ---- Utility: significance stars ----
p_stars <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
}

## ---- Utility: dual-format save (PDF for submission + PNG preview) ----
save_both <- function(p, base, w = 6, h = 5) {
  ggsave(paste0(base, ".pdf"), p, width = w, height = h, device = cairo_pdf)
  ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = 300)
  cat("  saved:", base, ".pdf / .png\n", sep = "")
}

# =============================================================================
# Part 1: Fig4a  SOCS3 donor-level logFC barplot
# Data source: SOCS3 rows from the 4 pseudobulk edgeR CSVs of script 28
# =============================================================================
cat("\n########## Part 1: Fig4a SOCS3 logFC barplot ##########\n")

socs3_files <- c(
  Infection_Control = "path28_blood_mono_pseudobulk_Infection_Control_vs_Healthy.csv",
  COPD              = "path28_blood_mono_pseudobulk_COPD_vs_Healthy.csv",
  Sepsis            = "path28_blood_mono_pseudobulk_Sepsis_vs_Healthy.csv",
  Sepsis_Pneumonia  = "path28_blood_mono_pseudobulk_Sepsis_Pneumonia_vs_Healthy.csv"
)

extract_socs3 <- function(f, grp) {
  fp <- file.path(out_dir, f)
  if (!file.exists(fp)) return(NULL)
  tab <- read.csv(fp, stringsAsFactors = FALSE)
  r <- tab[tab$gene == "SOCS3", , drop = FALSE]
  if (nrow(r) == 0) return(NULL)
  data.frame(group = grp, logFC = r$logFC[1], logCPM = r$logCPM[1],
             F = r$F[1], PValue = r$PValue[1], FDR = r$FDR[1],
             stringsAsFactors = FALSE)
}

socs3_list <- lapply(names(socs3_files), function(g) extract_socs3(socs3_files[[g]], g))
socs3_list <- socs3_list[!vapply(socs3_list, is.null, logical(1))]
socs3_df   <- do.call(rbind, socs3_list)

# Donor counts (aggregated from path29, for x-axis labels; only blood tissue=="blood"; omitted if unreadable)
donor_n <- NULL
if (file.exists(file.path(out_dir, "path29_dataset_donor_counts.csv"))) {
  d29 <- read.csv(file.path(out_dir, "path29_dataset_donor_counts.csv"), stringsAsFactors = FALSE)
  if (all(c("tissue", "group", "n_donors") %in% colnames(d29))) {
    d29b <- d29[d29$tissue == "blood", , drop = FALSE]
    donor_n <- tapply(d29b$n_donors, d29b$group, sum)
  }
}

# Group order (by severity) + labels + colors (red gradient)
socs3_df$group <- factor(socs3_df$group,
                         levels = c("Infection_Control", "COPD", "Sepsis", "Sepsis_Pneumonia"))
lab <- c(Infection_Control = "IC", COPD = "COPD", Sepsis = "Sepsis",
         Sepsis_Pneumonia = "Sepsis+PNA")
socs3_df$label <- lab[as.character(socs3_df$group)]
if (!is.null(donor_n)) {
  nn <- donor_n[as.character(socs3_df$group)]
  socs3_df$label <- paste0(socs3_df$label, "\n(n=", nn, ")")
}
socs3_df$stars <- p_stars(socs3_df$FDR)

# Save data before plotting
write.csv(socs3_df, file.path(out_dir, "Fig4a_SOCS3_logFC.csv"), row.names = FALSE)
cat("saved data: Fig4a_SOCS3_logFC.csv\n")
print(socs3_df)

# Barplot (logFC monotonically increasing, deepening red gradient)
p1 <- ggplot(socs3_df, aes(x = group, y = logFC, fill = group)) +
  geom_col(width = 0.65, color = "grey20", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", logFC)), vjust = -0.4, size = 3.2) +
  geom_text(aes(label = stars), vjust = 0.8, color = "white", fontface = "bold", size = 4) +
  scale_fill_manual(values = c("Infection_Control" = "#fcbba1", "COPD" = "#fb6a4a",
                               "Sepsis" = "#ef3b2c", "Sepsis_Pneumonia" = "#a50f15"),
                    guide = "none") +
  scale_x_discrete(labels = setNames(socs3_df$label, socs3_df$group)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = expression("SOCS3 log"[2]*"FC (vs Healthy)"),
       title = "SOCS3 — donor-level pseudobulk edgeR") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(face = "bold", size = 10))
save_both(p1, file.path(out_dir, "Fig4a_SOCS3_logFC_barplot"), w = 5, h = 5)

# =============================================================================
# Part 2: Fig5b  149 shared programs -> hdWGCNA module mapping
# Data source: path13_sc_149gene_module.csv (149 genes -> module/color)
# =============================================================================
cat("\n########## Part 2: Fig5b 149 -> hdWGCNA module ##########\n")

mod_file <- file.path(out_dir, "path13_sc_149gene_module.csv")
if (!file.exists(mod_file)) {
  cat("!! cannot find", mod_file, "; skipping Part 2\n")
} else {
  mod_df <- read.csv(mod_file, stringsAsFactors = FALSE)
  mod_df <- mod_df[mod_df$gene_name != "", ]
  # Genes per module
  mod_count <- as.data.frame(table(mod_df$module), stringsAsFactors = FALSE)
  names(mod_count) <- c("module", "n_genes")
  # Module colors (take the first-appearing color per module)
  mod_color <- mod_df[!duplicated(mod_df$module), c("module", "color")]
  mod_count <- merge(mod_count, mod_color, by = "module", all.x = TRUE)
  mod_count$color[is.na(mod_count$color)] <- "grey"
  # Gene lists (comma-joined)
  mod_genes <- aggregate(gene_name ~ module, data = mod_df,
                         FUN = function(x) paste(x, collapse = ","))
  mod_count <- merge(mod_count, mod_genes, by = "module", all.x = TRUE)
  # Order: non-grey first by descending gene count, grey (unassigned) last
  mod_count$is_grey <- mod_count$module == "grey"
  mod_count <- mod_count[order(mod_count$is_grey, -mod_count$n_genes), ]
  mod_count$module <- factor(mod_count$module, levels = mod_count$module)

  # Save data before plotting
  write.csv(mod_count, file.path(out_dir, "Fig5b_149_module_mapping.csv"), row.names = FALSE)
  cat("saved data: Fig5b_149_module_mapping.csv\n")
  print(mod_count[, c("module", "n_genes", "color")])

  # Highlight the interferon module (purple)
  mod_count$highlight <- ifelse(mod_count$module == "purple", "Interferon", "Other")

  p2 <- ggplot(mod_count, aes(x = module, y = n_genes, fill = color)) +
    geom_col(width = 0.7, color = "grey20", linewidth = 0.3) +
    geom_text(aes(label = n_genes), vjust = -0.4, size = 3.2) +
    geom_col(data = mod_count[mod_count$highlight == "Interferon", ],
             aes(x = module, y = n_genes), fill = NA, color = "black",
             linewidth = 0.8, width = 0.72) +
    scale_fill_identity() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "hdWGCNA module", y = "Number of 149 shared genes",
         title = "149 shared myeloid program → hdWGCNA modules",
         subtitle = "Enriched in purple (interferon) module") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
  save_both(p2, file.path(out_dir, "Fig5b_149_module_barplot"), w = 8, h = 5)
}

# =============================================================================
# Part 3: Fig6c  Three-cohort COPD signature ROC curves (predicting 28-day mortality)
# Data source: COPD inflammation signature scores from GSE65682 / E-MTAB-4451 / GSE95233
#   Priority: read existing scores.csv directly; if missing, rebuild from local raw
#   data (with tryCatch error tolerance)
# =============================================================================
cat("\n########## Part 3: Fig6c three-cohort ROC ##########\n")

copd_file <- file.path(out_dir, "path15_copd_signature_genes.txt")
copd_sig  <- if (file.exists(copd_file)) readLines(copd_file) else character(0)
copd_sig  <- trimws(copd_sig[copd_sig != ""])
cat("COPD signature gene count:", length(copd_sig), "\n")

# ---- Scoring function (z-score + mean, same convention as scripts 15/24) ----
score_geneset <- function(expr_mat, genes) {
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) < 2) return(rep(NA_real_, ncol(expr_mat)))
  sub <- expr_mat[genes, , drop = FALSE]
  z   <- t(scale(t(sub)))
  z[is.na(z)] <- 0
  colMeans(z)
}

# ---- Survival field detection (reuses scripts 15/24 find_surv_fields / to_status01) ----
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

# ---- Cohort 1: GSE95233 (existing scores.csv) ----
get_gse95233 <- function() {
  f <- file.path(out_dir, "path24_GSE95233_valid_scores.csv")
  if (!file.exists(f)) { cat("  [GSE95233] skipped: no path24_GSE95233_valid_scores.csv\n"); return(NULL) }
  d <- read.csv(f, stringsAsFactors = FALSE)
  if (!all(c("status", "copd") %in% colnames(d))) {
    cat("  [GSE95233] skipped: scores.csv lacks status/copd columns\n"); return(NULL)
  }
  cat("  [GSE95233] success:", nrow(d), "samples\n")
  data.frame(cohort = "GSE95233", status = d$status, score = d$copd,
             stringsAsFactors = FALSE)
}

# ---- Cohort 2: E-MTAB-4451 (read subdirectory, reuses script 16 logic; unzip if absent) ----
get_emtab4451 <- function() {
  if (!has_GEOquery) { cat("  [E-MTAB-4451] skipped: no GEOquery\n"); return(NULL) }
  emtab_dir <- file.path(out_dir, "E-MTAB-4451")
  expr_file <- file.path(emtab_dir, "Davenport_sepsis_Feb2016_normalised_106.txt")
  sdrf_file <- file.path(emtab_dir, "E-MTAB-4451.sdrf.txt")
  # Prefer the unzipped subdirectory from script 16; try unzipping the zip otherwise
  if (!file.exists(expr_file) || !file.exists(sdrf_file)) {
    zip_file <- file.path(out_dir, "E-MTAB-4451.zip")
    if (!file.exists(zip_file)) {
      cat("  [E-MTAB-4451] skipped: no E-MTAB-4451/ subdirectory and no zip\n")
      return(NULL)
    }
    if (!dir.exists(emtab_dir)) dir.create(emtab_dir)
    cat("  [E-MTAB-4451] unzipping zip ...\n")
    unzip(zip_file, exdir = emtab_dir)
  }
  if (!file.exists(expr_file) || !file.exists(sdrf_file)) {
    cat("  [E-MTAB-4451] skipped: files still missing after unzip\n"); return(NULL)
  }
  cat("  [E-MTAB-4451] reading SDRF + expression matrix ...\n")

  # SDRF survival
  sdrf <- read.table(sdrf_file, header = TRUE, sep = "\t", quote = "",
                     comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
  surv_col <- grep("28 day survival", colnames(sdrf), ignore.case = TRUE)
  if (length(surv_col) == 0) { cat("  [E-MTAB-4451] skipped: SDRF has no survival column\n"); return(NULL) }
  survival_raw <- as.character(sdrf[[surv_col[1]]])
  status <- ifelse(grepl("non.?survivor|dead|died|deceased", survival_raw, ignore.case = TRUE), 1, 0)
  sample_id <- as.character(sdrf[["Source Name"]])
  names(status) <- sample_id

  # Expression matrix (probe rows)
  expr_raw <- read.table(expr_file, header = TRUE, sep = "\t", row.names = 1,
                         check.names = FALSE, stringsAsFactors = FALSE)

  # GPL10558 probes -> symbols
  gpl <- tryCatch(getGEO("GPL10558"), error = function(e) NULL)
  if (is.null(gpl)) { cat("  [E-MTAB-4451] skipped: GPL10558 annotation download failed\n"); return(NULL) }
  gtab <- tryCatch(Table(gpl), error = function(e) NULL)
  if (is.null(gtab)) { cat("  [E-MTAB-4451] skipped: GPL10558 Table failed\n"); return(NULL) }
  sym_col <- intersect(c("Symbol", "ILMN_Gene", "Gene symbol", "GeneSymbol",
                         "gene_symbol", "GENE_SYMBOL", "Gene.Symbol"), colnames(gtab))
  if (length(sym_col) == 0) { cat("  [E-MTAB-4451] skipped: GPL10558 has no symbol column\n"); return(NULL) }
  probe2sym <- as.character(gtab[[sym_col[1]]])
  names(probe2sym) <- as.character(gtab$ID)
  probe2sym <- sub("\\|.*$", "", probe2sym)
  probe2sym <- sub("///.*$", "", probe2sym)
  probe2sym <- trimws(probe2sym)

  sym <- probe2sym[rownames(expr_raw)]
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

  common <- intersect(colnames(expr), sample_id)
  expr   <- expr[, common, drop = FALSE]
  status <- status[common]

  score <- score_geneset(expr, copd_sig)
  cat("  [E-MTAB-4451] success:", length(common), "samples, matched genes",
      sum(copd_sig %in% rownames(expr)), "\n")
  data.frame(cohort = "E-MTAB-4451", status = status, score = score,
             stringsAsFactors = FALSE)
}

# ---- Cohort 3: GSE65682 (getGEO online download, reuses script 15 logic) ----
get_gse65682 <- function() {
  if (!has_GEOquery) { cat("  [GSE65682] skipped: no GEOquery\n"); return(NULL) }
  # Prefer local series matrix; otherwise getGEO online download (script 15 style, no on-disk file dependency)
  matrix_file <- file.path(out_dir, "GSE65682_series_matrix.txt.gz")
  if (!file.exists(matrix_file)) matrix_file <- file.path(out_dir, "GSE65682_matrix.txt.gz")
  if (file.exists(matrix_file)) {
    gse <- tryCatch(getGEO(filename = matrix_file, getGPL = TRUE), error = function(e) NULL)
  } else {
    cat("  [GSE65682] getGEO online download (GSEMatrix) ...\n")
    gse <- tryCatch(getGEO("GSE65682", GSEMatrix = TRUE, getGPL = TRUE),
                    error = function(e) NULL)
    if (!is.null(gse)) gse <- gse[[1]]
  }
  if (is.null(gse)) { cat("  [GSE65682] skipped: getGEO download/parse failed\n"); return(NULL) }
  expr <- exprs(gse); pd <- pData(gse); fd <- fData(gse)

  # Probes -> symbols (reuses script 15 mapping logic)
  sym_col <- intersect(c("Gene Symbol", "GeneSymbol", "Symbol", "Gene.Symbol",
                         "gene_symbol", "GENE_SYMBOL"), colnames(fd))
  use_ga <- FALSE
  if (length(sym_col) == 0 && "gene_assignment" %in% colnames(fd)) {
    sym_col <- "gene_assignment"; use_ga <- TRUE
  }
  if (length(sym_col) == 0) { cat("  [GSE65682] skipped: no symbol column\n"); return(NULL) }
  sym <- as.character(fd[[sym_col[1]]])
  if (use_ga) {
    sym <- vapply(strsplit(sym, "//", fixed = TRUE),
                  function(x) if (length(x) >= 2) trimws(x[2]) else "", character(1))
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

  # Survival fields (fully reuses script 15's find_surv_fields + to_status01)
  sf <- find_surv_fields(pd)
  if (length(sf$status) == 0) {
    cat("  [GSE65682] skipped: status field not recognized; pData column names:",
        paste(head(colnames(pd), 8), collapse = ", "), "\n")
    return(NULL)
  }
  status_raw <- as.character(pd[[sf$status[1]]])
  status01   <- to_status01(status_raw)
  names(status01) <- rownames(pd)

  common <- intersect(colnames(expr), names(status01))
  expr   <- expr[, common, drop = FALSE]
  status01 <- status01[common]
  keep_s <- !is.na(status01)
  expr <- expr[, keep_s, drop = FALSE]; status01 <- status01[keep_s]

  score <- score_geneset(expr, copd_sig)
  cat("  [GSE65682] success:", sum(keep_s), "samples (dead=", sum(status01),
      "), matched genes", sum(copd_sig %in% rownames(expr)), "\n", sep = "")
  data.frame(cohort = "GSE65682", status = status01, score = score,
             stringsAsFactors = FALSE)
}

# ---- Collect three cohorts (each fault-tolerant, print reason on failure) ----
try_get <- function(tag, fun) {
  res <- tryCatch(fun(), error = function(e) {
    cat("  [", tag, "] error:", conditionMessage(e), "\n", sep = "")
    NULL
  })
  if (is.null(res)) cat("  [", tag, "] no scores obtained\n", sep = "")
  res
}
roc_list <- list()
roc_list$GSE95233    <- try_get("GSE95233", get_gse95233)
roc_list$EMTAB4451   <- try_get("E-MTAB-4451", get_emtab4451)
roc_list$GSE65682    <- try_get("GSE65682", get_gse65682)
roc_list <- roc_list[!vapply(roc_list, is.null, logical(1))]
cat("cohorts successfully collected:", length(roc_list), "\n")

if (length(roc_list) == 0) {
  cat("!! no cohort scores available; skipping Fig6c (run scripts 15/16/24 first to generate scores)\n")
} else {
  roc_scores <- do.call(rbind, roc_list)
  roc_scores <- roc_scores[!is.na(roc_scores$score), ]
  write.csv(roc_scores, file.path(out_dir, "Fig6c_ROC_scores.csv"), row.names = FALSE)
  cat("saved data: Fig6c_ROC_scores.csv (", nrow(roc_scores), " samples)\n", sep = "")

  # ---- ROC computation (pROC preferred; base-R manual fallback otherwise) ----
  calc_roc <- function(status, score) {
    if (length(unique(status)) < 2 || any(table(status) < 2)) return(NULL)
    if (has_pROC) {
      r <- pROC::roc(status, score, quiet = TRUE, direction = "auto")
      # Take coordinates
      df <- data.frame(fpr = 1 - r$specificities, tpr = r$sensitivities)
      df <- df[order(df$fpr, df$tpr), ]
      list(df = df, auc = as.numeric(r$auc))
    } else {
      # Mann-Whitney AUC + stepwise ROC
      n_pos <- sum(status == 1); n_neg <- sum(status == 0)
      auc <- (sum(rank(score)[status == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
      ord <- order(score, decreasing = TRUE)
      s <- status[ord]
      tpr <- cumsum(s == 1) / n_pos
      fpr <- cumsum(s == 0) / n_neg
      df <- data.frame(fpr = c(0, fpr, 1), tpr = c(0, tpr, 1))
      list(df = df, auc = auc)
    }
  }

  curve_list <- lapply(unique(roc_scores$cohort), function(cq) {
    sub <- roc_scores[roc_scores$cohort == cq, ]
    r <- calc_roc(sub$status, sub$score)
    if (is.null(r)) return(NULL)
    data.frame(cohort = cq, fpr = r$df$fpr, tpr = r$df$tpr,
               auc = r$auc, stringsAsFactors = FALSE)
  })
  curve_list <- curve_list[!vapply(curve_list, is.null, logical(1))]
  roc_curve <- do.call(rbind, curve_list)

  # Save ROC curve points before plotting
  write.csv(roc_curve, file.path(out_dir, "Fig6c_ROC_curve.csv"), row.names = FALSE)
  cat("saved data: Fig6c_ROC_curve.csv\n")

  # AUC labels
  auc_lab <- unique(roc_curve[, c("cohort", "auc")])
  auc_lab$label <- sprintf("%s  AUC=%.3f", auc_lab$cohort, auc_lab$auc)

  cohort_col <- c(GSE65682 = "#d62728", `E-MTAB-4451` = "#1f77b4", GSE95233 = "#2ca02c")

  p3 <- ggplot(roc_curve, aes(x = fpr, y = tpr, color = cohort)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = cohort_col) +
    coord_equal() +
    labs(x = "1 - Specificity (False positive rate)",
         y = "Sensitivity (True positive rate)",
         title = "COPD signature predicts 28-day mortality",
         color = "Cohort") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = c(0.75, 0.25))

  # Annotate AUC on the plot (per cohort)
  for (i in seq_len(nrow(auc_lab))) {
    p3 <- p3 + annotate("text", x = 0.6, y = 0.30 - 0.06 * (i - 1),
                        label = auc_lab$label[i], size = 3.4, hjust = 0,
                        color = cohort_col[auc_lab$cohort[i]])
  }
  save_both(p3, file.path(out_dir, "Fig6c_ROC"), w = 6, h = 6)
}

cat("\n35 three missing panels done\n")
cat("Data files: Fig4a_SOCS3_logFC.csv / Fig5b_149_module_mapping.csv / Fig6c_ROC_scores.csv / Fig6c_ROC_curve.csv\n")
