# =============================================================================
# 55_MR_two-sample_SOCS3_core_genes.R
# Two-sample Mendelian randomization (Two-sample MR):
#   Exposure = gene expression (eQTLGen whole-blood cis-eQTL)
#   Outcome = sepsis (Sepsis) + COPD
# Purpose: Upgrade SOCS3 (and 149 core genes / IFN module genes) from "observational association" to "genetic causality".
#   Expected direction: SOCS3 is a negative-feedback "brake" -> high SOCS3 expression -> low disease risk (OR < 1).
# Reference: YunShengXin MR approach (five methods + heterogeneity + pleiotropy + leave-one-out + Steiger directionality)
#       + TwoSampleMR standard workflow.
# -----------------------------------------------------------------------------
# Data sources (IEU OpenGWAS, https://gwas.mrcieu.ac.uk/):
#   Exposure eQTLGen whole-blood cis-eQTL : "eqtl-a-<ENSG>" (31,684 Europeans)
#   Outcome Sepsis  : "ieu-b-4980"  (UK Biobank sepsis, 486,484 individuals / 11,643 cases)
#                  alt "finn-b-AB1_SEPSIS" (FinnGen septicaemia)
#   Outcome COPD    : "finn-b-J10_COPD" (FinnGen COPD, replace if ID invalid)
# -----------------------------------------------------------------------------
# Prerequisites (server):
#   install.packages("TwoSampleMR"); install.packages("ieugwasr")
#   devtools::install_github("rondolab/MR-PRESSO")  # optional, skip if missing
#   Requires network access to IEU OpenGWAS API (clump uses IEU LD service)
# Run: Rscript 55_MR_two-sample_SOCS3_core_genes.R [OUT_DIR]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))))
if (length(script_dir) == 0) script_dir <- getwd()
OUT_DIR <- ifelse(length(args) >= 1 && nzchar(args[1]), args[1], script_dir)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("OUT_DIR:", OUT_DIR, "\n")

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

# =============================================================================
# 0. Environment check (IEU OpenGWAS API: requires ieugwasr >= 0.2.1 and an OPENGWAS_JWT token set)
# =============================================================================
igw_ver <- tryCatch(as.character(packageVersion("ieugwasr")), error = function(e) "not installed")
cat("ieugwasr version:", igw_ver, "\n")
if (identical(igw_ver, "not installed")) {
  stop("Missing ieugwasr package. Please run: install.packages('ieugwasr') or remotes::install_github('MRCIEU/ieugwasr')")
}
if (packageVersion("ieugwasr") < "0.2.0") {
  cat("!! Warning: ieugwasr version too old, extract_instruments will error 'unused argument (x_api_source...)'.\n",
      "  Please upgrade: remotes::install_github('MRCIEU/ieugwasr')\n")
}

jwt <- Sys.getenv("OPENGWAS_JWT")
token_file <- file.path(OUT_DIR, "opengwas_token.txt")
if (!file.exists(token_file)) token_file <- file.path(script_dir, "opengwas_token.txt")
if (!nzchar(jwt) && file.exists(token_file)) {
  jwt <- trimws(readLines(token_file, warn = FALSE)[1])
  Sys.setenv(OPENGWAS_JWT = jwt)
  cat("Token read from file:", token_file, "\n")
} else if (!nzchar(jwt)) {
  cat("!! Token file not found:", token_file, "\n")
  cat("   Please sync opengwas_token.txt to:", OUT_DIR, "\n")
  cat("   or set environment variable: Sys.setenv(OPENGWAS_JWT='your token')\n")
}
if (nzchar(jwt)) {
  cat("OPENGWAS_JWT token is in place\n")
} else {
  cat("!! OPENGWAS_JWT token not detected. Please first:\n",
      "  1) Visit https://api.opengwas.io/ log in with GitHub -> Generate a token\n",
      "  2) In R run Sys.setenv(OPENGWAS_JWT='your token') or write to ~/.Rprofile\n")
}

# =============================================================================
# 1. Configuration (adjustable parameters)
# =============================================================================

# ---- Exposure genes: gene symbol = ENSG ID ----
# SOCS3 is the main gene (core of this study); the rest are core genes of the IFN module / SOCS3-PPI network.
# ⚠️ Please verify ENSG IDs (queryable via bioMart or Ensembl), especially extended genes.
GENES <- c(
  SOCS3  = "ENSG00000164509",   # main gene, confirmed
  STAT1  = "ENSG00000115415",
  STAT3  = "ENSG00000168610",
  JAK1   = "ENSG00000162434",
  JAK2   = "ENSG00000096968",
  TYK2   = "ENSG00000105397",
  IL6ST  = "ENSG00000134352",   # gp130, SOCS3 inhibits its JAK recruitment
  ISG15  = "ENSG00000187608",
  MX1    = "ENSG00000157601",
  GBP1   = "ENSG00000117228",
  IFI44L = "ENSG00000137959"
)

# ---- Outcome GWAS (IEU OpenGWAS ID) ----
OUTCOMES <- c(
  sepsis = "ieu-b-4980",       # UK Biobank sepsis (11,643 cases)
  copd   = "finn-b-J10_COPD"   # FinnGen COPD (if ID invalid, replace with UK Biobank COPD)
)

# ---- Instrumental variable thresholds ----
EQTL_P   <- 5e-8      # cis-eQTL significance threshold (relax to 1e-5 when SOCS3 IVs are insufficient)
CLUMP_R2 <- 0.001     # LD clumping r2
CLUMP_KB <- 10000     # LD clumping window kb
F_STAT   <- 10        # weak-instrument F threshold

# ---- Output switches ----
SAVE_PLOTS  <- TRUE  # whether to generate plots (scatter/forest/funnel/leaveoneout)
DO_MRPRESSO <- TRUE  # whether to run MR-PRESSO (auto-skip if package not installed)
DO_STEIGER  <- TRUE  # whether to run Steiger directionality test

# =============================================================================
# 2. Utility functions
# =============================================================================

fstat <- function(b, se, eaf, n = 31684) {
  # Single-SNP F statistic (approximate), n = eQTLGen sample size
  ifelse(se > 0, (b / se)^2, NA_real_)
}

safe_extract_instruments <- function(gene_id, p1, r2, kb) {
  out <- tryCatch(
    extract_instruments(outcomes = gene_id, p1 = p1, clump = TRUE,
                        r2 = r2, kb = kb),
    error = function(e) {
      cat("  !! extract_instruments failed (", gene_id, "):",
          conditionMessage(e), "\n")
      cat("     Retrying extraction without clumping...\n")
      tryCatch(
        extract_instruments(outcomes = gene_id, p1 = p1, clump = FALSE),
        error = function(e2) {
          cat("  !! Without clumping also failed:", conditionMessage(e2), "\n")
          return(NULL)
        })
    })
  if (is.null(out) || nrow(out) == 0) return(NULL)
  out
}

safe_extract_outcome <- function(snps, out_id) {
  out <- tryCatch(
    extract_outcome_data(snps = snps, outcomes = out_id),
    error = function(e) {
      cat("  !! extract_outcome_data failed (", out_id, "):",
          conditionMessage(e), "\n")
      return(NULL)
    })
  out
}

# =============================================================================
# 3. Main loop: gene by gene x outcome by outcome
# =============================================================================

res_all <- list()       # mr() results
het_all <- list()
plei_all <- list()
loo_all <- list()
steiger_all <- list()
presso_all <- list()

for (gn in names(GENES)) {
  ensg <- GENES[[gn]]
  eqtl_id <- paste0("eqtl-a-", ensg)
  cat("\n========== Gene:", gn, "(", eqtl_id, ") ==========\n")

  expo <- safe_extract_instruments(eqtl_id, EQTL_P, CLUMP_R2, CLUMP_KB)
  if (is.null(expo)) { cat("  -> No usable instrumental variables, skipping\n"); next }
  cat("  Number of instrument SNPs:", nrow(expo), "\n")

  # Instrument strength F statistic
  expo$F <- fstat(expo$beta.exposure, expo$se.exposure, expo$eaf.exposure)
  weak <- sum(expo$F < F_STAT, na.rm = TRUE)
  cat("  F statistic range:", if (all(is.na(expo$F))) "NA" else
      paste0(round(min(expo$F, na.rm=TRUE),1), "~", round(max(expo$F, na.rm=TRUE),1)),
      "| weak instruments (F<", F_STAT, ") count:", weak, "\n")
  if (all(expo$F < F_STAT, na.rm = TRUE)) {
    cat("  !! All are weak instruments, skipping this gene\n"); next
  }

  for (oc in names(OUTCOMES)) {
    out_id <- OUTCOMES[[oc]]
    cat("  ---- Outcome:", oc, "(", out_id, ") ----\n")

    out <- safe_extract_outcome(expo$SNP, out_id)
    if (is.null(out) || nrow(out) == 0) {
      cat("    -> Outcome has no matching SNP, skipping\n"); next
    }

    dat <- tryCatch(
      harmonise_data(exposure_dat = expo, outcome_dat = out),
      error = function(e) {
        cat("    !! harmonise failed:", conditionMessage(e), "\n"); return(NULL)
      })
    if (is.null(dat) || nrow(dat) == 0) { cat("    -> No SNP after harmonise, skipping\n"); next }

    key <- paste(gn, oc, sep = "_")

    # ---- MR analysis: multiple SNPs use five methods; single SNP uses manual Wald ratio (to avoid mr() single-SNP bug) ----
    d_ok <- dat[dat$mr_keep, , drop = FALSE]
    res <- NULL
    if (nrow(d_ok) >= 2) {
      res <- tryCatch(
        mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                                "mr_weighted_median", "mr_simple_mode",
                                "mr_weighted_mode")),
        error = function(e) {
          cat("    !! mr() failed:", conditionMessage(e), "\n"); return(NULL)
        })
    } else if (nrow(d_ok) == 1) {
      x <- d_ok[1, ]
      if (!is.na(x$beta.exposure) && abs(x$beta.exposure) > 0) {
        b_w <- x$beta.outcome / x$beta.exposure
        se_w <- abs(x$se.outcome / x$beta.exposure)
        p_w <- 2 * pnorm(-abs(b_w / se_w))
        res <- data.frame(id.exposure = x$id.exposure, id.outcome = x$id.outcome,
                          exposure = x$exposure, outcome = x$outcome,
                          method = "Wald ratio", nsnp = 1L,
                          b = b_w, se = se_w, pval = p_w,
                          stringsAsFactors = FALSE)
        cat("    (Single SNP, manual Wald ratio: b=", round(b_w, 4),
            " p=", format(p_w, digits = 3), ")\n")
      } else {
        cat("    !! Exposure effect is 0/NA, cannot compute Wald ratio\n")
      }
    } else {
      cat("    !! No valid mr_keep SNP after harmonise\n")
    }
    if (!is.null(res) && nrow(res) > 0) {
      res$gene <- gn; res$outcome <- oc; res$n_snp <- nrow(d_ok)
      res_all[[key]] <- res
      main_method <- if ("Inverse variance weighted" %in% res$method) {
        "Inverse variance weighted"
      } else { res$method[1] }
      cat("    ", main_method, ": b=", round(res$b[res$method == main_method], 4),
          " p=", format(res$pval[res$method == main_method], digits = 3),
          " nSNP=", nrow(d_ok), "\n")
    }

    # ---- Heterogeneity (Cochran's Q) ----
    het <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
    if (!is.null(het) && nrow(het) > 0) { het$gene <- gn; het$outcome <- oc; het_all[[key]] <- het } else if (!is.null(het)) { cat("    (Heterogeneity: 0 rows, skipping)\n") }

    # ---- Pleiotropy (MR-Egger intercept) ----
    plei <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
    if (!is.null(plei) && nrow(plei) > 0) { plei$gene <- gn; plei$outcome <- oc; plei_all[[key]] <- plei }

    # ---- Leave-one-out ----
    loo <- tryCatch(mr_leaveoneout(dat), error = function(e) NULL)
    if (!is.null(loo) && nrow(loo) > 0) { loo$gene <- gn; loo$outcome <- oc; loo_all[[key]] <- loo } else if (!is.null(loo)) { cat("    (Leave-one-out: 0 rows, skipping)\n") }

    # ---- Steiger directionality ----
    if (DO_STEIGER) {
      st <- tryCatch(directionality_test(dat), error = function(e) NULL)
      if (!is.null(st) && nrow(st) > 0) { st$gene <- gn; st$outcome <- oc; steiger_all[[key]] <- st }
    }

    # ---- MR-PRESSO (optional, skip if package not installed) ----
    if (DO_MRPRESSO && nrow(dat) > 3) {
      pr <- tryCatch(
        run_mr_presso(dat, NbDistribution = 300),
        error = function(e) {
          cat("    (MR-PRESSO skipped:", conditionMessage(e), ")\n"); NULL
        })
      if (!is.null(pr)) {
        tryCatch({
          g <- pr[[1]]$`MR-PRESSO results`$`Global Test`
          df <- data.frame(gene = gn, outcome = oc,
                           RSSobs = g$RSSobs, Pvalue = g$Pvalue,
                           stringsAsFactors = FALSE)
          presso_all[[key]] <- df
        }, error = function(e) NULL)
      }
    }

    # ---- Plotting ----
    if (SAVE_PLOTS) {
      tag <- paste0(gn, "_", oc)
      tryCatch({
        p1 <- mr_scatter_plot(res, dat)[[1]]
        ggsave(file.path(OUT_DIR, paste0("path55_", tag, "_scatter.pdf")), p1, width=6, height=5)
        ggsave(file.path(OUT_DIR, paste0("path55_", tag, "_scatter.png")), p1, width=6, height=5, dpi=300)
      }, error = function(e) cat("    scatter plot failed:", conditionMessage(e), "\n"))
      tryCatch({
        p2 <- mr_forest_plot(mr_singlesnp(dat))[[1]]
        ggsave(file.path(OUT_DIR, paste0("path55_", tag, "_forest.pdf")), p2, width=6, height=4.5)
      }, error = function(e) cat("    forest plot failed:", conditionMessage(e), "\n"))
      tryCatch({
        p3 <- mr_funnel_plot(mr_singlesnp(dat))[[1]]
        ggsave(file.path(OUT_DIR, paste0("path55_", tag, "_funnel.pdf")), p3, width=6, height=5)
      }, error = function(e) cat("    funnel plot failed:", conditionMessage(e), "\n"))
      tryCatch({
        p4 <- mr_leaveoneout_plot(loo)[[1]]
        ggsave(file.path(OUT_DIR, paste0("path55_", tag, "_leaveoneout.pdf")), p4, width=6, height=4.5)
      }, error = function(e) cat("    leaveoneout plot failed:", conditionMessage(e), "\n"))
    }
  }
}

# =============================================================================
# 4. Result summary and saving
# =============================================================================

mr_df    <- if (length(res_all))    bind_rows(res_all)    else data.frame()
het_df   <- if (length(het_all))    bind_rows(het_all)    else data.frame()
plei_df  <- if (length(plei_all))   bind_rows(plei_all)   else data.frame()
loo_df   <- if (length(loo_all))    bind_rows(loo_all)    else data.frame()
stei_df  <- if (length(steiger_all)) bind_rows(steiger_all) else data.frame()
presso_df<- if (length(presso_all)) bind_rows(presso_all) else data.frame()

# Multiple-testing correction (Bonferroni by number of genes)
if (nrow(mr_df) > 0) {
  n_test <- length(GENES) * length(OUTCOMES)
  mr_df$bonf_threshold <- 0.05 / n_test
  mr_df$significant_bonf <- mr_df$pval < mr_df$bonf_threshold
}

write.csv(mr_df,    file.path(OUT_DIR, "path55_mr_results.csv"),    row.names = FALSE)
write.csv(het_df,   file.path(OUT_DIR, "path55_mr_heterogeneity.csv"), row.names = FALSE)
write.csv(plei_df,  file.path(OUT_DIR, "path55_mr_pleiotropy.csv"),   row.names = FALSE)
write.csv(stei_df,  file.path(OUT_DIR, "path55_mr_steiger.csv"),      row.names = FALSE)
if (nrow(presso_df)) write.csv(presso_df, file.path(OUT_DIR, "path55_mr_presso.csv"), row.names = FALSE)

# ---- Key summary (IVW results) ----
cat("\n\n========== IVW result summary ==========\n")
if (nrow(mr_df) > 0) {
  ivw <- mr_df[mr_df$method == "Inverse variance weighted", ]
  ivw <- ivw[order(ivw$pval), ]
  print(ivw[, c("gene","outcome","n_snp","b","se","pval","significant_bonf")],
        row.names = FALSE)
}

# ---- Pre-plot data (re-plottable locally) ----
saveRDS(
  list(genes = GENES, outcomes = OUTCOMES,
       mr = mr_df, heterogeneity = het_df, pleiotropy = plei_df,
       steiger = stei_df, presso = presso_df,
       params = list(eqtl_p = EQTL_P, clump_r2 = CLUMP_R2,
                     clump_kb = CLUMP_KB, f_stat = F_STAT)),
  file.path(OUT_DIR, "path55_mr_plotdata.rds")
)

cat("\nSaved: path55_mr_results.csv / heterogeneity / pleiotropy / steiger / presso + plotdata.rds\n")
cat("===== Script 55 complete =====\n")
