#!/usr/bin/env Rscript
# =============================================================================
# 68_permutation_overlap_test.R   (N = 10,000 version)
# Purpose: use a permutation test instead of the hypergeometric test to answer
#          "do the 149 shared programs significantly exceed random expectation"
# Background: a reviewer pointed out that the hypergeometric test's background
#       gene set is biased (detectable genes are not a uniform draw), and that
#       the binomial test (P=7.7e-40) rests on an invalid null hypothesis (50%
#       shared direction). The permutation test preserves the complete filtering
#       pipeline and donor structure, randomly permutes disease labels and
#       recomputes the overlap to get an empirical P value — the only approach
#       that sidesteps both the "background gene set" and "gene variability
#       bias" objections at once.
#
# Method:
#   1) Donor-level pseudobulk (blood monocytes Sepsis vs Healthy; lung myeloid COPD vs Control)
#   2) Real DEGs: donor-level edgeR (qlf), up = logFC>0.25 & FDR<0.05
#   3) Real overlap = blood up ∩ lung up (upward arm) / blood down ∩ lung down (downward arm)
#   4) Permute N times: randomly shuffle the "disease vs healthy" labels
#      (preserving donor structure and the filtering pipeline), rerun edgeR,
#      get up/down genes, compute the two-tissue overlap
#   5) Empirical P = (# permutations with overlap >= observed + 1) / (N + 1)
#
# ---------------------------------------------------------------------------
# Changes in v2 (10,000 perms) vs v1 (500 perms) — only for runnability +
# reproducibility, no change to the statistics:
#
#   (a) Reproducibility: v1 used set.seed(NULL) (reset from the clock), so
#       results were not reproducible and reviewers could not recompute.
#       v2 uses a fixed MASTER_SEED to pre-generate N per-permutation seeds.
#       Because each permutation sets its own seed, parallel / serial /
#       resumed runs give identical results.
#
#   (b) Exact caching: sample() is a permutation, group sizes are constant, so
#       filterByExpr's min.group.size is unchanged; TMM (calcNormFactors)
#       depends only on counts, not labels. Hence "gene filtering" and
#       "normalization factors" are identical across all permutations —
#       compute once. Note: this is exact equivalence, not an approximation.
#       Each permutation still reruns estimateDisp + glmQLFit + glmQLFTest
#       (the parts that truly depend on labels).
#
#   (c) Parallel + checkpoint/resume: 10,000 serial permutations could take
#       10+ hours. mclapply runs multi-core parallel with a checkpoint saved
#       every CHUNK=100 permutations; if the process is killed, rerunning
#       resumes automatically from the checkpoint. Core count is controlled
#       by the NPROC env var (default = min(detected logical cores, 16),
#       hard cap 128); at startup one parallel "wave" (ncores tasks) is
#       benchmarked to give an ETA and a linear extrapolation reference for
#       adding cores.
#       ⚠ Each worker needs ~1-2 GB RAM → keep NPROC ≤ available RAM (GB)/2;
#       if R uses multi-threaded BLAS (OpenBLAS/MKL), too many cores can
#       actually slow things down — trust the measured [calibration] output.
#
# Depends on: Seurat, edgeR, Matrix, parallel
# Run: Rscript 68_permutation_overlap_test.R
#       NPROC=16 Rscript 68_permutation_overlap_test.R    # specify core count
# =============================================================================

## ---- Auto-install dependencies ----
if (!requireNamespace("edgeR", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("edgeR", update = FALSE, ask = FALSE)
}
suppressPackageStartupMessages({
  library(Seurat); library(edgeR); library(Matrix); library(parallel)
})

## ---- Parameters ----
N_PERM      <- 10000      # number of permutations (reviewer-required magnitude; at 500 the P hits the 1/501 floor)
MASTER_SEED <- 20260905   # master seed — record it in the Methods to guarantee reproducibility
CHUNK       <- 100        # checkpoint granularity: save every 100 perms (a killed process loses at most ~100/ncores minutes)
LOGFCT      <- 0.25       # matches FindMarkers logfc.threshold
FDR_CUT     <- 0.05       # FDR threshold for up/down regulation

## Core count: the NPROC env var takes priority; otherwise default to
## min(detected logical cores, 16); hard cap 128.
## Memory reference: ~1-2 GB per worker (fork shares pseudobulk; mostly edgeR workspace).
##   Suggested NPROC ≤ available RAM (GB)/2. If this R uses multi-threaded BLAS
##   (OpenBLAS/MKL), adding too many cores can slow things down — judge linear
##   scaling from the measured [calibration] throughput below.
N_CORES <- suppressWarnings(as.integer(Sys.getenv("NPROC", unset = NA_character_)))
if (is.na(N_CORES)) N_CORES <- suppressWarnings(detectCores(logical = TRUE))
if (is.na(N_CORES) || N_CORES < 1L) N_CORES <- 1L
N_CORES <- max(1L, min(as.integer(N_CORES), 128L))

## ---- Single-thread BLAS/OMP: prevent forked workers' multi-threaded BLAS from oversubscribing cores ----
## Each forked worker inherits the parent R process's BLAS thread pool; if BLAS
## defaults to multi-threaded (common on big servers with OpenBLAS/MKL),
## NPROC workers × BLAS threads oversubscribes badly and slows down.
## If RhpcBLASctl is available, force single-thread per process (set in the
## parent, inherited by children after fork).
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
  try(RhpcBLASctl::omp_set_num_threads(1), silent = TRUE)
  cat("[BLAS] RhpcBLASctl available → BLAS/OMP set to single-thread\n")
} else {
  cat("[BLAS] RhpcBLASctl not available; if using multi-threaded BLAS, launch with:\n")
  cat("       OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NPROC=64 Rscript <this script>\n")
}

## ---- Paths ----
this_file <- commandArgs(trailingOnly = FALSE)
.f <- grep("--file=", this_file, value = TRUE)
if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1]) else .this_file <- "."
out_dir <- dirname(normalizePath(.this_file))
source(file.path(out_dir, "00_config.R"))

## ---- Donor field detection (reuses script 28's logic) ----
detect_donor <- function(seu) {
  meta <- seu@meta.data
  get <- function(f) if (f %in% colnames(meta)) as.character(meta[[f]]) else NULL
  patient <- get("patient"); sample <- get("sample")
  if (!is.null(patient)) {
    if (!is.null(sample)) {
      na_idx <- is.na(patient) | patient == ""
      patient[na_idx] <- sample[na_idx]
      return(patient)
    }
    return(patient)
  }
  for (cand in c("donor_id","donor","subject","subject_id","Subject_Identity","sample")) {
    v <- get(cand)
    if (!is.null(v)) return(v)
  }
  NULL
}

## ---- Donor-level pseudobulk construction ----
build_pb <- function(seu, keep, tag = "") {
  meta <- seu@meta.data
  donor_vec <- detect_donor(seu)
  if (is.null(donor_vec)) stop(sprintf("[%s] no donor-level field (patient/sample)", tag))
  if (length(keep) != nrow(meta))
    stop(sprintf("[%s] keep length %d != meta rows %d (internal error)", tag, length(keep), nrow(meta)))
  keep[is.na(keep)] <- FALSE
  if (sum(keep) == 0) stop(sprintf("[%s] keep matched no cells; check cell_type/group field values", tag))

  counts_all <- GetAssayData(seu, assay = "RNA", layer = "counts")
  miss <- setdiff(rownames(meta)[keep], colnames(counts_all))
  if (length(miss) > 0) {
    cat(sprintf("[%s] warning: %d/%d matched cells not in counts column names\n",
                tag, length(miss), sum(keep)))
    cat(sprintf("[%s]   missing barcode examples: %s\n", tag, paste(head(miss, 5), collapse = ", ")))
    cat(sprintf("[%s]   counts column name examples:   %s\n", tag, paste(head(colnames(counts_all), 5), collapse = ", ")))
    cat(sprintf("[%s]   meta row name examples:     %s\n", tag, paste(head(rownames(meta)[keep], 5), collapse = ", ")))
    if (length(miss) == sum(keep))
      stop(sprintf("[%s] counts column names completely mismatch meta row names: object may be Seurat v5 split layers or barcodes were rewritten", tag))
  }

  donor_all <- as.character(donor_vec)
  names(donor_all) <- rownames(meta)
  cells <- rownames(meta)[keep]
  cells <- cells[cells %in% colnames(counts_all)]
  meta <- meta[cells, , drop = FALSE]
  donor_vec <- donor_all[cells]
  counts <- counts_all[, cells, drop = FALSE]

  donors <- unique(donor_vec)
  donors <- donors[!is.na(donors) & donors != ""]
  cts_list <- lapply(donors, function(d) {
    cell_i <- names(donor_vec)[donor_vec == d]
    Matrix::rowSums(counts[, cell_i, drop = FALSE])
  })
  cts <- do.call(cbind, cts_list); colnames(cts) <- donors
  grp <- vapply(donors, function(d) unique(meta$group[donor_vec == d])[1], character(1))
  cat(sprintf("[%s] counts class=%s | cells into pseudobulk=%d | donors=%d\n",
              tag, class(counts_all)[1], length(cells), length(donors)))
  list(counts = cts, grp = grp)
}

## ---- Precompute: quantities invariant under label permutation (see header (b)) ----
prep_cache <- function(pb, ref, target) {
  grp_f <- factor(pb$grp, levels = c(ref, target))
  y0 <- DGEList(counts = pb$counts, group = grp_f)
  keep <- filterByExpr(y0, min.count = 10, min.total.count = 15)
  y0 <- y0[keep, , keep.lib.sizes = FALSE]
  y0 <- calcNormFactors(y0)               # DGEList method: returns the whole DGEList (not a vector)
  nf <- y0$samples$norm.factors           # pull the numeric vector back from the object (length = donor count)
  list(keep = keep, nf = nf, n_kept = sum(keep), levels = c(ref, target))
}

## ---- Donor-level DEG (using the cache; labels passed in by caller) ----
deg_arms <- function(pb, cache, grp_vec) {
  grp_f <- factor(grp_vec, levels = cache$levels)
  if (length(unique(grp_f)) < 2 || min(table(grp_f)) < 2) return(NULL)
  # norm.factors passed straight into the constructor (DGEList() accepts this arg);
  # lib.size defaults to recomputation from the filtered counts, matching
  # deg_full's keep.lib.sizes=FALSE path → exactly equivalent
  y <- DGEList(counts = pb$counts[cache$keep, , drop = FALSE],
               group = grp_f, norm.factors = cache$nf)
  design <- model.matrix(~ grp_f)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)           # target vs ref
  tab <- topTags(qlf, n = Inf)$table
  up   <- rownames(tab)[tab$logFC >  LOGFCT & tab$FDR < FDR_CUT]
  down <- rownames(tab)[tab$logFC < -LOGFCT & tab$FDR < FDR_CUT]
  list(up = up, down = down, ntested = nrow(tab))
}

## ---- Donor-level DEG (full recompute path, only for the startup self-check) ----
deg_full <- function(pb, ref, target) {
  grp_f <- factor(pb$grp, levels = c(ref, target))
  if (length(unique(grp_f)) < 2 || min(table(grp_f)) < 2) return(NULL)
  y <- DGEList(counts = pb$counts, group = grp_f)
  keep <- filterByExpr(y, min.count = 10, min.total.count = 15)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~ grp_f)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  tab <- topTags(qlf, n = Inf)$table
  up   <- rownames(tab)[tab$logFC >  LOGFCT & tab$FDR < FDR_CUT]
  down <- rownames(tab)[tab$logFC < -LOGFCT & tab$FDR < FDR_CUT]
  list(up = up, down = down, ntested = nrow(tab), nf = y$samples$norm.factors)
}

## ---- One permutation: explicit seed → independent of parallel/resume order, reproducible ----
perm_once <- function(seed, pb_blood, pb_lung, cache_b, cache_l) {
  set.seed(seed)
  b <- deg_arms(pb_blood, cache_b, sample(pb_blood$grp))
  l <- deg_arms(pb_lung,  cache_l, sample(pb_lung$grp))
  if (is.null(b) || is.null(l)) return(c(up = NA_real_, down = NA_real_))
  c(up   = length(intersect(b$up,   l$up)),
    down = length(intersect(b$down, l$down)))
}

## ============================================================================
## Part 1: read objects + build donor-level pseudobulk
## ============================================================================
cat("===== 68: permutation test (overlap significance) | N = ", N_PERM, " | cores = ", N_CORES,
    " | master seed = ", MASTER_SEED, " =====\n", sep = "")

blood_file <- file.path(out_dir, "path1_sepsis_copd_integrated.rds")
blood <- readRDS(blood_file)
b_keep <- blood$cell_type == "Monocyte" & blood$group %in% c("Healthy", "Sepsis")
pb_blood <- build_pb(blood, b_keep, tag = "blood-mono")
cat("blood monocyte donors: ", ncol(pb_blood$counts), " | groups: ", sep = "")
print(table(pb_blood$grp))

lung <- readRDS(file.path(out_dir, "path2_copd_lung.rds"))
if (!"group" %in% colnames(lung@meta.data)) lung$group <- lung$disease
ct_all <- unique(as.character(lung$cell_type))
cat("all lung object cell_type values: ", paste(ct_all, collapse = ", "), "\n", sep = "")
cat("all lung object group values:     ", paste(unique(as.character(lung$group)), collapse = ", "), "\n", sep = "")
myeloid_types <- grep("Macrophage|Monocyte|monocyte|macrophage|Myeloid|myeloid|DC|dendritic",
                      ct_all, value = TRUE)
if (length(myeloid_types) == 0)
  myeloid_types <- intersect(c("Myeloid", "Myeloid cells", "Myeloid"), ct_all)
cat("lung myeloid cell_type candidates: ", paste(myeloid_types, collapse = ", "), "\n", sep = "")
if (length(myeloid_types) == 0)
  stop("lung object cell_type matched no myeloid naming; adjust the grep pattern based on the output above")
l_keep <- lung$cell_type %in% myeloid_types & lung$group %in% c("Control", "COPD")
pb_lung <- build_pb(lung, l_keep, tag = "lung-myeloid")
cat("lung myeloid donors: ", ncol(pb_lung$counts), " | groups: ", sep = "")
print(table(pb_lung$grp))

## ============================================================================
## Part 2: real DEGs + real overlap
## ============================================================================
cache_b <- prep_cache(pb_blood, "Healthy", "Sepsis")
cache_l <- prep_cache(pb_lung,  "Control", "COPD")
cat("\n[cache] blood-side retained genes ", cache_b$n_kept, " | lung-side retained genes ", cache_l$n_kept,
    " (filterByExpr + TMM normalization factors computed once)\n", sep = "")

b_res <- deg_arms(pb_blood, cache_b, pb_blood$grp)
l_res <- deg_arms(pb_lung,  cache_l, pb_lung$grp)
cat("\nblood-side up (Sepsis vs Healthy): ", length(b_res$up),
    " | down: ", length(b_res$down), " (testing ", b_res$ntested, " genes)\n", sep = "")
cat("lung-side up (COPD vs Control): ", length(l_res$up),
    " | down: ", length(l_res$down), " (testing ", l_res$ntested, " genes)\n", sep = "")

## ---- Startup self-check: cached path vs full recompute path must agree gene by gene ----
## Purpose: empirically verify that "cached filterByExpr + TMM" is equivalent to
##       recomputing each time (as claimed in header (b)).
##       Runs once at startup (one full recompute each for blood+lung, ~1-2 extra
##       minutes); only run the 10,000 permutations if this passes.
cat("\n[self-check] comparing cached path vs full recompute path...\n", sep = "")
selftest <- function(tag, pb, cache, full) {
  cch <- deg_arms(pb, cache, pb$grp)         # real labels through the cached path
  up_same <- identical(sort(full$up),   sort(cch$up))
  dn_same <- identical(sort(full$down), sort(cch$down))
  nt_same <- full$ntested == cch$ntested
  nf_same <- isTRUE(all.equal(as.numeric(full$nf), as.numeric(cache$nf), tolerance = 1e-12))
  ok <- up_same && dn_same && nt_same && nf_same
  cat(sprintf("[self-check] %s: ntested match=%s | up-set match=%s | down-set match=%s | TMM factors match=%s -> %s\n",
              tag, nt_same, up_same, dn_same, nf_same,
              if (ok) "PASS" else "FAIL (do not continue; report this output)"))
  ok
}
b_full <- deg_full(pb_blood, "Healthy", "Sepsis")
l_full <- deg_full(pb_lung,  "Control", "COPD")
ok_b <- selftest("blood-mono", pb_blood, cache_b, b_full)
ok_l <- selftest("lung-myeloid", pb_lung,  cache_l, l_full)
if (!ok_b || !ok_l) stop("self-check failed: cached path and full recompute path disagree; check the filterByExpr/TMM caching logic")

real_up   <- length(intersect(b_res$up, l_res$up))
real_down <- length(intersect(b_res$down, l_res$down))
cat("\nreal overlap: up arm ", real_up, " genes | down arm ", real_down, " genes\n", sep = "")

## ============================================================================
## Part 3: permutation test (parallel + checkpoint/resume)
## ============================================================================
CKPT <- file.path(out_dir, "path68_perm10k_checkpoint.rds")

## Pre-generate per-permutation seeds (once globally, for reproducibility)
set.seed(MASTER_SEED)
perm_seeds <- sample.int(.Machine$integer.max - 1L, N_PERM)

perm_up <- rep(NA_real_, N_PERM); perm_down <- rep(NA_real_, N_PERM)
start_i <- 1L

## Resume: if a checkpoint exists with matching size, continue from it
if (file.exists(CKPT)) {
  ck <- try(readRDS(CKPT), silent = TRUE)
  if (!inherits(ck, "try-error") && !is.null(ck$perm_up) &&
      length(ck$perm_up) == N_PERM && identical(ck$seeds, perm_seeds)) {
    done <- sum(!is.na(ck$perm_up))
    if (done > 0) {
      perm_up <- ck$perm_up; perm_down <- ck$perm_down
      start_i <- done + 1L
      cat("\n[resume] checkpoint found: ", done, "/", N_PERM,
          " done, continuing from permutation ", start_i, "\n", sep = "")
    }
  } else {
    cat("\n[resume] checkpoint does not match current parameters (seeds/size); starting over\n")
  }
}

run_seeds <- function(seeds) {
  res <- if (N_CORES > 1L) {
    parallel::mclapply(seeds, perm_once,
                       pb_blood = pb_blood, pb_lung = pb_lung,
                       cache_b = cache_b, cache_l = cache_l,
                       mc.cores = N_CORES, mc.preschedule = TRUE)
  } else {
    lapply(seeds, perm_once,
           pb_blood = pb_blood, pb_lung = pb_lung,
           cache_b = cache_b, cache_l = cache_l)
  }
  pick <- function(z, nm) {
    if (is.null(z) || inherits(z, "try-error") || !(nm %in% names(z))) NA_real_ else as.numeric(z[[nm]])
  }
  list(up   = vapply(res, pick, numeric(1), nm = "up"),
       down = vapply(res, pick, numeric(1), nm = "down"))
}

if (start_i <= N_PERM) {
  ## Timing calibration: run one parallel "wave" (cal_n = ncores tasks, wall time ≈
  ## one permutation's cost) to measure throughput
  cal_n <- min(N_CORES, N_PERM)
  cat("\n[calibration] running one wave of ", cal_n, " permutations to measure throughput...\n", sep = "")
  t0 <- Sys.time()
  invisible(run_seeds(perm_seeds[seq_len(cal_n)]))
  t_wave <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  rem <- N_PERM - start_i + 1L
  rate_cal <- cal_n / t_wave                         # measured: perms/sec (cal_n cores concurrent)
  eta_cal_h <- rem / rate_cal / 3600
  cat(sprintf("[calibration] measured throughput on %d cores %.1f perms/min → %d remaining, ETA %.1f hours\n",
              cal_n, rate_cal * 60, rem, eta_cal_h))
  ## Linear extrapolation reference for adding cores (hyperthreading/memory/
  ## multi-threaded BLAS can make real scaling sublinear; order-of-magnitude only)
  hint <- function(c) rem / (rate_cal * c / cal_n) / 3600
  cat(sprintf("[calibration] scaling reference (linear extrapolation): NPROC=%d≈%.1fh | NPROC=%d≈%.1fh | NPROC=%d≈%.1fh\n",
              cal_n * 2, hint(cal_n * 2), cal_n * 4, hint(cal_n * 4), cal_n * 8, hint(cal_n * 8)))
  cat(sprintf("[calibration] to judge linearity: if one wave still takes ~%.0f s at NPROC=%d (throughput should double), adding cores helps;\n",
              t_wave, cal_n * 2))
  cat("         if the wave time rises noticeably (1.5-2×), you are bottlenecked by memory bandwidth / multi-threaded BLAS and more cores will not help.\n")

  cat("\n[permutation] starting; checkpoint saved every ", CHUNK, " permutations...\n", sep = "")
  t_start <- Sys.time()
  for (st in seq(start_i, N_PERM, by = CHUNK)) {
    en <- min(st + CHUNK - 1L, N_PERM)
    idx <- st:en
    r <- run_seeds(perm_seeds[idx])
    perm_up[idx]   <- r$up
    perm_down[idx] <- r$down
    saveRDS(list(perm_up = perm_up, perm_down = perm_down,
                 seeds = perm_seeds, n_perm = N_PERM,
                 master_seed = MASTER_SEED,
                 real_up = real_up, real_down = real_down), CKPT)
    el <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
    cat(sprintf("  %6d / %d  (%.0f minutes elapsed, ~%.0f minutes remaining)\n",
                en, N_PERM, el, el / en * (N_PERM - en)))
  }
}

## ============================================================================
## Part 4: results
## ============================================================================
emp_p <- function(obs, null) {
  ok <- !is.na(null)
  c(p        = (sum(null[ok] >= obs) + 1) / (sum(ok) + 1),
    n_ge     = sum(null[ok] >= obs),
    n_valid  = sum(ok),
    null_max = if (any(ok)) max(null[ok]) else NA_real_,
    null_mean= mean(null[ok]))
}
s_up   <- emp_p(real_up,   perm_up)
s_down <- emp_p(real_down, perm_down)

cat("\n===== permutation test results (N = ", sum(!is.na(perm_up)), ", master seed = ", MASTER_SEED, ") =====\n", sep = "")
cat(sprintf("up arm: observed overlap %d | null mean %.2f max %d | perms reaching/exceeding observed %d/%d | empirical P = %.3e\n",
            real_up, s_up["null_mean"], as.integer(s_up["null_max"]),
            as.integer(s_up["n_ge"]), as.integer(s_up["n_valid"]), s_up["p"]))
cat(sprintf("down arm: observed overlap %d | null mean %.2f max %d | perms reaching/exceeding observed %d/%d | empirical P = %.3e\n",
            real_down, s_down["null_mean"], as.integer(s_down["null_max"]),
            as.integer(s_down["n_ge"]), as.integer(s_down["n_valid"]), s_down["p"]))
cat(sprintf("\nNote: this design's P-value resolution floor = 1/(N+1) = %.2e. If n_ge = 0, P hits the floor; in that case report\n",
            1 / (sum(!is.na(perm_up)) + 1)))
cat("    \"P < 1/(N+1)\" or \"P = 1/(N+1) (0 of N permutations reached the observed overlap)\".\n")

## ---- Output ----
out <- data.frame(
  arm = c("upregulated", "downregulated"),
  observed_overlap  = c(real_up, real_down),
  permutation_mean  = c(s_up["null_mean"], s_down["null_mean"]),
  permutation_max   = c(s_up["null_max"],  s_down["null_max"]),
  n_perm_ge_observed= c(s_up["n_ge"],      s_down["n_ge"]),
  n_perm_valid      = c(s_up["n_valid"],   s_down["n_valid"]),
  empirical_p       = c(s_up["p"],         s_down["p"]),
  stringsAsFactors = FALSE
)
write.csv(out, file.path(out_dir, "path68_permutation_overlap.csv"), row.names = FALSE)
write.csv(out, file.path(out_dir, "path68_permutation_overlap_N10000.csv"), row.names = FALSE)
saveRDS(list(perm_up = perm_up, perm_down = perm_down, seeds = perm_seeds,
             n_perm = N_PERM, master_seed = MASTER_SEED,
             real_up = real_up, real_down = real_down),
        file.path(out_dir, "path68_permutation_null.rds"))
if (file.exists(CKPT)) file.remove(CKPT)
cat("\nOutput: path68_permutation_overlap.csv (+ _N10000.csv) / path68_permutation_null.rds\n")
cat("===== Script 68 done =====\n")
