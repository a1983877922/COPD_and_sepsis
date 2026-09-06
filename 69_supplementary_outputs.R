#!/usr/bin/env Rscript
# =============================================================================
# Script 69: supplementary outputs   (run after the 10,000 permutations are
# finished — to backfill numbers for reviewer R2)
# Outputs, in one go, the three small tables requested by the reviewers; no main
# analysis is changed:
#
#   (1) per_group_spearman.csv    donor-level Spearman rho + P for "SOCS3 vs the
#                                 65-gene IFN module" within each group of 2.5
#                                 (exact permutation when n <= 10, t approximation otherwise)
#                                 -> fixes the M3 wording "P<1e-4 in all five groups" (gives real values)
#
#   (2) donor_group_counts.csv    donor-level group counts in the blood object plus
#                                 "sample/group" details
#                                 -> explains Sepsis n=136 in 4.5 vs n=135 across the five groups (M1-3)
#
#   (3) observed_intersection_genes.csv
#                                 the 95 up + 55 down intersection genes under the
#                                 "true labels" of the 68 cached pipeline, with blood/lung
#                                 logFC and FDR -> rebuilds the ST4 rows (M1-2 / first review-4)
#
# Run: Rscript this script
# Requires: the same objects as script 68; edgeR/Seurat; ~5-10 min (reading objects + two edgeR runs)
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(edgeR); library(Matrix) })

this_file <- commandArgs(trailingOnly = FALSE)
.f <- grep("--file=", this_file, value = TRUE)
if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1]) else .this_file <- "."
out_dir <- dirname(normalizePath(.this_file))
source(file.path(out_dir, "00_config.R"))

LOGFCT <- 0.25; FDR_CUT <- 0.05

## ---- Reuse the donor detection / pseudobulk / caching logic of script 68 ----
detect_donor <- function(seu) {
  meta <- seu@meta.data
  get <- function(f) if (f %in% colnames(meta)) as.character(meta[[f]]) else NULL
  patient <- get("patient"); sample <- get("sample")
  if (!is.null(patient)) {
    if (!is.null(sample)) { na_idx <- is.na(patient) | patient == ""; patient[na_idx] <- sample[na_idx] }
    return(patient)
  }
  for (cand in c("donor_id","donor","subject","subject_id","Subject_Identity","sample")) {
    v <- get(cand); if (!is.null(v)) return(v)
  }
  NULL
}
build_pb <- function(seu, keep, tag = "") {
  meta <- seu@meta.data; donor_vec <- detect_donor(seu)
  if (is.null(donor_vec)) stop(sprintf("[%s] no donor field", tag))
  keep[is.na(keep)] <- FALSE
  counts_all <- GetAssayData(seu, assay = "RNA", layer = "counts")
  donor_all <- as.character(donor_vec); names(donor_all) <- rownames(meta)
  cells <- rownames(meta)[keep]; cells <- cells[cells %in% colnames(counts_all)]
  donor_vec <- donor_all[cells]
  meta <- meta[cells, , drop = FALSE]
  counts <- counts_all[, cells, drop = FALSE]
  donors <- unique(donor_vec); donors <- donors[!is.na(donors) & donors != ""]
  cts_list <- lapply(donors, function(d) Matrix::rowSums(counts[, names(donor_vec)[donor_vec == d], drop = FALSE]))
  cts <- do.call(cbind, cts_list); colnames(cts) <- donors
  grp <- vapply(donors, function(d) unique(meta$group[donor_vec == d])[1], character(1))
  list(counts = cts, grp = grp, donors = donors)
}
prep_cache <- function(pb, ref, target) {
  grp_f <- factor(pb$grp, levels = c(ref, target))
  y0 <- DGEList(counts = pb$counts, group = grp_f)
  keep <- filterByExpr(y0, min.count = 10, min.total.count = 15)
  y0 <- y0[keep, , keep.lib.sizes = FALSE]
  y0 <- calcNormFactors(y0)
  list(keep = keep, nf = y0$samples$norm.factors, levels = c(ref, target))
}
deg_arms <- function(pb, cache, grp_vec) {
  grp_f <- factor(grp_vec, levels = cache$levels)
  if (length(unique(grp_f)) < 2 || min(table(grp_f)) < 2) return(NULL)
  y <- DGEList(counts = pb$counts[cache$keep, , drop = FALSE], group = grp_f, norm.factors = cache$nf)
  design <- model.matrix(~ grp_f)
  y <- estimateDisp(y, design); fit <- glmQLFit(y, design)
  tab <- topTags(glmQLFTest(fit, coef = 2), n = Inf)$table
  list(tab = tab, up = rownames(tab)[tab$logFC > LOGFCT & tab$FDR < FDR_CUT],
       down = rownames(tab)[tab$logFC < -LOGFCT & tab$FDR < FDR_CUT])
}

## ---- Part A: read blood object + donor-level details (answers question 2) ----
blood <- readRDS(file.path(out_dir, "path1_sepsis_copd_integrated.rds"))
meta <- blood@meta.data
donor <- detect_donor(blood)
dd <- data.frame(donor = donor, group = as.character(meta$group),
                 stringsAsFactors = FALSE)
dn_cnt <- dd[!duplicated(paste(dd$donor, dd$group)), ]
write.csv(as.data.frame(table(dn_cnt$group)), file.path(out_dir, "69_donor_group_counts.csv"),
          row.names = FALSE)
cat("---- 69-2 donor-level group counts ----\n"); print(table(dn_cnt$group))
## If a second field is used for the five-group comparison (e.g. group5/five_group), list it as well
for (cand in c("group5", "five_group", "group_five", "Group")) {
  if (cand %in% colnames(meta)) {
    d5 <- data.frame(donor = donor, g5 = as.character(meta[[cand]]))
    d5 <- d5[!duplicated(paste(d5$donor, d5$g5)), ]
    cat("---- donor counts for group field", cand, "----\n"); print(table(d5$g5))
    write.csv(as.data.frame(table(d5$g5)), file.path(out_dir, paste0("69_donor_", cand, "_counts.csv")),
              row.names = FALSE)
  }
}

## ---- Part B: donor-level Spearman of SOCS3 vs the 65-gene IFN module within each group (question 1) ----
purple <- readLines(file.path(out_dir, "path13_sc_module_purple.txt"))
purple <- purple[nzchar(trimws(purple))]
purple <- intersect(purple, rownames(blood))
cat("\nGenes matched in the 65-gene module:", length(purple), "\n")
ex <- GetAssayData(blood, assay = "RNA", layer = "data")
mono <- which(as.character(blood$cell_type) == "Monocyte")
donor_mono <- detect_donor(blood)[mono]
## Donor level: mean SOCS3 vs the mean of the z-scored module 65 genes
donors_u <- unique(donor_mono); donors_u <- donors_u[!is.na(donors_u) & donors_u != ""]
score_mat <- Matrix::rowMeans(Matrix::t(ex[purple, mono, drop = FALSE]))
soc <- as.numeric(ex["SOCS3", mono])
grp_mono <- as.character(blood$group[mono])
get_donor_val <- function(vec, ids, d, fun = mean) {
  vapply(d, function(x) { ii <- which(ids == x); if (length(ii) == 0) NA_real_ else fun(vec[ii]) }, numeric(1))
}
zscore <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
df <- data.frame(donor = donor_mono, grp = grp_mono, socs3 = soc, mod = score_mat)
dflat <- aggregate(cbind(socs3, mod) ~ donor + grp, data = df, FUN = mean)
dflat$mod <- zscore(dflat$mod); dflat$socs3 <- zscore(dflat$socs3)
perm_p <- function(x, y, nperm = 20000) {
  obs <- suppressWarnings(cor(x, y, method = "spearman"))
  set.seed(1)
  cnt <- 0L
  for (k in 1:nperm) { cnt <- cnt + (suppressWarnings(cor(x, sample(y), method = "spearman")) >= obs) }
  (cnt + 1) / (nperm + 1)
}
res <- do.call(rbind, lapply(sort(unique(dflat$grp)), function(g) {
  sub <- dflat[dflat$grp == g, ]
  if (nrow(sub) < 3) return(data.frame(group = g, n = nrow(sub), rho = NA, p_two_sided = NA))
  rho <- suppressWarnings(cor(sub$socs3, sub$mod, method = "spearman"))
  p <- if (nrow(sub) <= 10) perm_p(sub$socs3, sub$mod) else
    suppressWarnings(cor.test(sub$socs3, sub$mod, method = "spearman")$p.value)
  data.frame(group = g, n = nrow(sub), rho = round(rho, 4), p_two_sided = signif(p, 3))
}))
write.csv(res, file.path(out_dir, "69_per_group_spearman.csv"), row.names = FALSE)
cat("\n---- 69-1 per-group Spearman (SOCS3 vs IFN module, donor level) ----\n"); print(res)

## ---- Part C: observed intersection gene table (question 3, rebuild ST4) ----
b_keep <- blood$cell_type == "Monocyte" & blood$group %in% c("Healthy", "Sepsis")
pb_blood <- build_pb(blood, b_keep, tag = "blood monocytes")
lung <- readRDS(file.path(out_dir, "path2_copd_lung.rds"))
if (!"group" %in% colnames(lung@meta.data)) lung$group <- lung$disease
ct_all <- unique(as.character(lung$cell_type))
myeloid_types <- grep("Macrophage|Monocyte|monocyte|macrophage|Myeloid|myeloid|DC|dendritic", ct_all, value = TRUE)
if (length(myeloid_types) == 0) myeloid_types <- intersect(c("Myeloid", "Myeloid cells"), ct_all)
l_keep <- lung$cell_type %in% myeloid_types & lung$group %in% c("Control", "COPD")
pb_lung <- build_pb(lung, l_keep, tag = "lung myeloid")

cache_b <- prep_cache(pb_blood, "Healthy", "Sepsis")
cache_l <- prep_cache(pb_lung,  "Control", "COPD")
b <- deg_arms(pb_blood, cache_b, pb_blood$grp)
l <- deg_arms(pb_lung,  cache_l, pb_lung$grp)
up <- intersect(b$up, l$up); dn <- intersect(b$down, l$down)
tb <- b$tab; tl <- l$tab
out_up <- data.frame(arm = "up_both", gene = up,
                     blood_logFC = tb[up, "logFC"], blood_FDR = tb[up, "FDR"],
                     lung_logFC  = tl[up, "logFC"], lung_FDR  = tl[up, "FDR"],
                     stringsAsFactors = FALSE)
out_dn <- data.frame(arm = "down_both", gene = dn,
                     blood_logFC = tb[dn, "logFC"], blood_FDR = tb[dn, "FDR"],
                     lung_logFC  = tl[dn, "logFC"], lung_FDR  = tl[dn, "FDR"],
                     stringsAsFactors = FALSE)
out_g <- rbind(out_up, out_dn)
write.csv(out_g, file.path(out_dir, "69_observed_intersection_genes.csv"), row.names = FALSE)
cat(sprintf("\n---- 69-3 observed intersection genes: up=%d down=%d (%d rows total), written to 69_observed_intersection_genes.csv ----\n",
            length(up), length(dn), nrow(out_g)))
cat("===== Script 69 done =====\n")
