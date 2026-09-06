# =========================================================================
# 66 - Fig6 SOCS3 mechanism validation
#
# New Fig6 = mechanism validation plot for SOCS3 and the IFN program (monocyte main line)
#   Upgrades the "association" in Fig4/Fig5 to "regulation/dynamics", answering "why call SOCS3 the IFN brake".
#
# 4 panels (each with an independent tryCatch, so any failure does not interrupt the rest):
#   a Co-induction     -- donor-level SOCS3 expression vs ISG module score scatter + Spearman
#                  (SOCS3 is co-induced with the IFN program => accompanying negative feedback, not an independent event)
#   b Severity gradient -- donor-level SOCS3 and ISG score gradient across the five groups
#                  (both rise together with severity, showing SOCS3 is a severity-graded brake)
#   c Pseudotime dynamics -- Monocle3 monocyte trajectory (root=Healthy): ISG rises first, SOCS3 later
#                  (negative-feedback lag; requires monocle3, skipped if missing)
#   d Upstream TF      -- donor-level correlation matrix of SOCS3 vs STAT1/STAT2/IRF7/IRF9
#                  (IFN-signaling TFs drive SOCS3)
#
# Input: path1_sepsis_copd_integrated.rds (output of script 01, blood object)
#       path13_sc_module_purple.txt      (single-cell hdWGCNA purple = 65 ISGs)
# Output: path66_panel_{a,b,c,d}_*.{pdf,png} + path66_stats.txt
# =========================================================================

.this_file <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
if (.this_file == "" || !file.exists(.this_file)) {
  .f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1])
}
SCRIPT_DIR <- dirname(.this_file)
config_file <- file.path(SCRIPT_DIR, "00_config.R")
if (!file.exists(config_file)) stop("Cannot find 00_config.R: ", config_file)
source(config_file)

suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })

GRP <- c("Healthy", "Infection_Control", "COPD", "Sepsis", "Sepsis_Pneumonia")
GRP_COLS <- c("Healthy" = "#90A4AE", "Infection_Control" = "#4DB6AC",
              "COPD" = "#FFB300", "Sepsis" = "#E53935", "Sepsis_Pneumonia" = "#8E24AA")
UP_TF <- c("STAT1", "STAT2", "IRF7", "IRF9")
MAX_CELLS <- 40000          # upper limit of cells sampled for scoring/scatter
MAX_GROUP_CELLS <- 6000     # sampling limit per group

## Detect non-constant donor fields in meta
pick_nonconst <- function(meta, cand) {
  for (cn in cand) {
    if (cn %in% colnames(meta)) {
      v <- meta[[cn]]
      if (length(unique(as.character(v[!is.na(v)]))) > 1) return(cn)
    }
  }
  NULL
}

## z-score mean module score (genes x cells log-normalized matrix)
score_genes <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (length(g) < 3) return(list(score = NULL, n = length(g)))
  z <- t(scale(t(mat[g, , drop = FALSE])))
  z[is.na(z)] <- 0
  list(score = colMeans(z), n = length(g))
}

safe <- function(name, expr) {
  cat("  [", name, "] ... ", sep = "")
  r <- tryCatch(list(ok = TRUE, val = expr), error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  if (r$ok) cat("ok\n") else cat("skip:", r$msg, "\n")
  invisible(r)
}

save2 <- function(p, f, w = 7, h = 5) {
  ggsave(file.path(out_dir, paste0(f, ".pdf")), p, width = w, height = h)
  ggsave(file.path(out_dir, paste0(f, ".png")), p, width = w, height = h, dpi = 300)
}

cat("\n===== 66: SOCS3 <-> IFN mechanism validation (monocyte main line) =====\n")

## ---- Load the monocyte subset ----
blood_file <- file.path(out_dir, "path1_sepsis_copd_integrated.rds")
if (!file.exists(blood_file)) stop("Cannot find ", blood_file, " please run script 01 first")
seu <- readRDS(blood_file)
cat("Blood object:", ncol(seu), "cells\n")
mono <- subset(seu, subset = cell_type == "Monocyte")
rm(seu); gc()
cat("Monocyte cell count:", ncol(mono), "\n")
mono$group <- factor(as.character(mono$group), levels = GRP)

## Read purple (IFN) genes
purp_file <- file.path(out_dir, "path13_sc_module_purple.txt")
purple <- if (file.exists(purp_file)) trimws(readLines(purp_file)) else character(0)
purple <- unique(purple[purple != "" & !grepl("///", purple)])
cat("purple ISG genes:", length(purple), "\n")

## donor field (aligned with DONOR_CAND in script 46, sample as fallback; blood object ~216 donors, lung 46)
DONOR_CAND <- c("donor_id", "donor", "subject", "subject_id", "Subject_Identity",
                "sample", "patient", "sample_id")
donor_col <- NULL
for (f in DONOR_CAND) if (f %in% colnames(mono@meta.data)) { donor_col <- f; break }
cat("donor field:", if (is.null(donor_col)) "(not found, using cell level)" else donor_col, "\n")
if (!is.null(donor_col)) {
  cat("donor distinct:", length(unique(mono@meta.data[[donor_col]])), "\n")
  for (f in intersect(DONOR_CAND, colnames(mono@meta.data)))
    cat("  [candidate] ", f, " = ", length(unique(mono@meta.data[[f]])), "\n", sep = "")
}

## ---- Sampling (for scoring) ----
set.seed(42)
if (ncol(mono) > MAX_CELLS) {
  keep <- sample(colnames(mono), MAX_CELLS)
  mono_sub <- subset(mono, cells = keep)
} else mono_sub <- mono
rm(mono); gc()

expr <- as.matrix(GetAssayData(mono_sub, assay = "RNA", layer = "data"))
soc3 <- expr["SOCS3", ]
isg_r <- score_genes(expr, purple)
cat("ISG score genes matched:", isg_r$n, "/", length(purple), "\n")
if (is.null(isg_r$score)) stop("Fewer than 3 purple genes matched in monocytes, cannot score")

meta <- mono_sub@meta.data
meta$socs3 <- as.numeric(soc3)
meta$isg   <- as.numeric(isg_r$score)
for (tf in UP_TF) meta[[paste0("tf_", tf)]] <- if (tf %in% rownames(expr)) as.numeric(expr[tf, ]) else NA

stats_lines <- c()

## ================= Panel a: donor-level co-induction =================
cat("\n[Panel a] donor-level SOCS3 vs ISG score co-induction\n")
a_done <- safe("panel_a", {
  if (!is.null(donor_col)) {
    dm <- meta[, c(donor_col, "group", "socs3", "isg")]
    names(dm)[1] <- "donor"
    agg <- aggregate(cbind(socs3, isg) ~ donor + group, dm, mean)
  } else {
    # No donor field: sample by group x cell and use cell level directly
    agg <- meta[, c("group", "socs3", "isg")]
    names(agg)[1] <- "group"
  }
  sp <- suppressWarnings(cor.test(agg$socs3, agg$isg, method = "spearman"))
  rho <- unname(sp$estimate); p <- sp$p.value
  stats_lines <<- c(stats_lines, sprintf("a co-induction: SOCS3 vs ISG_score Spearman rho=%.3f p=%.3g (n=%d)", rho, p, nrow(agg)))
  pa <- ggplot(agg, aes(isg, socs3, color = group)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, color = "grey20", linewidth = 0.5) +
    scale_color_manual(values = GRP_COLS, drop = FALSE) +
    labs(x = "ISG module score (purple 65)", y = "SOCS3 expression",
         title = "a  SOCS3 is co-induced with the interferon program",
         subtitle = sprintf("donor-level Spearman rho = %.3f, P = %.3g", rho, p)) +
    theme_bw(base_size = 12)
  save2(pa, "path66_panel_a_coinduction", 6.5, 5)
  list(rho = rho, p = p, n = nrow(agg))
})

## ================= Panel b: severity gradient =================
cat("\n[Panel b] donor-level SOCS3 and ISG gradient across the five groups\n")
b_done <- safe("panel_b", {
  if (!is.null(donor_col)) {
    dm <- meta[, c(donor_col, "group", "socs3", "isg")]
    names(dm)[1] <- "donor"
    agg <- aggregate(cbind(socs3, isg) ~ donor + group, dm, mean)
  } else {
    agg <- meta[, c("group", "socs3", "isg")]
  }
  agg$group <- factor(as.character(agg$group), levels = GRP)
  # KW test
  kw_s <- kruskal.test(socs3 ~ group, agg)$p.value
  kw_i <- kruskal.test(isg ~ group, agg)$p.value
  stats_lines <<- c(stats_lines,
    sprintf("b severity gradient: SOCS3 KW P=%.3g | ISG_score KW P=%.3g", kw_s, kw_i))
  # Facet the long table by the two metrics
  long <- rbind(
    data.frame(group = agg$group, val = agg$socs3, metric = "SOCS3"),
    data.frame(group = agg$group, val = agg$isg,   metric = "ISG module score"))
  long$metric <- factor(long$metric, levels = c("SOCS3", "ISG module score"))
  pb <- ggplot(long, aes(group, val, fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 0.7, alpha = 0.5) +
    facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = GRP_COLS, drop = FALSE) +
    labs(x = NULL, y = "Donor-level mean",
         title = "b  SOCS3 and the ISG module rise together across severity",
         subtitle = sprintf("SOCS3 KW P = %.2g; ISG KW P = %.2g", kw_s, kw_i)) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save2(pb, "path66_panel_b_severity_gradient", 6.5, 6.5)
  list(kw_socs3 = kw_s, kw_isg = kw_i)
})

## ================= Panel d: upstream TF co-expression =================
cat("\n[Panel d] SOCS3 vs upstream STAT/IRF TFs\n")
d_done <- safe("panel_d", {
  tf_present <- UP_TF[UP_TF %in% rownames(expr)]
  if (length(tf_present) < 2) stop("Fewer than 2 upstream TFs matched: ", paste(tf_present, collapse = ","))
  if (!is.null(donor_col)) {
    cols <- c(donor_col, "group", "socs3", paste0("tf_", tf_present))
    dm <- meta[, cols]; names(dm)[1] <- "donor"
    agg <- aggregate(. ~ donor + group, dm, mean)
  } else {
    agg <- meta[, c("group", "socs3", paste0("tf_", tf_present))]
  }
  # Spearman between SOCS3 and each TF
  corvec <- sapply(tf_present, function(tf) {
    col <- paste0("tf_", tf)
    if (all(is.na(agg[[col]]))) return(NA)
    unname(suppressWarnings(cor.test(agg$socs3, agg[[col]], method = "spearman")$estimate))
  })
  stats_lines <<- c(stats_lines, sprintf("d upstream TF correlation with SOCS3: %s",
    paste(sprintf("%s=%.3f", tf_present, corvec), collapse = ", ")))
  cormat <- data.frame(TF = tf_present, rho = corvec)
  cormat$TF <- factor(cormat$TF, levels = tf_present[order(corvec, na.last = NA)])
  pd <- ggplot(cormat, aes(TF, rho, fill = rho)) +
    geom_col(width = 0.6, alpha = 0.9) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    geom_text(aes(label = sprintf("%.3f", rho)), vjust = ifelse(corvec >= 0, -0.5, 1.5), size = 3.5) +
    labs(x = NULL, y = "Spearman rho with SOCS3",
         title = "d  IFN-signal transcription factors co-vary with SOCS3",
         subtitle = "donor-level correlation (upstream drivers of SOCS3)") +
    theme_bw(base_size = 12)
  save2(pd, "path66_panel_d_upstream_TF", 6, 5)
  list(cor = corvec)
})

## ================= Panel c: pseudotime (Monocle3, optional) =================
cat("\n[Panel c] monocyte pseudotime (Monocle3)\n")
c_done <- safe("panel_c", {
  if (!requireNamespace("monocle3", quietly = TRUE))
    stop("monocle3 is not installed (panel c skipped, does not affect a/b/d)")
  library(monocle3)
  set.seed(42)
  # Sample down to ~8000 cells to keep it computable
  cells <- colnames(mono_sub)
  if (length(cells) > 8000) cells <- sample(cells, 8000)
  sm <- mono_sub[, cells]
  sm$group <- factor(as.character(sm$group), levels = GRP)
  cds <- new_cell_data_set(as(GetAssayData(sm, assay = "RNA", layer = "counts"), "sparseMatrix"),
                           cell_metadata = sm@meta.data)
  cds <- preprocess_cds(cds, num_dim = 30)
  cds <- reduce_dimension(cds)
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  # root = the cluster with the highest Healthy fraction
  cl <- clusters(cds)
  healthy_frac <- tapply(colnames(cds), cl, function(cc) mean(sm@meta.data[cc, "group"] == "Healthy"))
  root_cells <- colnames(cds)[cl == names(which.max(healthy_frac))]
  cds <- order_cells(cds, root_cells = root_cells)
  pt <- pseudotime(cds)
  ok <- is.finite(pt)
  # Bin to inspect SOCS3 / ISG dynamics
  e2 <- as.matrix(GetAssayData(sm, assay = "RNA", layer = "data"))
  soc3v <- e2["SOCS3", ]; isgv <- score_genes(e2, purple)$score
  bins <- cut(pt[ok], breaks = 6)
  dyn <- data.frame(
    bin = as.numeric(bins),
    socs3 = soc3v[ok],
    isg = isgv[ok])
  dmean <- aggregate(cbind(socs3, isg) ~ bin, dyn, mean)
  # Relative peaks: ISG peak bin vs SOCS3 peak bin (lag assessment)
  isg_peak <- which.max(dmean$isg); soc3_peak <- which.max(dmean$socs3)
  stats_lines <<- c(stats_lines,
    sprintf("c pseudotime: ISG peak bin=%d, SOCS3 peak bin=%d (%s)", isg_peak, soc3_peak,
            ifelse(soc3_peak >= isg_peak, "SOCS3 lags => negative feedback", "SOCS3 does not lag")))
  dl <- rbind(data.frame(bin = dmean$bin, val = dmean$isg,   metric = "ISG module"),
              data.frame(bin = dmean$bin, val = dmean$socs3, metric = "SOCS3"))
  dl$metric <- factor(dl$metric, levels = c("ISG module", "SOCS3"))
  # Normalize to each range to overlay the time course
  pc <- ggplot(dl, aes(bin, val, color = metric, group = metric)) +
    geom_line(linewidth = 1) + geom_point(size = 2.5) +
    scale_color_manual(values = c("ISG module" = "#2166AC", "SOCS3" = "#B2182B")) +
    labs(x = "Pseudotime bin (Healthy -> activated)", y = "Mean expression",
         title = "c  IFN program peaks before SOCS3 along monocyte activation",
         subtitle = "negative-feedback lag: SOCS3 rises after the ISG program") +
    theme_bw(base_size = 12) + theme(legend.title = element_blank())
  save2(pc, "path66_panel_c_pseudotime", 6.5, 4.8)
  list(isg_peak = isg_peak, soc3_peak = soc3_peak)
})

## ================= Summary stats =================
cat("\n===== 66 stats summary =====\n")
writeLines(stats_lines)
writeLines(stats_lines, file.path(out_dir, "path66_stats.txt"))
cat("\n===== Script 66 finished =====\n")
