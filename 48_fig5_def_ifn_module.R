# =========================================================================
# 48 - Fig5 D/E/F: interferon module
#
# Purpose: panels D / E / F of Fig5
#   D  Membership overlap of the IFN module between two orthogonal methods
#      (sc hdWGCNA purple 65 ∩ bulk WGCNA magenta 103)
#      -> 34 genes, expected 0.33, hypergeometric P ~ 1.1e-62 (recomputed and printed in this script)
#   E  purple module activity (module eigengene ME) at donor level x five groups
#      -> whether this interferon module activity rises with severity
#   F  Coupling of SOCS3 with purple ME (donor-level Spearman; cell level added)
#      -> echoes Fig4: whether the brake (SOCS3) and the target (IFN module) co-occur
#
# Input:
#   D : path13_sc_gene_module.csv / path13_sc_module_purple.txt (single-cell modules)
#       path12_bulk_gene_module.csv (bulk modules)
#   E/F: path13_sc_hdwgcna.rds (preferred, contains meta+MEs, monocytes only) or
#        path1_sepsis_copd_integrated.rds (fallback, subset monocytes, MEs aligned from CSV)
#
# Output (Fig5D_* / Fig5E_* / Fig5F_*, PDF+PNG):
#   Fig5D_module_overlap.pdf/.png          membership alignment strip + overlap statistics
#   Fig5D_overlap_genes.csv                 shared gene list
#   Fig5E_purple_ME_by_group.pdf/.png       donor-level purple ME x five groups
#   Fig5F_SOCS3_vs_purpleME.pdf/.png        donor-level SOCS3 vs purple ME scatter
#   path48_fig5def_stats.txt                statistics summary
#   path48_fig5def_plotdata.rds             pre-plot data (for local re-rendering)
#
# Server: Rscript this script
# Requires: ggplot2/dplyr/tidyr; run 12_WGCNA.R and the single-cell hdWGCNA script (13) first
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

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})
suppressPackageStartupMessages(library(Seurat))   # for donor-level aggregation (only meta/assays are used after loading the object)
has_hdwgcna <- requireNamespace("hdWGCNA", quietly = TRUE)

## Colors for the five groups (consistent with Fig2-F/G and Fig4)
GROUP_LEVELS5 <- c("Healthy", "Infection_Control", "COPD", "Sepsis", "Sepsis_Pneumonia")
GROUP_COLS5   <- c(Healthy = "#3C5488", Infection_Control = "#91D1C2",
                   COPD = "#00A087", Sepsis = "#E64B35", Sepsis_Pneumonia = "#F39B7F")

## Donor field detection (patient first, sample as fallback) [same as 28b/31]
detect_donor <- function(meta) {
  get <- function(f) if (f %in% colnames(meta)) as.character(meta[[f]]) else NULL
  patient <- get("patient"); sample <- get("sample")
  if (!is.null(patient)) {
    if (!is.null(sample)) {
      na_idx <- is.na(patient) | patient == ""
      patient[na_idx] <- sample[na_idx]
    }
    return(list(field = "patient(+sample fallback)", vec = patient))
  }
  for (cand in c("donor_id","donor","subject","subject_id","sample")) {
    v <- get(cand)
    if (!is.null(v)) return(list(field = cand, vec = v))
  }
  NULL
}

cat("\n==============================================================\n")
cat("Fig5 D/E/F: cross-method overlap and activity of the interferon module\n")
cat("==============================================================\n")

stats_lines <- c("===== Fig5 D/E/F statistics summary =====")

## ===========================================================================
## D - module membership overlap between the two methods
## ===========================================================================
cat("\n########## D: sc purple ∩ bulk magenta ##########\n")

# Single-cell purple members (module CSV preferred, txt fallback)
purple <- character(0)
sc_mod <- file.path(out_dir, "path13_sc_gene_module.csv")
if (file.exists(sc_mod)) {
  m <- read.csv(sc_mod, check.names = FALSE, stringsAsFactors = FALSE)
  gcol <- intersect(c("gene_name", "gene", "gene_id"), colnames(m))[1]
  if (!is.na(gcol) && "module" %in% colnames(m)) {
    purple <- as.character(m[[gcol]][m$module == "purple"])
  }
}
if (length(purple) == 0) {
  f_txt <- file.path(out_dir, "path13_sc_module_purple.txt")
  if (file.exists(f_txt)) purple <- trimws(readLines(f_txt))
}
purple <- unique(purple[purple != "" & !grepl("///", purple)])
cat("sc hdWGCNA purple gene count:", length(purple), "\n")

# bulk magenta members
bulk_mod <- file.path(out_dir, "path12_bulk_gene_module.csv")
if (!file.exists(bulk_mod)) stop("Missing path12_bulk_gene_module.csv, please run 12_WGCNA.R first")
bm <- read.csv(bulk_mod, check.names = FALSE, stringsAsFactors = FALSE)
gcol_b <- intersect(c("gene", "gene_name", "gene_id"), colnames(bm))[1]
if (is.na(gcol_b) || !"module" %in% colnames(bm)) stop("bulk module CSV column names do not match")
magenta <- unique(as.character(bm[[gcol_b]][bm$module == "magenta"]))
magenta <- magenta[magenta != "" & !is.na(magenta) & !grepl("///", magenta)]
cat("bulk WGCNA magenta gene count:", length(magenta), "\n")

ov <- intersect(purple, magenta)
only_p <- setdiff(purple, ov); only_m <- setdiff(magenta, ov)
cat("overlapping genes:", length(ov), " (expected 0.33)\n")

# Hypergeometric test: background ~20000 expressed genes
BG <- 20000L
p_hyper <- phyper(length(ov) - 1L, length(purple), BG - length(purple),
                  length(magenta), lower.tail = FALSE)
cat(sprintf("Hypergeometric P(X>=%d) = %.3e\n", length(ov), p_hyper))
stats_lines <- c(stats_lines,
  sprintf("D overlap: purple(%d) ∩ magenta(%d) = %d; expected=%.2f; hypergeom P=%.2e",
          length(purple), length(magenta), length(ov),
          length(purple) * length(magenta) / BG, p_hyper))

# Write shared genes sorted alphabetically
ov_sorted <- sort(ov)
writeLines(ov_sorted, file.path(out_dir, "Fig5D_overlap_genes.txt"))
write.csv(data.frame(gene = ov_sorted, module = "purple_and_magenta"),
          file.path(out_dir, "Fig5D_overlap_genes.csv"), row.names = FALSE)

# ---- Membership alignment strip plot (three rows: segment ribbon / single-cell row / bulk row) ----
if (length(ov) > 0) {
  union_all <- c(ov_sorted, sort(only_p), sort(only_m))
  ord <- factor(union_all, levels = union_all)

  seg_col <- c(rep("#7B3294", length(ov_sorted)),
               rep("#A06CD5", length(only_p)),
               rep("#E377C2", length(only_m)))
  ylevs <- c("single-cell hdWGCNA (purple)", "bulk WGCNA (magenta)", "segment")
  seg_df <- data.frame(gene = ord,
                       y = factor("segment", levels = ylevs),
                       fill = seg_col)
  seg_lab <- data.frame(
    xmid = c((1 + length(ov_sorted)) / 2,
             length(ov_sorted) + (1 + length(only_p)) / 2,
             length(ov_sorted) + length(only_p) + (1 + length(only_m)) / 2),
    label = paste0(c("Overlap (", "purple-only (", "magenta-only ("),
                   c(length(ov_sorted), length(only_p), length(only_m)), ")"))

  d_tile <- rbind(
    data.frame(gene = ord, method = "single-cell hdWGCNA (purple)",
               present = as.integer(union_all %in% purple)),
    data.frame(gene = ord, method = "bulk WGCNA (magenta)",
               present = as.integer(union_all %in% magenta))
  )
  p_d <- ggplot(d_tile, aes(gene, method)) +
    geom_tile(data = seg_df, aes(x = gene, y = y), fill = seg_col,
              width = 0.92, height = 0.55) +
    annotate("text", x = seg_lab$xmid, y = "segment", label = seg_lab$label,
             size = 2.7, colour = "white", fontface = "bold") +
    geom_tile(aes(fill = method, alpha = present),
              width = 0.92, height = 0.72) +
    scale_fill_manual(values = c("single-cell hdWGCNA (purple)" = "#A06CD5",
                                 "bulk WGCNA (magenta)" = "#E377C2"),
                      name = NULL) +
    scale_alpha_continuous(range = c(0, 1), guide = "none") +
    scale_y_discrete(limits = ylevs) +
    labs(title = "IFN module membership: single-cell hdWGCNA vs bulk WGCNA",
         subtitle = sprintf("purple (n=%d)  ∩  magenta (n=%d)  =  %d shared genes;  expected 0.33;  hypergeometric P = %.1e",
                            length(purple), length(magenta), length(ov), p_hyper),
         caption = "Full shared-gene list: Fig5D_overlap_genes.txt / .csv (e.g. ISG15, MX1, STAT1, IRF7, IFIT1-3, GBP1-5, OAS2/3, IFI44L, IFI44, IFI6, UBE2L6)",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          plot.subtitle = element_text(size = 8.5, colour = "grey30"),
          plot.caption = element_text(size = 7, colour = "grey45"))
  ggsave(file.path(out_dir, "Fig5D_module_overlap.pdf"), p_d,
         width = 11, height = 3.2)
  ggsave(file.path(out_dir, "Fig5D_module_overlap.png"), p_d,
         width = 11, height = 3.2, dpi = 300, device = grDevices::png)
  cat("Saved: Fig5D_module_overlap.{pdf,png}\n")
} else {
  cat("!! No overlapping genes, panel D skipped\n")
}

## ===========================================================================
## E/F - load the monocyte object (hdWGCNA rds preferred), extract purple ME + meta
## ===========================================================================
cat("\n########## E/F: purple ME activity ##########\n")

load_mes <- function() {
  hd <- file.path(out_dir, "path13_sc_hdwgcna.rds")
  if (file.exists(hd)) {
    cat("Reading path13_sc_hdwgcna.rds\n")
    seu <- readRDS(hd)
    meta <- seu@meta.data
    # meta must contain cell_type/group; otherwise align via CSV
    if (!all(c("cell_type") %in% colnames(meta))) {
      cat("hdWGCNA rds meta lacks cell_type, switching to path1 + CSV MEs\n")
      return(NULL)
    }
    mes <- NULL
    if (has_hdwgcna) {
      mes <- tryCatch(as.data.frame(hdWGCNA::GetMEs(seu, harmonized = TRUE)),
                      error = function(e) NULL)
    }
    if (is.null(mes) || nrow(mes) == 0) {
      csv_mes <- file.path(out_dir, "path13_sc_MEs.csv")
      if (file.exists(csv_mes)) {
        mes <- read.csv(csv_mes, row.names = 1, check.names = FALSE)
        cat("MEs read from path13_sc_MEs.csv:", nrow(mes), "cells\n")
      }
    }
    if (!is.null(mes)) {
      rownames(mes) <- gsub("^\"|\"$", "", rownames(mes))
      return(list(seu = seu, meta = meta, mes = mes))
    }
    return(NULL)
  }
  NULL
}

obj <- load_mes()
if (is.null(obj)) {
  cat("!! No usable hdWGCNA object/CSV MEs found. Please run the single-cell hdWGCNA script (13) first (saving path13_sc_hdwgcna.rds and path13_sc_MEs.csv)\n")
} else {
  meta <- obj$meta; mes <- obj$mes
  cat("object cell count:", nrow(meta), "| MEs row count:", nrow(mes), "\n")

  # Align MEs rownames with meta rownames
  rownames(meta) <- gsub("^\"|\"$", "", rownames(meta))
  common <- intersect(rownames(mes), rownames(meta))
  cat("MEs aligned directly with meta:", length(common), "cells\n")
  if (length(common) < 1000 && !is.null(obj$seu)) {
    cn <- colnames(obj$seu); names(cn) <- NULL
    common2 <- intersect(rownames(mes), cn)
    if (length(common2) > length(common)) {
      cat("Switching to colnames alignment:", length(common2), "\n")
      rownames(meta) <- cn
      common <- common2
    }
  }
  if (length(common) < 1000) {
    cat("!! MEs cannot be aligned with object cells, E/F skipped (check path13_sc_MEs.csv rownames)\n")
  } else {
    if (!"purple" %in% colnames(mes)) {
      cat("!! MEs have no purple column, available modules:", paste(colnames(mes), collapse = ", "), "\n")
    } else {
      d_ef <- data.frame(cell = common,
                         purple_ME = as.numeric(mes[common, "purple"]),
                         group = as.character(meta[common, "group"]),
                         cell_type = as.character(meta[common, "cell_type"]),
                         stringsAsFactors = FALSE)
      d_ef$group[is.na(d_ef$group) | d_ef$group == ""] <- "Healthy"
      d_ef$group[!d_ef$group %in% GROUP_LEVELS5] <- "Healthy"   # fallback (Pediatric, etc.)
      d_ef <- d_ef[d_ef$cell_type %in% c("Monocyte", "DC"), ]
      cat("Cells analyzed (Monocyte/DC):", nrow(d_ef), "\n")

      don <- detect_donor(meta)
      if (!is.null(don)) {
        d_ef$donor <- don$vec[match(d_ef$cell, rownames(meta))]
        cat("Donor field:", don$field, " donor count:", length(unique(d_ef$donor)), "\n")
      } else {
        d_ef$donor <- d_ef$cell
        cat("!! No donor field, donor level falls back to cells\n")
      }

      # SOCS3 expression (log-normalized), from the object data layer (if the object holds expression)
      socs3_expr <- NULL
      if (!is.null(obj$seu)) {
        tryCatch({
          dat <- tryCatch(GetAssayData(obj$seu, assay = "RNA", layer = "data"),
                          error = function(e) GetAssayData(obj$seu, assay = "RNA", slot = "data"))
          if ("SOCS3" %in% rownames(dat)) {
            socs3_expr <- as.numeric(dat["SOCS3", common])
            names(socs3_expr) <- common
          }
        }, error = function(e) cat("Failed to get SOCS3 expression:", conditionMessage(e), "\n"))
      }
      if (is.null(socs3_expr)) {
        cat("!! Object has no SOCS3 expression, panel F can only plot ME-based alternatives; path1 fallback skipped\n")
      }

      ## ---------- E: donor-level purple ME x five groups ----------
      d_don <- d_ef %>%
        filter(!is.na(purple_ME), !is.na(donor), donor != "") %>%
        group_by(donor, group) %>%
        summarise(purple_ME = mean(purple_ME), n_cell = n(), .groups = "drop") %>%
        filter(group %in% GROUP_LEVELS5)
      d_don$group <- factor(d_don$group, levels = GROUP_LEVELS5)
      cat("Donor-level purple ME summary: donor count =", nrow(d_don), "\n")
      print(table(d_don$group))

      if (nrow(d_don) > 0) {
        kw <- suppressWarnings(kruskal.test(purple_ME ~ group, d_don))
        wt <- lapply(GROUP_LEVELS5[-1], function(g) {
          x <- d_don$purple_ME[d_don$group == "Healthy"]
          y <- d_don$purple_ME[d_don$group == g]
          if (length(x) >= 3 && length(y) >= 3) {
            p <- suppressWarnings(wilcox.test(x, y)$p.value)
            data.frame(comparison = paste0(g, " vs Healthy"), p = p)
          } else NULL
        })
        wt_df <- do.call(rbind, wt[!vapply(wt, is.null, logical(1))])
        cat(sprintf("Kruskal-Wallis P = %.3e\n", kw$p.value))
        stats_lines <- c(stats_lines,
          sprintf("E purple ME: KW P=%.3e (n donors per group: %s)",
                  kw$p.value,
                  paste(paste0(names(table(d_don$group)), "=", as.integer(table(d_don$group))),
                        collapse = ", ")))
        if (!is.null(wt_df) && nrow(wt_df) > 0) {
          for (i in seq_len(nrow(wt_df)))
            stats_lines <- c(stats_lines, sprintf("   %s: Wilcoxon P=%.3e",
                                                  wt_df$comparison[i], wt_df$p[i]))
        }
        p_e <- ggplot(d_don, aes(group, purple_ME)) +
          geom_boxplot(aes(fill = group), outlier.shape = NA, width = 0.55,
                       alpha = 0.75, linewidth = 0.35) +
          geom_jitter(width = 0.12, size = 1.1, alpha = 0.5, colour = "grey25") +
          scale_fill_manual(values = GROUP_COLS5, guide = "none") +
          labs(title = "IFN module activity (purple ME) by group",
               subtitle = sprintf("Donor-level mean module eigengene; Kruskal-Wallis P = %.2e", kw$p.value),
               x = NULL, y = "purple module eigengene (donor mean)") +
          theme_bw(base_size = 11) +
          theme(axis.text.x = element_text(angle = 40, hjust = 1),
                panel.grid.minor = element_blank(),
                plot.subtitle = element_text(size = 8.5, colour = "grey30"))
        ggsave(file.path(out_dir, "Fig5E_purple_ME_by_group.pdf"), p_e,
               width = 5.6, height = 4.6)
        ggsave(file.path(out_dir, "Fig5E_purple_ME_by_group.png"), p_e,
               width = 5.6, height = 4.6, dpi = 300, device = grDevices::png)
        cat("Saved: Fig5E_purple_ME_by_group.{pdf,png}\n")
      }

      ## ---------- F: SOCS3 x purple ME (donor level) ----------
      if (!is.null(socs3_expr)) {
        d_ef$socs3 <- socs3_expr[d_ef$cell]
        d_don2 <- d_ef %>%
          filter(!is.na(purple_ME), !is.na(socs3), !is.na(donor), donor != "") %>%
          group_by(donor, group) %>%
          summarise(purple_ME = mean(purple_ME),
                    socs3 = mean(socs3),
                    socs3_pct = mean(socs3 > 0) * 100,
                    .groups = "drop") %>%
          filter(group %in% GROUP_LEVELS5)
        d_don2$group <- factor(d_don2$group, levels = GROUP_LEVELS5)
        rr <- suppressWarnings(cor.test(d_don2$socs3, d_don2$purple_ME, method = "spearman"))
        cat(sprintf("Donor-level SOCS3 ~ purple ME: Spearman rho=%.3f, P=%.3e (n=%d)\n",
                    rr$estimate, rr$p.value, nrow(d_don2)))
        stats_lines <- c(stats_lines,
          sprintf("F donor-level SOCS3 vs purple ME: rho=%.3f, P=%.3e (n donors=%d)",
                  rr$estimate, rr$p.value, nrow(d_don2)))

        p_f <- ggplot(d_don2, aes(purple_ME, socs3)) +
          geom_point(aes(colour = group), size = 2, alpha = 0.75) +
          geom_smooth(method = "lm", se = FALSE, colour = "grey40",
                      linewidth = 0.5, linetype = "dashed") +
          scale_colour_manual(values = GROUP_COLS5, name = NULL) +
          labs(title = "SOCS3 couples to IFN module activity",
               subtitle = sprintf("Donor-level mean; Spearman rho = %.3f, P = %.2e (n = %d)",
                                  rr$estimate, rr$p.value, nrow(d_don2)),
               x = "purple module eigengene (donor mean)",
               y = "SOCS3 expression (donor mean, log)") +
          theme_bw(base_size = 11) +
          theme(legend.position = "bottom",
                panel.grid.minor = element_blank(),
                plot.subtitle = element_text(size = 8.5, colour = "grey30"))
        ggsave(file.path(out_dir, "Fig5F_SOCS3_vs_purpleME.pdf"), p_f,
               width = 5.6, height = 5.2)
        ggsave(file.path(out_dir, "Fig5F_SOCS3_vs_purpleME.png"), p_f,
               width = 5.6, height = 5.2, dpi = 300, device = grDevices::png)
        cat("Saved: Fig5F_SOCS3_vs_purpleME.{pdf,png}\n")
      }

      ## ---------- Pre-plot data ----------
      plotdata <- list(
        d = list(purple = purple, magenta = magenta, overlap = ov,
                 p_hyper = p_hyper,
                 tile_df = if (length(ov) > 0) d_tile else NULL),
        e = list(donor_df = if (exists("d_don")) d_don else NULL,
                 kw_p = if (exists("kw")) kw$p.value else NULL),
        f = list(donor_df = if (exists("d_don2")) d_don2 else NULL,
                 rho = if (exists("rr")) rr$estimate else NULL,
                 p = if (exists("rr")) rr$p.value else NULL),
        group_levels = GROUP_LEVELS5, group_cols = GROUP_COLS5,
        stats = stats_lines,
        ggplots = list(d = if (length(ov) > 0) p_d else NULL,
                       e = if (exists("p_e")) p_e else NULL,
                       f = if (exists("p_f")) p_f else NULL)
      )
      saveRDS(plotdata, file.path(out_dir, "path48_fig5def_plotdata.rds"))
      sz <- round(file.info(file.path(out_dir, "path48_fig5def_plotdata.rds"))$size / 1e6, 2)
      cat("Saved path48_fig5def_plotdata.rds (", sz, " MB )\n")
    }
  }
}

writeLines(stats_lines, file.path(out_dir, "path48_fig5def_stats.txt"))
cat("\n----- statistics summary -----\n"); cat(paste(stats_lines, collapse = "\n"), "\n")
cat("\n===== Script 48 finished =====\n")
