# =============================================================================
# 62 FigS blood/lung QC annotation  (v2: panel-level fault tolerance + DotPlot large-object hardening)
# Purpose: add a supplementary figure — a single-cell "QC + annotation" two-sided version
#          (similar to Fig2):
#   blood, 6 datasets (cp3_annotated.rds) / lung GSE136831 (path2_copd_lung.rds)
#   each side: QC violin + celltype/dataset/disease UMAP + marker gene DotPlot
# Hardening:
#   - safe() wraps every panel; a single failure prints [skip] and continues instead of
#     crashing the whole script
#   - DotPlot first removes NA groupings and drops empty groups; if it still fails it
#     automatically downsamples by cluster and retries
#     (>100k cells: Seurat DotPlot has a known "$<- id NA / 0 row" problem)
# Outputs: path62_{blood,lung}_*.{pdf,png} + path62_meta_summary.txt
# Run: Rscript 62_blood_lung_QC_figures.R [OUT_DIR]
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
if (length(args) >= 1 && nzchar(args[1])) script_dir <- args[1]
out_dir <- if (length(args) >= 2 && nzchar(args[2])) args[2] else script_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

safe <- function(lbl, code) tryCatch({ code }, error = function(e) {
  cat("  [skip]", lbl, "-", conditionMessage(e), "\n"); NULL
})
save2 <- function(p, base, w = 8, h = 6) {
  if (is.null(p)) return(invisible(NULL))
  ggplot2::ggsave(file.path(out_dir, paste0(base, ".pdf")), p, width = w, height = h)
  ggplot2::ggsave(file.path(out_dir, paste0(base, ".png")), p, width = w, height = h, dpi = 300)
  cat("  [ok]", base, "\n")
}
pick <- function(seu, cand) {
  hit <- cand[cand %in% colnames(seu@meta.data)]
  if (length(hit) > 0) hit[1] else NULL
}
## take the first "non-constant" candidate field (avoid mistakenly picking an all-identical column such as sample)
pick_nonconst <- function(seu, cand) {
  cand <- cand[cand %in% colnames(seu@meta.data)]
  for (cc in cand) if (length(unique(seu@meta.data[[cc]])) > 1) return(cc)
  if (length(cand) > 0) cand[1] else NULL
}
## make sure the object has a "umap" reduction (for DimPlot); otherwise use a candidate reduction name or run RunUMAP automatically
ensure_umap <- function(seu, tag) {
  reds <- Reductions(seu)
  cat("  [", tag, "] Reductions:", paste(reds, collapse = ", "), "\n", sep = "")
  if ("umap" %in% reds) return(seu)
  cand_red <- intersect(c("umap","UMAP","tsne","pca","harmony"), reds)
  if (length(cand_red) > 0) {
    seu[["umap"]] <- seu[[cand_red[1]]]
    cat("  [", tag, "] using the existing reduction as umap:", cand_red[1], "\n", sep = "")
    return(seu)
  }
  src <- intersect(c("harmony","pca"), reds)
  if (length(src) > 0) {
    cat("  [", tag, "] no 2D reduction, trying RunUMAP(reduction=", src[1], ")\n", sep = "")
    seu <- tryCatch(RunUMAP(seu, reduction = src[1],
                            dims = 1:min(30, length(Embeddings(seu, src[1])[1, ])),
                            verbose = FALSE), error = function(e) seu)
    return(seu)
  }
  ## no reduction at all (e.g. the lung object was only processed with NormalizeData) -> run the full dimensionality reduction pipeline automatically
  cat("  [", tag, "] no dimensionality reduction at all, auto Normalize→VarFeatures→Scale→PCA→UMAP\n", sep = "")
  seu <- tryCatch({
    DefaultAssay(seu) <- "RNA"
    if (!"data" %in% Layers(seu)) seu <- NormalizeData(seu, verbose = FALSE)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
    seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 30, verbose = FALSE)
    seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
    cat("  [", tag, "] automatic dimensionality reduction done\n", sep = "")
    seu
  }, error = function(e) {
    cat("  [", tag, "] automatic dimensionality reduction failed:", conditionMessage(e), "\n", sep = "")
    seu
  })
  seu
}
qc_panel <- function(seu, grp) {
  have <- intersect(c("nFeature_RNA","nCount_RNA","percent.mt"),
                    colnames(seu@meta.data))
  if (length(have) < 2 || is.null(grp)) return(NULL)
  df <- as.data.frame(seu@meta.data[, c(grp, have)])
  colnames(df)[1] <- "grp"; df$grp <- factor(df$grp)
  ylabs <- c(nFeature_RNA = "Genes", nCount_RNA = "UMIs", percent.mt = "Mito %")
  plist <- lapply(have, function(m) {
    ggplot(df, aes(grp, .data[[m]])) +
      geom_violin(fill = "#4DBBD5", alpha = 0.6, scale = "width") +
      labs(x = NULL, y = ylabs[[m]]) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  })
  if (length(plist) == 1) plist[[1]] else Reduce(`+`, plist)
}
umap_p <- function(seu, grp, label = FALSE) {
  DimPlot(seu, reduction = "umap", group.by = grp, label = label, label.size = 3) +
    theme(legend.text = element_text(size = 8))
}
dot_p <- function(seu, markers, grp) {
  ft <- intersect(unique(unlist(markers)), rownames(seu))
  if (length(ft) < 6 || is.null(grp)) return(NULL)
  v <- as.character(seu[[grp]][[1]])
  keep <- !is.na(v)
  seu <- seu[, keep]
  v <- factor(v[keep])
  seu[[grp]] <- v
  if (length(levels(v)) < 2) return(NULL)
  Idents(seu) <- grp
  ## try directly first; large objects often raise the "$<- id NA/0row" error -> then downsample and retry
  p <- tryCatch(DotPlot(seu, features = ft, cols = c("grey85", "#B2182B")) +
                  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)),
                error = function(e) NULL)
  if (is.null(p) && ncol(seu) > 1e5) {
    cat("    (DotPlot failed on the large object, sampling <=30000 per cluster and retrying)\n")
    set.seed(2026)
    cells <- unlist(lapply(levels(v), function(lv) {
      ids <- WhichCells(seu, idents = lv)
      if (length(ids) > 30000) sample(ids, 30000) else ids
    }))
    seu2 <- subset(seu, cells = cells)
    Idents(seu2) <- grp
    p <- tryCatch(DotPlot(seu2, features = ft, cols = c("grey85", "#B2182B")) +
                    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)),
                  error = function(e) NULL)
    if (!is.null(p)) attr(p, "downsampled") <- TRUE
  }
  p
}
MAJOR_MARKERS <- list(
  Monocyte = c("CD14","LYZ","S100A8","FCN1"),
  Macrophage = c("CD68","C1QA","C1QB","CSF1R"),
  DC = c("ITGAX","CD1C","CLEC9A"),
  pDC = c("IL3RA","CLEC4C"),
  T_cell = c("CD3D","CD3E","IL7R"),
  NK = c("NKG7","GNLY","KLRD1"),
  B_cell = c("MS4A1","CD79A","CD19"),
  Plasma = c("MZB1","SDC1","IGHG1"),
  Neutrophil = c("FCGR3B","CSF3R","S100A8"),
  Mast = c("TPSAB1","TPSB2","CPA3"),
  RBC = c("HBB","HBA1","AHSP"),
  Platelet = c("PPBP","PF4"),
  AT2 = c("SFTPC","SFTPB","NAPSA"),
  AT1 = c("AGER","CAV1","PDPN"),
  Endothelial = c("PECAM1","VWF","CLDN5"),
  Fibroblast = c("COL1A1","DCN","LUM"),
  Epithelial = c("EPCAM","KRT18","KRT19"),
  SmoothMuscle = c("ACTA2","MYH11","TAGLN")
)

cat("===== 62 v2: blood/lung QC/annotation supplementary figure =====\n")
summ <- character()

## ---- Blood ----
blood_file <- file.path(out_dir, "cp3_annotated.rds")
if (!file.exists(blood_file)) blood_file <- file.path(out_dir, "path1_sepsis_copd_integrated.rds")
bl <- NULL
if (file.exists(blood_file)) {
  cat("Loading blood object:", blood_file, "\n")
  bl <- readRDS(blood_file)
  ct <- pick(bl, c("cell_type","celltype","annotation"))
  gp <- pick(bl, c("condition","group","clinical_group"))
  ds <- pick(bl, c("dataset","Dataset","orig.ident"))
  cat("Blood fields: cell_type=", ct, " group=", gp, " dataset=", ds,
      "| cells=", ncol(bl), "\n", sep = "")
  safe("blood umap check", {
    reds <- Reductions(bl)
    if (!("umap" %in% reds)) { um <- grep("umap", reds, ignore.case = TRUE, value = TRUE)
      if (length(um) > 0) bl[["umap"]] <- bl[[um[1]]] }
  })
  if (!is.null(ct)) {
    safe("blood celltype umap", save2(umap_p(bl, ct, label = TRUE),
                                      "path62_blood_celltype_umap", 9, 7))
    safe("blood dotplot", save2(dot_p(bl, MAJOR_MARKERS, ct), "path62_blood_dotplot", 9, 7))
  }
  if (!is.null(gp)) safe("blood group umap", save2(umap_p(bl, gp), "path62_blood_group_umap"))
  if (!is.null(ds)) safe("blood dataset umap", save2(umap_p(bl, ds), "path62_blood_dataset_umap"))
  safe("blood qc", save2(qc_panel(bl, if (!is.null(ds)) ds else gp), "path62_blood_qc", 10, 5))
  summ <- c(summ, sprintf("BLOOD: cells=%d", ncol(bl)))
  if (!is.null(ds)) summ <- c(summ, sprintf("  per dataset: %s",
        paste(sprintf("%s=%d", names(table(bl[[ds]])), table(bl[[ds]])), collapse = "; ")))
} else cat("!! blood object not found\n")

## ---- Lung ----
lung_file <- file.path(out_dir, "path2_copd_lung.rds")
lu <- NULL
if (file.exists(lung_file)) {
  cat("Loading lung object:", lung_file, "\n")
  lu <- readRDS(lung_file)
  ct <- pick(lu, c("cell_type","Manuscript_Identity","celltype","annotation"))
  dis <- pick(lu, c("disease","Disease_Identity","group"))
  ## sample = 46 donors (GSE136831 donor ID), highest priority; orig.ident=62 is a different identifier
  dn <- pick_nonconst(lu, c("sample","donor","donor_id","patient","subject","geo_accession","orig.ident"))
  cat("Lung fields: cell_type=", ct, " disease=", dis, " donor=", dn,
      "| cells=", ncol(lu), "\n", sep = "")
  if (!is.null(dn)) cat("  lung donor distinct:", length(unique(lu[[dn]])), "\n")
  lu <- ensure_umap(lu, "lung")
  if (!is.null(ct)) {
    if ("umap" %in% Reductions(lu)) {
      safe("lung celltype umap", save2(umap_p(lu, ct, label = FALSE), "path62_lung_celltype_umap", 9, 7))
    } else cat("  [skip] lung celltype umap (no umap)\n")
    safe("lung dotplot", save2(dot_p(lu, MAJOR_MARKERS, ct), "path62_lung_dotplot", 10, 8))
  }
  cc <- pick(lu, c("cell_category"))
  if (!is.null(cc) && "umap" %in% Reductions(lu)) {
    safe("lung cell_category umap", save2(umap_p(lu, cc, label = TRUE), "path62_lung_cell_category_umap", 9, 7))
  }
  if (!is.null(dis)) {
    if ("umap" %in% Reductions(lu)) safe("lung disease umap", save2(umap_p(lu, dis), "path62_lung_disease_umap"))
    else cat("  [skip] lung disease umap (no umap)\n")
  }
  safe("lung qc", save2(qc_panel(lu, if (!is.null(dis)) dis else dn), "path62_lung_qc", 8, 5))
  summ <- c(summ, sprintf("LUNG: cells=%d", ncol(lu)))
  if (!is.null(dis)) summ <- c(summ, sprintf("  per disease: %s",
        paste(sprintf("%s=%d", names(table(lu[[dis]])), table(lu[[dis]])), collapse = "; ")))
  if (!is.null(dn)) summ <- c(summ, sprintf("  donors=%d", length(unique(lu[[dn]]))))
} else cat("!! lung object not found\n")

writeLines(summ, file.path(out_dir, "path62_meta_summary.txt"))
cat(paste(summ, collapse = "\n"), "\n")
cat("===== script 62 v2 done =====\n")
