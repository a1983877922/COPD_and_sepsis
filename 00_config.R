###############################################################################
# Sepsis + COPD comorbidity analysis — shared configuration & utility functions
# File: 00_config.R
#
# Purpose: shared by 01_assemble_blood_integrated.R and 02_build_lung_cross_tissue.R.
#
# Data sources:
#   Sepsis (Sepsis)   — reuses the integrated objects from the 20260726 task:
#     merged_seu_qc.rds    : raw counts after QC (unnormalized, for re-integration in path 1)
#     sepsis_integrated.rds: integrated + annotated cell_type (for path 2 cross-tissue comparison)
#   COPD (COPD)       — under the "慢阻肺COPD/" directory:
#     GSE249584 : peripheral blood PBMC, 10x 3' v3.1 (7 HC + 8 COPD), for path 1 blood-blood merging
#     GSE136831 : lung tissue (Adams 2020 IPF/COPD/Control atlas), for path 2 cross-tissue
#
# Analysis design (two paths):
#   Path 1 (blood-blood direct merge): sepsis PBMC + COPD PBMC → Harmony integration → four-group comorbidity analysis
#   Path 2 (lung-blood cross-tissue)  : sepsis blood PBMC vs COPD lung tissue → independent annotation + cross-tissue comparison
#
# Follows 20260726/integration_script.R's:
#   - Seurat v5 + Harmony integration strategy
#   - MAD (median absolute deviation, nmads=3) QC
#   - configurable parameter section + classic marker gene annotation + CellChat
###############################################################################

## =========================================================================
## 0. Path configuration
## =========================================================================

# =====================================================================
# ★★★ On the server, only this line needs changing (raw data root dir RAW_DIR) ★★★
#   Server:   "/media/desk16/ysx5991/sepsis_copd/00rawdata"
#   Local test: "D:/biosoft/others/脓毒症"
#   Note: RAW_DIR must keep the following two subdirectories:
#         RAW_DIR/脓毒症sepsis/    -> 5 sepsis datasets
#         RAW_DIR/慢阻肺copd/      -> COPD data (GSE249584 / GSE136831)
RAW_DIR <- "/media/desk16/ysx5991/sepsis_copd/00rawdata"
# =====================================================================
sepsis_dir <- file.path(RAW_DIR, "脓毒症sepsis")
copd_dir   <- file.path(RAW_DIR, "慢阻肺copd")

# Output directory = the directory containing this script (located automatically
# by scripts 01~05 and passed in via SCRIPT_DIR)
if (exists("SCRIPT_DIR") && nzchar(SCRIPT_DIR) && dir.exists(SCRIPT_DIR)) {
  out_dir <- SCRIPT_DIR
} else {
  out_dir <- getwd()
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## =========================================================================
## 0.1 Dependencies (load order: plyr must be loaded before dplyr, otherwise
##     it masks summarise/n/arrange)
## =========================================================================

library(Seurat)
library(harmony)
library(Matrix)
library(plyr)
library(dplyr)
library(stringr)
library(data.table)
library(reshape2)
library(ggplot2)
library(patchwork)

# future parallel memory ceiling: with 1.17M cells, the globals passed by
# FindNeighbors/FindMarkers (FUN/query/index) reach ~700MiB, which exceeds
# future's default 500MiB limit and errors out.
# Server has 512GB RAM, so 8GB is plenty. To enable multi-core parallel
# acceleration, uncomment the next line:
# library(future); plan(multisession, workers = 16)
options(future.globals.maxSize = 8 * 1024^3)   # 8GB

## =========================================================================
## 0.2 Configurable parameters
## =========================================================================

# Harmony parameters
harmony_theta  <- 2      # degree of batch mixing
harmony_lambda <- 1      # degree of biological structure preservation

# Normalization parameters
varGeneNum   <- 5000   # number of highly variable genes (plenty of 512GB RAM, raised from 3000 to 5000 for finer integration)
computePCs   <- 50
pcadim       <- 1:30
vars_to_regress <- c("percent.mt")

# Clustering parameters
resolution <- 0.5

# UMAP parameters
umap_nneighbors <- 30

# Marker Gene parameters
markerMethod  <- "wilcox"
markerThresh  <- 0.25
markerMinPct  <- 0.1
markerOnlyPos <- TRUE
markerBase    <- 2

# Downsampling parameters
# Server with 512GB RAM: no limit at all (Inf), keeping all ~1.17M cells
# (peak memory ~100-120GB, well within 512GB).
max_cells_per_sample <- Inf   # global per-sample downsampling cap (Inf = no limit)
gse279452_max_cells  <- Inf   # separate per-sample cap for GSE279452 (Inf = no limit)

# Whether to include GSE279452 (133 samples, ~890k cells, the bulk of the data)
#   With 512GB server RAM it can be fully included; set to TRUE
include_gse279452 <- TRUE

# CellChat per-group downsampling cap
#   computeCommunProb is a permutation-based computation with complexity
#   roughly O(n_cells x 100 permutations x number of ligand-receptor pairs).
#   Without downsampling, the Sepsis group has ~660k cells and would run for
#   days; it must be reduced to tens of thousands.
#   CellChat officially recommends downsampling large datasets to 20-50k;
#   30000 is used here (communication network built independently per group, robust enough).
cellchat_max_cells <- 30000

## =========================================================================
## 0.3 MAD QC functions (same as 20260726/integration_script.R)
## =========================================================================
## For each dataset separately, compute the median and MAD of QC metrics;
## cells deviating from the median by more than nmads MADs are treated as
## outliers and removed.

qc_mad <- function(vals, nmads = 3, lower_only = FALSE, upper_only = FALSE,
                   hard_lower = -Inf, hard_upper = Inf,
                   log_transform = FALSE, label = "metric") {
  vals_q <- if (log_transform) log1p(vals) else vals
  med    <- median(vals_q, na.rm = TRUE)
  mad_v  <- mad(vals_q, na.rm = TRUE)
  lower  <- med - nmads * mad_v
  upper  <- med + nmads * mad_v
  if (log_transform) { lower <- expm1(lower); upper <- expm1(upper) }
  if (mad_v == 0)    { lower <- hard_lower; upper <- hard_upper }
  lower <- max(lower, hard_lower)
  upper <- min(upper, hard_upper)
  keep <- vals >= lower & vals <= upper
  if (lower_only) keep <- vals <= upper
  if (upper_only) keep <- vals >= lower
  med_orig <- if (log_transform) median(vals, na.rm = TRUE) else med
  cat(sprintf("  %-14s: median=%.1f, keep %.1f%%  (range: %.1f ~ %.1f)\n",
              label, med_orig, 100 * mean(keep, na.rm = TRUE), lower, upper))
  return(keep)
}

# Run MAD QC on a single Seurat object (thresholds computed per dataset), returns the filtered object
mad_qc <- function(seu, nmads = 3) {
  seu[["percent.mt"]]  <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]\\d")
  hb_genes <- grep("^HB[^PS]", rownames(seu), value = TRUE)
  seu[["percent.hb"]] <- if (length(hb_genes) > 0) {
    PercentageFeatureSet(seu, features = hb_genes)
  } else 0

  qc_results <- list()
  for (ds in unique(seu$dataset)) {
    cells_idx <- which(seu$dataset == ds)
    keep_umi  <- qc_mad(seu$nCount_RNA[cells_idx],  nmads = nmads,
                        log_transform = TRUE, hard_lower = 500, label = "nCount_RNA")
    keep_gene <- qc_mad(seu$nFeature_RNA[cells_idx], nmads = nmads,
                        log_transform = TRUE, hard_lower = 100, label = "nFeature_RNA")
    keep_mt   <- qc_mad(seu$percent.mt[cells_idx],   nmads = nmads,
                        lower_only = TRUE, hard_upper = 30, label = "percent.mt")
    keep_all  <- keep_umi & keep_gene & keep_mt
    qc_results[[ds]] <- cells_idx[keep_all]
    cat(sprintf("  [%s] kept %d / %d cells (%.1f%%)\n",
                ds, sum(keep_all), length(keep_all), 100 * mean(keep_all)))
  }
  seu[, unlist(qc_results)]
}

## =========================================================================
## 0.4 COPD GSE249584 condition mapping (peripheral blood PBMC)
## =========================================================================
## Filenames contain only patient IDs (e.g. 20_00264_LI_SING), no COPD/Control labels.
## Labels were verified against GEO GSM titles (eutils esummary): 15 samples = 7 Control + 8 COPD.
## Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249584
copd_condition_map <- c(
  # Healthy controls (7 cases)
  "20_00262" = "Control",
  "20_00263" = "Control",
  "20_00264" = "Control",
  "20_00266" = "Control",
  "20_00267" = "Control",
  "20_00193" = "Control",
  "20_00198" = "Control",
  # COPD (8 cases)
  "20_00265" = "COPD",
  "20_00269" = "COPD",
  "20_00270" = "COPD",
  "20_00271" = "COPD",
  "20_00194" = "COPD",
  "20_00195" = "COPD",
  "20_00196" = "COPD",
  "20_00197" = "COPD"
)

## =========================================================================
## 0.5 Read GSE249584 (COPD peripheral blood PBMC, standard 10x mtx)
## =========================================================================
## Filenames: GSM7950677_20_00269_LI_SING_matrix.mtx.gz
##           GSM7950677_20_00269_LI_SING_barcodes.tsv.gz
##           GSM7950677_20_00269_LI_SING_features.tsv.gz
## Same format as the sepsis 10x datasets; reuses the read_10x_from_gsm logic.

read_copd_10x <- function(gse_dir, gse_id = "GSE249584",
                          condition_map = copd_condition_map,
                          max_cells = NULL) {
  files <- list.files(gse_dir, pattern = "_matrix\\.mtx\\.gz$")
  seu_list <- list()
  for (f in files) {
    gsm_prefix <- sub("_matrix\\.mtx\\.gz$", "", f)
    sample_dir <- file.path(gse_dir, gsm_prefix)
    dir.create(sample_dir, showWarnings = FALSE)

    file.copy(file.path(gse_dir, paste0(gsm_prefix, "_barcodes.tsv.gz")),
              file.path(sample_dir, "barcodes.tsv.gz"), overwrite = TRUE)
    file.copy(file.path(gse_dir, paste0(gsm_prefix, "_features.tsv.gz")),
              file.path(sample_dir, "features.tsv.gz"), overwrite = TRUE)
    file.copy(file.path(gse_dir, f),
              file.path(sample_dir, "matrix.mtx.gz"), overwrite = TRUE)

    # Sample name: "20_00269_LI_SING" -> patient ID "20_00269"
    sample_name <- sub("^GSM\\d+_", "", gsm_prefix)
    patient_id  <- sub("_LI_SING$", "", sample_name)

    mat <- Read10X(sample_dir)
    if (is.list(mat)) mat <- mat[["Gene Expression"]]

    if (!is.null(max_cells) && ncol(mat) > max_cells) {
      set.seed(42)
      mat <- mat[, sample(ncol(mat), max_cells)]
    }

    seu_obj <- CreateSeuratObject(mat, project = gse_id,
                                  min.cells = 3, min.features = 200)
    rm(mat); gc(verbose = FALSE)

    seu_obj$sample    <- sample_name
    seu_obj$patient   <- patient_id
    seu_obj$dataset   <- gse_id
    seu_obj$condition <- unname(condition_map[patient_id])

    seu_list[[gsm_prefix]] <- seu_obj
    cat(gse_id, sample_name, ":", ncol(seu_obj), "cells, condition =",
        seu_obj$condition[1], "\n")
  }
  seu_list
}

## =========================================================================
## 0.6 Read GSE136831 (COPD lung tissue, sparse matrix + metadata)
## =========================================================================
## Files:
##   GSE136831_RawCounts_Sparse.mtx.gz                 : sparse matrix (genes x cells)
##   GSE136831_AllCells.cellBarcodes.txt.gz            : cell barcodes (Subject_barcode)
##   GSE136831_AllCells.GeneIDs.txt.gz                 : 2 columns (Ensembl_GeneID, HGNC symbol)
##   GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz : metadata
## Key metadata columns:
##   CellBarcode_Identity   = "001C_AAACCTGCATCGGGTC" (matches barcode)
##   Disease_Identity       = Control / COPD / IPF
##   Manuscript_Identity    = pre-annotated cell type (ncMonocyte, Macrophage_Alveolar, ...)
##   Subject_Identity       = "001C" (patient/donor)
## Note: this atlas contains IPF samples; the comorbidity study keeps only
## Control + COPD (IPF excluded).

read_gse136831 <- function(gse_dir, gse_id = "GSE136831",
                           keep_disease = c("Control", "COPD")) {
  mtx_file  <- file.path(gse_dir, "GSE136831_RawCounts_Sparse.mtx.gz")
  bc_file   <- file.path(gse_dir, "GSE136831_AllCells.cellBarcodes.txt.gz")
  gene_file <- file.path(gse_dir, "GSE136831_AllCells.GeneIDs.txt.gz")
  meta_file <- file.path(gse_dir, "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz")

  cat("Reading GSE136831 sparse matrix (2.15GB, may take a few minutes)...\n")
  mat <- readMM(gzfile(mtx_file))

  barcodes <- readLines(gzfile(bc_file))
  genes    <- read.delim(gzfile(gene_file), header = TRUE, stringsAsFactors = FALSE)
  # Column 2 of GeneIDs is the HGNC gene symbol; fall back to Ensembl ID if missing
  gene_sym <- genes[[2]]
  gene_sym[is.na(gene_sym) | gene_sym == ""] <- genes[[1]][is.na(gene_sym) | gene_sym == ""]
  gene_sym <- make.unique(gene_sym)

  # Dimension sanity check
  stopifnot(nrow(mat) == length(gene_sym),
            ncol(mat) == length(barcodes))
  rownames(mat) <- gene_sym
  colnames(mat) <- barcodes

  # Metadata (column 1 is the barcode, used as row names)
  meta <- read.delim(gzfile(meta_file), header = TRUE,
                     stringsAsFactors = FALSE, check.names = FALSE)
  rownames(meta) <- meta$CellBarcode_Identity

  # Keep only Control + COPD (IPF excluded)
  keep_cells <- rownames(meta)[meta$Disease_Identity %in% keep_disease]
  meta <- meta[keep_cells, , drop = FALSE]
  mat  <- mat[, keep_cells, drop = FALSE]

  cat(sprintf("GSE136831 after filtering: %d cells (%s)\n", ncol(mat),
              paste(names(table(meta$Disease_Identity)), collapse = ", ")))

  seu <- CreateSeuratObject(mat, project = gse_id,
                            min.cells = 3, min.features = 200)
  rm(mat); gc(verbose = FALSE)

  # Built-in annotation + disease labels written directly
  seu$dataset       <- gse_id
  seu$sample        <- meta[colnames(seu), "Subject_Identity"]
  seu$disease       <- meta[colnames(seu), "Disease_Identity"]
  seu$condition     <- ifelse(seu$disease == "COPD", "COPD", "Healthy")
  seu$cell_type     <- meta[colnames(seu), "Manuscript_Identity"]
  seu$cell_category <- meta[colnames(seu), "CellType_Category"]

  cat("GSE136831 cell types:", paste(names(table(seu$cell_type)), collapse = ", "), "\n")
  seu
}

## =========================================================================
## 0.6b Read raw sepsis data (5 core PBMC datasets)
## =========================================================================
## Re-read from raw data (does not depend on 20260726/merged_seu_qc.rds):
##   - GSE167363 / GSE175453 / GSE279452 : standard 10x mtx
##   - GSE151263                          : processed UMI txt
##   - SCP548                             : CSV expression matrix + metadata
## Reading logic follows 20260726/integration_script.R's read_10x_from_gsm /
## read_umi_txt / read_scp548.

# Standard 10x mtx reading (GSE167363 / GSE175453 / GSE279452)
read_10x_from_gsm <- function(gse_id, condition_map = NULL, max_cells = NULL) {
  gse_dir <- file.path(sepsis_dir, gse_id)
  if (!dir.exists(gse_dir)) stop("Directory not found: ", gse_dir, ". Please check RAW_DIR and the dataset directory name (case-sensitive).")
  files <- list.files(gse_dir, pattern = "[._]matrix\\.mtx\\.gz$")
  if (length(files) == 0) stop("No *matrix.mtx.gz files found in ", gse_dir, ".")
  seu_list <- list()
  for (f in files) {
    gsm_prefix <- sub("[._]matrix\\.mtx\\.gz$", "", f)
    sample_dir <- file.path(gse_dir, gsm_prefix)
    dir.create(sample_dir, showWarnings = FALSE)

    bc_src <- file.path(gse_dir, paste0(gsm_prefix, "_barcodes.tsv.gz"))
    if (!file.exists(bc_src)) bc_src <- file.path(gse_dir, paste0(gsm_prefix, ".barcodes.tsv.gz"))
    file.copy(bc_src, file.path(sample_dir, "barcodes.tsv.gz"), overwrite = TRUE)
    ft_src <- file.path(gse_dir, paste0(gsm_prefix, "_features.tsv.gz"))
    if (!file.exists(ft_src)) ft_src <- file.path(gse_dir, paste0(gsm_prefix, ".features.tsv.gz"))
    file.copy(ft_src, file.path(sample_dir, "features.tsv.gz"), overwrite = TRUE)
    file.copy(file.path(gse_dir, f), file.path(sample_dir, "matrix.mtx.gz"), overwrite = TRUE)

    sample_name <- sub("^GSM\\d+_", "", gsm_prefix)
    mat <- Read10X(sample_dir)
    if (is.list(mat)) mat <- mat[["Gene Expression"]]   # keep gene expression only, skip ADT

    if (!is.null(max_cells) && ncol(mat) > max_cells) {
      set.seed(42)
      mat <- mat[, sample(ncol(mat), max_cells)]
    }
    seu_obj <- CreateSeuratObject(mat, project = gse_id,
                                  min.cells = 3, min.features = 200)
    rm(mat); gc(verbose = FALSE)

    seu_obj$sample <- sample_name
    seu_obj$dataset <- gse_id
    if (!is.null(condition_map)) {
      seu_obj$condition <- condition_map(sample_name)
    } else {
      seu_obj$condition <- ifelse(grepl("HC|Healthy|H\\d", sample_name, ignore.case = TRUE),
                                  "Healthy", "Sepsis")
    }
    seu_list[[gsm_prefix]] <- seu_obj
    cat(gse_id, sample_name, ":", ncol(seu_obj), "cells\n")
  }
  seu_list
}

# GSE151263: processed UMI txt (gene x cell text matrix)
read_umi_txt_151263 <- function() {
  gse_dir <- file.path(sepsis_dir, "GSE151263")
  if (!dir.exists(gse_dir)) stop("Directory not found: ", gse_dir, ". Please check RAW_DIR and the dataset directory name (case-sensitive).")
  files <- list.files(gse_dir, pattern = "_processed_UMI\\.txt\\.gz$")
  if (length(files) == 0) stop("No *_processed_UMI.txt.gz files found in ", gse_dir, ".")
  seu_list <- list()
  for (f in files) {
    gsm_prefix <- sub("_processed_UMI\\.txt\\.gz$", "", f)
    sample_name <- sub("^GSM\\d+_", "", gsm_prefix)

    mat <- fread(file.path(gse_dir, f), header = TRUE, data.table = FALSE)
    rn <- make.unique(mat[, 1])       # deduplicate gene names (e.g. 1-Mar, 2-Mar)
    mat <- mat[, -1]
    mat <- as(as.matrix(mat), "dgCMatrix")
    rownames(mat) <- rn

    seu_obj <- CreateSeuratObject(mat, project = "GSE151263",
                                  min.cells = 3, min.features = 200)
    rm(mat); gc(verbose = FALSE)
    seu_obj$sample <- sample_name
    seu_obj$dataset <- "GSE151263"
    seu_obj$condition <- ifelse(grepl("ARDS", sample_name), "ARDS_and_Sepsis", "Sepsis")
    seu_list[[gsm_prefix]] <- seu_obj
    cat("GSE151263", sample_name, ":", ncol(seu_obj), "cells\n")
  }
  seu_list
}

# SCP548: CSV expression matrix + metadata (Reyes et al. 2020)
read_scp548 <- function() {
  scp_dir <- file.path(sepsis_dir, "SCP548")
  mtx_file <- file.path(scp_dir, "expression", "scp_gex_matrix_raw.csv.gz")
  if (!file.exists(mtx_file)) stop("SCP548 expression matrix not found: ", mtx_file, ". Please check that the file was uploaded completely.")
  cat("Reading SCP548 expression matrix (chunked reading)...\n")

  header_dt <- fread(mtx_file, nrow = 0, data.table = FALSE)
  cell_names <- colnames(header_dt)[-1]
  n_cells <- length(cell_names)
  gene_col <- fread(mtx_file, select = 1, header = TRUE, data.table = FALSE)
  gene_names <- make.unique(gene_col[, 1])
  cat(sprintf("  Dimensions: %d genes x %d cells\n", length(gene_names), n_cells))

  chunk_size <- 10000
  n_chunks <- ceiling(n_cells / chunk_size)
  mat_list <- vector("list", n_chunks)
  for (i in seq_len(n_chunks)) {
    start <- (i - 1) * chunk_size + 1
    end <- min(i * chunk_size, n_cells)
    col_idx <- c(1, (start + 1):(end + 1))
    chunk <- fread(mtx_file, select = col_idx, header = TRUE, data.table = FALSE)
    chunk_mat <- as(as.matrix(chunk[, -1]), "dgCMatrix")
    colnames(chunk_mat) <- cell_names[start:end]
    mat_list[[i]] <- chunk_mat
    cat(sprintf("  Chunk %d/%d: cells %d-%d\n", i, n_chunks, start, end))
    rm(chunk); gc()
  }
  mat <- do.call(cbind, mat_list)
  rm(mat_list); gc()
  rownames(mat) <- gene_names

  meta_lines <- readLines(file.path(scp_dir, "metadata", "scp_meta_updated.txt"))
  meta_header <- strsplit(meta_lines[1], "\t")[[1]]
  meta <- read.delim(textConnection(meta_lines[-c(1, 2)]), header = FALSE,
                     stringsAsFactors = FALSE, sep = "\t")
  colnames(meta) <- trimws(meta_header)

  seu <- CreateSeuratObject(mat, project = "SCP548",
                            min.cells = 3, min.features = 200)
  seu$dataset <- "SCP548"
  meta_idx <- match(colnames(seu), meta$NAME)
  seu$sample <- unname(meta$donor_id[meta_idx])
  # SCP548 ships author annotations (Cell_Type, values like Mono/DC/T/NK/B/Megakaryocyte
  # abbreviations). Save to the separate scp548_cell_type column for reference only; do NOT
  # occupy the cell_type column — otherwise after merge the other datasets would have
  # cell_type = NA, and the naming mismatch with blood_markers would make downstream
  # "already annotated?" checks wrongly assume annotation is done and skip annotate_by_markers.
  seu$scp548_cell_type <- unname(meta$Cell_Type[meta_idx])
  cohort_map <- function(cohort) {
    if (is.na(cohort)) return(NA)
    if (cohort == "Control") return("Healthy")
    if (cohort %in% c("ICU-NoSEP", "Leuk-UTI")) return("Infection_Control")
    return("Sepsis")
  }
  seu$condition <- unname(sapply(meta$Cohort[meta_idx], cohort_map))
  rm(mat, gene_col); gc()
  cat("SCP548:", ncol(seu), "cells\n")
  list(SCP548 = seu)
}

# Read the 5 sepsis PBMC datasets -> per-sample downsampling -> merge -> MAD QC -> adult subset
read_sepsis_pbmc <- function(max_cells = max_cells_per_sample) {
  cat("\n=== Reading 5 sepsis PBMC datasets (raw data) ===\n")

  seu_167363 <- read_10x_from_gsm("GSE167363", function(s) {
    if (grepl("^HC", s)) "Healthy" else "Sepsis"
  })
  seu_175453 <- read_10x_from_gsm("GSE175453", function(s) {
    if (grepl("^HC|^H\\d", s)) "Healthy" else "Sepsis"
  })
  seu_151263 <- read_umi_txt_151263()
  seu_scp548 <- read_scp548()

  all_seu <- c(seu_167363, seu_175453, seu_151263, seu_scp548)
  rm(seu_167363, seu_175453, seu_151263, seu_scp548)

  # GSE279452 has 133 samples (12GB); uses the separate gse279452_max_cells cap.
  # When include_gse279452=FALSE it is skipped (local testing); set TRUE on the
  # server to include the full data.
  if (include_gse279452) {
    cat("  Reading GSE279452 (133 samples, separate cap ", gse279452_max_cells, ")...\n")
    seu_279452 <- read_10x_from_gsm("GSE279452", function(s) {
      if (grepl("^Abd-PS|^Res-PS", s)) "Pediatric_Sepsis"
      else if (grepl("^PHC", s)) "Pediatric_Healthy"
      else if (grepl("^PIC", s)) "Pediatric_Infection_Control"
      else if (grepl("^HC", s)) "Healthy"
      else if (grepl("^IC", s)) "Infection_Control"
      else "Sepsis"
    }, max_cells = gse279452_max_cells)
    all_seu <- c(all_seu, seu_279452)
    rm(seu_279452)
  } else {
    cat("  [Skipping GSE279452] include_gse279452 = FALSE (local testing)\n")
  }
  gc()

  # Per-sample downsampling (max_cells_per_sample)
  cat("\n--- Per-sample downsampling (max = ", max_cells, ") ---\n")
  set.seed(42)
  for (i in seq_along(all_seu)) {
    obj <- all_seu[[i]]
    meta_df <- obj@meta.data
    cells_keep <- unlist(lapply(unique(meta_df$sample), function(s) {
      cells_s <- rownames(meta_df)[meta_df$sample == s]
      if (length(cells_s) > max_cells) sample(cells_s, max_cells) else cells_s
    }))
    rm(meta_df)
    if (length(cells_keep) < ncol(obj)) {
      counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
      counts_mat <- counts_mat[, cells_keep, drop = FALSE]
      meta_new <- obj@meta.data[cells_keep, , drop = FALSE]
      rm(obj); gc()
      all_seu[[i]] <- CreateSeuratObject(counts = counts_mat,
                                         meta.data = meta_new,
                                         min.cells = 0, min.features = 0,
                                         project = names(all_seu)[i])
      rm(counts_mat, meta_new); gc()
    }
    cat(sprintf("  %s: %d cells\n", names(all_seu)[i], ncol(all_seu[[i]])))
  }

  # Merge
  cat("\n=== Merging 5 datasets ===\n")
  total_cells <- sum(sapply(all_seu, ncol))
  cat("Total cells before merging:", total_cells, "\n")
  if (total_cells > 600000) {
    cat("[Hint] Total cell count is large; if memory is insufficient, lower gse279452_max_cells or max_cells_per_sample further\n")
  }
  gc()
  # Batch merging: all_seu holds ~162 per-sample Seurat objects; a single merge would
  # hold all 162 objects plus the new object simultaneously, doubling peak memory and
  # causing std::bad_alloc.
  # Instead, merge 20 objects per batch and immediately JoinLayers into a single layer,
  # greatly reducing peak memory.
  batch_size <- 20
  n_obj <- length(all_seu)
  n_batches <- ceiling(n_obj / batch_size)
  cat("  Batch-merging", n_obj, "sample objects (", batch_size, "per batch)...\n")
  batch_list <- list()
  for (bi in seq_len(n_batches)) {
    idx <- ((bi - 1) * batch_size + 1):min(bi * batch_size, n_obj)
    b <- all_seu[idx]
    if (length(b) == 1) {
      m <- b[[1]]
    } else {
      m <- merge(b[[1]], y = b[-1], add.cell.ids = names(b),
                 project = "Sepsis_PBMC")
    }
    m <- JoinLayers(m, assay = "RNA")
    batch_list[[bi]] <- m
    rm(b, m); gc()
    cat("    Batch", bi, "/", n_batches, ":", ncol(batch_list[[bi]]), "cells\n")
  }
  rm(all_seu); gc()
  # Merge batches (each batch is already a single counts layer, low memory overhead)
  if (length(batch_list) == 1) {
    merged <- batch_list[[1]]
  } else {
    merged <- merge(batch_list[[1]], y = batch_list[-1], project = "Sepsis_PBMC")
  }
  merged <- JoinLayers(merged, assay = "RNA")
  rm(batch_list); gc()
  cat("After merging:", ncol(merged), "cells\n")
  print(table(merged$dataset))

  # MAD QC
  cat("\n=== MAD QC ===\n")
  merged <- mad_qc(merged, nmads = 3)
  cat("After QC:", ncol(merged), "cells\n")

  # Keep only the adult subset (GSE279452 pediatric samples excluded)
  merged <- subset(merged, subset = !(condition %in%
    c("Pediatric_Sepsis", "Pediatric_Healthy", "Pediatric_Infection_Control")))
  cat("After adult subset:", ncol(merged), "cells\n")
  print(table(merged$dataset))
  print(table(merged$condition))
  merged
}

## =========================================================================
## 0.7 Cell type marker genes (peripheral blood / whole blood)
## =========================================================================

blood_markers <- list(
  "T_cell"     = c("CD3D", "CD3E", "CD4", "CD8A", "IL7R"),
  "B_cell"     = c("MS4A1", "CD79A", "CD79B", "IGHM", "IGHD"),
  "NK"         = c("GNLY", "NKG7", "KLRF1", "KLRD1"),
  "Monocyte"   = c("CD14", "CD68", "LYZ", "S100A8", "S100A9"),
  "DC"         = c("FCER1A", "CST3", "CLEC9A", "ITGAX"),
  "Neutrophil" = c("FCGR3B", "CSF3R", "CXCR2"),
  "Platelet"   = c("PPBP", "PF4"),
  "RBC"        = c("HBB", "HBA1", "ALAS2"),
  "Mast"       = c("CPA3", "MS4A2", "TPSAB1")
)

## =========================================================================
## 0.8 Automatic annotation function (classic marker gene ModuleScore, same as 20260726 Step 8)
## =========================================================================

annotate_by_markers <- function(seu, marker_list = blood_markers,
                                cluster_col = "seurat_clusters",
                                col_name = "cell_type") {
  # Match genes against rownames(seu) (after Seurat v5 merge the layers are split;
  # GetAssayData(layer="data") errors with "multiple layers"; taking the union via
  # rownames is sufficient)
  avail_genes <- rownames(seu)
  scored_ct <- character(0)   # record cell types successfully scored
  for (ct in names(marker_list)) {
    genes <- intersect(marker_list[[ct]], avail_genes)
    if (length(genes) < 2) {
      cat(sprintf("  [Skipping %s] too few usable marker genes (%d)\n", ct, length(genes)))
      next
    }
    ok <- tryCatch({
      # AddModuleScore writes scores to the "<ct>1" column of meta.data (auto numeric suffix)
      seu <- AddModuleScore(seu, features = list(genes), name = ct)
      # Copy to the score_<ct> column for downstream use (cannot index with [,1], that returns a Seurat object)
      seu[[paste0("score_", ct)]] <- seu@meta.data[[paste0(ct, "1")]]
      TRUE
    }, error = function(e) {
      cat(sprintf("  [Skipping %s] AddModuleScore error: %s\n", ct, conditionMessage(e)))
      FALSE
    })
    if (ok) scored_ct <- c(scored_ct, ct)
  }
  # Assign clusters only for cell types that were successfully scored; skipped types
  # have missing scores and must not get 0 (otherwise which.max would wrongly pick
  # them in clusters that "match nothing")
  score_cols <- paste0("score_", scored_ct)
  score_mat <- sapply(score_cols, function(col) {
    if (col %in% colnames(seu@meta.data)) seu@meta.data[[col]] else NA_real_
  })
  cluster_scores <- as.data.frame(score_mat) %>%
    cbind(cluster = seu@meta.data[[cluster_col]]) %>%
    group_by(cluster) %>%
    summarise(across(starts_with("score_"), \(x) mean(x, na.rm = TRUE)),
              .groups = "drop")

  cluster_celltype <- as.data.frame(cluster_scores[, -1])
  cluster_assignment <- character(nrow(cluster_scores))
  for (i in seq_len(nrow(cluster_celltype))) {
    max_col <- which.max(cluster_celltype[i, ])
    cluster_assignment[i] <- sub("score_", "", colnames(cluster_celltype)[max_col])
    cat(sprintf("  Cluster %s -> %s\n", cluster_scores$cluster[i],
                cluster_assignment[i]))
  }
  seu[[col_name]] <- plyr::mapvalues(seu@meta.data[[cluster_col]],
                                     from = cluster_scores$cluster,
                                     to = cluster_assignment)
  seu
}

## =========================================================================
## 0.8a DEG wrapper + low-quality cluster identification
## =========================================================================
## Note: FindMarkers' wilcox test does not support latent.vars (passing it only
## warns and is ignored), so percent.mt regression correction must use MAST.
## Two strategies are provided:
##   use_mast_regression = TRUE and MAST installed → MAST + latent.vars="percent.mt" (rigorous but slow)
##   use_mast_regression = FALSE (default)         → wilcox + post-hoc removal of mitochondrial/ATP genes (pragmatic)
use_mast_regression <- FALSE   # local default FALSE (fast); on the server can set TRUE (requires MAST)

# Mitochondrial/ATP gene pattern (for filtering DEG results, removes false positives
# caused by percent.mt differences between groups)
mt_gene_pattern <- "^(MT-|ATP5|ATP6|COX[0-9]|NDUF|ND[1-6]|CYTB|UQCR|SDH[A-D]|MRPS|MRPL)"

# Long non-coding RNA / low-information gene pattern (for low-quality cluster identification)
lnc_gene_pattern <- "^(AC[0-9]|AL[0-9]|AP[0-9]|RP11|RP1[0-9]|CTD-|LINC|XX|GS[0-9]|HCG|RN7S|RNU|RNA5S|SCARN|SNORD|SNORA)"

run_deg <- function(seu, ident.1, ident.2, min.pct = 0.1, logfc.threshold = 0.25) {
  if (use_mast_regression && requireNamespace("MAST", quietly = TRUE)) {
    cat("  [DEG] MAST + percent.mt regression correction\n")
    res <- FindMarkers(seu, ident.1 = ident.1, ident.2 = ident.2,
                       only.pos = FALSE, min.pct = min.pct,
                       logfc.threshold = logfc.threshold,
                       test.use = "MAST", latent.vars = "percent.mt")
  } else {
    cat("  [DEG] wilcox (post-hoc removal of mitochondrial/ATP genes)\n")
    res <- FindMarkers(seu, ident.1 = ident.1, ident.2 = ident.2,
                       only.pos = FALSE, min.pct = min.pct,
                       logfc.threshold = logfc.threshold,
                       test.use = "wilcox")
    drop <- grepl(mt_gene_pattern, rownames(res))
    if (any(drop)) {
      cat("    Removed mitochondrial/ATP genes:", sum(drop), "\n")
      res <- res[!drop, , drop = FALSE]
    }
  }
  res
}

# Low-quality cluster identification: fraction of lncRNA/low-information genes among
# the cluster's top 10 markers >= threshold
identify_low_quality_clusters <- function(all_markers, top_n = 10,
                                          lnc_frac_thresh = 0.5) {
  low_q <- character(0)
  for (cl in unique(all_markers$cluster)) {
    sub <- all_markers[all_markers$cluster == cl, ]
    top <- head(sub[order(sub$avg_log2FC, decreasing = TRUE), ], top_n)
    lnc_frac <- mean(grepl(lnc_gene_pattern, top$gene))
    if (lnc_frac >= lnc_frac_thresh) {
      low_q <- c(low_q, as.character(cl))
      cat(sprintf("  [Low quality] Cluster %s: lncRNA fraction %.2f, top: %s\n",
                  cl, lnc_frac, paste(head(top$gene, 3), collapse = ", ")))
    }
  }
  low_q
}

## =========================================================================
## 0.8b Integrated integration + annotation function
##      (normalize → HVG → PCA → Harmony → clustering → UMAP → annotation)
## =========================================================================
## Reused by path 1 / path 2; performs full integration and annotation on
## "already merged + QC'd" objects.

integrate_and_annotate <- function(seu, batch_col = "dataset",
                                   marker_list = blood_markers,
                                   col_name = "cell_type") {
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst",
                              nfeatures = varGeneNum, verbose = FALSE)
  seu <- ScaleData(seu, features = VariableFeatures(seu),
                   vars.to.regress = vars_to_regress, verbose = FALSE)
  seu <- RunPCA(seu, features = VariableFeatures(seu),
                npcs = computePCs, verbose = FALSE)
  seu <- RunHarmony(seu, group.by.vars = batch_col,
                    theta = harmony_theta, lambda = harmony_lambda, verbose = FALSE)
  seu <- FindNeighbors(seu, reduction = "harmony", dims = pcadim)
  seu <- FindClusters(seu, resolution = resolution, verbose = FALSE)
  seu <- RunUMAP(seu, dims = pcadim, n.neighbors = umap_nneighbors,
                 reduction = "harmony", reduction.name = "umap", verbose = FALSE)
  seu <- annotate_by_markers(seu, marker_list = marker_list, col_name = col_name)
  seu
}

cat("00_config.R loaded\n")
cat("  Output directory:", out_dir, "\n")
cat("  Sepsis data:", sepsis_dir, "\n")
cat("  COPD data:", copd_dir, "\n")
