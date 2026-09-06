#!/usr/bin/env Rscript
# =============================================================================
# 29_dataset_donor_summary.R
# Purpose: output donor counts + cell counts per single-cell dataset, for use in
#       Table 1 (dataset summary table). Reviewers want to see "donor counts per
#       cohort" at first glance; it is not enough to report only COPD blood 7/8
#       and COPD lung 17/15.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

this_file <- commandArgs(trailingOnly = FALSE)
.f <- grep("--file=", this_file, value = TRUE)
if (length(.f) > 0) .this_file <- sub("^--file=", "", .f[1]) else .this_file <- "."
out_dir <- dirname(normalizePath(.this_file))
source(file.path(out_dir, "00_config.R"))

## ---- Donor field detection (patient first, sample as fallback) ----
detect_donor <- function(seu) {
  meta <- seu@meta.data
  get <- function(f) if (f %in% colnames(meta)) as.character(meta[[f]]) else NULL
  patient <- get("patient"); sample <- get("sample")
  if (!is.null(patient)) {
    if (!is.null(sample)) {
      na_idx <- is.na(patient) | patient == ""
      patient[na_idx] <- sample[na_idx]
    }
    return(patient)
  }
  for (cand in c("donor_id","donor","subject","subject_id","Subject_Identity","sample")) {
    v <- get(cand); if (!is.null(v)) return(v)
  }
  NULL
}

## ---- Blood object: donor count + cell count by dataset x group ----
cat("########## Blood object (dataset x group) ##########\n")
blood_file <- file.path(out_dir, "path1_sepsis_copd_integrated.rds")
if (!file.exists(blood_file)) stop("Cannot find path1_sepsis_copd_integrated.rds")
blood <- readRDS(blood_file)
donor <- detect_donor(blood)
meta <- blood@meta.data

ds <- sort(unique(meta$dataset))
grp <- c("Healthy","Infection_Control","COPD","Sepsis","Sepsis_Pneumonia")

out <- list()
for (d in ds) {
  for (g in grp) {
    idx <- meta$dataset == d & !is.na(meta$group) & meta$group == g
    if (sum(idx) == 0) next
    nd <- length(unique(donor[idx]))
    nc <- sum(idx)
    out[[length(out) + 1]] <- data.frame(
      dataset = d, tissue = "blood", group = g,
      n_donors = nd, n_cells = nc, stringsAsFactors = FALSE)
  }
}
blood_tab <- do.call(rbind, out)

# Total donor count per dataset (deduplicated)
tot_donor <- sapply(ds, function(d) length(unique(donor[meta$dataset == d])))
tot_cell  <- sapply(ds, function(d) sum(meta$dataset == d))
cat("Total donors/cells per dataset in the blood object:\n")
for (d in ds) cat("  ", d, ": ", tot_donor[d], " donors / ", tot_cell[d], " cells\n", sep = "")
cat("\nBlood object dataset x group details:\n")
print(blood_tab)

## ---- Lung object: GSE136831 disease x donor count ----
cat("\n########## Lung object GSE136831 ##########\n")
lung_seu <- tryCatch(
  read_gse136831(file.path(copd_dir, "GSE136831"), gse_id = "GSE136831"),
  error = function(e) { cat("Failed to read lung object:", conditionMessage(e), "\n"); NULL })
lung_tab <- NULL
if (!is.null(lung_seu)) {
  ldonor <- detect_donor(lung_seu)
  lmeta <- lung_seu@meta.data
  for (g in sort(unique(lmeta$disease))) {
    idx <- !is.na(lmeta$disease) & lmeta$disease == g
    lung_tab <- rbind(lung_tab, data.frame(
      dataset = "GSE136831", tissue = "lung", group = g,
      n_donors = length(unique(ldonor[idx])), n_cells = sum(idx),
      stringsAsFactors = FALSE))
  }
  cat("Lung object disease x donor count:\n"); print(lung_tab)
}

## ---- Save summary ----
final <- rbind(blood_tab, lung_tab)
write.csv(final, file.path(out_dir, "path29_dataset_donor_counts.csv"), row.names = FALSE)
cat("\nSaved: path29_dataset_donor_counts.csv\n")
