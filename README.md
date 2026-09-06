# COPD_and_sepsis

Custom analysis code for the study:

> **A shared myeloid program links COPD to sepsis, and nominates SOCS3 as a candidate regulator**
> (Manuscript under review at *Respiratory Research*)

The repository contains the numbered R (and a few Python) analysis scripts used to
generate all main and supplementary figures and the supplementary tables of the
manuscript. The pipeline is organised so that running the scripts in ascending
numerical order reproduces the analyses on the public datasets listed below.

---

## Contents

- `00_config.R` — global configuration: data directories, colours
  (`GROUP_COLS`), cohort maps, download/caching switches, output directory.
- `01_*`–`02_*` — single-cell data assembly: six blood PBMC datasets
  (GSE151263, GSE167363, GSE175453, GSE249584, GSE279452, SCP548) and the lung
  atlas GSE136831; quality control, adult-only filtering, Seurat v5 + Harmony
  integration, and the cross-tissue (lung myeloid) object.
- `03_*`–`13_*` — subclustering, annotation, enrichment, cell–cell
  communication, WGCNA / single-cell hdWGCNA and related analyses.
- `14_*`–`16_*`, `23_*`, `24_*` — bulk validation and external prognosis
  cohorts (GSE66099, GSE57148, GSE248493, GSE65682, E-MTAB-4451, GSE95233).
- `17_*`–`22_*`, `25_*`–`27_*` — complementary mechanistic / trajectory /
  severity analyses (some reflect earlier lines of the project that were not
  carried into the final manuscript).
- `28_*`–`36_*` — donor-level pseudobulk differential expression, severity
  gradient, meta-analysis, leave-one-dataset-out, ROC reconstruction and
  supplementary analyses.
- `37_*`–`66_*` — main and supplementary **figure rendering** scripts
  (Fig. 1–6, Fig. S1–S9), including the unified style helper scripts.
- `67_*` — hypergeometric background recomputation.
- `68_permutation_overlap_test.R` — **donor-label permutation test** of the cross-tissue overlap
  (10,000 permutations; reproducible via `MASTER_SEED`; parallel with
  `NPROC`; checkpoint/resume built in).
- `69_supplementary_outputs.R`, `69b_intersection_gene_table.R`, `70_candidate_survival_checks.R` — auxiliary outputs requested during review (per-group
  Spearman statistics, observed intersection gene lists, single-gene/composite
  survival checks).
- `make_GPL13667_symbol_map.py`, `62*.py` — small helpers used by specific
  cohorts (probe-to-symbol maps) or for vector figure assembly.

Output files written by the scripts (PDF/PNG figures, CSV/`rds` tables) are
**not** committed; they are produced into the directory given by the scripts /
`00_config.R` (`OUT_DIR`), which defaults to the script directory.

---

## Software requirements

- R ≥ 4.3 (developed under R 4.4.3), Bioconductor ≥ 3.19
- Key packages:
  `Seurat` (v5), `Harmony`, `edgeR`, `Matrix`, `WGCNA`, `hdWGCNA`,
  `CellChat` (v2), `Monocle3`, `clusterProfiler`, `GEOquery`, `survival`,
  `glmnet`, `pROC`, `metafor` (or `rmeta`), `ggplot2`, `dplyr`, `tidyr`,
  `ggrepel`, `patchwork`, `RColorBrewer`
- Python 3 (helpers): `pymupdf`/`fitz`, `Pillow` for vector figure assembly

Install missing Bioconductor packages with
`BiocManager::install("edgeR")` and CRAN packages with
`install.packages(...)`.

---

## Input data

All datasets are public:

| Dataset | Type | GEO / ArrayExpress |
|---|---|---|
| GSE151263 | Blood PBMC scRNA-seq | GEO GSE151263 |
| GSE167363 | Blood PBMC scRNA-seq | GEO GSE167363 |
| GSE175453 | Blood PBMC scRNA-seq | GEO GSE175453 |
| GSE249584 | Blood PBMC scRNA-seq (COPD) | GEO GSE249584 |
| GSE279452 | Blood PBMC scRNA-seq (sepsis atlas) | GEO GSE279452 |
| SCP548 | Blood PBMC scRNA-seq (paediatric) | SCP (CZ Biohub) SCP548 |
| GSE136831 | Lung tissue atlas (COPD vs control) | GEO GSE136831 |
| GSE66099 | Bulk whole-blood (sepsis) | GEO GSE66099 |
| GSE57148 | Bulk lung (COPD) | GEO GSE57148 |
| GSE248493 | Bulk whole-blood RNA-seq (COPD signature) | GEO GSE248493 |
| GSE65682 | Whole-blood microarray (MARS; 28-day survival) | GEO GSE65682 |
| E-MTAB-4451 | Whole-blood leukocyte RNA (28-day survival) | ArrayExpress E-MTAB-4451 |
| GSE95233 | Whole-blood microarray (28-day survival) | GEO GSE95233 |

Intermediate Seurat objects (e.g. the integrated blood object and the lung
object used by scripts `28_*`–`70_*`) are large and are **not** included;
`00_config.R` contains the expected paths. Scripts that can fetch data online
fall back to `GEOquery` downloads when local `*_series_matrix.txt.gz` files are
not present (see `00_config.R`).

---

## Running

Run scripts from within the repository directory (or pass an explicit output
directory as the first argument where supported). For example:

```bash
# configuration paths must be adapted to your local data layout first
Rscript 00_config.R                 # config is sourced by all other scripts
Rscript 28_pseudobulk_donor_deg.R
# the 10,000-permutation test can take hours; run in parallel and it
# checkpoints automatically:
NPROC=64 Rscript 68_permutation_overlap_test.R
```

Notes:

- File names contain non-ASCII (Chinese) characters and are UTF-8 encoded;
  they work on Linux/macOS and on Windows R ≥ 4.2 with UTF-8 enabled.
- Some figure scripts assume the working directory layout of the original
  analysis server (an `01script/` folder with a shared `00_config.R`).
- The permutation script (`68_permutation_overlap_test.R`) uses a fixed `MASTER_SEED` so results are
  reproducible; it can be resumed from its checkpoint file if interrupted.

---

## Output → manuscript mapping

- Main figures: `Fig1`–`Fig6` (assembled panels from `38_*`–`66_*` outputs).
- Supplementary figures: `FigS1`–`FigS9` (incl. vector assembly with the
  `62*.py` helpers).
- Supplementary tables ST1–ST9 (`Supplementary_Tables.xlsx` + one CSV per
  sheet): assembled from the CSV outputs of the corresponding scripts
  (e.g. `28_*`, `30_*`, `31_*`, `68_permutation_overlap_test.R`–`70_*`).
- Key statistics (hypergeometric/permutation P values, donor-level LODO,
  meta-analysis HR/CI, per-group Spearman) are printed by the respective
  scripts to `path*.txt` / `path*.csv` outputs.

---

## Citation

If you use or adapt this code, please cite the manuscript (details will be
updated once published).
