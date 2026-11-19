# HMGA1 and the immune microenvironment of IDH-wildtype glioblastoma

Reproducible analysis pipeline for:

> Siraj Ud Din. *HMGA1 expression is associated with CD276 (B7-H3) and macrophage-related signatures, but not with reduced immune content, in IDH-wildtype glioblastoma: a three-cohort and single-cell analysis.* (manuscript under submission)

## Overview
Confounder-adjusted analysis of *HMGA1* in three IDH-wildtype primary glioblastoma cohorts (TCGA n = 160; CGGA_325 n = 74; CGGA_693 n = 109):

- partial Spearman correlations under prespecified adjustment models (M0–M6; post hoc M7)
- fixed-effect meta-analysis (Fisher z, I²) with DerSimonian–Laird and Hartung–Knapp–Sidik–Jonkman (HKSJ) sensitivity estimates
- Benjamini–Hochberg correction within hypothesis families
- Cox models with Schoenfeld tests, time-varying sensitivity analysis and numbers-at-risk tables
- single-cell localization in Neftel et al. 2019 (GSE131928, adult Smart-seq2 tumors), including an expression-based chromosome 7/10 copy-number proxy for malignant cells
- quality control for duplicated patients between the two CGGA cohorts

## Scripts
| Script | Content |
|---|---|
| `scripts/00_functions.R` | settings, gene sets, adjustment models, statistical and plotting helpers (sourced by 02 and 03) |
| `scripts/01_prepare_data.R` | run once: downloads TCGA (GDC), builds the cohorts, prepares CGGA files, computes ESTIMATE and MCP-counter scores |
| `scripts/02_bulk_analyses.R` | cohorts, duplicate check, survival (Figure 1, Fig. S1), immune composition (Figure 2), CCL2/TAM and STING/IFN (Figure 3), checkpoints and subtype (Figure 4, Fig. S7), all-grade analysis (Fig. S6), cohort characteristics |
| `scripts/03_single_cell_and_tables.R` | TISCH2 (Fig. S2), Neftel single-cell analysis and copy-number proxy (Figure 5, Figs. S3–S5, S8), myeloid-restricted TAM signatures and bootstrap comparison, Supplementary Tables workbook |

## Data (public; not redistributed here)
| Data | Source | Place in |
|---|---|---|
| TCGA GBM expression (STAR, TPM) and clinical; GDC Data Release 45.0 | GDC via TCGAbiolinks (downloaded by `01_prepare_data.R`) | `share/` |
| CGGA mRNAseq_325 / mRNAseq_693 (RSEM, FPKM) and clinical, all grades (4 files, manual download) | http://www.cgga.org.cn | `share/raw/cgga_all_grades/` |
| Neftel Smart-seq2 processed TPM (GSM3828672) and cell annotation | GEO GSE131928 | `share/raw/` |
| TISCH2 cell-type means (HMGA1, HAVCR2, CCL2, CD274) | http://tisch.compbio.cn | `share/` |

Gene sets included in `share/`: Wang et al. 2017 (`wang2017_subtype_signatures.csv`), Müller et al. 2017 (`muller2017_tam_signatures.csv`), Neftel et al. 2019 (`neftel2019_modules.csv`). The MSigDB Hallmark interferon-α set is loaded with `msigdbr`.

## Requirements
R 4.6.0 with survival 3.8-9, survminer 0.5.2, ggplot2 4.0.3, patchwork 1.3.2, data.table 1.18.4, readxl 1.5.0, writexl, GSVA 2.6.3, msigdbr 26.1.1, R.utils, BiocManager, org.Hs.eg.db 3.23.1; for data preparation TCGAbiolinks 2.40.0, SummarizedExperiment, dplyr, stringr, estimate 1.0.13, immunedeconv 2.1.4 and MCPcounter 1.2.0.

## Run
```r
setwd("path/to/project")                            # folder containing share/ and scripts/
source("scripts/01_prepare_data.R")                 # once; does nothing if share/ is complete
source("scripts/02_bulk_analyses.R")                # writes results/log_02_bulk.txt
source("scripts/03_single_cell_and_tables.R")       # writes results/log_03_single_cell.txt
```
Re-run selected sections only, for example:
```r
ONLY <- c("survival"); source("scripts/02_bulk_analyses.R")
ONLY <- c("tables");   source("scripts/03_single_cell_and_tables.R")
```

Outputs: `results/` (tables, logs, `Supplementary_Tables_S1-S5.xlsx`) and `figures/` (PDF and PNG).

## License
MIT License (see `LICENSE`).
