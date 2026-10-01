# =====================================================================
# 00a_download_and_prepare_data.R
# Downloads and prepares the bulk RNA-seq input data for the HMGA1
# IDH-wildtype glioblastoma pipeline and writes the files in "share/".
#
# Run ONCE, from the project folder, BEFORE 00b_microenvironment_scores.R
# and run_all.R.
#
#   Part A  Download TCGA-GBM RNA-seq (GDC, STAR - Counts) + clinical
#   Part B  Build the TCGA IDH-wildtype primary GBM cohort
#   Part C  Check the CGGA files (manual download, see below)
#   Part D  Export compact input files to share/
#
# Requirements: R >= 4.3, internet access, ~3 GB free disk space.
# =====================================================================

# ---------------------------------------------------------------------
# 0. Packages
# ---------------------------------------------------------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("TCGAbiolinks", "SummarizedExperiment"))
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
for (p in c("dplyr", "stringr"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(TCGAbiolinks); library(SummarizedExperiment)
  library(dplyr); library(stringr)
})

RAW_DIR  <- "Data/raw"          # GDC download cache
PROC_DIR <- "Data/processed"
CGGA_DIR <- "share/raw/cgga_all_grades"   # place the four CGGA files here (Part C)
SHARE    <- "share"
for (d in c(RAW_DIR, PROC_DIR, CGGA_DIR, SHARE)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Numbers obtained in this study (GDC data release used for the
# manuscript). If a later GDC release changes the data, the script
# warns but continues; report the release number printed below.
EXPECTED <- c(patients = 289, primary = 185, idh_wt = 162, final = 160, normals = 5)
check_n <- function(label, observed) {
  exp <- EXPECTED[[label]]
  if (observed == exp) cat(sprintf("  ✓ %-9s n = %d (as expected)\n", label, observed))
  else warning(sprintf("%s: n = %d, expected %d. The GDC data may have changed since the analysis.",
                       label, observed, exp), call. = FALSE)
}

# ---------------------------------------------------------------------
# Part A. Download and prepare TCGA-GBM RNA-seq
# ---------------------------------------------------------------------
gbm_file <- file.path(PROC_DIR, "gbm_data.rds")
if (file.exists(gbm_file)) {
  cat("Part A: found", gbm_file, "- skipping download\n")
  gbm_data <- readRDS(gbm_file)
} else {
  cat("Part A: querying GDC ...\n")
  info <- tryCatch(getGDCInfo(), error = function(e) NULL)
  if (!is.null(info)) cat("  GDC data release:", info$data_release, "\n")

  query <- GDCquery(
    project       = "TCGA-GBM",
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type   = c("Primary Tumor", "Solid Tissue Normal"))
  print(table(getResults(query)$sample_type))

  GDCdownload(query, method = "api", files.per.chunk = 10, directory = RAW_DIR)
  gbm_data <- GDCprepare(query, directory = RAW_DIR, summarizedExperiment = TRUE)
  saveRDS(gbm_data, gbm_file)
  if (!is.null(info)) writeLines(paste("GDC data release:", info$data_release),
                                 file.path(PROC_DIR, "GDC_release.txt"))
}
stopifnot("tpm_unstrand" %in% assayNames(gbm_data))
cat("  Samples:", ncol(gbm_data), "| genes:", nrow(gbm_data), "\n")

# ---------------------------------------------------------------------
# Part B. Build the TCGA IDH-wildtype primary glioblastoma cohort
# ---------------------------------------------------------------------
cat("\nPart B: building the IDH-wildtype primary GBM cohort\n")
clin <- as.data.frame(colData(gbm_data))
needed <- c("barcode", "patient", "sample_type", "classification_of_tumor", "paper_IDH.status",
            "vital_status", "days_to_death", "days_to_last_follow_up")
miss <- setdiff(needed, colnames(clin))
if (length(miss)) stop("Missing clinical columns: ", paste(miss, collapse = ", "))

# B1. One sample per patient: prefer Primary Tumor, then the newer
#     GDC re-processing plate (plate code starting with "A96").
clin_unique <- clin %>%
  mutate(plate_code = str_extract(barcode, "(?<=-)[A-Z0-9]{4}(?=-\\d{2}$)")) %>%
  arrange(patient,
          desc(sample_type == "Primary Tumor"),
          desc(str_starts(plate_code, "A96"))) %>%
  distinct(patient, .keep_all = TRUE)
check_n("patients", nrow(clin_unique))

# B2. Primary (not recurrent) tumors
gbm_primary <- clin_unique %>%
  filter(sample_type == "Primary Tumor", classification_of_tumor == "primary")
check_n("primary", nrow(gbm_primary))

# B3. IDH-wildtype (curated pan-glioma annotation, Ceccarelli et al. 2016)
gbm_idhwt <- gbm_primary %>% filter(paper_IDH.status == "WT")
check_n("idh_wt", nrow(gbm_idhwt))

# B4. Complete survival information
cohort <- gbm_idhwt %>%
  filter(!is.na(vital_status), vital_status %in% c("Alive", "Dead"),
         !(vital_status == "Dead"  & is.na(days_to_death)),
         !(vital_status == "Alive" & is.na(days_to_last_follow_up)))
check_n("final", nrow(cohort))

# ---------------------------------------------------------------------
# Part C. CGGA files (manual download)
# ---------------------------------------------------------------------
# From http://www.cgga.org.cn  (Download -> mRNAseq_325 and mRNAseq_693)
# download and UNZIP these four files into share/raw/cgga_all_grades/ :
#   CGGA.mRNAseq_325.RSEM-genes.20200506.txt
#   CGGA.mRNAseq_325_clinical.20200506.txt
#   CGGA.mRNAseq_693.RSEM-genes.20200506.txt
#   CGGA.mRNAseq_693_clinical.20200506.txt
cat("\nPart C: checking CGGA files\n")
cgga_files <- file.path(CGGA_DIR, c(
  "CGGA.mRNAseq_325.RSEM-genes.20200506.txt", "CGGA.mRNAseq_325_clinical.20200506.txt",
  "CGGA.mRNAseq_693.RSEM-genes.20200506.txt", "CGGA.mRNAseq_693_clinical.20200506.txt"))
absent <- cgga_files[!file.exists(cgga_files)]
if (length(absent)) stop("Download and unzip the CGGA files first (see Part C comments). Missing:\n  ",
                         paste(absent, collapse = "\n  "))
cat("  ✓ all four CGGA files present\n")

# ---------------------------------------------------------------------
# Part D. Export compact input files to share/
# ---------------------------------------------------------------------
cat("\nPart D: exporting to", SHARE, "\n")
normal_bc <- colnames(gbm_data)[colData(gbm_data)$sample_type == "Solid Tissue Normal"]
check_n("normals", length(normal_bc))
keep_bc <- c(rownames(cohort), normal_bc)
stopifnot(all(keep_bc %in% colnames(gbm_data)))

gi <- as.data.frame(rowData(gbm_data))
gi <- gi[, intersect(c("gene_id", "gene_name", "gene_type"), colnames(gi))]
saveRDS(gi, file.path(SHARE, "tcga_gene_info.rds"), compress = "xz")

tpm <- round(assay(gbm_data, "tpm_unstrand")[, keep_bc], 3)
tpm <- tpm[rowSums(tpm) > 0, ]
saveRDS(tpm, file.path(SHARE, "tcga_tpm.rds"), compress = "xz")


saveRDS(cohort, file.path(SHARE, "tcga_cohort_clinical.rds"), compress = "xz")
saveRDS(as.data.frame(colData(gbm_data))[normal_bc, ], file.path(SHARE, "tcga_normals_clinical.rds"))

for (id in c("325", "693")) {
  expr <- read.delim(file.path(CGGA_DIR, paste0("CGGA.mRNAseq_", id, ".RSEM-genes.20200506.txt")), check.names = FALSE)
  clin_c <- read.delim(file.path(CGGA_DIR, paste0("CGGA.mRNAseq_", id, "_clinical.20200506.txt")), check.names = FALSE)
  colnames(clin_c) <- sub(" \\(.*\\)$", "", colnames(clin_c))
  gbm_ids <- clin_c$CGGA_ID[trimws(clin_c$PRS_type) == "Primary" & trimws(clin_c$Histology) == "GBM"]
  expr_sub <- expr[, c(colnames(expr)[1], intersect(gbm_ids, colnames(expr)))]
  saveRDS(expr_sub, file.path(SHARE, paste0("cgga", id, "_expr_gbm.rds")), compress = "xz")
  saveRDS(clin_c,   file.path(SHARE, paste0("cgga", id, "_clinical.rds")))
  cat(sprintf("  ✓ CGGA_%s: %d primary GBM samples (expected: %s)\n",
              id, ncol(expr_sub) - 1, ifelse(id == "325", "85", "140")))
}

writeLines(capture.output(sessionInfo()), file.path(SHARE, "sessionInfo_data_preparation.txt"))
cat("\n✓ Data preparation complete.\n",
    "Next: (1) source('scripts/00b_microenvironment_scores.R')  (2) place the TISCH2 CSVs and GSE131928 files in share/ (see README)  (3) source('scripts/run_all.R')\n")
