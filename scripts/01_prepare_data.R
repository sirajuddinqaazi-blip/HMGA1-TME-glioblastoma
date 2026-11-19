# =====================================================================
# 01_prepare_data.R   (run ONCE; merges the former 00a and 00b scripts)
# Part 1  Download TCGA-GBM RNA-seq (GDC, STAR - Counts) + clinical,
#         build the TCGA IDH-wildtype primary GBM cohort, check and
#         export the CGGA files (manual download) to share/
# Part 2  ESTIMATE immune/stromal scores and MCP-counter populations
#         -> share/scores_TCGA.csv, scores_CGGA325.csv, scores_CGGA693.csv
#
# If share/ is already complete the script does nothing. To rebuild
# everything, run:  FORCE <- TRUE; source("scripts/01_prepare_data.R")
# The analysis code is unchanged from 00a/00b.
# Requirements: R >= 4.3, internet access, ~3 GB free disk space;
# estimate 1.0.13, immunedeconv 2.1.4, MCPcounter 1.2.0
# =====================================================================
if (!exists("FORCE")) FORCE <- FALSE
needed_files <- file.path("share", c("tcga_tpm.rds", "tcga_gene_info.rds", "tcga_cohort_clinical.rds",
                                      "tcga_normals_clinical.rds", "cgga325_expr_gbm.rds", "cgga693_expr_gbm.rds",
                                      "cgga325_clinical.rds", "cgga693_clinical.rds",
                                      "scores_TCGA.csv", "scores_CGGA325.csv", "scores_CGGA693.csv"))

prepare_part1 <- function() {
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
  cat("\n✓ Part 1 complete.\n")
}

prepare_part2 <- function() {
  library(estimate)
  library(immunedeconv)

  dir.create("share/tmp_estimate", recursive = TRUE, showWarnings = FALSE)

  # ---------- helpers ----------
  collapse_by_symbol <- function(mat, symbols) {
    keep <- !is.na(symbols) & symbols != ""
    mat <- mat[keep, , drop = FALSE]; symbols <- symbols[keep]
    mean_expr <- rowMeans(mat)
    idx <- tapply(seq_len(nrow(mat)), symbols, function(i) i[which.max(mean_expr[i])])
    out <- mat[unlist(idx), , drop = FALSE]
    rownames(out) <- names(idx)
    out
  }

  run_estimate <- function(mat, label) {
    in_f  <- file.path("share/tmp_estimate", paste0(label, "_input.txt"))
    cg_f  <- file.path("share/tmp_estimate", paste0(label, "_common.gct"))
    out_f <- file.path("share/tmp_estimate", paste0(label, "_estimate.gct"))
    df <- data.frame(NAME = rownames(mat), log2(mat + 1), check.names = FALSE)
    write.table(df, in_f, sep = "\t", quote = FALSE, row.names = FALSE)
    filterCommonGenes(input.f = in_f, output.f = cg_f, id = "GeneSymbol")
    estimateScore(input.ds = cg_f, output.ds = out_f, platform = "illumina")
    sc <- read.table(out_f, skip = 2, header = TRUE, sep = "\t", check.names = FALSE)
    rownames(sc) <- sc$NAME
    sc <- t(sc[, -(1:2)])
    data.frame(sample = colnames(mat),        # keep original IDs (read.table may alter them)
               ESTIMATE_stromal = as.numeric(sc[, "StromalScore"]),
               ESTIMATE_immune  = as.numeric(sc[, "ImmuneScore"]),
               ESTIMATE_score   = as.numeric(sc[, "ESTIMATEScore"]),
               check.names = FALSE)
  }

  run_mcp <- function(mat) {
    res <- immunedeconv::deconvolute(mat, "mcp_counter")
    m <- as.data.frame(t(as.matrix(res[, -1])))
    colnames(m) <- paste0("MCP_", gsub("[^A-Za-z0-9]+", "_", res$cell_type))
    m$sample <- rownames(m)
    m
  }

  process <- function(mat, label) {
    cat("\n=== ", label, ": ", nrow(mat), " genes x ", ncol(mat), " samples ===\n", sep = "")
    est <- run_estimate(mat, label)
    mcp <- run_mcp(mat)
    out <- merge(est, mcp, by = "sample")
    stopifnot(nrow(out) == ncol(mat))
    write.csv(out, paste0("share/scores_", label, ".csv"), row.names = FALSE)
    cat("✓ Saved share/scores_", label, ".csv\n", sep = "")
    print(summary(out$ESTIMATE_immune))
  }

  # ---------- TCGA ----------
  tpm    <- readRDS("share/tcga_tpm.rds")
  gi     <- readRDS("share/tcga_gene_info.rds")
  cohort <- readRDS("share/tcga_cohort_clinical.rds")
  sym    <- gi$gene_name[match(rownames(tpm), rownames(gi))]
  tcga   <- collapse_by_symbol(tpm[, rownames(cohort)], sym)
  process(tcga, "TCGA")

  # ---------- CGGA ----------
  for (id in c("325", "693")) {
    ex <- readRDS(paste0("share/cgga", id, "_expr_gbm.rds"))
    mat <- as.matrix(ex[, -1]); mode(mat) <- "numeric"
    mat <- collapse_by_symbol(mat, ex[[1]])
    process(mat, paste0("CGGA", id))
  }

  cat("\n✓ Part 2 complete: share/scores_TCGA.csv, share/scores_CGGA325.csv, share/scores_CGGA693.csv\n")
}

if (all(file.exists(needed_files)) && !FORCE) {
  message("share/ is already complete -> nothing to do. (Set FORCE <- TRUE to rebuild.)")
} else {
  prepare_part1()
  prepare_part2()
  cat("\nNext: place the TISCH2 CSVs and the GSE131928 files in share/ and share/raw/ (see README), then run\n",
      "  source('scripts/02_bulk_analyses.R'); source('scripts/03_single_cell_and_tables.R')\n")
}
