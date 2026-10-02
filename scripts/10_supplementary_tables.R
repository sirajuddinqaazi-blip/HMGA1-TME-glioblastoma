# =====================================================================
# 10_supplementary_tables.R
# Assembles Supplementary Tables S1-S5 from the pipeline result files
# into one Excel workbook: results/Supplementary_Tables_S1-S5.xlsx
# Read-only with respect to all analyses (no numbers are recomputed).
# =====================================================================
source("scripts/00_config.R", local = TRUE)
if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")

rd <- function(f) { p <- file.path(RES_DIR, f); if (!file.exists(p)) stop("Missing result file: ", p, " (run run_all.R first)"); fread(p) }
rnd <- function(dt) { dt <- as.data.frame(dt)
  for (c in names(dt)) if (is.numeric(dt[[c]])) dt[[c]] <- signif(dt[[c]], 4)
  dt }

# ---- S1: all association results ----
s1 <- rbind(rd("T3_H4_immune_composition.csv"), rd("T4_H1_H2_mechanism_axes.csv"), rd("T5_H3_checkpoints.csv"),
            rd("T5b_HMGA1_vs_subtype_scores.csv"), rd("T5d_CD276_sensitivity.csv"),
            if (file.exists(file.path(RES_DIR, "T11_tam_restricted.csv"))) rd("T11_tam_restricted.csv"), fill = TRUE)
s1$model_description <- unname(c(MODEL_LABELS, M3b_TCGA_labels = "+ prolif + ESTIMATE + subtype labels (TCGA)")[s1$model])
s1$post_hoc <- s1$model %in% c("M7_vascular_stromal") | grepl("^TAM_", s1$outcome)
setcolorder(s1, c("family", "outcome", "model", "model_description", "post_hoc"))

# ---- S2: survival ----
s2a <- rd("T1_HMGA1_proliferation_age.csv")
s2b <- rd("T2_HMGA1_cox.csv")

# ---- S3: single cell ----
s3a <- rd("T7a_neftel_cell_counts.csv"); setnames(s3a, 1, "tumor")
s3b <- rd("T7b_neftel_HMGA1_by_class.csv")
s3c <- rd("T7e_neftel_targets_by_class.csv")
s3d <- rd("T7d_neftel_within_malignant.csv")
s3e <- rd("T7f_neftel_malignant_pseudobulk.csv")
s3f <- rd("T6_TISCH2_Neftel_Smartseq2.csv")
s3g <- if (file.exists(file.path(RES_DIR, "T7g_muller_gene_specificity.csv"))) rd("T7g_muller_gene_specificity.csv") else data.table()

# ---- S4: all-grade analysis ----
s4 <- rd("T8_all_grades_reconciliation.csv")

# ---- S5: gene sets ----
gs <- list()
add <- function(set, genes, source) gs[[length(gs) + 1]] <<- data.table(gene_set = set, gene = genes, source = source)
add("Proliferation score", PROLIF_GENES, "This study (10 canonical proliferation genes)")
add("He et al. ISGs", SIGNATURES$He_ISG, "He et al., Nat Commun 2025")
add("Checkpoint genes", CHECKPOINTS, "This study (prespecified, H3)")
for (f in c(WANG_FILE, MULLER_FILE, NEFTEL_FILE)) if (file.exists(f)) {
  w <- fread(f); src <- c("Wang et al., Cancer Cell 2017 (Table S1)", "Mueller et al., Genome Biol 2017 (Table S4)",
                          "Neftel et al., Cell 2019 (Table S2)")[match(f, c(WANG_FILE, MULLER_FILE, NEFTEL_FILE))]
  add(paste(basename(tools::file_path_sans_ext(f)), w$set, sep = ": "), w$gene, src)
}
if (requireNamespace("msigdbr", quietly = TRUE)) {
  h <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
                error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
  add("HALLMARK_INTERFERON_ALPHA_RESPONSE", sort(unique(h$gene_symbol[h$gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE"])),
      paste0("Liberzon et al., Cell Syst 2015; msigdbr ", packageVersion("msigdbr")))
}
add("Neftel cell-class markers: myeloid", c("CD14", "AIF1", "FCER1G", "FCGR3A", "TYROBP", "CSF1R"), "Neftel et al., Cell 2019 (STAR Methods)")
add("Neftel cell-class markers: T cell", c("CD2", "CD3D", "CD3E", "CD3G"), "Neftel et al., Cell 2019 (STAR Methods)")
add("Neftel cell-class markers: oligodendrocyte", c("MBP", "TF", "PLP1", "MAG", "MOG", "CLDN11"), "Neftel et al., Cell 2019 (STAR Methods)")
add("Endothelial markers", c("PECAM1", "VWF", "CLDN5"), "This study")
s5 <- rbindlist(gs)

readme <- data.frame(Sheet = c("S1_associations", "S2a_proliferation_age", "S2b_cox", "S3a_cell_counts", "S3b_HMGA1_by_class",
                               "S3c_targets_by_class", "S3d_within_malignant", "S3e_pseudobulk", "S3f_TISCH2", "S3g_TAM_gene_specificity", "S4_all_grades", "S5_gene_sets"),
  Content = c(
    "Partial Spearman correlations of HMGA1 with all outcomes: per-cohort r, P and n; pooled fixed-effect r (95% CI), P, Benjamini-Hochberg q (within family and model), I2, direction consistency, robustness flag and random-effects (DerSimonian-Laird) sensitivity estimates (r_RE, 95% CI, P). Models M0-M7 as defined in Methods.",
    "Spearman correlations of HMGA1 with proliferation score, MKI67 and age (per cohort and pooled).",
    "Cox models: HR per 1 SD log2 HMGA1, 95% CI, P, Schoenfeld tests, sensitivity models and pooled estimates.",
    "Neftel Smart-seq2 adult tumors: cells per tumor and class (Neftel markers, cutoff 4). Proportions reflect CD45 sorting.",
    "HMGA1 in malignant vs myeloid/T cells, paired within tumors, at cutoffs 3, 4 and 5.",
    "Mean log2(TPM/10+1) and fraction detected of HMGA1 and target genes by cell class.",
    "Within-malignant-cell correlations with HMGA1 per target, pooled across tumors (unadjusted and cell-cycle/complexity adjusted).",
    "Tumor-level malignant-cell pseudo-bulk correlations (n = 20 tumors).",
    "TISCH2 cell-type means for the Neftel Smart-seq2 dataset (independent cross-check).",
    "Mueller TAM signature genes: mean expression and detection in myeloid vs malignant cells; myeloid-restricted flag (prespecified criterion).",
    "All-grade CGGA analysis (post hoc): HR per 1 SD HMGA1 and MKI67 (positive control) across sequential models; within-molecular-group estimates.",
    "All gene sets used, with sources."))

out <- file.path(RES_DIR, "Supplementary_Tables_S1-S5.xlsx")
writexl::write_xlsx(list(README = readme, S1_associations = rnd(s1), S2a_proliferation_age = rnd(s2a), S2b_cox = rnd(s2b),
                         S3a_cell_counts = rnd(s3a), S3b_HMGA1_by_class = rnd(s3b), S3c_targets_by_class = rnd(s3c),
                         S3d_within_malignant = rnd(s3d), S3e_pseudobulk = rnd(s3e), S3f_TISCH2 = rnd(s3f), S3g_TAM_gene_specificity = rnd(s3g),
                         S4_all_grades = rnd(s4), S5_gene_sets = as.data.frame(s5)), out)
say("Saved %s (S1: %d rows; S5: %d genes)", out, nrow(s1), nrow(s5))
