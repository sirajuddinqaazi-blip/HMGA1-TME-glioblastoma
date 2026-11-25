# =====================================================================
# 04_cptac_proteomics.R
# Protein-level analysis:
#   1 cptac  HMGA1 vs CD276 / CCL2 / HAVCR2 / STING1 PROTEIN in CPTAC-GBM
#            (Wang et al., Cancer Cell 2021; TMT global proteome), IDH-wt only;
#            optional RNA-protein concordance in the same tumors
#   2 hpa    Human Protein Atlas glioma immunohistochemistry summary
#
# Manual downloads:
#   share/raw/cptac_gbm/  (LinkedOmics, https://www.linkedomics.org/data_download/CPTAC-GBM/)
#     HS_CPTAC_GBM_proteome_mssm_per_gene.cct                 gene-level proteome (log2 ratio)
#     HS_CPTAC_GBM_rnaseq_fpkm_uq_log2.cct                    RNA-seq (optional; RNA-protein check)
#     somaticwrapper_all_cases_filtered.v4.0.20200430.maf.gz  somatic mutations (fallback for IDH status)
#   share/raw/
#     NIHMS1665743-supplement-2.xlsx   Wang et al. 2021 Table S1 (IDH1 status; primary source)
#     pathology.tsv                    Human Protein Atlas cancer tissue data (optional)
#
# Usage:  source("scripts/04_cptac_proteomics.R")
# Log: results/log_04_cptac.txt
# =====================================================================
source("scripts/00_functions.R")

CPTAC_DIR <- file.path(RAW_DIR, "cptac_gbm")
PROT_TARGETS <- c("CD276", "CCL2", "HAVCR2", "STING1", "CD274")
MIN_NONMISSING <- 0.5     # a protein is analysed only if quantified in >= 50% of tumors

# read a genes x samples table (first column = gene symbol), return numeric matrix
read_matrix <- function(pattern) {
  f <- list.files(CPTAC_DIR, pattern = pattern, ignore.case = TRUE, full.names = TRUE)
  if (!length(f)) return(NULL)
  if (length(f) > 1) say("  several files match '%s'; using %s", pattern, basename(f[1]))
  x <- fread(f[1], check.names = FALSE)
  m <- as.matrix(x[, -1, with = FALSE]); storage.mode(m) <- "numeric"
  rownames(m) <- as.character(x[[1]])
  m <- m[!is.na(rownames(m)) & !duplicated(rownames(m)), , drop = FALSE]
  say("  %s: %d genes x %d samples", basename(f[1]), nrow(m), ncol(m))
  m
}

sec_cptac <- function() {
  P <- read_matrix("_proteome_.*per_gene")   # gene-level proteome (not phosphoproteome)
  if (is.null(P)) stop("No proteome file in ", CPTAC_DIR, " (expected HS_CPTAC_GBM_proteome_mssm_per_gene.cct).")
  say("Value range %.1f to %.1f (expected: log2 ratios centred near 0)", min(P, na.rm = TRUE), max(P, na.rm = TRUE))

  # ---- IDH-wildtype only ----
  # Primary source: Wang et al. 2021 Table S1 (sheet additional_annotations, IDH1_alter_status).
  # WHO 2021: only hotspot (R132) mutations define IDH-mutant tumors; the non-hotspot R222C tumor
  # (not G-CIMP) is kept as IDH-wildtype. Fallback: MAF file from LinkedOmics.
  base_id <- function(x) sub("[_.].*$", "", x)   # e.g. C3L-00104_T -> C3L-00104
  ts1 <- file.path(RAW_DIR, "NIHMS1665743-supplement-2.xlsx")
  mf <- list.files(CPTAC_DIR, pattern = "\\.maf(\\.gz)?$", full.names = TRUE)
  if (file.exists(ts1)) {
    aa <- as.data.frame(read_excel(ts1, sheet = "additional_annotations"))
    mut <- aa$case[aa$IDH1_alter_status %in% "Missense"]
    kept <- aa$case[aa$IDH1_alter_status %in% "Missense_not_hotspot"]
    drop <- base_id(colnames(P)) %in% mut
    say("Wang 2021 Table S1: IDH1 hotspot-mutant tumors excluded: %d (%s); non-hotspot IDH1 kept as IDH-wt: %s",
        sum(drop), paste(colnames(P)[drop], collapse = ", "), paste(kept, collapse = ", "))
    P <- P[, !drop, drop = FALSE]
  } else if (length(mf)) {
    if (grepl("\\.gz$", mf[1]) && !requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
    maf <- fread(mf[1], select = c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode"))
    hit <- maf[Hugo_Symbol %in% c("IDH1", "IDH2") & Variant_Classification == "Missense_Mutation"]
    mut <- unique(base_id(hit$Tumor_Sample_Barcode))
    drop <- base_id(colnames(P)) %in% mut
    say("MAF: %d tumors with mutation data; IDH1/IDH2 missense-mutant tumors excluded: %d (%s)",
        length(unique(maf$Tumor_Sample_Barcode)), sum(drop), paste(colnames(P)[drop], collapse = ", "))
    if (length(mut) && !any(drop)) say("WARNING: IDH-mutant IDs did not match proteome column names - check: %s vs %s",
                                       paste(head(mut, 3), collapse = ","), paste(head(colnames(P), 3), collapse = ","))
    P <- P[, !drop, drop = FALSE]
  } else say("WARNING: no .maf file -> IDH-mutant tumors NOT excluded (state this in Methods).")
  say("Tumors analysed: %d", ncol(P))

  get <- function(g) { s <- find_symbol(g, rownames(P)); if (is.na(s)) rep(NA_real_, ncol(P)) else P[s, ] }
  if (all(is.na(get("HMGA1")))) stop("HMGA1 not quantified in the proteome file.")
  d <- data.frame(sample = colnames(P), HMGA1 = get("HMGA1"),
                  sapply(PROT_TARGETS, get), check.names = FALSE)

  # protein proliferation score: z-mean of proliferation proteins quantified in >= 80% of tumors
  pg <- intersect(PROLIF_GENES, rownames(P))
  pg <- pg[rowMeans(!is.na(P[pg, , drop = FALSE])) >= 0.8]
  say("Proliferation proteins used (%d): %s", length(pg), paste(pg, collapse = ", "))
  d$prolif <- if (length(pg) >= 3) rowMeans(scale(t(P[pg, , drop = FALSE])), na.rm = TRUE) else NA

  cat("\n---- HMGA1 protein vs target proteins (partial Spearman) ----\n")
  res <- rbindlist(lapply(PROT_TARGETS, function(y) {
    fr <- mean(!is.na(d[[y]]))
    if (fr < MIN_NONMISSING) {
      say("%-7s quantified in %3.0f%% of tumors -> not analysed", y, 100 * fr)
      return(data.table(target = y, frac_quantified = fr, note = "not quantified"))
    }
    rbindlist(lapply(list("none" = character(0), "+ proliferation" = "prolif"), function(Z) {
      if (length(Z) && all(is.na(d$prolif))) return(NULL)
      v <- pspear("HMGA1", y, Z, d)
      data.table(target = y, frac_quantified = fr, adjustment = if (length(Z)) "+ proliferation" else "none",
                 n = v[["n"]], rho = v[["r"]], p = v[["p"]])
    }))
  }), fill = TRUE)
  for (i in seq_len(nrow(res))) if (!is.na(res$rho[i]))
    say("%-7s %-16s n = %2d  rho = %5.2f  P = %s", res$target[i], res$adjustment[i], res$n[i], res$rho[i], fmt_p(res$p[i]))
  say("NOTE: one cohort, nominal P values; protein-level check of the bulk RNA findings.")
  fwrite(res, file.path(RES_DIR, "T12_cptac_protein.csv"))

  # ---- RNA-protein concordance (optional) ----
  R <- read_matrix("rnaseq")
  if (!is.null(R)) {
    common <- intersect(colnames(P), colnames(R))
    say("\nRNA-protein concordance in %d shared tumors:", length(common))
    if (max(R, na.rm = TRUE) > 50) R <- log2(R + 1)
    conc <- rbindlist(lapply(c("HMGA1", PROT_TARGETS), function(g) {
      sp <- find_symbol(g, rownames(P)); sr <- find_symbol(g, rownames(R))
      if (is.na(sp) || is.na(sr)) return(NULL)
      ok <- !is.na(P[sp, common]) & !is.na(R[sr, common]); if (sum(ok) < 10) return(NULL)
      ct <- cor.test(P[sp, common][ok], R[sr, common][ok], method = "spearman", exact = FALSE)
      data.table(gene = g, n = sum(ok), rho_RNA_protein = unname(ct$estimate), p = ct$p.value)
    }))
    print(conc, digits = 2)
    fwrite(conc, file.path(RES_DIR, "T12b_cptac_rna_protein.csv"))
  }

  # ---- Figure S9 ----
  shown <- intersect(c("CD276", "CCL2"), res$target[!is.na(res$rho)])
  if (length(shown)) {
    pd <- rbindlist(lapply(shown, function(y) data.table(target = pretty_lab(y), HMGA1 = d$HMGA1, y = d[[y]])))
    lab <- res[target %in% shown & adjustment == "none"]; lab$target <- pretty_lab(lab$target)
    p <- ggplot(pd[!is.na(y) & !is.na(HMGA1)], aes(HMGA1, y)) + geom_point(size = 1, alpha = 0.7) +
      facet_wrap(~target, scales = "free_y") +
      geom_text(data = lab, aes(x = -Inf, y = Inf, label = sprintf("rho = %.2f, P = %s", rho, fmt_p(p))),
                hjust = -0.1, vjust = 1.4, size = 3) +
      labs(x = "HMGA1 protein (log2 ratio)", y = "Protein (log2 ratio)",
           title = "CPTAC-GBM, IDH-wildtype tumors: protein level") + theme_pub
    save_fig(p, "FigureS9_CPTAC_protein", 4 + 3 * length(shown), 4)
    cat("OK Figure S9 saved\n")
  }
}

sec_hpa <- function() {
  f <- file.path(RAW_DIR, "pathology.tsv")
  if (!file.exists(f)) { say("pathology.tsv not found -> HPA section skipped (use the HPA web pages instead)."); return(invisible(NULL)) }
  h <- fread(f)
  gcol <- grep("^Gene name$", names(h), value = TRUE)[1]; ccol <- grep("^Cancer$", names(h), value = TRUE)[1]
  lv <- intersect(c("High", "Medium", "Low", "Not detected"), names(h))
  out <- h[get(gcol) %in% c("HMGA1", "CD276", "CCL2", "HAVCR2") & grepl("glioma", get(ccol), ignore.case = TRUE),
           c(gcol, ccol, lv), with = FALSE]
  cat("\n---- Human Protein Atlas: glioma IHC (number of patients per staining level) ----\n"); print(out)
  fwrite(out, file.path(RES_DIR, "T12c_HPA_glioma_IHC.csv"))
}

SECTIONS <- list(cptac = sec_cptac, hpa = sec_hpa)
if (exists("ONLY")) { SECTIONS <- SECTIONS[names(SECTIONS) %in% ONLY]; rm(ONLY) }
run_logged(file.path(RES_DIR, "log_04_cptac.txt"), SECTIONS)
