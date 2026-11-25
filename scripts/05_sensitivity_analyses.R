# =====================================================================
# 05_sensitivity_analyses.R
# Post hoc sensitivity analyses:
#   1 collinearity  covariate correlations and variance-inflation diagnostics
#                   (bulk cohorts, models M1-M5)
#   2 cptac         CPTAC sensitivity analyses:
#                   (a) excluding the non-hotspot IDH1 R222C tumor (n = 92)
#                   (b) adjustment for proliferation + ESTIMATE immune/stromal scores
#                       computed from the CPTAC RNA data (same model as M2)
#   3 complete_case included vs excluded patients of the multivariable Cox models
# Requires: results/cohorts_hmga1.rds (02_bulk_analyses.R), the CPTAC files used by
#           04_cptac_proteomics.R, share/raw/NIHMS1665743-supplement-2.xlsx, package estimate
# Usage:  source("scripts/05_sensitivity_analyses.R")
# Log: results/log_05_sensitivity_analyses.txt
# =====================================================================
source("scripts/00_functions.R")

sec_collinearity <- function() {
  cohorts <- readRDS(file.path(RES_DIR, "cohorts_hmga1.rds"))$cohorts
  COV <- c("prolif", "ESTIMATE_immune", "ESTIMATE_stromal", "MCP_Monocyte", "sub_MES", "sub_PN", "sub_CL")
  r2 <- function(y, X) { if (!ncol(X)) return(0); summary(lm(y ~ ., data = as.data.frame(X)))$r.squared }

  cat("\n---- Spearman correlations between covariates (per cohort) ----\n")
  cm <- rbindlist(lapply(names(cohorts), function(k) {
    d <- cohorts[[k]]; v <- intersect(COV, names(d))
    m <- cor(d[, v], method = "spearman", use = "pairwise.complete.obs")
    cat("\n", k, "\n"); print(round(m, 2))
    data.table(cohort = k, var1 = rep(v, each = length(v)), var2 = rep(v, length(v)), rho = as.vector(m))[var1 < var2]
  }))
  fwrite(cm, file.path(RES_DIR, "T13a_covariate_correlations.csv"))

  # Diagnostics on the rank scale (the scale of the partial Spearman correlation):
  #  VIF_max         largest VIF among the covariates of the model
  #  HMGA1_R2        share of HMGA1 rank variance explained by the covariates
  #  HMGA1_remaining 1 - HMGA1_R2 = share of HMGA1 variation left for the partial correlation
  cat("\n---- Variance inflation and HMGA1 variance remaining after adjustment ----\n")
  vt <- rbindlist(lapply(names(cohorts), function(k) {
    d0 <- cohorts[[k]]
    rbindlist(lapply(names(MODELS)[-1], function(mn) {
      Z <- MODELS[[mn]]
      if (!all(Z %in% names(d0))) return(NULL)
      d <- na.omit(d0[, c("log2HMGA1", Z)]); if (nrow(d) < 20) return(NULL)
      R <- as.data.frame(lapply(d, rank))
      vif <- if (length(Z) > 1) sapply(Z, function(z) 1 / (1 - r2(R[[z]], R[, setdiff(Z, z), drop = FALSE]))) else c(1)
      h <- r2(R$log2HMGA1, R[, Z, drop = FALSE])
      data.table(cohort = k, model = mn, n = nrow(d), k = length(Z), VIF_max = max(vif),
                 VIF_max_covariate = if (length(Z) > 1) names(which.max(vif)) else Z,
                 HMGA1_R2 = h, HMGA1_remaining = 1 - h)
    }))
  }))
  print(vt, digits = 2)
  say("Rule of thumb: VIF > 5-10 indicates problematic collinearity among covariates.")
  fwrite(vt, file.path(RES_DIR, "T13b_vif_models.csv"))
}

sec_complete_case <- function() {
  cohorts <- readRDS(file.path(RES_DIR, "cohorts_hmga1.rds"))$cohorts
  covs <- list(TCGA = c("age", "KPS", "MGMT"), CGGA_325 = c("age", "MGMT", "radio", "chemo"),
               CGGA_693 = c("age", "MGMT", "radio", "chemo"))
  cat("\n---- Patients included vs excluded from the multivariable Cox models (patients with survival data) ----\n")
  out <- rbindlist(lapply(names(cohorts), function(k) {
    d <- cohorts[[k]]; d <- d[!is.na(d$OS_months) & !is.na(d$event), ]
    inc <- complete.cases(d[, c("log2HMGA1", covs[[k]])])
    miss <- sapply(covs[[k]], function(v) sum(is.na(d[[v]])))
    say("%-9s included %d, excluded %d | missing values per covariate: %s", k, sum(inc), sum(!inc),
        paste(names(miss), miss, sep = " ", collapse = ", "))
    if (sum(!inc) < 3) return(data.table(cohort = k, n_included = sum(inc), n_excluded = sum(!inc)))
    w <- function(v) wilcox.test(d[[v]][inc], d[[v]][!inc], exact = FALSE)$p.value
    mg <- !is.na(d$MGMT)
    lr <- survdiff(Surv(OS_months, event) ~ inc, data = d)
    km <- summary(survfit(Surv(OS_months, event) ~ inc, data = d))$table[, "median"]
    data.table(cohort = k, n_included = sum(inc), n_excluded = sum(!inc),
               age_included = median(d$age[inc], na.rm = TRUE), age_excluded = median(d$age[!inc], na.rm = TRUE), age_p = w("age"),
               HMGA1_included = median(d$log2HMGA1[inc]), HMGA1_excluded = median(d$log2HMGA1[!inc]), HMGA1_p = w("log2HMGA1"),
               MGMT_meth_included = mean(d$MGMT[inc & mg] == "Methylated"), MGMT_meth_excluded = mean(d$MGMT[!inc & mg] == "Methylated"),
               MGMT_known_excluded = sum(!inc & mg),
               median_OS_included = unname(km[2]), median_OS_excluded = unname(km[1]),
               OS_logrank_p = pchisq(lr$chisq, 1, lower.tail = FALSE))
  }), fill = TRUE)
  print(out, digits = 2)
  fwrite(out, file.path(RES_DIR, "T13d_complete_case_comparison.csv"))
}

sec_cptac <- function() {
  CPTAC_DIR <- file.path(RAW_DIR, "cptac_gbm")
  rd <- function(f) { x <- fread(file.path(CPTAC_DIR, f)); m <- as.matrix(x[, -1, with = FALSE])
                      storage.mode(m) <- "numeric"; rownames(m) <- as.character(x[[1]]); m[!duplicated(rownames(m)), ] }
  P <- rd("HS_CPTAC_GBM_proteome_mssm_per_gene.cct")
  R <- rd("HS_CPTAC_GBM_rnaseq_fpkm_uq_log2.cct")

  aa <- as.data.frame(read_excel(file.path(RAW_DIR, "NIHMS1665743-supplement-2.xlsx"), sheet = "additional_annotations"))
  hot <- aa$case[aa$IDH1_alter_status %in% "Missense"]
  r222 <- aa$case[aa$IDH1_alter_status %in% "Missense_not_hotspot"]
  # IDH2: report coding mutations in the MAF (if present)
  mf <- list.files(CPTAC_DIR, pattern = "\\.maf(\\.gz)?$", full.names = TRUE)
  if (length(mf)) {
    maf <- fread(mf[1], select = c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode"))
    idh2 <- maf[Hugo_Symbol == "IDH2" & !Variant_Classification %in% c("Silent", "Intron", "5'Flank", "3'Flank", "5'UTR", "3'UTR")]
    say("IDH2 coding mutations in MAF: %d", nrow(idh2))
  }
  # Methylation annotation: G-CIMP flag vs IDH1 status and methylation classes
  aa$IDH1 <- ifelse(aa$IDH1_alter_status %in% "Missense", "R132", ifelse(aa$IDH1_alter_status %in% "Missense_not_hotspot", "R222C", "wt"))
  cat("\nG-CIMP flag by IDH1 status:\n"); print(table(gcimp = aa$is_gcimp, IDH1 = aa$IDH1, useNA = "ifany"))
  cat("\nMethylation classes of G-CIMP-flagged and IDH1-mutant tumors:\n")
  print(aa[aa$is_gcimp %in% "TRUE" | aa$IDH1 != "wt", c("case", "IDH1", "is_gcimp", "dna_methyl", "dna_methyl_consensus", "dna_methyl_mnp")], row.names = FALSE)

  keep <- setdiff(intersect(colnames(P), colnames(R)), hot)
  say("Tumors with proteome and RNA, IDH1 R132 excluded: %d", length(keep))

  suppressPackageStartupMessages(library(estimate))   # attaches the common_genes data used by filterCommonGenes
  # ESTIMATE from CPTAC RNA (values are already log2 FPKM-UQ: passed as they are)
  dir.create("share/tmp_estimate", recursive = TRUE, showWarnings = FALSE)
  Rk <- R[, keep]; Rk <- Rk[rowSums(is.na(Rk)) == 0, ]
  in_f <- "share/tmp_estimate/CPTAC_input.txt"; cg_f <- "share/tmp_estimate/CPTAC_common.gct"; out_f <- "share/tmp_estimate/CPTAC_estimate.gct"
  write.table(data.frame(NAME = rownames(Rk), Rk, check.names = FALSE), in_f, sep = "\t", quote = FALSE, row.names = FALSE)
  filterCommonGenes(input.f = in_f, output.f = cg_f, id = "GeneSymbol")
  estimateScore(input.ds = cg_f, output.ds = out_f, platform = "illumina")
  sc <- read.table(out_f, skip = 2, header = TRUE, sep = "\t", check.names = FALSE)
  rownames(sc) <- sc$NAME; sc <- t(sc[, -(1:2)])
  est <- data.frame(sample = colnames(Rk), ESTIMATE_stromal = as.numeric(sc[, "StromalScore"]),
                    ESTIMATE_immune = as.numeric(sc[, "ImmuneScore"]))

  get <- function(M, g) { s <- find_symbol(g, rownames(M)); if (is.na(s)) rep(NA_real_, length(keep)) else M[s, keep] }
  pz <- function(M, genes, minfrac) { g <- intersect(genes, rownames(M)); g <- g[rowMeans(!is.na(M[g, keep, drop = FALSE])) >= minfrac]
                                      rowMeans(scale(t(M[g, keep, drop = FALSE])), na.rm = TRUE) }
  TG <- c("CD276", "CCL2", "HAVCR2", "STING1", "CD274")
  d <- data.frame(sample = keep, prot_prolif = pz(P, PROLIF_GENES, 0.8), rna_prolif = pz(R, PROLIF_GENES, 1),
                  est[match(keep, est$sample), c("ESTIMATE_immune", "ESTIMATE_stromal")],
                  p_HMGA1 = get(P, "HMGA1"), r_HMGA1 = get(R, "HMGA1"), check.names = FALSE)
  for (g in TG) { d[[paste0("p_", g)]] <- get(P, g); d[[paste0("r_", g)]] <- get(R, g) }
  say("ESTIMATE immune vs RNA HMGA1 (Spearman): %.2f", cor(d$ESTIMATE_immune, d$r_HMGA1, method = "spearman"))

  sets <- list("All IDH-wt (n incl. R222C)" = rep(TRUE, nrow(d)), "Excluding R222C" = !d$sample %in% r222)
  res <- rbindlist(lapply(names(sets), function(sn) {
    dd <- d[sets[[sn]], ]
    rbindlist(lapply(TG, function(g) rbindlist(list(
      { v <- pspear("r_HMGA1", paste0("r_", g), character(0), dd); data.table(level = "RNA", model = "unadjusted", n = v[["n"]], rho = v[["r"]], p = v[["p"]]) },
      { v <- pspear("r_HMGA1", paste0("r_", g), c("rna_prolif", "ESTIMATE_immune", "ESTIMATE_stromal"), dd)
        data.table(level = "RNA", model = "+ prolif + ESTIMATE", n = v[["n"]], rho = v[["r"]], p = v[["p"]]) },
      if (mean(!is.na(dd[[paste0("p_", g)]])) >= 0.5) {
        v <- pspear("p_HMGA1", paste0("p_", g), character(0), dd); data.table(level = "protein", model = "unadjusted", n = v[["n"]], rho = v[["r"]], p = v[["p"]]) },
      if (mean(!is.na(dd[[paste0("p_", g)]])) >= 0.5) {
        v <- pspear("p_HMGA1", paste0("p_", g), c("prot_prolif", "ESTIMATE_immune", "ESTIMATE_stromal"), dd)
        data.table(level = "protein", model = "+ prolif + ESTIMATE", n = v[["n"]], rho = v[["r"]], p = v[["p"]]) }
    ))[, gene := g]))[, set := sn]
  }))
  setcolorder(res, c("set", "gene", "level", "model"))
  cat("\n---- CPTAC: HMGA1 vs targets, sensitivity analyses ----\n")
  for (i in seq_len(nrow(res))) say("%-26s %-7s %-8s %-20s n = %2d  rho = %5.2f  P = %s",
                                    res$set[i], res$gene[i], res$level[i], res$model[i], res$n[i], res$rho[i], fmt_p(res$p[i]))
  fwrite(res, file.path(RES_DIR, "T13c_cptac_sensitivity.csv"))
}

SECTIONS <- list(collinearity = sec_collinearity, cptac = sec_cptac, complete_case = sec_complete_case)
if (exists("ONLY")) { SECTIONS <- SECTIONS[names(SECTIONS) %in% ONLY]; rm(ONLY) }
run_logged(file.path(RES_DIR, "log_05_sensitivity_analyses.txt"), SECTIONS)
