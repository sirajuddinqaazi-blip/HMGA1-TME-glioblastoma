# =====================================================================
# 11_tam_restricted_signatures.R  (post hoc; Supplementary Table S1)
# The Mueller et al. microglial and blood-derived macrophage signatures
# were derived within sorted TAMs. In bulk tissue, genes also expressed
# by malignant cells can make the scores non-specific. This script keeps
# only genes that are myeloid-restricted in the Neftel single-cell data
# (criterion prespecified in 07_neftel_single_cell.R), re-scores the
# signatures in each bulk cohort (mean of within-cohort z-scores, as in
# Mueller et al.) and repeats the H1 association analysis.
# Requires: results/cohorts_hmga1.rds (01) and T7g_muller_gene_specificity.csv (07)
# =====================================================================
source("scripts/00_config.R", local = TRUE)

run_script11 <- function() {
  obj  <- readRDS(file.path(RES_DIR, "cohorts_hmga1.rds"))
  spec_f <- file.path(RES_DIR, "T7g_muller_gene_specificity.csv")
  if (is.null(obj$muller_expr) || !length(obj$muller_expr) || !file.exists(spec_f)) {
    say("Mueller gene data or specificity table not available -> script 11 skipped"); return(invisible(NULL))
  }
  spec <- fread(spec_f)
  cohorts <- obj$cohorts
  MIN_GENES <- 3

  cat("\n---- Myeloid-restricted Mueller signatures ----\n")
  keep <- list()
  for (st in c("microglia", "macrophage")) {
    g <- spec[set == st & myeloid_restricted %in% TRUE, gene]
    say("%-10s restricted genes (%d): %s", st, length(g), paste(g, collapse = ", "))
    if (length(g) >= MIN_GENES) keep[[st]] <- g else say("  -> fewer than %d genes: restricted %s signature not estimable", MIN_GENES, st)
  }
  if (!length(keep)) return(invisible(NULL))

  for (k in names(cohorts)) {
    m <- obj$muller_expr[[k]]
    for (st in names(keep)) {
      g <- intersect(keep[[st]], rownames(m)); g <- g[apply(m[g, , drop = FALSE], 1, var) > 0]
      cohorts[[k]][[paste0("TAM_", st, "_restricted")]] <- if (length(g) >= MIN_GENES)
        colMeans(t(scale(t(m[g, , drop = FALSE])))) else NA_real_
      say("  [%s] %s restricted score uses %d genes", k, st, length(g))
    }
  }

  outs <- paste0("TAM_", names(keep), "_restricted")
  tab <- assoc_table(cohorts, "log2HMGA1", outs,
                     MODELS[c("M0_unadjusted", "M1_prolif", "M2_prolif_ESTIMATE", "M4_TCGA_ABSOLUTE", "M5_composition")],
                     "H1_TAM_restricted")
  cat("\n---- H1 (post hoc): HMGA1 vs myeloid-restricted TAM-ontogeny signatures ----\n")
  print_assoc(tab[order(tab$outcome, tab$model), ])
  fwrite(tab, file.path(RES_DIR, "T11_tam_restricted.csv"))
  cat("OK T11 saved\n")
}
run_script11()
