# =====================================================================
# 03_immune_composition.R  (Figure 2, H4)
# Does the pan-cancer "immune-excluded" association (Shahzadi 2025)
# hold in IDH-wt GBM once proliferation and purity are accounted for?
# Family H4: ESTIMATE immune/stromal + 10 MCP-counter populations.
# =====================================================================
source("scripts/00_config.R", local = TRUE)
cohorts <- readRDS(file.path(RES_DIR, "cohorts_hmga1.rds"))$cohorts

H4_OUT <- c("ESTIMATE_immune", "ESTIMATE_stromal",
            "MCP_T_cell", "MCP_T_cell_CD8_", "MCP_cytotoxicity_score", "MCP_NK_cell", "MCP_B_cell",
            "MCP_Monocyte", "MCP_Myeloid_dendritic_cell", "MCP_Neutrophil",
            "MCP_Endothelial_cell", "MCP_Cancer_associated_fibroblast")

tab <- assoc_table(cohorts, "log2HMGA1", H4_OUT, MODELS, "H4_immune_composition")
cat("\n---- H4: HMGA1 vs immune/stromal composition (per-cohort r: TCGA CGGA_325 CGGA_693) ----\n")
print_assoc(tab[order(match(tab$outcome, H4_OUT), tab$model), ])
fwrite(tab, file.path(RES_DIR, "T3_H4_immune_composition.csv"))

# Direct test of the pan-cancer claim: raw vs proliferation-adjusted ESTIMATE immune
e0 <- tab[tab$outcome == "ESTIMATE_immune" & tab$model == "M0_unadjusted", ]
e1 <- tab[tab$outcome == "ESTIMATE_immune" & tab$model == "M1_prolif", ]
cat("\n---- H4 summary ----\n")
say("ESTIMATE immune: unadjusted pooled r = %.2f (q = %s); + proliferation r = %.2f (q = %s)",
    e0$r_pooled, fmt_p(e0$q_pooled), e1$r_pooled, fmt_p(e1$q_pooled))
say("Immune exclusion (negative, robust association) supported? %s",
    ifelse(any(tab$outcome == "ESTIMATE_immune" & tab$r_pooled < 0 & tab$robust), "YES", "NO"))

# ---------------- Figure 2 ----------------
pA <- forest_models(tab[tab$model %in% c("M0_unadjusted", "M1_prolif", "M2_prolif_ESTIMATE", "M5_composition"), ],
                    "A  HMGA1 and immune/stromal composition (pooled, 3 cohorts)", H4_OUT)
cd <- rbindlist(lapply(c("M0_unadjusted", "M1_prolif", "M4_TCGA_ABSOLUTE"), function(mn) {
  t <- tab[tab$outcome == "ESTIMATE_immune" & tab$model == mn, ]
  rbindlist(lapply(c("TCGA", "CGGA_325", "CGGA_693"), function(k) {
    r <- t[[paste0("r_", k)]]; n <- t[[paste0("n_", k)]]
    if (is.null(r) || is.na(r)) return(NULL)
    se <- 1 / sqrt(n - 3 - t$k)
    data.table(model = MODEL_LABELS[mn], cohort = k, r = r, lo = tanh(atanh(r) - 1.96 * se), hi = tanh(atanh(r) + 1.96 * se))
  }))
}))
cd$model <- factor(cd$model, levels = MODEL_LABELS)
pB <- forest_dodge(cd, "cohort", "model", "r", "lo", "hi", c("TCGA", "CGGA_325", "CGGA_693"), COL_MODELS,
                   "Partial Spearman r, HMGA1 vs ESTIMATE immune", "B  Per cohort: ESTIMATE immune score") +
  guides(colour = guide_legend(ncol = 1))
save_fig(pA + pB + plot_layout(widths = c(1.6, 1)), "Figure2_immune_composition_H4", 13, 7)
cat("OK Figure 2 saved\n")
