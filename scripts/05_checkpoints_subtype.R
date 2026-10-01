# =====================================================================
# 05_checkpoints_subtype.R  (Figure 4, H3)
# H3: HMGA1 - TIM-3 (HAVCR2) independent of mesenchymal subtype
# (Ahmady et al. 2025: HAVCR2 highest in mesenchymal GBM).
# =====================================================================
source("scripts/00_config.R", local = TRUE)
cohorts <- readRDS(file.path(RES_DIR, "cohorts_hmga1.rds"))$cohorts
has_sub <- "sub_MES" %in% names(cohorts$TCGA)
if (!has_sub) say("WARNING: no Wang 2017 subtype scores -> M3 skipped; only TCGA published-label sensitivity is run.")

# TCGA sensitivity: published transcriptome subtype labels as dummies (Neural = reference)
tc <- cohorts$TCGA
tc$lab_ME <- as.numeric(tc$subtype == "Mesenchymal")
tc$lab_CL <- as.numeric(tc$subtype == "Classical")
tc$lab_PN <- as.numeric(tc$subtype == "Proneural")
cohorts$TCGA <- tc
MODELS_H3 <- c(MODELS, list(M3b_TCGA_labels = c("prolif", "ESTIMATE_immune", "ESTIMATE_stromal",
                                                "lab_ME", "lab_CL", "lab_PN")))
MODEL_LABELS["M3b_TCGA_labels"] <- "+ prolif + ESTIMATE + subtype labels (TCGA)"
COL_MODELS[MODEL_LABELS[["M3b_TCGA_labels"]]] <- "#B8860B"

tab <- assoc_table(cohorts, "log2HMGA1", CHECKPOINTS, MODELS_H3, "H3_checkpoints")
cat("\n---- H3: HMGA1 vs checkpoint genes (per-cohort r: TCGA CGGA_325 CGGA_693) ----\n")
print_assoc(tab[order(match(tab$outcome, CHECKPOINTS), tab$model), ])
fwrite(tab, file.path(RES_DIR, "T5_H3_checkpoints.csv"))

cat("\n---- HMGA1 and subtype ----\n")
kw <- kruskal.test(log2HMGA1 ~ subtype, data = tc)
say("TCGA HMGA1 by published subtype: %s; Kruskal-Wallis P = %s",
    paste(sprintf("%s %.2f (n=%d)", levels(tc$subtype), tapply(tc$log2HMGA1, tc$subtype, median),
                  as.integer(table(tc$subtype))), collapse = ", "), fmt_p(kw$p.value))
sub_tab <- NULL
if (has_sub) {
  sub_tab <- assoc_table(cohorts, "log2HMGA1", c("sub_MES", "sub_PN", "sub_CL"), MODELS[c("M0_unadjusted", "M1_prolif")],
                         "HMGA1_vs_subtype_scores")
  print_assoc(sub_tab)
  fwrite(sub_tab, file.path(RES_DIR, "T5b_HMGA1_vs_subtype_scores.csv"))
}
# HAVCR2 within mesenchymal vs non-mesenchymal TCGA tumors (small n; exploratory)
for (g in list(c("Mesenchymal"), c("Classical", "Neural", "Proneural"))) {
  dd <- tc[tc$subtype %in% g, ]
  v <- pspear("log2HMGA1", "HAVCR2", "prolif", dd)
  say("TCGA %-30s HMGA1-HAVCR2 (+prolif) r = %.2f, n = %d, P = %s",
      paste(g, collapse = "/"), v[["r"]], v[["n"]], fmt_p(v[["p"]]))
}
h <- tab[tab$outcome == "HAVCR2", c("model", "r_pooled", "q_pooled", "I2", "robust")]
cat("\nHAVCR2 across models:\n"); print(h, row.names = FALSE, digits = 2)

# ---------------- B7-H3 (CD276): vascular / stromal sensitivity ----------------
# CD276 is expressed by malignant AND endothelial cells (single-cell data), so its association is
# additionally adjusted for endothelial, fibroblast-like and monocytic content (post hoc sensitivity).
cat("\n---- CD276 (B7-H3) vascular/stromal sensitivity (post hoc) ----\n")
b7 <- assoc_table(cohorts, "log2HMGA1", "CD276",
                  list(M7_vascular_stromal = c("prolif", "ESTIMATE_immune", "ESTIMATE_stromal", "MCP_Monocyte",
                                               "MCP_Endothelial_cell", "MCP_Cancer_associated_fibroblast")), "CD276_sensitivity")
print_assoc(b7)
fwrite(b7, file.path(RES_DIR, "T5d_CD276_sensitivity.csv"))

# ---------------- Figure 4 ----------------
keep <- c("M0_unadjusted", "M1_prolif", "M2_prolif_ESTIMATE", "M3_plus_subtype", "M5_composition")
pA <- forest_models(tab[tab$model %in% keep, ], "A  Checkpoint genes (pooled, 3 cohorts)", CHECKPOINTS)
per_cohort <- function(tb, models) {
  tb <- tb[tb$model %in% models, ]
  d <- rbindlist(lapply(seq_len(nrow(tb)), function(i) rbindlist(lapply(c("TCGA", "CGGA_325", "CGGA_693"), function(k) {
    r <- tb[[paste0("r_", k)]][i]; n <- tb[[paste0("n_", k)]][i]
    if (is.null(r) || is.na(r)) return(NULL)
    se <- 1 / sqrt(n - 3 - tb$k[i])
    data.table(model = MODEL_LABELS[tb$model[i]], cohort = k, r = r, lo = tanh(atanh(r) - 1.96 * se), hi = tanh(atanh(r) + 1.96 * se))
  }))))
  d$model <- factor(d$model, levels = MODEL_LABELS); d
}
b7d <- per_cohort(rbind(tab[tab$outcome == "CD276", ], b7, fill = TRUE),
                  c("M0_unadjusted", "M2_prolif_ESTIMATE", "M3_plus_subtype", "M5_composition", "M7_vascular_stromal"))
pB <- forest_dodge(b7d, "cohort", "model", "r", "lo", "hi", c("TCGA", "CGGA_325", "CGGA_693"), COL_MODELS,
                   "Partial Spearman r, HMGA1 vs CD276", "B  HMGA1-CD276 (B7-H3) per cohort") +
  guides(colour = guide_legend(ncol = 1))
hd <- per_cohort(tab[tab$outcome == "HAVCR2", ], c(keep, "M3b_TCGA_labels"))
pC <- forest_dodge(hd, "cohort", "model", "r", "lo", "hi", c("TCGA", "CGGA_325", "CGGA_693"), COL_MODELS,
                   "Partial Spearman r, HMGA1 vs HAVCR2", "C  HMGA1-HAVCR2 (TIM-3) per cohort") +
  guides(colour = guide_legend(ncol = 1))
save_fig(pA + (pB / pC) + plot_layout(widths = c(1.3, 1)), "Figure4_checkpoints_H3", 13, 11)
pS <- ggplot(tc[!is.na(tc$subtype), ], aes(subtype, log2HMGA1)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) + geom_jitter(width = 0.15, size = 0.7, alpha = 0.6) +
  labs(x = NULL, y = "HMGA1 log2(TPM + 1)", title = sprintf("TCGA IDH-wt GBM: HMGA1 by transcriptional subtype (P = %s)", fmt_p(kw$p.value))) +
  theme_pub
save_fig(pS, "FigureS7_TCGA_HMGA1_by_subtype", 5.5, 4.5)
cat("OK Figure 4 and S7 saved\n")
