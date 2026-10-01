# =====================================================================
# 04_mechanism_axes.R  (Figure 3, H1 and H2)
# H1: HMGA1 - CCL2 - TAM axis (Chen et al. 2022, HCC)
# H2: HMGA1 - STING1 / interferon axis (He et al. 2025, ESCC; predicted inverse)
# =====================================================================
source("scripts/00_config.R", local = TRUE)
cohorts <- readRDS(file.path(RES_DIR, "cohorts_hmga1.rds"))$cohorts

H1_OUT <- c("CCL2", "MCP_Monocyte", "TAM_microglia", "TAM_macrophage")
H2_OUT <- c("STING1", "CXCL10", "CCL5", "He_ISG", "IFNa_hallmark")

# H2 composition sensitivity: in brain, STING1 is expressed largely by
# myeloid and endothelial cells, so both are adjusted for (prespecified).
MODELS_H2 <- c(MODELS, list(M6_STING_composition =
  c("prolif", "ESTIMATE_immune", "ESTIMATE_stromal", "MCP_Monocyte", "MCP_Endothelial_cell")))

t1 <- assoc_table(cohorts, "log2HMGA1", H1_OUT, MODELS,    "H1_CCL2_TAM")
t2 <- assoc_table(cohorts, "log2HMGA1", H2_OUT, MODELS_H2, "H2_STING_IFN")

cat("\n---- H1: HMGA1 vs CCL2 / TAM (per-cohort r: TCGA CGGA_325 CGGA_693) ----\n")
print_assoc(t1[order(match(t1$outcome, H1_OUT), t1$model), ])
cat("\n---- H2: HMGA1 vs STING / interferon (prediction: NEGATIVE) ----\n")
print_assoc(t2[order(match(t2$outcome, H2_OUT), t2$model), ])
fwrite(rbind(t1, t2, fill = TRUE), file.path(RES_DIR, "T4_H1_H2_mechanism_axes.csv"))

cat("\n---- H1/H2 summary (fully adjusted model M2) ----\n")
for (tb in list(t1, t2)) for (y in unique(tb$outcome)) {
  t <- tb[tb$outcome == y & tb$model == "M2_prolif_ESTIMATE", ]
  if (nrow(t)) say("%-15s M2 pooled r = %5.2f  q = %-7s robust = %s", y, t$r_pooled, fmt_p(t$q_pooled), t$robust)
}
s <- t2[t2$outcome == "STING1" & t2$model == "M6_STING_composition", ]
if (nrow(s)) say("STING1 after myeloid + endothelial adjustment: r = %.2f (q = %s). He et al. predicted r < 0.",
                 s$r_pooled, fmt_p(s$q_pooled))

# ---------------- Figure 3 ----------------
keep <- c("M0_unadjusted", "M1_prolif", "M2_prolif_ESTIMATE", "M5_composition", "M6_STING_composition")
pA <- forest_models(t1[t1$model %in% keep, ], "A  H1: CCL2 and TAMs", H1_OUT)
pB <- forest_models(t2[t2$model %in% keep, ], "B  H2: STING1 and interferon genes", H2_OUT)

# C: CCL2 residual scatter after proliferation + ESTIMATE (ranks)
resid_rank <- function(d, v, Z) resid(lm(rank(d[[v]]) ~ ., data = as.data.frame(lapply(d[, Z], rank))))
sd <- rbindlist(lapply(cohorts, function(d) {
  d <- d[complete.cases(d[, c("log2HMGA1", "CCL2", MODELS$M2_prolif_ESTIMATE)]), ]
  data.table(cohort = d$cohort[1], x = resid_rank(d, "log2HMGA1", MODELS$M2_prolif_ESTIMATE),
             y = resid_rank(d, "CCL2", MODELS$M2_prolif_ESTIMATE))
}))
sd$cohort <- factor(sd$cohort, levels = names(cohorts))
pC <- ggplot(sd, aes(x, y)) + geom_point(size = 0.8, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, colour = "#C0392B", linewidth = 0.6) + facet_wrap(~cohort, scales = "free") +
  labs(x = "HMGA1 rank residual", y = "CCL2 rank residual",
       title = "C  HMGA1 vs CCL2 after proliferation + ESTIMATE adjustment") + theme_pub
lims <- unname(MODEL_LABELS[keep])
ab <- (pA + pB) + plot_layout(guides = "collect") &
  scale_colour_manual(values = COL_MODELS[lims], limits = lims, drop = FALSE, name = NULL) &
  theme(legend.position = "bottom")
save_fig(ab / pC +
           plot_layout(heights = c(1.3, 1)), "Figure3_CCL2_STING_H1_H2", 13, 10)
cat("OK Figure 3 saved\n")
