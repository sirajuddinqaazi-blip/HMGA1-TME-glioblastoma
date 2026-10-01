# =====================================================================
# 02_expression_survival.R  (Figure 1)
# HMGA1 in tumor vs normal brain (TCGA), coupling to proliferation,
# and prognostic value (Cox, HR per 1 SD of log2 HMGA1; KM median split).
# =====================================================================
source("scripts/00_config.R", local = TRUE)
obj <- readRDS(file.path(RES_DIR, "cohorts_hmga1.rds")); cohorts <- obj$cohorts

# ---------------- A. tumor vs normal brain (TCGA) ----------------
nt <- rbind(data.frame(group = "Normal brain", HMGA1 = log2(obj$tcga_normals$HMGA1 + 1)),
            data.frame(group = "IDH-wt GBM",   HMGA1 = cohorts$TCGA$log2HMGA1))
nt$group <- factor(nt$group, levels = c("Normal brain", "IDH-wt GBM"))
wt <- wilcox.test(HMGA1 ~ group, data = nt)
cat("\n---- A. HMGA1 tumor vs normal brain (TCGA) ----\n")
say("median log2(TPM+1): normal %.2f (n=%d) vs GBM %.2f (n=%d); Wilcoxon P = %s (normal n is small: descriptive only)",
    median(nt$HMGA1[nt$group == "Normal brain"]), sum(nt$group == "Normal brain"),
    median(nt$HMGA1[nt$group == "IDH-wt GBM"]), sum(nt$group == "IDH-wt GBM"), fmt_p(wt$p.value))

# ---------------- B. coupling to proliferation, age ----------------
cat("\n---- B. HMGA1 vs proliferation / MKI67 / age (Spearman) ----\n")
coup <- rbindlist(lapply(c("prolif", "MKI67", "age"), function(v) {
  per <- rbindlist(lapply(names(cohorts), function(k) {
    x <- pspear("log2HMGA1", v, character(0), cohorts[[k]])
    data.table(variable = v, cohort = k, n = x[["n"]], rho = x[["r"]], p = x[["p"]])
  }))
  m <- meta_r(per$rho, per$n)
  rbind(per, data.table(variable = v, cohort = "Pooled", n = sum(per$n), rho = m[["r"]], p = m[["p"]]), fill = TRUE)[
    , I2 := c(rep(NA, nrow(per)), m[["I2"]])]
}))
for (i in seq_len(nrow(coup))) say("%-7s %-9s n=%3d rho = %5.2f  P = %s%s", coup$variable[i], coup$cohort[i],
                                   coup$n[i], coup$rho[i], fmt_p(coup$p[i]),
                                   ifelse(is.na(coup$I2[i]), "", sprintf("  I2 = %.0f%%", 100 * coup$I2[i])))
fwrite(coup, file.path(RES_DIR, "T1_HMGA1_proliferation_age.csv"))

# ---------------- C. Cox models ----------------
drop_constant <- function(dd, covars) covars[vapply(covars, function(v) length(unique(na.omit(dd[[v]]))) > 1, logical(1))]

cox_fit <- function(d, covars, label, tmax = Inf) {
  d$age10 <- d$age / 10
  if (!is.null(d$radio)) d$radio <- factor(d$radio)
  if (!is.null(d$chemo)) d$chemo <- factor(d$chemo)
  vars <- c("OS_months", "event", "log2HMGA1", covars)
  dd <- d[complete.cases(d[, vars]), ]
  if (is.finite(tmax)) { dd$event[dd$OS_months > tmax] <- 0; dd$OS_months <- pmin(dd$OS_months, tmax) }
  covars <- drop_constant(dd, covars)
  dd$HMGA1_z <- as.numeric(scale(dd$log2HMGA1))                      # HR per 1 SD within analysed patients
  f <- as.formula(paste("Surv(OS_months, event) ~", paste(c("HMGA1_z", covars), collapse = " + ")))
  fit <- coxph(f, data = dd)
  s <- summary(fit)
  z <- tryCatch(cox.zph(fit)$table, error = function(e) NULL)
  data.table(model = label, n = nrow(dd), events = sum(dd$event),
             HR = s$conf.int["HMGA1_z", 1], lo = s$conf.int["HMGA1_z", 3], hi = s$conf.int["HMGA1_z", 4],
             p = s$coefficients["HMGA1_z", 5], logHR = coef(fit)[["HMGA1_z"]],
             se = s$coefficients["HMGA1_z", 3],
             zph_HMGA1_p = if (!is.null(z)) z["HMGA1_z", "p"] else NA,
             zph_global_p = if (!is.null(z)) z["GLOBAL", "p"] else NA,
             covariates = paste(covars, collapse = "+"))
}
covs <- list(TCGA = c("age10", "KPS", "MGMT"),
             CGGA_325 = c("age10", "MGMT", "radio", "chemo"),
             CGGA_693 = c("age10", "MGMT", "radio", "chemo"))

cat("\n---- C. Cox regression, HR per 1 SD log2 HMGA1 ----\n")
for (k in c("CGGA_325", "CGGA_693")) {
  say("  %s radio coding: %s | chemo coding: %s", k,
      paste(names(table(cohorts[[k]]$radio, useNA = "ifany")), table(cohorts[[k]]$radio, useNA = "ifany"), sep = ":", collapse = " "),
      paste(names(table(cohorts[[k]]$chemo, useNA = "ifany")), table(cohorts[[k]]$chemo, useNA = "ifany"), sep = ":", collapse = " "))
}
cox <- rbindlist(lapply(names(cohorts), function(k) {
  d <- cohorts[[k]]
  rbind(cox_fit(d, character(0), "Univariable"),
        cox_fit(d, covs[[k]], "Multivariable"),
        cox_fit(d, covs[[k]], "Multivariable, 24-month truncation", tmax = 24),
        cox_fit(d, c("age10", "MGMT"), "Common (age + MGMT)"))[, cohort := k]
}))
pool_hr <- function(tab) {
  w <- 1 / tab$se^2; b <- sum(w * tab$logHR) / sum(w); se <- sqrt(1 / sum(w))
  Q <- sum(w * (tab$logHR - b)^2); df <- nrow(tab) - 1
  data.table(model = tab$model[1], cohort = "Pooled", n = sum(tab$n), events = sum(tab$events),
             HR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)),
             I2 = ifelse(Q > 0, max(0, (Q - df) / Q), 0))
}
cox_pooled <- rbindlist(lapply(split(cox, cox$model), pool_hr))
cox_all <- rbind(cox, cox_pooled, fill = TRUE)
for (i in seq_len(nrow(cox_all))) {
  t <- cox_all[i]
  say("%-9s %-36s n=%3d ev=%3d HR = %.2f [%.2f-%.2f] P = %-7s%s%s", t$cohort, t$model, t$n, t$events,
      t$HR, t$lo, t$hi, fmt_p(t$p),
      ifelse(is.na(t$zph_HMGA1_p), "", sprintf(" | PH test HMGA1 P = %s, global P = %s", fmt_p(t$zph_HMGA1_p), fmt_p(t$zph_global_p))),
      ifelse(is.na(t$I2), "", sprintf(" | I2 = %.0f%%", 100 * t$I2)))
}
viol <- cox[model == "Multivariable" & (zph_HMGA1_p < 0.05 | zph_global_p < 0.05)]
if (nrow(viol)) say("  NOTE: PH assumption questionable in: %s -> see 24-month truncation sensitivity row",
                    paste(viol$cohort, collapse = ", "))
fwrite(cox_all, file.path(RES_DIR, "T2_HMGA1_cox.csv"))

# ---------------- D. Kaplan-Meier (median split) ----------------
km_plot <- function(d, k) {
  d <- d[!is.na(d$OS_months) & !is.na(d$event), ]
  fit <- survfit(Surv(OS_months, event) ~ HMGA1_group, data = d)
  lr  <- survdiff(Surv(OS_months, event) ~ HMGA1_group, data = d)
  p   <- pchisq(lr$chisq, 1, lower.tail = FALSE)
  say("  KM %-9s log-rank P = %s (n = %d)", k, fmt_p(p), nrow(d))
  g <- ggsurvplot(fit, data = d, palette = unname(COL_LOWHIGH), legend.labs = c("Low", "High"),
                  legend.title = "HMGA1", censor.size = 2, ggtheme = theme_pub)$plot
  g + annotate("text", x = max(d$OS_months) * 0.65, y = 0.9, size = 3,
               label = paste0("log-rank P = ", fmt_p(p))) +
    labs(title = k, x = "Overall survival (months)", y = "Survival probability")
}
cat("\n---- D. Kaplan-Meier ----\n")
kms <- lapply(names(cohorts), function(k) km_plot(cohorts[[k]], k))

# ---------------- Figure 1 ----------------
pA <- ggplot(nt, aes(group, HMGA1)) + geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) +
  labs(x = NULL, y = "HMGA1 log2(TPM + 1)", title = "A  TCGA: normal brain vs IDH-wt GBM") + theme_pub
bd <- rbindlist(lapply(cohorts, function(d) d[, c("cohort", "log2HMGA1", "prolif")]))
bd$cohort <- factor(bd$cohort, levels = names(cohorts))
lab <- coup[variable == "prolif" & cohort != "Pooled"]
lab$cohort <- factor(lab$cohort, levels = names(cohorts))
pB <- ggplot(bd, aes(prolif, log2HMGA1)) + geom_point(size = 0.8, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, colour = "#C0392B", linewidth = 0.6) +
  facet_wrap(~cohort, scales = "free") +
  geom_text(data = lab, aes(x = -Inf, y = Inf, label = sprintf("rho = %.2f", rho)), hjust = -0.1, vjust = 1.3, size = 3) +
  labs(x = "Proliferation score", y = "HMGA1 log2(expr + 1)", title = "A  HMGA1 vs proliferation") + theme_pub
fd <- cox_all[model %in% c("Univariable", "Multivariable")]
pC <- forest_dodge(fd, "cohort", "model", "HR", "lo", "hi", c(names(cohorts), "Pooled"),
                   c(Univariable = "grey50", Multivariable = "#C0392B"),
                   "HR per 1 SD HMGA1 (95% CI, log scale)", "B  Cox regression", ref = 1, log_x = TRUE)
kms[[1]] <- kms[[1]] + labs(title = "C  TCGA")
fig1 <- pB / pC / wrap_plots(kms, nrow = 1) + plot_layout(heights = c(1, 0.8, 1))
save_fig(pA, "FigureS1_TCGA_normal_vs_tumor", 4, 4)   # descriptive only (n = 5 normals)
save_fig(fig1, "Figure1_expression_survival", 11, 11)
cat("OK Figure 1 saved\n")
