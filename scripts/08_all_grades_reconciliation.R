# =====================================================================
# 08_all_grades_reconciliation.R  (Figure 6)
# Links this study to Din et al. 2025 (J Neurooncol; grades 1-4):
#   Is HMGA1 prognostic ACROSS glioma grades (reproducing 2025) and does
#   that effect disappear once grade and IDH/1p19q status are accounted
#   for (explaining why it is not prognostic within IDH-wt GBM)?
# Data: CGGA mRNAseq_325 and mRNAseq_693, ALL grades, primary tumors.
# All Cox models use the SAME patients (complete cases for every
# covariate), so attenuation cannot come from a changing sample.
# =====================================================================
source("scripts/00_config.R", local = TRUE)
ALLG_DIR <- file.path(RAW_DIR, "cgga_all_grades")

read_cgga_all <- function(id) {
  fe <- list.files(ALLG_DIR, pattern = paste0("mRNAseq_", id, "\\.RSEM-genes.*\\.txt$"), full.names = TRUE)
  fc <- list.files(ALLG_DIR, pattern = paste0("mRNAseq_", id, "_clinical.*\\.txt$"), full.names = TRUE)
  if (length(fe) != 1 || length(fc) != 1) stop("CGGA_", id, ": expression/clinical file not found in ", ALLG_DIR)
  ex <- fread(fe[1])
  hm <- ex[ex[[1]] == "HMGA1"]
  if (nrow(hm) != 1) stop("HMGA1 not found (or duplicated) in ", basename(fe[1]))
  mk <- ex[ex[[1]] == "MKI67"]                                        # positive-control gene
  cl <- as.data.frame(fread(fc[1]))
  names(cl) <- trimws(sub("\\s*\\(.*$", "", names(cl)))              # strip unit text in headers
  for (cc in names(cl)) if (is.character(cl[[cc]])) cl[[cc]] <- trimws(cl[[cc]])
  cod <- names(cl)[grepl("1p19q", names(cl))][1]
  d <- data.frame(sample = cl$CGGA_ID, prs = cl$PRS_type, histology = cl$Histology, grade_txt = cl$Grade,
                  age = suppressWarnings(as.numeric(cl$Age)),
                  OS_months = suppressWarnings(as.numeric(cl$OS)) / 30.44,
                  event = suppressWarnings(as.numeric(cl$Censor)),
                  IDH = cl$IDH_mutation_status, codel = if (!is.na(cod)) cl[[cod]] else NA)
  d$HMGA1 <- as.numeric(unlist(hm[, -1, with = FALSE]))[match(d$sample, names(hm)[-1])]
  d$MKI67 <- if (nrow(mk) == 1) as.numeric(unlist(mk[, -1, with = FALSE]))[match(d$sample, names(mk)[-1])] else NA
  d$grade <- c("WHO II" = 2, "WHO III" = 3, "WHO IV" = 4)[d$grade_txt]
  d$group <- ifelse(d$IDH == "Wildtype", "IDH-wt",
             ifelse(d$IDH == "Mutant" & d$codel == "Codel", "IDH-mut, 1p/19q-codel",
             ifelse(d$IDH == "Mutant" & d$codel == "Non-codel", "IDH-mut, non-codel", NA)))
  d$group <- factor(d$group, levels = c("IDH-mut, 1p/19q-codel", "IDH-mut, non-codel", "IDH-wt"))
  d$cohort <- paste0("CGGA_", id)
  say("CGGA_%s: %d samples in file | PRS types: %s | grades: %s", id, nrow(d),
      paste(names(table(d$prs)), table(d$prs), sep = "=", collapse = " "),
      paste(names(table(d$grade_txt)), table(d$grade_txt), sep = "=", collapse = " "))
  d <- d[d$prs %in% "Primary" & !is.na(d$grade) & !is.na(d$HMGA1), ]
  d$log2HMGA1 <- log2(d$HMGA1 + 1); d$log2MKI67 <- log2(d$MKI67 + 1)
  d
}
all_g <- list(CGGA_325 = read_cgga_all("325"), CGGA_693 = read_cgga_all("693"))

# ---------------- A. HMGA1 by grade and molecular group ----------------
cat("\n---- A. HMGA1 by WHO grade and molecular group (primary gliomas) ----\n")
for (k in names(all_g)) {
  d <- all_g[[k]]
  rg <- cor.test(d$log2HMGA1, d$grade, method = "spearman", exact = FALSE)
  kw <- kruskal.test(log2HMGA1 ~ group, data = d)
  say("%s n = %d | median log2(FPKM+1) by grade: %s | Spearman with grade rho = %.2f, P = %s",
      k, nrow(d), paste(sprintf("G%d %.2f (n=%d)", 2:4, tapply(d$log2HMGA1, factor(d$grade, 2:4), median),
                                as.integer(table(factor(d$grade, 2:4)))), collapse = ", "), rg$estimate, fmt_p(rg$p.value))
  say("        by molecular group: %s | Kruskal-Wallis P = %s",
      paste(sprintf("%s %.2f (n=%d)", levels(d$group), tapply(d$log2HMGA1, d$group, median),
                    as.integer(table(d$group))), collapse = ", "), fmt_p(kw$p.value))
}

# ---------------- B. Cox models on a common patient set ----------------
cat("\n---- B. Cox: HR per 1 SD log2 HMGA1, same patients in every model ----\n")
MODELS_ALLG <- list("Unadjusted" = character(0), "+ age" = "age10",
                    "+ age + grade" = c("age10", "grade_f"),
                    "+ age + grade + IDH/1p19q" = c("age10", "grade_f", "group"))
cox_allg <- rbindlist(lapply(names(all_g), function(k) {
  d <- all_g[[k]]
  d$age10 <- d$age / 10; d$grade_f <- factor(d$grade)
  d <- d[complete.cases(d[, c("OS_months", "event", "log2HMGA1", "age10", "grade_f", "group")]) & d$OS_months > 0, ]
  d$HMGA1_z <- as.numeric(scale(d$log2HMGA1))
  say("%s common analysis set: n = %d, events = %d", k, nrow(d), sum(d$event))
  rbindlist(lapply(names(MODELS_ALLG), function(mn) {
    f <- as.formula(paste("Surv(OS_months, event) ~", paste(c("HMGA1_z", MODELS_ALLG[[mn]]), collapse = " + ")))
    s <- summary(coxph(f, data = d))
    data.table(cohort = k, model = mn, n = nrow(d), events = sum(d$event),
               HR = s$conf.int["HMGA1_z", 1], lo = s$conf.int["HMGA1_z", 3], hi = s$conf.int["HMGA1_z", 4],
               p = s$coefficients["HMGA1_z", 5], logHR = s$coefficients["HMGA1_z", 1], se = s$coefficients["HMGA1_z", 3])
  }))
}))
pool <- rbindlist(lapply(split(cox_allg, cox_allg$model), function(t) {
  w <- 1 / t$se^2; b <- sum(w * t$logHR) / sum(w); se <- sqrt(1 / sum(w)); Q <- sum(w * (t$logHR - b)^2)
  data.table(cohort = "Pooled", model = t$model[1], n = sum(t$n), events = sum(t$events), HR = exp(b),
             lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)),
             I2 = ifelse(Q > 0, max(0, (Q - (nrow(t) - 1)) / Q), 0))
}))
cox_allg <- rbind(cox_allg, pool, fill = TRUE)
cox_allg$model <- factor(cox_allg$model, levels = names(MODELS_ALLG))
setorder(cox_allg, cohort, model)
for (i in seq_len(nrow(cox_allg))) {
  t <- cox_allg[i]
  say("%-9s %-28s n=%4d ev=%4d HR = %.2f [%.2f-%.2f] P = %s%s", t$cohort, t$model, t$n, t$events, t$HR, t$lo, t$hi,
      fmt_p(t$p), ifelse(is.na(t$I2), "", sprintf(" | I2 = %.0f%%", 100 * t$I2)))
}

# ---------------- B2. Positive control: MKI67 in the same patients and models ----------------
# If MKI67 (a grade-linked gene with known prognostic value across gliomas) behaves as expected
# in a cohort, a null HMGA1 result there is informative rather than a data artefact.
cat("\n---- B2. Positive control: MKI67 (same patients, same models) ----\n")
mk_rows <- list()
for (k in names(all_g)) {
  d <- all_g[[k]]; d$age10 <- d$age / 10; d$grade_f <- factor(d$grade)
  d <- d[complete.cases(d[, c("OS_months", "event", "log2HMGA1", "log2MKI67", "age10", "grade_f", "group")]) & d$OS_months > 0, ]
  rc <- cor.test(d$log2HMGA1, d$log2MKI67, method = "spearman", exact = FALSE)
  rg <- cor.test(d$log2MKI67, d$grade, method = "spearman", exact = FALSE)
  say("%s n = %d | HMGA1-MKI67 rho = %.2f | MKI67-grade rho = %.2f", k, nrow(d), rc$estimate, rg$estimate)
  d$MKI67_z <- as.numeric(scale(d$log2MKI67))
  for (mn in names(MODELS_ALLG)) {
    s <- summary(coxph(as.formula(paste("Surv(OS_months, event) ~", paste(c("MKI67_z", MODELS_ALLG[[mn]]), collapse = " + "))), data = d))
    say("   MKI67 %-28s HR = %.2f [%.2f-%.2f] P = %s", mn, s$conf.int["MKI67_z", 1], s$conf.int["MKI67_z", 3],
        s$conf.int["MKI67_z", 4], fmt_p(s$coefficients["MKI67_z", 5]))
    mk_rows[[length(mk_rows) + 1]] <- data.table(cohort = k, model = mn, HR = s$conf.int["MKI67_z", 1],
                                                 lo = s$conf.int["MKI67_z", 3], hi = s$conf.int["MKI67_z", 4],
                                                 p = s$coefficients["MKI67_z", 5])
  }
  s <- summary(coxph(Surv(OS_months, event) ~ grade_f, data = d))
  say("   grade alone: HR WHO III vs II = %.2f, WHO IV vs II = %.2f", s$conf.int[1, 1], s$conf.int[2, 1])
}

# ---------------- C. Within each molecular group ----------------
cat("\n---- C. Within molecular groups (+ age + grade), pooled over both cohorts ----\n")
within <- rbindlist(lapply(levels(all_g[[1]]$group), function(g) {
  per <- rbindlist(lapply(names(all_g), function(k) {
    d <- all_g[[k]]; d$age10 <- d$age / 10; d$grade_f <- factor(d$grade)
    d <- d[d$group %in% g & complete.cases(d[, c("OS_months", "event", "log2HMGA1", "age10", "grade_f")]) & d$OS_months > 0, ]
    if (nrow(d) < 20 || sum(d$event) < 10) return(NULL)
    d$HMGA1_z <- as.numeric(scale(d$log2HMGA1))
    cv <- c("age10", if (length(unique(d$grade_f)) > 1) "grade_f")
    s <- summary(coxph(as.formula(paste("Surv(OS_months, event) ~", paste(c("HMGA1_z", cv), collapse = " + "))), data = d))
    data.table(cohort = k, n = nrow(d), events = sum(d$event), logHR = s$coefficients["HMGA1_z", 1],
               se = s$coefficients["HMGA1_z", 3])
  }))
  if (is.null(per) || nrow(per) == 0) return(NULL)
  w <- 1 / per$se^2; b <- sum(w * per$logHR) / sum(w); se <- sqrt(1 / sum(w)); Q <- sum(w * (per$logHR - b)^2)
  data.table(group = g, n = sum(per$n), events = sum(per$events), HR = exp(b), lo = exp(b - 1.96 * se),
             hi = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)),
             I2 = ifelse(Q > 0 && nrow(per) > 1, max(0, (Q - (nrow(per) - 1)) / Q), NA))
}))
for (i in seq_len(nrow(within))) {
  t <- within[i]
  say("%-24s n=%4d ev=%4d HR = %.2f [%.2f-%.2f] P = %s", t$group, t$n, t$events, t$HR, t$lo, t$hi, fmt_p(t$p))
}
mk <- rbindlist(mk_rows)
fwrite(rbind(cox_allg[, gene := "HMGA1"], mk[, gene := "MKI67 (positive control)"], within, fill = TRUE),
       file.path(RES_DIR, "T8_all_grades_reconciliation.csv"))

# ---------------- Supplementary Figure S6 ----------------
bd <- rbindlist(lapply(all_g, function(d) d[, c("cohort", "grade_txt", "log2HMGA1")]))
bd$grade_txt <- factor(bd$grade_txt, levels = c("WHO II", "WHO III", "WHO IV"))
pA <- ggplot(bd, aes(grade_txt, log2HMGA1)) + geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, size = 0.4, alpha = 0.35) + facet_wrap(~cohort) +
  labs(x = NULL, y = "HMGA1 log2(FPKM + 1)", title = "A  HMGA1 by WHO grade (primary gliomas)") + theme_pub
fd <- rbind(cox_allg[cohort != "Pooled", .(cat = paste(cohort, "- HMGA1"), model = as.character(model), HR, lo, hi)],
            mk[, .(cat = paste(cohort, "- MKI67 (control)"), model, HR, lo, hi)])
fd$model <- factor(fd$model, levels = names(MODELS_ALLG))
pB <- forest_dodge(fd, "cat", "model", "HR", "lo", "hi",
                   c("CGGA_325 - HMGA1", "CGGA_325 - MKI67 (control)", "CGGA_693 - HMGA1", "CGGA_693 - MKI67 (control)"),
                   c("Unadjusted" = "grey50", "+ age" = "#3B6FB6", "+ age + grade" = "#E67E22",
                     "+ age + grade + IDH/1p19q" = "#C0392B"),
                   "HR per 1 SD (95% CI, log scale)", "B  Prognostic value across grades, per cohort",
                   ref = 1, log_x = TRUE) + guides(colour = guide_legend(ncol = 1))
save_fig(pA + pB + plot_layout(widths = c(1.1, 1)), "FigureS6_all_grades_reconciliation", 13, 5.5)
cat("OK Figure S6 saved\n")
