# =====================================================================
# 09_manuscript_numbers.R
# Prints every manuscript number that is not already in run_log.txt:
# Table 1 (cohort characteristics), HMGA1 median (IQR), ABSOLUTE n,
# median OS (95% CI) and median follow-up (reverse Kaplan-Meier).
# Read-only: changes no results. Output: results/T9_manuscript_numbers.txt
# =====================================================================
source("scripts/00_config.R", local = TRUE)
cohorts <- readRDS(file.path(RES_DIR, "cohorts_hmga1.rds"))$cohorts
out <- file.path(RES_DIR, "T9_manuscript_numbers.txt")
sink(out, split = TRUE)

q3 <- function(x, d = 1) { x <- x[!is.na(x)]; q <- quantile(x, c(.5, .25, .75)); sprintf(paste0("%.", d, "f (%.", d, "f-%.", d, "f)"), q[1], q[2], q[3]) }
npct <- function(x, level) { m <- sum(is.na(x)); x <- x[!is.na(x)]
  sprintf("%d (%.0f%%); missing %d", sum(x == level), 100 * mean(x == level), m) }

for (k in names(cohorts)) {
  d <- cohorts[[k]]
  cat("\n=====", k, "=====\n")
  say("Patients, n: %d", nrow(d))
  say("Age, years, median (IQR): %s", q3(d$age, 0))
  say("Male sex: %s", npct(as.character(d$sex), "male"))
  say("MGMT promoter methylated: %s", npct(as.character(d$MGMT), "Methylated"))
  if (k == "TCGA") {
    say("KPS, median (IQR): %s; missing %d", q3(d$KPS, 0), sum(is.na(d$KPS)))
    say("ABSOLUTE purity available, n: %d", sum(!is.na(d$ABSOLUTE_purity)))
    say("Transcriptional subtype available, n: %d", sum(!is.na(d$subtype)))
  } else {
    say("Radiotherapy: %s", npct(as.character(d$radio), "1"))
    say("Chemotherapy: %s", npct(as.character(d$chemo), "1"))
  }
  s <- d[!is.na(d$OS_months) & !is.na(d$event), ]
  fit <- survfit(Surv(OS_months, event) ~ 1, data = s)
  tb <- summary(fit)$table
  rev <- survfit(Surv(OS_months, 1 - event) ~ 1, data = s)
  say("Patients with survival data: %d; deaths: %d", nrow(s), sum(s$event))
  say("Median OS, months (95%% CI): %.1f (%.1f-%.1f)", tb["median"], tb["0.95LCL"], tb["0.95UCL"])
  say("Median follow-up, months (reverse KM): %.1f", summary(rev)$table["median"])
  say("HMGA1 expression, median (IQR) [%s]: %s", ifelse(k == "TCGA", "TPM", "FPKM"), q3(d$HMGA1, 1))
}
sink()
cat("Saved", out, "\n")
