# =====================================================================
# run_all.R  -  runs the full HMGA1 pipeline and writes results/run_log.txt
# Usage in RStudio:
#   setwd("D:/Siraj_HMGA1")
#   source("scripts/run_all.R")
# To re-run only some steps, set STEPS before sourcing, e.g.
#   STEPS <- c("03","04"); source("scripts/run_all.R")
# =====================================================================
if (!dir.exists("share")) stop("Run setwd('D:/Siraj_HMGA1') first.")
dir.create("results", showWarnings = FALSE)

ALL <- c("01_cohorts.R", "02_expression_survival.R", "03_immune_composition.R",
         "04_mechanism_axes.R", "05_checkpoints_subtype.R", "06_single_cell_TISCH2.R",
         "07_neftel_single_cell.R", "08_all_grades_reconciliation.R",
         "11_tam_restricted_signatures.R")
if (!exists("STEPS")) STEPS <- substr(ALL, 1, 2)
run <- ALL[substr(ALL, 1, 2) %in% STEPS]

log_file <- if (length(run) == length(ALL)) "results/run_log.txt" else
  paste0("results/run_log_steps_", paste(STEPS, collapse = "_"), ".txt")
con <- file(log_file, open = "wt")
sink(con, split = TRUE)
cat("HMGA1 pipeline run:", format(Sys.time()), "\n")
cat("Working directory:", getwd(), "\n")

ok <- TRUE
for (s in run) {
  cat("\n\n=====================================================================\n")
  cat("  ", s, "\n")
  cat("=====================================================================\n")
  t0 <- Sys.time()
  res <- tryCatch({
    withCallingHandlers(source(file.path("scripts", s), local = new.env()),
      warning = function(w) { cat("  [warning]", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
    TRUE
  }, error = function(e) { cat("\n*** ERROR in", s, ":", conditionMessage(e), "\n"); FALSE })
  cat(sprintf("  (%s finished in %.1f min)\n", s, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (!res) { ok <- FALSE; break }
}

cat("\n\n---- sessionInfo ----\n"); print(sessionInfo())
cat("\nRun", ifelse(ok, "COMPLETED", "STOPPED WITH ERROR"), "at", format(Sys.time()), "\n")
sink(); close(con)
rm(STEPS)
message("Log written to ", log_file)
