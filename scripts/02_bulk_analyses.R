# =====================================================================
# 02_bulk_analyses.R
# All bulk RNA-seq analyses (former scripts 01, 02, 03, 04, 05, 08, 09):
#   1 cohorts      three IDH-wildtype primary GBM cohorts (TCGA 160, CGGA_325 74,
#                  CGGA_693 109), signature scores; NEW: CGGA duplicate-sample check
#   2 survival     proliferation coupling, Cox, Kaplan-Meier (Figure 1, Fig. S1);
#                  NEW: full proportional-hazards report, time-varying HMGA1 effect,
#                  numbers-at-risk tables, HKSJ pooled HR
#   3 immune       H4 immune/stromal composition (Figure 2)
#   4 mechanism    H1 CCL2/TAM and H2 STING/IFN (Figure 3); FIXED: no junk TRUE row
#   5 checkpoints  H3 checkpoints, CD276 sensitivity, subtype (Figure 4, Fig. S7)
#   6 all_grades   CGGA grades 2-4 with MKI67 control (post hoc; Fig. S6)
#   7 numbers      cohort characteristics, complete-case n, GDC release
# Every analysis is otherwise identical to the previous scripts, so all earlier
# numbers reproduce. Pooled correlations now also carry HKSJ 95% CIs.
#
# Usage (project folder = the one containing share/):
#   source("scripts/02_bulk_analyses.R")
# To re-run only some sections:  ONLY <- c("survival"); source("scripts/02_bulk_analyses.R")
# Log: results/log_02_bulk.txt
# =====================================================================
source("scripts/00_functions.R")

sec_cohorts <- function() {
  num <- function(x) suppressWarnings(as.numeric(as.character(x)))

  # ------------------------- TCGA -------------------------
  tpm <- readRDS(file.path(DATA_DIR, "tcga_tpm.rds"))
  gi  <- readRDS(file.path(DATA_DIR, "tcga_gene_info.rds"))
  cl  <- readRDS(file.path(DATA_DIR, "tcga_cohort_clinical.rds"))
  nc  <- readRDS(file.path(DATA_DIR, "tcga_normals_clinical.rds"))
  sym <- gi$gene_name[match(rownames(tpm), rownames(gi))]
  tcga_mat <- collapse_symbols(as.matrix(tpm), sym)   # highest-mean Ensembl ID per symbol
  rm(tpm); invisible(gc())

  tcga <- data.frame(sample = rownames(cl), pull_genes(tcga_mat, GENES_NEEDED, rownames(cl), "TCGA"),
                     check.names = FALSE)
  os_days <- ifelse(cl$vital_status == "Dead", cl$days_to_death, cl$days_to_last_follow_up)
  tcga$OS_months <- os_days / 30.44
  tcga$event     <- as.integer(cl$vital_status == "Dead")
  tcga$age       <- num(cl$age_at_index)
  tcga$sex       <- factor(tolower(cl$paper_Gender), levels = c("female", "male"))
  tcga$MGMT      <- factor(cl$paper_MGMT.promoter.status, levels = c("Methylated", "Unmethylated"))
  tcga$KPS       <- num(cl$paper_Karnofsky.Performance.Score)
  tcga$subtype   <- factor(cl$paper_Transcriptome.Subtype, levels = c("CL", "ME", "NE", "PN"),
                           labels = c("Classical", "Mesenchymal", "Neural", "Proneural"))
  tcga$ABSOLUTE_purity <- num(cl$paper_ABSOLUTE.purity)
  tcga$cohort <- "TCGA"

  normals <- data.frame(sample = rownames(nc),
                        HMGA1 = as.numeric(tcga_mat[find_symbol("HMGA1", rownames(tcga_mat)), rownames(nc)]))

  # ------------------------- CGGA -------------------------
  cgga_mats <- list()
  build_cgga <- function(id) {
    ex <- readRDS(file.path(DATA_DIR, paste0("cgga", id, "_expr_gbm.rds")))
    c2 <- readRDS(file.path(DATA_DIR, paste0("cgga", id, "_clinical.rds")))
    for (cc in names(c2)) if (is.character(c2[[cc]])) c2[[cc]] <- trimws(c2[[cc]])
    keep <- c2$PRS_type == "Primary" & c2$Histology == "GBM" & c2$Grade == "WHO IV" &
            c2$IDH_mutation_status %in% "Wildtype"
    ids <- intersect(c2$CGGA_ID[keep], colnames(ex))
    m <- as.matrix(ex[, -1]); storage.mode(m) <- "numeric"
    m <- collapse_symbols(m, as.character(ex[[1]]))
    cgga_mats[[paste0("CGGA_", id)]] <<- m[, ids]
    d <- data.frame(sample = ids, pull_genes(m, GENES_NEEDED, ids, paste0("CGGA_", id)), check.names = FALSE)
    c2 <- c2[match(ids, c2$CGGA_ID), ]
    d$OS_months <- num(c2$OS) / 30.44
    d$event     <- num(c2$Censor)             # CGGA: 1 = dead, 0 = alive
    d$age       <- num(c2$Age)
    d$sex       <- factor(tolower(c2$Gender), levels = c("female", "male"))
    d$MGMT      <- factor(c2$MGMTp_methylation_status, levels = c("methylated", "un-methylated"),
                          labels = c("Methylated", "Unmethylated"))
    d$radio     <- c2$Radio_status; d$chemo <- c2$Chemo_status
    d$cohort    <- paste0("CGGA_", id)
    d
  }
  cohorts <- list(TCGA = tcga, CGGA_325 = build_cgga("325"), CGGA_693 = build_cgga("693"))
  expr_mats <- c(list(TCGA = tcga_mat[, rownames(cl)]), cgga_mats)
  rm(tcga_mat); invisible(gc())

  # ------------------- scores and derived variables -------------------
  score_files <- c(TCGA = "scores_TCGA.csv", CGGA_325 = "scores_CGGA325.csv", CGGA_693 = "scores_CGGA693.csv")
  for (k in names(cohorts)) {
    d  <- cohorts[[k]]
    sc <- read.csv(file.path(DATA_DIR, score_files[[k]]), check.names = FALSE)
    d  <- merge(d, sc, by = "sample", all.x = TRUE, sort = FALSE)
    if (any(is.na(d$ESTIMATE_immune))) stop("Missing ESTIMATE scores in ", k)
    d$prolif       <- rowMeans(scale(log2(d[, PROLIF_GENES] + 1)))   # 10-gene proliferation score
    d$log2HMGA1    <- log2(d$HMGA1 + 1)
    d$HMGA1_group  <- factor(ifelse(d$HMGA1 >= median(d$HMGA1), "High", "Low"), levels = c("Low", "High"))
    for (s in names(SIGNATURES)) {
      ss <- sig_score(d, SIGNATURES[[s]]); d[[s]] <- ss$score
      say("  [%s] %-15s uses %d/%d genes", k, s, ss$n_used, length(SIGNATURES[[s]]))
    }
    cohorts[[k]] <- d[match(cohorts[[k]]$sample, d$sample), ]
  }

  # ------------------- gene-set scores by ssGSEA (full transcriptome) -------------------
  read_sets <- function(f) {   # long (columns 'set','gene' or 'subtype','gene') or wide (one column per set)
    w <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
    key <- intersect(c("set", "subtype"), names(w))
    if (length(key) && "gene" %in% names(w)) return(split(trimws(w$gene), trimws(w[[key[1]]])))
    lapply(w, function(v) { v <- trimws(v); v[!is.na(v) & v != ""] })
  }
  run_ssgsea <- function(gsets, prefix) {
    if (!requireNamespace("GSVA", quietly = TRUE))
      stop("Package GSVA needed: install.packages('BiocManager'); BiocManager::install('GSVA')")
    for (k in names(cohorts)) {
      m  <- log2(expr_mats[[k]][, cohorts[[k]]$sample] + 1)
      m  <- m[apply(m, 1, var) > 0, ]
      gs <- lapply(gsets, intersect, rownames(m))
      if (any(lengths(gs) < 5)) stop(prefix, " gene sets with < 5 genes found in ", k, ": ", paste(names(gs)[lengths(gs) < 5], collapse = ","))
      say("  [%s] %s genes found: %s", k, prefix,
          paste(sprintf("%s %d/%d", names(gsets), lengths(gs), lengths(gsets)), collapse = ", "))
      es <- suppressMessages(GSVA::gsva(GSVA::ssgseaParam(m, gs), verbose = FALSE))
      for (n in names(gsets)) cohorts[[k]][[paste0(prefix, n)]] <<- as.numeric(es[n, ])
    }
  }

  # Wang et al. 2017 subtype signatures (MES, PN, CL) -> sub_MES, sub_PN, sub_CL
  if (file.exists(WANG_FILE)) {
    gsets <- read_sets(WANG_FILE); names(gsets) <- toupper(names(gsets))
    if (!all(c("MES", "PN", "CL") %in% names(gsets))) stop("wang2017 file must contain MES, PN and CL gene sets")
    run_ssgsea(gsets[c("MES", "PN", "CL")], "sub_")
    tt <- cohorts$TCGA
    say("  Check: TCGA mean MES score by published subtype: %s",
        paste(sprintf("%s=%.2f", levels(tt$subtype), tapply(tt$sub_MES, tt$subtype, mean)), collapse = ", "))
  } else say("  WARNING: %s not found -> subtype scores NOT computed; model M3 will be skipped.", WANG_FILE)

  # MSigDB Hallmark interferon-alpha response (Liberzon et al. 2015) -> IFNa_hallmark
  if (requireNamespace("msigdbr", quietly = TRUE)) {
    h <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
                  error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
    ifa <- unique(h$gene_symbol[h$gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE"])
    say("  Hallmark IFN-alpha response: %d genes (msigdbr %s)", length(ifa), as.character(packageVersion("msigdbr")))
    run_ssgsea(list(hallmark = ifa), "IFNa_")
  } else say("  WARNING: package msigdbr not installed -> Hallmark IFN-alpha score skipped (install.packages('msigdbr')).")

  # Mueller et al. 2017 TAM ontogeny signatures -> TAM_microglia, TAM_macrophage
  # Scored as in Mueller et al.: average of z-scored genes of each signature.
  muller_expr <- list()
  if (file.exists(MULLER_FILE)) {
    gsets <- read_sets(MULLER_FILE); names(gsets) <- tolower(names(gsets))
    if (!all(c("microglia", "macrophage") %in% names(gsets))) stop("muller2017 file must contain microglia and macrophage sets")
    for (k in names(cohorts)) {
      m <- log2(expr_mats[[k]][, cohorts[[k]]$sample] + 1)
      muller_expr[[k]] <- m[intersect(unlist(gsets), rownames(m)), , drop = FALSE]   # used by script 11
      for (n in c("microglia", "macrophage")) {
        g <- intersect(gsets[[n]], rownames(m)); g <- g[apply(m[g, , drop = FALSE], 1, var) > 0]
        cohorts[[k]][[paste0("TAM_", n)]] <- colMeans(t(scale(t(m[g, , drop = FALSE]))))
        say("  [%s] Mueller %s signature: %d/%d genes", k, n, length(g), length(gsets[[n]]))
      }
    }
  } else say("  NOTE: %s not found -> microglia/blood-derived TAM split skipped.", MULLER_FILE)


  # ------------------- QC: duplicate patients between and within the CGGA cohorts -------------------
  # Genome-wide expression correlation. Raw correlations between different glioblastomas are high
  # (shared tissue programme), so genes are also centred within each cohort: a duplicated tumor then
  # appears as an isolated outlier pair (centred r far above all other pairs).
  cat("\n---- QC: duplicate samples between CGGA_325 and CGGA_693 ----\n")
  say("Shared CGGA IDs between cohorts: %d", length(intersect(colnames(expr_mats$CGGA_325), colnames(expr_mats$CGGA_693))))
  la <- log2(expr_mats$CGGA_325 + 1); lb <- log2(expr_mats$CGGA_693 + 1)
  g  <- intersect(rownames(la), rownames(lb))
  g  <- g[rowMeans(la[g, ]) > 1 & rowMeans(lb[g, ]) > 1]
  v  <- (apply(la[g, ], 1, var) + apply(lb[g, ], 1, var)) / 2
  g  <- names(sort(v, decreasing = TRUE))[seq_len(min(5000, length(g)))]
  ac <- la[g, ] - rowMeans(la[g, ]); bc <- lb[g, ] - rowMeans(lb[g, ])
  r_raw <- cor(la[g, ], lb[g, ]); r_cen <- cor(ac, bc)
  w_a <- cor(ac); w_b <- cor(bc); diag(w_a) <- NA; diag(w_b) <- NA
  say("Genes used: %d | between cohorts: raw r median %.3f, max %.3f; centred r median %.3f, 99.9th pct %.3f, max %.3f",
      length(g), median(r_raw), max(r_raw), median(r_cen), quantile(r_cen, 0.999), max(r_cen))
  say("Within CGGA_325: max centred r = %.3f | within CGGA_693: max centred r = %.3f",
      max(w_a, na.rm = TRUE), max(w_b, na.rm = TRUE))
  top <- order(r_cen, decreasing = TRUE)[1:5]
  ij  <- arrayInd(top, dim(r_cen))
  cl_a <- cohorts$CGGA_325; cl_b <- cohorts$CGGA_693
  dupt <- data.table(sample_325 = rownames(r_cen)[ij[, 1]], sample_693 = colnames(r_cen)[ij[, 2]],
                     r_centred = r_cen[top], r_raw = r_raw[top],
                     age_325 = cl_a$age[match(rownames(r_cen)[ij[, 1]], cl_a$sample)],
                     age_693 = cl_b$age[match(colnames(r_cen)[ij[, 2]], cl_b$sample)],
                     sex_325 = cl_a$sex[match(rownames(r_cen)[ij[, 1]], cl_a$sample)],
                     sex_693 = cl_b$sex[match(colnames(r_cen)[ij[, 2]], cl_b$sample)],
                     OS_325 = round(cl_a$OS_months[match(rownames(r_cen)[ij[, 1]], cl_a$sample)], 1),
                     OS_693 = round(cl_b$OS_months[match(colnames(r_cen)[ij[, 2]], cl_b$sample)], 1))
  cat("Five most similar cross-cohort pairs (check age/sex/OS for a match):\n"); print(dupt)
  say("Interpretation: a duplicate is an isolated outlier (centred r far above the 99.9th percentile, e.g. > 0.6) with matching age and sex.")
  fwrite(dupt[, comparison := "CGGA_325 vs CGGA_693, top 5 pairs"][], file.path(RES_DIR, "T0_cgga_duplicate_check.csv"))

  saveRDS(list(cohorts = cohorts, tcga_normals = normals, muller_expr = muller_expr), file.path(RES_DIR, "cohorts_hmga1.rds"))

  cat("\n---- Cohort summary ----\n")
  for (k in names(cohorts)) {
    d <- cohorts[[k]]
    say("%-9s n = %3d | survival n = %3d | events = %3d | HMGA1 median = %.2f",
        k, nrow(d), sum(!is.na(d$OS_months) & !is.na(d$event)), sum(d$event, na.rm = TRUE), median(d$HMGA1))
  }
  say("TCGA normal brain n = %d", nrow(normals))
  say("Checkpoint: cohort sizes must be 160 / 74 / 109")
  cat("OK results/cohorts_hmga1.rds saved\n")
}

sec_survival <- function() {
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
    # time-varying HMGA1 effect: HMGA1_z x log(t + 1); P tests whether the HR changes with time
    ft <- tryCatch(coxph(update(f, . ~ . + tt(HMGA1_z)), data = dd, tt = function(x, t, ...) x * log(t + 1)),
                   error = function(e) NULL)
    tvc <- if (!is.null(ft)) coef(ft)[c("HMGA1_z", "tt(HMGA1_z)")] else c(NA, NA)
    tvp <- if (!is.null(ft)) summary(ft)$coefficients["tt(HMGA1_z)", "Pr(>|z|)"] else NA
    data.table(model = label, n = nrow(dd), events = sum(dd$event),
               HR = s$conf.int["HMGA1_z", 1], lo = s$conf.int["HMGA1_z", 3], hi = s$conf.int["HMGA1_z", 4],
               p = s$coefficients["HMGA1_z", 5], logHR = coef(fit)[["HMGA1_z"]],
               se = s$coefficients["HMGA1_z", 3],
               zph_HMGA1_p = if (!is.null(z)) z["HMGA1_z", "p"] else NA,
               zph_global_p = if (!is.null(z)) z["GLOBAL", "p"] else NA,
               tvc_p = tvp, HR_6m = exp(tvc[1] + tvc[2] * log(7)), HR_12m = exp(tvc[1] + tvc[2] * log(13)),
               HR_24m = exp(tvc[1] + tvc[2] * log(25)),
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
    m <- meta_generic(tab$logHR, tab$se)
    data.table(model = tab$model[1], cohort = "Pooled", n = sum(tab$n), events = sum(tab$events),
               HR = exp(m[["est"]]), lo = exp(m[["lo"]]), hi = exp(m[["hi"]]), p = m[["p"]], I2 = m[["I2"]],
               HR_RE = exp(m[["est_RE"]]), lo_HK = exp(m[["lo_HK"]]), hi_HK = exp(m[["hi_HK"]]), p_HK = m[["p_HK"]])
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
  cat("\n---- Proportional hazards: ALL models with a Schoenfeld P < 0.05 ----\n")
  viol <- cox[zph_HMGA1_p < 0.05 | zph_global_p < 0.05]
  for (i in seq_len(nrow(viol))) say("  %-9s %-36s HMGA1-specific P = %s%s | global P = %s%s | time-varying HMGA1 x log(t) P = %s",
    viol$cohort[i], viol$model[i], fmt_p(viol$zph_HMGA1_p[i]), ifelse(viol$zph_HMGA1_p[i] < 0.05, " (violated)", ""),
    fmt_p(viol$zph_global_p[i]), ifelse(viol$zph_global_p[i] < 0.05, " (violated)", ""), fmt_p(viol$tvc_p[i]))
  say("  HMGA1-specific violations: %d of %d models; global violations: %d of %d models",
      sum(cox$zph_HMGA1_p < 0.05, na.rm = TRUE), nrow(cox), sum(cox$zph_global_p < 0.05, na.rm = TRUE), nrow(cox))
  cat("Time-varying HMGA1 HR (per SD) at 6 / 12 / 24 months, models with an HMGA1-specific violation:\n")
  for (i in which(cox$zph_HMGA1_p < 0.05)) say("  %-9s %-36s %.2f / %.2f / %.2f", cox$cohort[i], cox$model[i],
                                                cox$HR_6m[i], cox$HR_12m[i], cox$HR_24m[i])
  pm <- cox_pooled[model == "Multivariable"]
  say("Pooled multivariable HR %.2f (95%% CI %.2f-%.2f; HKSJ %.2f-%.2f): effects above HR %.2f per SD are excluded at the 95%% level",
      pm$HR, pm$lo, pm$hi, pm$lo_HK, pm$hi_HK, pm$hi)
  fwrite(cox_all, file.path(RES_DIR, "T2_HMGA1_cox.csv"))

  # ---------------- D. Kaplan-Meier (median split) ----------------
  km_plot <- function(d, k) {
    d <- d[!is.na(d$OS_months) & !is.na(d$event), ]
    fit <- survfit(Surv(OS_months, event) ~ HMGA1_group, data = d)
    lr  <- survdiff(Surv(OS_months, event) ~ HMGA1_group, data = d)
    p   <- pchisq(lr$chisq, 1, lower.tail = FALSE)
    say("  KM %-9s log-rank P = %s (n = %d)", k, fmt_p(p), nrow(d))
    brk <- if (max(d$OS_months) > 100) 24 else 12
    sv <- ggsurvplot(fit, data = d, palette = unname(COL_LOWHIGH), legend.labs = c("Low", "High"),
                     legend.title = "HMGA1", censor.size = 2, ggtheme = theme_pub, break.time.by = brk,
                     risk.table = TRUE, risk.table.title = "Number at risk", fontsize = 2.8,
                     tables.theme = theme_cleantable())
    g <- sv$plot + annotate("text", x = max(d$OS_months) * 0.65, y = 0.9, size = 3,
                            label = paste0("log-rank P = ", fmt_p(p))) +
      labs(title = k, x = "Overall survival (months)", y = "Survival probability")
    tb <- sv$table + labs(x = NULL, y = NULL) + theme(plot.margin = margin(5, 15, 5, 5))
    g / tb + plot_layout(heights = c(3, 1))
  }
  cat("\n---- D. Kaplan-Meier ----\n")
  kms <- lapply(names(cohorts), function(k) km_plot(cohorts[[k]], k))

  # ---------------- Figure 1 ----------------
  pA <- ggplot(nt, aes(group, HMGA1)) + geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) +
    labs(x = NULL, y = "HMGA1 log2(TPM + 1)", title = "TCGA: normal brain vs IDH-wt GBM") + theme_pub
  bd <- rbindlist(lapply(cohorts, function(d) d[, c("cohort", "log2HMGA1", "prolif")]))
  bd$cohort <- factor(bd$cohort, levels = names(cohorts))
  lab <- coup[variable == "prolif" & cohort != "Pooled"]
  lab$cohort <- factor(lab$cohort, levels = names(cohorts))
  pB <- ggplot(bd, aes(prolif, log2HMGA1)) + geom_point(size = 0.8, alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, colour = "#C0392B", linewidth = 0.6) +
    facet_wrap(~cohort, scales = "free") +
    geom_text(data = lab, aes(x = -Inf, y = Inf, label = sprintf("rho = %.2f", rho)), hjust = -0.1, vjust = 1.3, size = 3) +
    labs(x = "Proliferation score", y = "HMGA1 log2(expr + 1)", title = "a  HMGA1 vs proliferation") + theme_pub
  fd <- cox_all[model %in% c("Univariable", "Multivariable")]
  pC <- forest_dodge(fd, "cohort", "model", "HR", "lo", "hi", c(names(cohorts), "Pooled"),
                     c(Univariable = "grey50", Multivariable = "#C0392B"),
                     "HR per 1 SD HMGA1 (95% CI, log scale)", "b  Cox regression", ref = 1, log_x = TRUE)
  kms[[1]][[1]] <- kms[[1]][[1]] + labs(title = "c  TCGA")
  fig1 <- pB / pC / wrap_plots(kms, nrow = 1) + plot_layout(heights = c(1, 0.8, 1.35))
  save_fig(pA, "FigureS1_TCGA_normal_vs_tumor", 4, 4)   # descriptive only (n = 5 normals)
  save_fig(fig1, "Figure1_expression_survival", 11, 12.5)
  cat("OK Figure 1 saved\n")
}

sec_immune <- function() {
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
                      "a  HMGA1 and immune/stromal composition (pooled, 3 cohorts)", H4_OUT)
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
                     "Partial Spearman r, HMGA1 vs ESTIMATE immune", "b  Per cohort: ESTIMATE immune score") +
    guides(colour = guide_legend(ncol = 1))
  save_fig(pA + pB + plot_layout(widths = c(1.6, 1)), "Figure2_immune_composition_H4", 13, 7)
  cat("OK Figure 2 saved\n")
}

sec_mechanism <- function() {
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
  fwrite(rbindlist(list(t1, t2), fill = TRUE),   # rbind() on data frames would add a row of TRUE
         file.path(RES_DIR, "T4_H1_H2_mechanism_axes.csv"))

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
  pA <- forest_models(t1[t1$model %in% keep, ], "a  H1: CCL2 and TAMs", H1_OUT)
  pB <- forest_models(t2[t2$model %in% keep, ], "b  H2: STING1 and interferon genes", H2_OUT)

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
         title = "c  HMGA1 vs CCL2 after proliferation + ESTIMATE adjustment") + theme_pub
  lims <- unname(MODEL_LABELS[keep])
  # One shared legend: panel a has none; panel b carries all five models and both symbols.
  # Settings are added to each panel with "+" (not "&"), so they cannot override each other.
  pA <- suppressMessages(pA + scale_colour_manual(values = COL_MODELS[lims], limits = lims, drop = FALSE,
                                                  name = NULL, guide = "none") +
                           guides(shape = "none", colour = "none") + theme(legend.position = "none"))
  pB <- suppressMessages(pB + scale_colour_manual(values = COL_MODELS[lims], limits = lims, drop = FALSE, name = NULL) +
                           guides(shape = guide_legend(override.aes = list(colour = "black"), order = 1),
                                  colour = guide_legend(ncol = 3, order = 2)) +
                           theme(legend.position = "bottom", legend.box = "vertical"))
  ab <- (pA + pB) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")   # position only; panel a has no guides, so nothing is re-enabled
  save_fig(ab / pC +
             plot_layout(heights = c(1.3, 1)), "Figure3_CCL2_STING_H1_H2", 13, 10)
  cat("OK Figure 3 saved\n")
}

sec_checkpoints <- function() {
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
  pA <- forest_models(tab[tab$model %in% keep, ], "a  Checkpoint genes (pooled, 3 cohorts)", CHECKPOINTS)
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
  b7d <- per_cohort(as.data.frame(rbindlist(list(tab[tab$outcome == "CD276", ], b7), fill = TRUE)),
                    c("M0_unadjusted", "M2_prolif_ESTIMATE", "M3_plus_subtype", "M5_composition", "M7_vascular_stromal"))
  pB <- forest_dodge(b7d, "cohort", "model", "r", "lo", "hi", c("TCGA", "CGGA_325", "CGGA_693"), COL_MODELS,
                     "Partial Spearman r, HMGA1 vs CD276", "b  HMGA1-CD276 (B7-H3) per cohort") +
    guides(colour = guide_legend(ncol = 1))
  hd <- per_cohort(tab[tab$outcome == "HAVCR2", ], c(keep, "M3b_TCGA_labels"))
  pC <- forest_dodge(hd, "cohort", "model", "r", "lo", "hi", c("TCGA", "CGGA_325", "CGGA_693"), COL_MODELS,
                     "Partial Spearman r, HMGA1 vs HAVCR2", "c  HMGA1-HAVCR2 (TIM-3) per cohort") +
    guides(colour = guide_legend(ncol = 1))
  save_fig(pA + (pB / pC) + plot_layout(widths = c(1.3, 1)), "Figure4_checkpoints_H3", 13, 11)
  pS <- ggplot(tc[!is.na(tc$subtype), ], aes(subtype, log2HMGA1)) +
    geom_boxplot(outlier.shape = NA, width = 0.55) + geom_jitter(width = 0.15, size = 0.7, alpha = 0.6) +
    labs(x = NULL, y = "HMGA1 log2(TPM + 1)", title = sprintf("TCGA IDH-wt GBM: HMGA1 by transcriptional subtype (P = %s)", fmt_p(kw$p.value))) +
    theme_pub
  save_fig(pS, "FigureS7_TCGA_HMGA1_by_subtype", 5.5, 4.5)
  cat("OK Figure 4 and S7 saved\n")
}

sec_all_grades <- function() {
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
    labs(x = NULL, y = "HMGA1 log2(FPKM + 1)", title = "a  HMGA1 by WHO grade (primary gliomas)") + theme_pub
  fd <- rbind(cox_allg[cohort != "Pooled", .(cat = paste(cohort, "- HMGA1"), model = as.character(model), HR, lo, hi)],
              mk[, .(cat = paste(cohort, "- MKI67 (control)"), model, HR, lo, hi)])
  fd$model <- factor(fd$model, levels = names(MODELS_ALLG))
  pB <- forest_dodge(fd, "cat", "model", "HR", "lo", "hi",
                     c("CGGA_325 - HMGA1", "CGGA_325 - MKI67 (control)", "CGGA_693 - HMGA1", "CGGA_693 - MKI67 (control)"),
                     c("Unadjusted" = "grey50", "+ age" = "#3B6FB6", "+ age + grade" = "#E67E22",
                       "+ age + grade + IDH/1p19q" = "#C0392B"),
                     "HR per 1 SD (95% CI, log scale)", "b  Prognostic value across grades, per cohort",
                     ref = 1, log_x = TRUE) + guides(colour = guide_legend(ncol = 1))
  save_fig(pA + pB + plot_layout(widths = c(1.1, 1)), "FigureS6_all_grades_reconciliation", 13, 5.5)
  cat("OK Figure S6 saved\n")
}

sec_numbers <- function() {
  cohorts <- readRDS(file.path(RES_DIR, "cohorts_hmga1.rds"))$cohorts
  out <- file.path(RES_DIR, "T9_manuscript_numbers.txt")
  sink(out, split = TRUE)
  on.exit(sink(), add = TRUE)
  gr <- "Data/processed/GDC_release.txt"
  say("GDC data release used: %s", if (file.exists(gr)) paste(readLines(gr), collapse = " ") else "file Data/processed/GDC_release.txt not found - look up the release on the GDC download date")

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
    cc <- list(TCGA = c("age", "KPS", "MGMT"), CGGA_325 = c("age", "MGMT", "radio", "chemo"), CGGA_693 = c("age", "MGMT", "radio", "chemo"))[[k]]
    say("Complete cases for the multivariable Cox model (%s): %d", paste(cc, collapse = ", "),
        sum(complete.cases(s[, c("log2HMGA1", cc)])))
    say("Median OS, months (95%% CI): %.1f (%.1f-%.1f)", tb["median"], tb["0.95LCL"], tb["0.95UCL"])
    say("Median follow-up, months (reverse KM): %.1f", summary(rev)$table["median"])
    say("HMGA1 expression, median (IQR) [%s]: %s", ifelse(k == "TCGA", "TPM", "FPKM"), q3(d$HMGA1, 1))
  }
  cat("Saved", out, "\n")
}

SECTIONS <- list(cohorts = sec_cohorts, survival = sec_survival, immune = sec_immune, mechanism = sec_mechanism,
                 checkpoints = sec_checkpoints, all_grades = sec_all_grades, numbers = sec_numbers)
if (exists("ONLY")) { SECTIONS <- SECTIONS[names(SECTIONS) %in% ONLY]; rm(ONLY) }
run_logged(file.path(RES_DIR, "log_02_bulk.txt"), SECTIONS)
