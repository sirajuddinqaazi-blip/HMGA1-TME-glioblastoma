# =====================================================================
# 01_cohorts.R
# Builds the three IDH-wildtype primary GBM cohorts:
# TCGA (n=160), CGGA_325 (n=74), CGGA_693 (n=109).
# Adds HMGA1 targets, signature scores and (optional) Wang 2017 subtype
# scores by ssGSEA.
# Output: results/cohorts_hmga1.rds
# =====================================================================
source("scripts/00_config.R", local = TRUE)
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
