# =====================================================================
# 03_single_cell_and_tables.R
# Single-cell analyses, post hoc TAM signatures and the supplementary workbook
# (former scripts 06, 07, 11, 10):
#   1 tisch2         TISCH2 summaries of the SAME Neftel dataset (Fig. S2; not independent)
#   2 neftel         Neftel Smart-seq2 adult tumors (Figure 5, Figs. S3, S4, S5);
#                    NEW: chr7/chr10 copy-number check of malignant cells with sensitivity
#                    analysis (Fig. S8); FIXED: Fig. 5a medians over the paired tumors only
#   3 tam_restricted myeloid-restricted Mueller signatures (post hoc);
#                    NEW: bootstrap test of macrophage vs microglial difference
#   4 tables         Supplementary Tables workbook (new sheets S3h, S3i, S3j, QC)
# Requires results/cohorts_hmga1.rds from 02_bulk_analyses.R.
# First run after this update re-reads the Neftel TPM file once (chr7/chr10 genes
# added to the cache; several minutes).
#
# Usage:  source("scripts/03_single_cell_and_tables.R")
# To re-run only some sections:  ONLY <- c("tables"); source("scripts/03_single_cell_and_tables.R")
# Log: results/log_03_single_cell.txt
# =====================================================================
source("scripts/00_functions.R")

sec_tisch2 <- function() {
  TISCH_GENES <- c("HMGA1", "HAVCR2", "CCL2", "CD274")
  CONFIRMED   <- c("CD8Tex", "Mono/Macro", "AC-like", "MES-like", "NPC-like", "OPC-like", "Oligodendrocyte")

  read_tisch <- function(g) {
    f <- file.path(DATA_DIR, paste0("TISCH_", g, "_heatmap.csv"))
    if (!file.exists(f)) { say("  missing %s", f); return(NULL) }
    x <- read.csv(f, fileEncoding = "UTF-8-BOM", check.names = FALSE, stringsAsFactors = FALSE)
    x <- x[, 1:3]; names(x) <- c("celltype", "dataset", "value")   # "Category", "... (y)", "... (value)"
    x$gene <- g; x$value <- as.numeric(x$value); x
  }
  ti <- rbindlist(lapply(TISCH_GENES, read_tisch))

  cat("\n---- TISCH2 datasets containing GSE131928 ----\n")
  print(unique(ti$dataset[grepl("131928", ti$dataset)]))
  ss <- ti[grepl("131928", dataset) & grepl("smart", dataset, ignore.case = TRUE)]
  if (nrow(ss) == 0) stop("No GSE131928 Smart-seq2 dataset found; check dataset names printed above.")
  cat("Cell types listed for the Smart-seq2 dataset:\n"); print(unique(ss$celltype))

  ss$ct <- NA_character_
  for (c0 in CONFIRMED) ss$ct[grepl(c0, ss$celltype, fixed = TRUE)] <- c0
  ss <- ss[!is.na(ct)]
  found <- unique(ss$ct)
  if (length(found) < length(CONFIRMED))
    say("  WARNING: confirmed cell types not matched: %s (edit CONFIRMED to the exact names printed above)",
        paste(setdiff(CONFIRMED, found), collapse = ", "))

  wide <- dcast(ss, gene ~ ct, value.var = "value", fun.aggregate = mean)
  cat("\nMean log(TPM/10+1), Neftel Smart-seq2, confirmed cell types only:\n")
  print(wide, digits = 2)
  fwrite(wide, file.path(RES_DIR, "T6_TISCH2_Neftel_Smartseq2.csv"))

  ss$ct   <- factor(ss$ct, levels = CONFIRMED)
  ss$gene <- factor(ss$gene, levels = rev(TISCH_GENES))
  p <- ggplot(ss, aes(ct, gene, fill = value)) + geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
    scale_fill_gradient(low = "white", high = "#C0392B", name = "log(TPM/10+1)") +
    labs(x = NULL, y = NULL, title = "TISCH2: Neftel et al. 2019 Smart-seq2 (GSE131928), confirmed cell types") +
    theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_fig(p, "FigureS2_TISCH2_heatmap", 8, 3.5)
  cat("OK Figure S2 saved (TISCH2 summaries of the SAME Neftel dataset; secondary reference, not independent validation)\n")
}

sec_neftel <- function() {
  GZ    <- file.path(RAW_DIR, "GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv.gz")
  TSV   <- sub("\\.gz$", "", GZ)
  XLSX  <- file.path(RAW_DIR, "GSE131928_single_cells_tumor_name_and_adult_or_peidatric.xlsx")
  CACHE <- file.path(RES_DIR, "neftel_subset.rds")

  # Cell classification follows Neftel et al. 2019 (STAR Methods): mean log2(TPM/10+1) of published
  # marker sets, "classified to each of these cell types by scores above 4".
  THR      <- 4          # Neftel et al. 2019 cutoff
  THR_SENS <- c(3, 5)    # sensitivity thresholds
  MIN_CELLS_CLASS <- 10  # min cells of a class per tumor for comparisons
  MIN_MAL  <- 50         # min malignant cells per tumor for within-malignant correlations
  MIN_DET  <- 0.10       # both genes must be detected in >= 10% of the tumor's malignant cells

  MARKERS <- list(       # Neftel et al. 2019 marker sets (+ endothelial, added here)
    myeloid     = c("CD14", "AIF1", "FCER1G", "FCGR3A", "TYROBP", "CSF1R"),
    tcell       = c("CD2", "CD3D", "CD3E", "CD3G"),
    oligo       = c("MBP", "TF", "PLP1", "MAG", "MOG", "CLDN11"),
    endothelial = c("PECAM1", "VWF", "CLDN5")
  )
  TARGETS <- c("HMGA1", "CCL2", "STING1", "HAVCR2", "CD274", "CD276")

  neftel_mod <- NULL
  if (file.exists(NEFTEL_FILE)) {
    w <- read.csv(NEFTEL_FILE, check.names = FALSE, stringsAsFactors = FALSE)
    neftel_mod <- if (all(c("set", "gene") %in% names(w))) split(trimws(w$gene), trimws(w$set)) else
      lapply(w, function(v) { v <- trimws(v); v[!is.na(v) & v != ""] })
  }

  # ---------------- metadata: adult Smart-seq2 cells ----------------
  x0 <- read_excel(XLSX, col_names = FALSE, .name_repair = "minimal")
  h  <- which(x0[[1]] == "Sample name")[1]
  meta <- as.data.frame(x0[(h + 1):nrow(x0), ]); names(meta) <- as.character(unlist(x0[h, ]))
  meta <- meta[!is.na(meta[["Sample name"]]), ]
  meta <- data.frame(cell = meta[["Sample name"]], file = meta[["processed data file"]],
                     tumor = meta[["tumour name"]], age_group = meta[["adult/pediatric"]])
  cat("\n---- Metadata ----\n")
  print(table(platform = ifelse(grepl("Smartseq2", meta$file), "Smart-seq2", "10X"), meta$age_group))
  meta <- meta[grepl("Smartseq2", meta$file) & meta$age_group == "adult", ]
  say("Adult Smart-seq2: %d cells from %d tumors", nrow(meta), length(unique(meta$tumor)))

  # ---------------- expression subset (cached) ----------------
  muller_sets <- NULL
  if (file.exists(MULLER_FILE)) { w <- read.csv(MULLER_FILE, stringsAsFactors = FALSE); muller_sets <- split(trimws(w$gene), trimws(w$set)) }
  # chromosome 7 and 10 genes for the copy-number check (malignant cells of IDH-wt GBM typically
  # carry chr7 gain / chr10 loss); gene-to-chromosome map from org.Hs.eg.db
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db", ask = FALSE, update = FALSE)
  s2c <- merge(AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egSYMBOL), AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egCHR), by = "gene_id")
  s2c <- s2c[!s2c$symbol %in% s2c$symbol[duplicated(s2c$symbol) & !duplicated(s2c[, c("symbol", "chromosome")])], ]
  CHR_GENES <- list(chr7 = unique(s2c$symbol[s2c$chromosome == "7"]), chr10 = unique(s2c$symbol[s2c$chromosome == "10"]))
  want <- unique(c(TARGETS, unlist(MARKERS), PROLIF_GENES, unlist(neftel_mod), unlist(muller_sets), unlist(CHR_GENES)))
  rebuild <- !file.exists(CACHE) || !all(want %in% readRDS(CACHE)$want)
  if (rebuild) {
    if (!file.exists(TSV)) {
      if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
      say("Decompressing %s (one time only) ...", basename(GZ))
      R.utils::gunzip(GZ, destname = TSV, remove = FALSE)
    }
    genes_all <- fread(TSV, select = 1)[[1]]
    cols_all  <- names(fread(TSV, nrows = 0))
    sym  <- setNames(sapply(want, find_symbol, available = genes_all), want)
    if (any(is.na(sym[TARGETS]))) stop("Target genes missing in TPM file: ", paste(TARGETS[is.na(sym[TARGETS])], collapse = ","))
    sym <- sym[!is.na(sym)]
    idx <- match(sym, genes_all)
    cells <- intersect(meta$cell, cols_all)
    say("Cells matched between Excel and TPM file: %d of %d", length(cells), nrow(meta))
    chunks <- split(cells, ceiling(seq_along(cells) / 1000))
    E <- NULL; ngenes <- NULL
    for (i in seq_along(chunks)) {
      say("  reading chunk %d/%d", i, length(chunks))
      dt <- fread(TSV, select = c(cols_all[1], chunks[[i]]), showProgress = FALSE)
      ngenes <- c(ngenes, vapply(dt[, -1], function(v) sum(v > 0), 1))
      m <- as.matrix(dt[idx, -1]); rownames(m) <- names(sym)
      E <- cbind(E, m); rm(dt); invisible(gc())
    }
    saveRDS(list(E = E, ngenes = ngenes, want = want), CACHE)
    say("Cached %d genes x %d cells in %s (the .tsv can now be deleted to save disk space)", nrow(E), ncol(E), CACHE)
  }
  cc <- readRDS(CACHE); E <- cc$E; ngenes <- cc$ngenes
  meta <- meta[match(colnames(E), meta$cell), ]

  mx <- max(E)
  say("Max value in matrix = %.1f; %s", mx,
      ifelse(mx > 100, "treated as TPM -> log2(TPM/10 + 1)",
             "file is already log2(TPM/10 + 1) (Neftel processed format) -> used as is"))
  L <- if (mx > 100) log2(E / 10 + 1) else E

  # ---------------- cell classification ----------------
  msc <- function(g) { g <- intersect(g, rownames(L)); colMeans(L[g, , drop = FALSE]) }
  S <- data.frame(myeloid = msc(MARKERS$myeloid), tcell = msc(MARKERS$tcell),
                  oligo = msc(MARKERS$oligo), endo = msc(MARKERS$endothelial))
  classify <- function(S, thr) {
    cls <- ifelse(S$tcell >= thr & S$tcell > S$myeloid, "T cell",
           ifelse(S$myeloid >= thr, "Myeloid",
           ifelse(S$oligo >= thr, "Oligodendrocyte",
           ifelse(S$endo >= thr, "Endothelial", "Malignant"))))
    factor(cls, levels = c("Malignant", "Myeloid", "T cell", "Oligodendrocyte", "Endothelial"))
  }
  meta$class <- classify(S, THR)
  cat("\n---- Cell classes (Neftel et al. 2019 markers, cutoff 4) ----\n")
  ct <- table(meta$tumor, meta$class); print(ct); print(colSums(ct))
  fwrite(as.data.frame.matrix(ct), file.path(RES_DIR, "T7a_neftel_cell_counts.csv"), row.names = TRUE)
  say("NOTE: immune proportions reflect CD45 FACS sorting (MGH143 immune-depleted), not tumor composition.")

  prolif_cell <- colMeans(t(scale(t(L[intersect(PROLIF_GENES, rownames(L)), ]))), na.rm = TRUE)
  meta$cycling <- prolif_cell; meta$ngenes <- ngenes[meta$cell]
  meta$HMGA1 <- L["HMGA1", ]

  # ---------------- (a) HMGA1 by cell class, paired within tumor ----------------
  by_class <- function(cls_vec, label) {
    d <- data.table(tumor = meta$tumor, class = cls_vec, HMGA1 = meta$HMGA1, det = meta$HMGA1 > 0,
                    tpm = 10 * (2^meta$HMGA1 - 1))
    pt <- d[, .(n = .N, mean_log = mean(HMGA1), det = mean(det), mean_tpm = mean(tpm)), by = .(tumor, class)]
    pt <- pt[n >= MIN_CELLS_CLASS]
    out <- NULL
    for (cmp in c("Myeloid", "T cell")) {
      w <- dcast(pt[class %in% c("Malignant", cmp)], tumor ~ class, value.var = c("mean_log", "mean_tpm"))
      a <- w[[paste0("mean_log_Malignant")]]; b <- w[[paste0("mean_log_", cmp)]]
      ok <- !is.na(a) & !is.na(b)
      if (sum(ok) < 3) { say("  [%s] Malignant vs %s: fewer than 3 tumors with both classes", label, cmp); next }
      wt <- wilcox.test(a[ok], b[ok], paired = TRUE, exact = FALSE)
      fc <- w[[paste0("mean_tpm_Malignant")]][ok] / w[[paste0("mean_tpm_", cmp)]][ok]
      say("  [%s] Malignant vs %-7s: %2d tumors | median of per-tumor mean log2(TPM/10+1) %.2f vs %.2f | median TPM ratio %.2f | paired Wilcoxon P = %s",
          label, cmp, sum(ok), median(a[ok]), median(b[ok]), median(fc), fmt_p(wt$p.value))
      out <- rbind(out, data.table(threshold = label, comparison = paste("Malignant vs", cmp), n_tumors = sum(ok),
                                   median_malignant = median(a[ok]), median_other = median(b[ok]),
                                   median_TPM_ratio = median(fc), p = wt$p.value))
    }
    list(per_tumor = pt, test = out)
  }
  cat("\n---- (a) HMGA1 by cell class (per-tumor means, paired) ----\n")
  A <- by_class(meta$class, "thr=4")
  sens <- lapply(THR_SENS, function(t) by_class(classify(S, t), paste0("thr=", t)))
  fwrite(rbind(A$test, rbindlist(lapply(sens, `[[`, "test"))), file.path(RES_DIR, "T7b_neftel_HMGA1_by_class.csv"))
  fwrite(A$per_tumor, file.path(RES_DIR, "T7c_neftel_HMGA1_per_tumor_class.csv"))
  cat("\nPer-class summary over tumors (thr=4):\n")
  print(A$per_tumor[, .(tumors = .N, median_mean_log = median(mean_log), median_detection = median(det)), by = class])

  # ---------------- (b) within-malignant correlations ----------------
  cat("\n---- (b) Within malignant cells: HMGA1 vs targets (per tumor, then Fisher-z pooled) ----\n")
  mal <- meta$class == "Malignant"
  make_mdat <- function(mask) {
    Md <- data.frame(tumor = meta$tumor[mask], HMGA1 = L["HMGA1", mask], cycling = meta$cycling[mask],
                     ngenes = meta$ngenes[mask], t(L[setdiff(TARGETS, "HMGA1"), mask, drop = FALSE]), check.names = FALSE)
    if (!is.null(neftel_mod)) {
      z <- t(scale(t(L[, mask])))
      for (m in names(neftel_mod)) {
        g <- setdiff(intersect(neftel_mod[[m]], rownames(z)), "HMGA1");   # HMGA1 excluded: avoids circularity
        if (length(g) < 5) { say("  Neftel module %-5s: only %d genes found -> skipped", m, length(g)); next }
        Md[[paste0("mod_", m)]] <- colMeans(z[g, , drop = FALSE], na.rm = TRUE)
        say("  Neftel module %-5s genes found: %d/%d", m, length(g), length(neftel_mod[[m]]))
      }
      mes <- grep("^mod_MES", names(Md), value = TRUE)
      if (length(mes)) Md$MES_like <- rowMeans(Md[, mes, drop = FALSE])
    }
    Md
  }
  within_mal <- function(Md) {
    targets_b <- c(setdiff(TARGETS, "HMGA1"), if ("MES_like" %in% names(Md)) "MES_like")
    wm <- rbindlist(lapply(targets_b, function(y) {
      rbindlist(lapply(list(c(), c("cycling", "ngenes")), function(Z) {
        per <- rbindlist(lapply(split(Md, Md$tumor), function(d) {
          if (nrow(d) < MIN_MAL) return(NULL)
          if (y != "MES_like" && (mean(d$HMGA1 > 0) < MIN_DET || mean(d[[y]] > 0) < MIN_DET)) return(NULL)
          v <- pspear("HMGA1", y, Z, d); data.table(tumor = d$tumor[1], n = v[["n"]], r = v[["r"]], p = v[["p"]])
        }))
        if (nrow(per) < 3) return(data.table(target = y, adjustment = ifelse(length(Z), "cycle + complexity", "none"),
                                             n_tumors = nrow(per)))
        m <- meta_r(per$r, per$n, length(Z))
        data.table(target = y, adjustment = ifelse(length(Z), "cycle + complexity", "none"), n_tumors = nrow(per),
                   n_cells = sum(per$n), r_pooled = m[["r"]], lo = m[["lo"]], hi = m[["hi"]], p = m[["p"]], I2 = m[["I2"]],
                   median_r = median(per$r), frac_positive = mean(per$r > 0))
      }), fill = TRUE)
    }), fill = TRUE)
    for (cn in c("n_cells", "r_pooled", "lo", "hi", "p", "I2", "median_r", "frac_positive"))
      if (!cn %in% names(wm)) wm[, (cn) := NA_real_]
    wm
  }
  Mdat <- make_mdat(mal)
  targets_b <- c(setdiff(TARGETS, "HMGA1"), if ("MES_like" %in% names(Mdat)) "MES_like")
  wm <- within_mal(Mdat)
  for (i in seq_len(nrow(wm))) {
    t <- wm[i]
    if (is.na(t$r_pooled)) { say("%-9s %-19s only %d tumors pass detection filters -> not estimable", t$target, t$adjustment, t$n_tumors); next }
    say("%-9s %-19s tumors = %2d cells = %4d | pooled r = %5.2f [%5.2f, %5.2f] P = %-7s I2 = %3.0f%% | median r = %5.2f, positive in %3.0f%% of tumors",
        t$target, t$adjustment, t$n_tumors, t$n_cells, t$r_pooled, t$lo, t$hi, fmt_p(t$p), 100 * t$I2, t$median_r, 100 * t$frac_positive)
  }
  say("NOTE: single-cell P values are nominal (not corrected for multiple testing); exploratory.")
  fwrite(wm, file.path(RES_DIR, "T7d_neftel_within_malignant.csv"))

  # ---------------- (a2) where are the target genes expressed? ----------------
  cat("\n---- (a2) Target genes by cell class (all adult cells) ----\n")
  tg <- rbindlist(lapply(TARGETS, function(g)
    data.table(gene = g, class = meta$class, v = L[g, ])[, .(cells = .N, mean = mean(v), detected = mean(v > 0)), by = .(gene, class)]))
  cat("Mean log2(TPM/10+1):\n");   print(dcast(tg, gene ~ class, value.var = "mean"), digits = 2)
  cat("Fraction of cells with expression > 0:\n"); print(dcast(tg, gene ~ class, value.var = "detected"), digits = 2)
  fwrite(tg, file.path(RES_DIR, "T7e_neftel_targets_by_class.csv"))

  # ---------------- (a3) cell-type specificity of the Mueller TAM signature genes ----------------
  # Mueller et al. derived their signatures WITHIN sorted TAMs; in bulk tissue some genes are also
  # expressed by malignant or other cells. Criterion for a myeloid-restricted gene:
  # mean log2(TPM/10+1) in myeloid cells >= 1 unit above malignant cells AND detection in myeloid
  # cells >= 2x detection in malignant cells. (Defined after the primary analysis: POST HOC.)
  if (!is.null(muller_sets)) {
    cat("\n---- (a3) Mueller signature genes: myeloid vs malignant expression ----\n")
    spec <- rbindlist(lapply(names(muller_sets), function(st) rbindlist(lapply(muller_sets[[st]], function(g) {
      if (!g %in% rownames(L)) return(data.table(set = st, gene = g, in_data = FALSE))
      v <- L[g, ]; my <- meta$class == "Myeloid"; ma <- meta$class == "Malignant"
      data.table(set = st, gene = g, in_data = TRUE, mean_myeloid = mean(v[my]), mean_malignant = mean(v[ma]),
                 det_myeloid = mean(v[my] > 0), det_malignant = mean(v[ma] > 0))
    }))), fill = TRUE)
    spec[, myeloid_restricted := in_data & (mean_myeloid - mean_malignant >= 1) & (det_myeloid >= 2 * det_malignant)]
    for (st in names(muller_sets)) {
      x <- spec[set == st]
      say("  %-10s %2d/%2d genes myeloid-restricted; not restricted: %s", st, sum(x$myeloid_restricted, na.rm = TRUE),
          nrow(x), paste(x$gene[!x$myeloid_restricted %in% TRUE], collapse = ", "))
    }
    fwrite(spec, file.path(RES_DIR, "T7g_muller_gene_specificity.csv"))
  }


  # ---------------- (c) copy-number check of the "malignant" class (sensitivity) ----------------
  # Malignant cells were defined by exclusion. Most IDH-wt GBMs carry chr7 gain and chr10 loss, so an
  # expression-based chr7-minus-chr10 score (genes centred on non-malignant reference cells, as in
  # inferCNV, without positional smoothing) is used to check the classification. Reference cells:
  # myeloid, T cells and oligodendrocytes. Cells above the 95th percentile of the reference
  # distribution are called CNV-supported. Tumors without +7/-10 will show few CNV-supported cells.
  cat("\n---- (c) chr7 gain / chr10 loss score: check of the malignant classification ----\n")
  ref <- meta$class %in% c("Myeloid", "T cell", "Oligodendrocyte")
  gch <- lapply(CHR_GENES, function(g) { g <- intersect(g, rownames(L)); g[rowMeans(L[g, , drop = FALSE]) >= 1] })
  say("Expressed genes used (mean log2(TPM/10+1) >= 1): chr7 %d, chr10 %d; reference cells %d", length(gch$chr7), length(gch$chr10), sum(ref))
  if (min(lengths(gch)) < 30) stop("Too few expressed chr7/chr10 genes for the copy-number check; send the log to Claude.")
  X <- L[c(gch$chr7, gch$chr10), , drop = FALSE]
  X <- X - rowMeans(X[, ref, drop = FALSE]); X <- pmin(pmax(X, -3), 3)
  meta$cnv <- colMeans(X[gch$chr7, , drop = FALSE]) - colMeans(X[gch$chr10, , drop = FALSE])
  thr_cnv <- quantile(meta$cnv[ref], 0.95)
  meta$cnv_pos <- meta$cnv > thr_cnv
  say("CNV threshold (95th percentile of reference cells) = %.3f", thr_cnv)
  print(as.data.table(meta)[, .(cells = .N, median_score = median(cnv), frac_CNV_supported = mean(cnv_pos)), by = class])
  cnv_t <- as.data.table(meta)[class == "Malignant", .(malignant_cells = .N, median_score = median(cnv),
                                                       frac_CNV_supported = mean(cnv_pos)), by = tumor][order(tumor)]
  cat("Per tumor (malignant-class cells):\n"); print(cnv_t, digits = 2)
  fwrite(cnv_t, file.path(RES_DIR, "T7h_cnv_by_tumor.csv"))
  say("Malignant-class cells CNV-supported: %d of %d (%.0f%%)", sum(mal & meta$cnv_pos), sum(mal), 100 * mean(meta$cnv_pos[mal]))

  # Sensitivity: repeat the key single-cell results in CNV-supported malignant cells only
  malc <- mal & meta$cnv_pos
  sens_det <- rbindlist(lapply(TARGETS, function(g) rbindlist(lapply(list(
    "Malignant (all, by exclusion)" = mal, "Malignant, CNV-supported" = malc, "Malignant class, not CNV-supported" = mal & !meta$cnv_pos,
    "Myeloid" = meta$class == "Myeloid", "T cell" = meta$class == "T cell"), function(mk)
      data.table(cells = sum(mk), mean = mean(L[g, mk]), detected = mean(L[g, mk] > 0))), idcol = "group")[, gene := g]))
  setcolorder(sens_det, c("gene", "group"))
  cat("Detection by group (all vs CNV-supported malignant cells):\n")
  print(dcast(sens_det, gene ~ group, value.var = "detected"), digits = 2)
  cls2 <- as.character(meta$class); cls2[cls2 == "Malignant" & !meta$cnv_pos] <- "Unconfirmed"
  A_cnv <- by_class(factor(cls2), "thr=4, CNV-supported malignant")
  wm_cnv <- within_mal(make_mdat(malc))
  cat("Within CNV-supported malignant cells:\n"); print(wm_cnv[!is.na(r_pooled), .(target, adjustment, n_tumors, n_cells, r_pooled, lo, hi, p)], digits = 2)
  fwrite(rbind(sens_det[, .(analysis = "detection", gene, group, cells, mean, detected)],
               A_cnv$test[, .(analysis = "HMGA1 paired malignant vs other", gene = "HMGA1", group = comparison, n_tumors,
                              median_malignant, median_other, median_TPM_ratio, p)],
               wm_cnv[, .(analysis = "within CNV-supported malignant", gene = target, group = adjustment, n_tumors, n_cells,
                          r_pooled, lo, hi, p)], fill = TRUE),
         file.path(RES_DIR, "T7i_cnv_sensitivity.csv"))
  hc <- data.frame(cnv = meta$cnv, group = ifelse(ref, "Reference (myeloid, T, oligodendrocyte)", "Malignant class"))
  pS8a <- ggplot(hc, aes(cnv, fill = group)) + geom_histogram(bins = 80, alpha = 0.6, position = "identity") +
    geom_vline(xintercept = thr_cnv, linetype = 2, colour = "#C0392B") + scale_y_sqrt() +
    scale_fill_manual(values = c("#7F8C8D", "#3B6FB6"), name = NULL) +
    labs(x = "chr7 minus chr10 expression score", y = "Cells (sqrt scale)", title = "a  Copy-number score by cell group") + theme_pub
  pS8b <- ggplot(as.data.frame(meta[mal, ]), aes(tumor, cnv)) + geom_boxplot(outlier.size = 0.3) +
    geom_hline(yintercept = thr_cnv, linetype = 2, colour = "#C0392B") +
    labs(x = NULL, y = "chr7 minus chr10 score", title = "b  Malignant-class cells per tumor") +
    theme_pub + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7))
  save_fig(pS8a / pS8b, "FigureS8_cnv_check", 9, 8)

  # ---------------- (b2) between tumors, malignant cells only (pseudo-bulk) ----------------
  cat("\n---- (b2) Between tumors, malignant-cell pseudo-bulk (n = tumors; exploratory) ----\n")
  pb <- as.data.table(Mdat)[, c(list(n_cells = .N), lapply(.SD, mean)), by = tumor,
                            .SDcols = setdiff(names(Mdat), "tumor")]
  pb <- as.data.frame(pb[n_cells >= MIN_MAL])
  pbt <- rbindlist(lapply(targets_b, function(y) {
    if (sd(pb[[y]]) == 0) return(data.table(target = y, note = "not expressed in malignant cells"))
    rbindlist(lapply(list(character(0), "cycling"), function(Z) {
      v <- pspear("HMGA1", y, Z, pb)
      data.table(target = y, adjustment = ifelse(length(Z), "cycle", "none"), n_tumors = v[["n"]], rho = v[["r"]], p = v[["p"]])
    }))
  }), fill = TRUE)
  for (i in seq_len(nrow(pbt))) {
    t <- pbt[i]
    if (!is.null(t$note) && !is.na(t$note)) { say("%-9s %s", t$target, t$note); next }
    say("%-9s %-5s n = %2d tumors | rho = %5.2f  P = %s", t$target, t$adjustment, t$n_tumors, t$rho, fmt_p(t$p))
  }
  fwrite(pbt, file.path(RES_DIR, "T7f_neftel_malignant_pseudobulk.csv"))

  # ---------------- Figure 5 ----------------
  pt <- A$per_tumor[class %in% c("Malignant", "Myeloid")]
  paired <- intersect(pt[class == "Malignant", tumor], pt[class == "Myeloid", tumor])
  pt[, paired := tumor %in% paired]
  pA <- ggplot(pt, aes(class, mean_log, group = tumor)) + geom_line(data = pt[paired == TRUE], colour = "grey60") +
    geom_point(aes(colour = paired), size = 1.2) +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey75"), breaks = c("TRUE", "FALSE"),
                        labels = c(`TRUE` = sprintf("paired tumors (n = %d)", length(paired)), `FALSE` = "malignant cells only"), name = NULL) +
    stat_summary(data = pt[paired == TRUE], aes(group = 1), fun = median, geom = "crossbar", width = 0.4, colour = "#C0392B") +
    labs(x = NULL, y = "Mean HMGA1 log2(TPM/10 + 1) per tumor", title = "a  HMGA1, paired within tumors") +
    theme_pub + theme(legend.position = "bottom", legend.direction = "vertical")
  bd <- tg[class %in% c("Malignant", "Myeloid", "T cell", "Oligodendrocyte", "Endothelial")]
  bd$gene <- factor(bd$gene, levels = rev(TARGETS))
  pB <- ggplot(bd, aes(class, gene)) + geom_point(aes(size = detected, colour = mean)) +
    scale_size_continuous(range = c(0.5, 7), name = "Fraction\ndetected") +
    scale_colour_gradient(low = "grey85", high = "#C0392B", name = "Mean\nlog2(TPM/10+1)") +
    labs(x = NULL, y = NULL, title = "b  Target genes by cell class") +
    theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")
  rc <- pbt[target == "CCL2" & adjustment == "none"]
  pC <- ggplot(pb, aes(HMGA1, CCL2)) + geom_point() +   # no OLS line: the statistic is rank-based (Spearman)
    labs(x = "Mean HMGA1 in malignant cells", y = "Mean CCL2 in malignant cells",
         title = if (nrow(rc)) sprintf("Malignant-cell pseudo-bulk per tumor (n = %d): rho = %.2f, P = %s", rc$n_tumors, rc$rho, fmt_p(rc$p))
                 else "Malignant-cell pseudo-bulk per tumor") + theme_pub
  wf <- wm[!is.na(r_pooled)]; wf$target <- pretty_lab(wf$target)
  if (nrow(wf)) save_fig(forest_dodge(wf, "target", "adjustment", "r_pooled", "lo", "hi", NULL,
                     c(none = "grey50", "cycle + complexity" = "#C0392B"),
                     "Pooled within-tumor Spearman r with HMGA1", "Within malignant cells (cell level)"),
                     "FigureS4_within_malignant_cells", 7, 4)
  cd <- as.data.frame(ct); names(cd) <- c("tumor", "class", "n"); cd <- cd[cd$class %in% names(which(colSums(ct) > 0)), ]
  pD <- ggplot(cd, aes(tumor, n, fill = class)) + geom_col() +
    scale_fill_manual(values = c(Malignant = "#7F8C8D", Myeloid = "#C0392B", "T cell" = "#3B6FB6",
                                 Oligodendrocyte = "#27AE60", Endothelial = "#8E44AD"), name = NULL) +
    labs(x = NULL, y = "Cells", title = "c  Cells per tumor (proportions set by CD45 sorting)") +
    theme_pub + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7))
  save_fig((pA + pB + plot_layout(widths = c(0.8, 1.4))) / pD + plot_layout(heights = c(1, 0.8)),
           "Figure5_single_cell_H5", 13, 9)
  save_fig(pC, "FigureS5_neftel_pseudobulk_CCL2", 5, 4.5)

  hd <- rbind(data.frame(score = "myeloid", v = S$myeloid), data.frame(score = "T cell", v = S$tcell),
              data.frame(score = "oligodendrocyte", v = S$oligo), data.frame(score = "endothelial", v = S$endo))
  pS <- ggplot(hd, aes(v)) + geom_histogram(bins = 60) + geom_vline(xintercept = THR, colour = "#C0392B") +
    geom_vline(xintercept = THR_SENS, linetype = 2, colour = "#C0392B") +
    facet_wrap(~score, scales = "free_y") + scale_y_sqrt() +
    labs(x = "Marker score, mean log2(TPM/10 + 1)", y = "Cells (sqrt scale)",
         title = "Neftel et al. marker scores (solid: cutoff 4; dashed: sensitivity 3 and 5)") + theme_pub
  save_fig(pS, "FigureS3_neftel_marker_scores", 9, 5)
  cat("OK Figure 5 and Figs. S3, S4, S5, S8 saved\n")
}

sec_tam_restricted <- function() {
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

    # ---- Formal test: is HMGA1 more strongly associated with the macrophage than the microglial signature? ----
    # Difference of partial Spearman correlations (macrophage - microglia) within each cohort, SE by
    # bootstrap over patients (2,000 resamples), pooled by inverse variance (fixed, DL, HKSJ).
    cat("\n---- Macrophage vs microglial signature: difference in partial correlation with HMGA1 ----\n")
    set.seed(2026); B <- 2000
    pairs <- list(restricted = c("TAM_macrophage_restricted", "TAM_microglia_restricted"),
                  unrestricted = c("TAM_macrophage", "TAM_microglia"))
    dres <- list()
    for (pn in names(pairs)) for (mn in c("M0_unadjusted", "M2_prolif_ESTIMATE")) {
      Z <- MODELS[[mn]]; a <- pairs[[pn]][1]; b <- pairs[[pn]][2]
      per <- rbindlist(lapply(names(cohorts), function(k) {
        d <- cohorts[[k]]
        if (!all(c(a, b) %in% names(d)) || all(is.na(d[[a]])) || all(is.na(d[[b]]))) return(NULL)
        d <- as.data.frame(d)[complete.cases(as.data.frame(d)[, c("log2HMGA1", a, b, Z)]), ]
        dfun <- function(x) pspear("log2HMGA1", a, Z, x)[["r"]] - pspear("log2HMGA1", b, Z, x)[["r"]]
        est <- dfun(d)
        bs <- replicate(B, dfun(d[sample(nrow(d), replace = TRUE), ]))
        data.table(signatures = pn, model = mn, cohort = k, n = nrow(d), diff = est, se = sd(bs),
                   lo = unname(quantile(bs, 0.025)), hi = unname(quantile(bs, 0.975)))
      }))
      m <- meta_generic(per$diff, per$se)
      pooled <- data.table(signatures = pn, model = mn, cohort = "Pooled", n = sum(per$n), diff = m[["est"]],
                           lo = m[["lo"]], hi = m[["hi"]], p = m[["p"]], I2 = m[["I2"]], diff_RE = m[["est_RE"]],
                           lo_HK = m[["lo_HK"]], hi_HK = m[["hi_HK"]], p_HK = m[["p_HK"]])
      dres[[length(dres) + 1]] <- rbind(per, pooled, fill = TRUE)
      say("%-12s %-20s per cohort diff: %s | pooled %.2f (%.2f to %.2f), P = %s, I2 = %.0f%% | HKSJ %.2f to %.2f, P = %s",
          pn, mn, paste(sprintf("%s %.2f", per$cohort, per$diff), collapse = ", "), m[["est"]], m[["lo"]], m[["hi"]],
          fmt_p(m[["p"]]), 100 * m[["I2"]], m[["lo_HK"]], m[["hi_HK"]], fmt_p(m[["p_HK"]]))
    }
    fwrite(rbindlist(dres, fill = TRUE), file.path(RES_DIR, "T11b_tam_difference.csv"))
  }
  run_script11()
}

sec_tables <- function() {
  if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")

  rd <- function(f) { p <- file.path(RES_DIR, f); if (!file.exists(p)) stop("Missing result file: ", p, " (run run_all.R first)"); fread(p) }
  rnd <- function(dt) { dt <- as.data.frame(dt)
    for (c in names(dt)) if (is.numeric(dt[[c]])) dt[[c]] <- signif(dt[[c]], 4)
    dt }

  # ---- S1: all association results ----
  s1 <- rbind(rd("T3_H4_immune_composition.csv"), rd("T4_H1_H2_mechanism_axes.csv"), rd("T5_H3_checkpoints.csv"),
              rd("T5b_HMGA1_vs_subtype_scores.csv"), rd("T5d_CD276_sensitivity.csv"),
              if (file.exists(file.path(RES_DIR, "T11_tam_restricted.csv"))) rd("T11_tam_restricted.csv"), fill = TRUE)
  s1$model_description <- unname(c(MODEL_LABELS, M3b_TCGA_labels = "+ prolif + ESTIMATE + subtype labels (TCGA)")[s1$model])
  s1$post_hoc <- s1$model %in% c("M7_vascular_stromal") | grepl("^TAM_", s1$outcome)
  s1 <- s1[!is.na(model) & model %in% c(names(MODEL_LABELS), "M3b_TCGA_labels")]   # safety: drop malformed rows
  s1$model_description <- unname(c(MODEL_LABELS, M3b_TCGA_labels = "+ prolif + ESTIMATE + subtype labels (TCGA)")[s1$model])
  setcolorder(s1, c("family", "outcome", "model", "model_description", "post_hoc"))

  # ---- S2: survival ----
  s2a <- rd("T1_HMGA1_proliferation_age.csv")
  s2b <- rd("T2_HMGA1_cox.csv")

  # ---- S3: single cell ----
  s3a <- rd("T7a_neftel_cell_counts.csv"); setnames(s3a, 1, "tumor")
  s3b <- rd("T7b_neftel_HMGA1_by_class.csv")
  s3c <- rd("T7e_neftel_targets_by_class.csv")
  s3d <- rd("T7d_neftel_within_malignant.csv")
  s3e <- rd("T7f_neftel_malignant_pseudobulk.csv")
  s3f <- rd("T6_TISCH2_Neftel_Smartseq2.csv")
  s3g <- if (file.exists(file.path(RES_DIR, "T7g_muller_gene_specificity.csv"))) rd("T7g_muller_gene_specificity.csv") else data.table()

  # ---- S4: all-grade analysis ----
  s4 <- rd("T8_all_grades_reconciliation.csv")
  opt <- function(f) if (file.exists(file.path(RES_DIR, f))) rd(f) else data.table(note = paste(f, "not produced"))
  s3h <- opt("T7h_cnv_by_tumor.csv"); s3i <- opt("T7i_cnv_sensitivity.csv"); s3j <- opt("T11b_tam_difference.csv")
  qc  <- opt("T0_cgga_duplicate_check.csv")

  # ---- S5: gene sets ----
  gs <- list()
  add <- function(set, genes, source) gs[[length(gs) + 1]] <<- data.table(gene_set = set, gene = genes, source = source)
  add("Proliferation score", PROLIF_GENES, "This study (10 canonical proliferation genes)")
  add("He et al. ISGs", SIGNATURES$He_ISG, "He et al., Nat Commun 2025")
  add("Checkpoint genes", CHECKPOINTS, "This study (prespecified, H3)")
  for (f in c(WANG_FILE, MULLER_FILE, NEFTEL_FILE)) if (file.exists(f)) {
    w <- fread(f); src <- c("Wang et al., Cancer Cell 2017 (Table S1)", "Mueller et al., Genome Biol 2017 (Table S4)",
                            "Neftel et al., Cell 2019 (Table S2)")[match(f, c(WANG_FILE, MULLER_FILE, NEFTEL_FILE))]
    add(paste(basename(tools::file_path_sans_ext(f)), w$set, sep = ": "), w$gene, src)
  }
  if (requireNamespace("msigdbr", quietly = TRUE)) {
    h <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
                  error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
    add("HALLMARK_INTERFERON_ALPHA_RESPONSE", sort(unique(h$gene_symbol[h$gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE"])),
        paste0("Liberzon et al., Cell Syst 2015; msigdbr ", packageVersion("msigdbr")))
  }
  add("Neftel cell-class markers: myeloid", c("CD14", "AIF1", "FCER1G", "FCGR3A", "TYROBP", "CSF1R"), "Neftel et al., Cell 2019 (STAR Methods)")
  add("Neftel cell-class markers: T cell", c("CD2", "CD3D", "CD3E", "CD3G"), "Neftel et al., Cell 2019 (STAR Methods)")
  add("Neftel cell-class markers: oligodendrocyte", c("MBP", "TF", "PLP1", "MAG", "MOG", "CLDN11"), "Neftel et al., Cell 2019 (STAR Methods)")
  add("Endothelial markers", c("PECAM1", "VWF", "CLDN5"), "This study")
  s5 <- rbindlist(gs)

  readme <- data.frame(Sheet = c("S1_associations", "S2a_proliferation_age", "S2b_cox", "S3a_cell_counts", "S3b_HMGA1_by_class",
                                 "S3c_targets_by_class", "S3d_within_malignant", "S3e_pseudobulk", "S3f_TISCH2", "S3g_TAM_gene_specificity",
                                 "S3h_CNV_by_tumor", "S3i_CNV_sensitivity", "S3j_TAM_difference", "S4_all_grades", "S5_gene_sets", "QC_CGGA_duplicates"),
    Content = c(
      "Partial Spearman correlations of HMGA1 with all outcomes: per-cohort r, P and n; pooled fixed-effect r (95% CI), P, Benjamini-Hochberg q (within family and model), I2, direction consistency, robustness flag (q < 0.05, I2 < 50%, same direction in all cohorts), DerSimonian-Laird random-effects r (r_RE) with its 95% CI, and Hartung-Knapp-Sidik-Jonkman 95% CI and P (lo_HK, hi_HK, p_HK). Models M0-M7 as defined in Methods; post_hoc marks post hoc analyses.",
      "Spearman correlations of HMGA1 with proliferation score, MKI67 and age (per cohort and pooled).",
      "Cox models: HR per 1 SD log2 HMGA1, 95% CI, P; Schoenfeld tests for HMGA1 (zph_HMGA1_p) and for the whole model (zph_global_p); time-varying sensitivity (HMGA1 x log(t+1): tvc_p, HR at 6/12/24 months); complete-case n; pooled fixed-effect, random-effects and HKSJ estimates.",
      "Neftel Smart-seq2 adult tumors: cells per tumor and class (Neftel markers, cutoff 4). No cell reached the endothelial cutoff. Proportions reflect CD45 sorting, not tumor composition.",
      "HMGA1 in malignant vs myeloid/T cells, paired within tumors with >= 10 cells of both classes, at cutoffs 3, 4 and 5.",
      "Mean log2(TPM/10+1) and fraction detected of HMGA1 and target genes by cell class.",
      "Within-malignant-cell correlations with HMGA1 per target, pooled across tumors (unadjusted and cell-cycle/complexity adjusted). Nominal P values (exploratory).",
      "Tumor-level malignant-cell pseudo-bulk correlations (n = 20 tumors). Nominal P values (exploratory).",
      "TISCH2 cell-type means for the same Neftel Smart-seq2 dataset (secondary reference; not independent validation).",
      "Mueller TAM signature genes: mean expression and detection in myeloid vs malignant cells; myeloid-restricted flag (post hoc criterion).",
      "Copy-number check: chr7-minus-chr10 expression score of malignant-class cells per tumor and fraction above the reference threshold.",
      "Key single-cell results repeated in CNV-supported malignant cells: detection by group, paired HMGA1 comparison, within-malignant correlations.",
      "Difference in partial correlation of HMGA1 with the blood-derived macrophage vs microglial signature (bootstrap SE, 2,000 resamples), per cohort and pooled (fixed, random-effects, HKSJ); post hoc.",
      "All-grade CGGA analysis (post hoc): HR per 1 SD HMGA1 and MKI67 (positive control) across sequential models; within-molecular-group estimates.",
      "All gene sets used, with sources.",
      "Quality control: five most similar CGGA_325 / CGGA_693 sample pairs by genome-wide expression (raw and within-cohort centred correlation) with age, sex and overall survival."))

  out <- file.path(RES_DIR, "Supplementary_Tables_S1-S5.xlsx")
  writexl::write_xlsx(list(README = readme, S1_associations = rnd(s1), S2a_proliferation_age = rnd(s2a), S2b_cox = rnd(s2b),
                           S3a_cell_counts = rnd(s3a), S3b_HMGA1_by_class = rnd(s3b), S3c_targets_by_class = rnd(s3c),
                           S3d_within_malignant = rnd(s3d), S3e_pseudobulk = rnd(s3e), S3f_TISCH2 = rnd(s3f), S3g_TAM_gene_specificity = rnd(s3g),
                           S3h_CNV_by_tumor = rnd(s3h), S3i_CNV_sensitivity = rnd(s3i), S3j_TAM_difference = rnd(s3j),
                           S4_all_grades = rnd(s4), S5_gene_sets = as.data.frame(s5), QC_CGGA_duplicates = rnd(qc)), out)
  say("Junk rows check: %d rows in S1 with an unknown model (must be 0)", sum(!s1$model %in% c(names(MODEL_LABELS), "M3b_TCGA_labels")))
  say("Saved %s (S1: %d rows; S5: %d genes)", out, nrow(s1), nrow(s5))
}

SECTIONS <- list(tisch2 = sec_tisch2, neftel = sec_neftel, tam_restricted = sec_tam_restricted, tables = sec_tables)
if (exists("ONLY")) { SECTIONS <- SECTIONS[names(SECTIONS) %in% ONLY]; rm(ONLY) }
run_logged(file.path(RES_DIR, "log_03_single_cell.txt"), SECTIONS)
