# =====================================================================
# 07_neftel_single_cell.R  (Figure 5, H5)
# Neftel et al. 2019 Smart-seq2 (GSE131928), ADULT tumors only.
#
# Design note: Smart-seq2 cells were FACS-sorted into CD45+ and CD45-
# plates (and MGH143 was CD45-depleted), so the immune FRACTION per
# tumor reflects sorting, not biology. This script therefore does NOT
# correlate HMGA1 with immune fraction. Instead it tests:
#   (a) HMGA1 level by cell class, paired within tumors (H5)
#   (b) HMGA1 vs CCL2 / STING1 / HAVCR2 / CD274 WITHIN malignant cells,
#       adjusted for cell cycle and cell complexity (composition-free)
#   (c) optional: HMGA1 vs MES-like state score (needs neftel2019_modules.csv)
# First run: decompresses the TPM file once and caches the needed genes
# (takes several minutes). Later runs use results/neftel_subset.rds.
# =====================================================================
source("scripts/00_config.R", local = TRUE)

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
want <- unique(c(TARGETS, unlist(MARKERS), PROLIF_GENES, unlist(neftel_mod), unlist(muller_sets)))
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
Mdat <- data.frame(tumor = meta$tumor[mal], HMGA1 = L["HMGA1", mal], cycling = meta$cycling[mal],
                   ngenes = meta$ngenes[mal], t(L[setdiff(TARGETS, "HMGA1"), mal, drop = FALSE]), check.names = FALSE)
if (!is.null(neftel_mod)) {
  z <- t(scale(t(L[, mal])))
  for (m in names(neftel_mod)) {
    g <- setdiff(intersect(neftel_mod[[m]], rownames(z)), "HMGA1");   # HMGA1 excluded: avoids circularity
    if (length(g) < 5) { say("  Neftel module %-5s: only %d genes found -> skipped", m, length(g)); next }
    Mdat[[paste0("mod_", m)]] <- colMeans(z[g, , drop = FALSE], na.rm = TRUE)
    say("  Neftel module %-5s genes found: %d/%d", m, length(g), length(neftel_mod[[m]]))
  }
  mes <- grep("^mod_MES", names(Mdat), value = TRUE)
  if (length(mes)) Mdat$MES_like <- rowMeans(Mdat[, mes, drop = FALSE])
}
targets_b <- c(setdiff(TARGETS, "HMGA1"), if ("MES_like" %in% names(Mdat)) "MES_like")
wm <- rbindlist(lapply(targets_b, function(y) {
  rbindlist(lapply(list(c(), c("cycling", "ngenes")), function(Z) {
    per <- rbindlist(lapply(split(Mdat, Mdat$tumor), function(d) {
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
for (i in seq_len(nrow(wm))) {
  t <- wm[i]
  if (is.na(t$r_pooled)) { say("%-9s %-19s only %d tumors pass detection filters -> not estimable", t$target, t$adjustment, t$n_tumors); next }
  say("%-9s %-19s tumors = %2d cells = %4d | pooled r = %5.2f [%5.2f, %5.2f] P = %-7s I2 = %3.0f%% | median r = %5.2f, positive in %3.0f%% of tumors",
      t$target, t$adjustment, t$n_tumors, t$n_cells, t$r_pooled, t$lo, t$hi, fmt_p(t$p), 100 * t$I2, t$median_r, 100 * t$frac_positive)
}
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
# expressed by malignant or other cells. Prespecified criterion for a myeloid-restricted gene:
# mean log2(TPM/10+1) in myeloid cells >= 1 unit above malignant cells AND detection in myeloid
# cells >= 2x detection in malignant cells.
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
pA <- ggplot(pt, aes(class, mean_log, group = tumor)) + geom_line(colour = "grey70") + geom_point(size = 1.2) +
  stat_summary(aes(group = 1), fun = median, geom = "crossbar", width = 0.4, colour = "#C0392B") +
  labs(x = NULL, y = "Mean HMGA1 log2(TPM/10 + 1) per tumor", title = "A  HMGA1, paired within tumors") + theme_pub
bd <- tg[class %in% c("Malignant", "Myeloid", "T cell", "Oligodendrocyte", "Endothelial")]
bd$gene <- factor(bd$gene, levels = rev(TARGETS))
pB <- ggplot(bd, aes(class, gene)) + geom_point(aes(size = detected, colour = mean)) +
  scale_size_continuous(range = c(0.5, 7), name = "Fraction\ndetected") +
  scale_colour_gradient(low = "grey85", high = "#C0392B", name = "Mean\nlog2(TPM/10+1)") +
  labs(x = NULL, y = NULL, title = "B  Target genes by cell class") +
  theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")
rc <- pbt[target == "CCL2" & adjustment == "none"]
pC <- ggplot(pb, aes(HMGA1, CCL2)) + geom_point() +
  geom_smooth(method = "lm", se = FALSE, colour = "#C0392B", linewidth = 0.6, formula = y ~ x) +
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
  labs(x = NULL, y = "Cells", title = "C  Cells per tumor (proportions set by CD45 sorting)") +
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
cat("OK Figure 5 and S2 saved\n")
