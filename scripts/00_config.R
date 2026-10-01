# =====================================================================
# 00_config.R
# HMGA1 and the immune microenvironment of IDH-wildtype glioblastoma
# Shared settings, gene sets, models and helper functions.
# Sourced by every script. Working directory MUST be D:/Siraj_HMGA1
# =====================================================================

suppressPackageStartupMessages({
  library(survival); library(survminer)
  library(ggplot2); library(patchwork)
  library(data.table); library(readxl)
})

if (!dir.exists("share")) stop("Folder 'share' not found. Run setwd('D:/Siraj_HMGA1') first.")
set.seed(2026)
DATA_DIR <- "share"
RAW_DIR  <- file.path(DATA_DIR, "raw")
RES_DIR  <- "results"
FIG_DIR  <- "figures"
dir.create(RES_DIR, showWarnings = FALSE); dir.create(FIG_DIR, showWarnings = FALSE)

# ---------------------------------------------------------------------
# Gene sets
# ---------------------------------------------------------------------
PROLIF_GENES <- c("MKI67","TOP2A","CCNB1","CDK1","BUB1","CENPF","AURKA","PCNA","MCM2","TPX2")
CHECKPOINTS  <- c("HAVCR2","PDCD1","CD274","CTLA4","LAG3","TIGIT","CD276")

# Published gene sets only (every set must be traceable to a citation)
#   He_ISG      : ISGs measured by He et al. 2025 (Nat Commun 16:5098) downstream of
#                 HMGA1-STING (CXCL10 and CCL5 are tested as separate outcomes)
#   Hallmark IFN-alpha response (Liberzon et al. 2015, Cell Systems) -> loaded with msigdbr in 01
#   TAM abundance: MCP-counter monocytic lineage (Becht et al. 2016) -> scores files
#   Microglia vs blood-derived TAMs: Mueller et al. 2017 (Genome Biol 18:234) -> optional file
SIGNATURES <- list(He_ISG = c("IFIT1", "IFIT2", "IFIT3"))

# Genes that must exist in every cohort (script stops if missing)
REQUIRED_GENES <- unique(c("HMGA1","CCL2","STING1","CXCL10","CCL5", CHECKPOINTS, PROLIF_GENES))
GENES_NEEDED   <- unique(c(REQUIRED_GENES, unlist(SIGNATURES)))

# Official symbol -> older symbols used by some datasets (CGGA, Neftel use TMEM173)
ALIASES <- list(STING1 = c("STING1","TMEM173"))

# Optional user-supplied signature files (analyses are skipped if absent)
WANG_FILE   <- file.path(DATA_DIR, "wang2017_subtype_signatures.csv")   # Wang et al. 2017, Table S1
NEFTEL_FILE <- file.path(DATA_DIR, "neftel2019_modules.csv")            # Neftel et al. 2019, Table S2
MULLER_FILE <- file.path(DATA_DIR, "muller2017_tam_signatures.csv")      # Mueller et al. 2017; columns microglia, macrophage

# ---------------------------------------------------------------------
# Prespecified adjustment models (Decided; handover Sect. 5)
# M4 is TCGA-only (ABSOLUTE purity exists only in TCGA).
# A model is skipped for an outcome that is itself one of its covariates.
# ---------------------------------------------------------------------
MODELS <- list(
  M0_unadjusted      = character(0),
  M1_prolif          = "prolif",
  M2_prolif_ESTIMATE = c("prolif","ESTIMATE_immune","ESTIMATE_stromal"),
  M3_plus_subtype    = c("prolif","ESTIMATE_immune","ESTIMATE_stromal","sub_MES","sub_PN","sub_CL"),
  M4_TCGA_ABSOLUTE   = c("prolif","ABSOLUTE_purity"),
  M5_composition     = c("prolif","ESTIMATE_immune","ESTIMATE_stromal","MCP_Monocyte")
)
MODEL_LABELS <- c(M0_unadjusted = "Unadjusted", M1_prolif = "+ proliferation",
                  M2_prolif_ESTIMATE = "+ prolif + ESTIMATE", M3_plus_subtype = "+ prolif + ESTIMATE + subtype",
                  M4_TCGA_ABSOLUTE = "+ prolif + ABSOLUTE (TCGA)", M5_composition = "+ prolif + ESTIMATE + monocytic",
                  M6_STING_composition = "+ prolif + ESTIMATE + monocytic + endothelial",
                  M7_vascular_stromal = "+ prolif + ESTIMATE + monocytic + endothelial + fibroblast-like")

# ---------------------------------------------------------------------
# Core statistical helper functions
# ---------------------------------------------------------------------
# Partial Spearman: rank all variables, regress covariates Z out of x and y,
# correlate residuals; t-test with n - 2 - k degrees of freedom.
pspear <- function(x, y, Z = character(0), data) {
  dd <- na.omit(data[, c(x, y, Z), drop = FALSE])
  r  <- as.data.frame(lapply(dd, rank))
  n  <- nrow(dd); k <- length(Z)
  if (k == 0) {
    rr <- cor(r[[x]], r[[y]])
  } else {
    rx <- resid(lm(r[[x]] ~ ., data = r[, Z, drop = FALSE]))
    ry <- resid(lm(r[[y]] ~ ., data = r[, Z, drop = FALSE]))
    rr <- cor(rx, ry)
  }
  tt <- rr * sqrt((n - 2 - k) / (1 - rr^2))
  c(n = n, r = rr, p = 2 * pt(-abs(tt), n - 2 - k))
}

# Fixed-effect meta-analysis of correlations (Fisher z)
meta_r <- function(r, n, k = 0) {
  z <- atanh(r); w <- n - 3 - k
  zp <- sum(w * z) / sum(w); se <- 1 / sqrt(sum(w))
  Q  <- sum(w * (z - zp)^2); df <- length(r) - 1
  c(r = tanh(zp), lo = tanh(zp - 1.96 * se), hi = tanh(zp + 1.96 * se),
    p = 2 * pnorm(-abs(zp / se)), Q_p = pchisq(Q, df, lower.tail = FALSE),
    I2 = ifelse(Q > 0, max(0, (Q - df) / Q), 0))
}

# ---------------------------------------------------------------------
# New helpers
# ---------------------------------------------------------------------
say <- function(...) cat(sprintf(...), "\n", sep = "")
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.0001, "<0.0001", formatC(p, format = "g", digits = 2)))

# resolve a gene symbol (or its alias) among available symbols
find_symbol <- function(g, available) {
  cands <- unique(c(g, ALIASES[[g]]))
  hit <- cands[cands %in% available]
  if (length(hit)) hit[1] else NA_character_
}

# collapse a gene x sample matrix to one row per symbol (highest mean expression)
collapse_symbols <- function(mat, sym) {
  ok  <- !is.na(sym) & sym != ""
  mat <- mat[ok, , drop = FALSE]; sym <- sym[ok]
  o   <- order(rowMeans(mat), decreasing = TRUE)
  mat <- mat[o, , drop = FALSE]; sym <- sym[o]
  keep <- !duplicated(sym)
  mat <- mat[keep, , drop = FALSE]; rownames(mat) <- sym[keep]
  mat
}

# pull genes (canonical names) for given samples; required genes must exist
pull_genes <- function(mat, genes, samples, label) {
  out <- sapply(genes, function(g) {
    s <- find_symbol(g, rownames(mat))
    if (is.na(s)) {
      if (g %in% REQUIRED_GENES) stop("Required gene missing in ", label, ": ", g)
      return(rep(NA_real_, length(samples)))
    }
    as.numeric(mat[s, samples])
  })
  miss <- genes[colSums(!is.na(out)) == 0]
  if (length(miss)) say("  [%s] optional genes not found: %s", label, paste(miss, collapse = ", "))
  out
}

# mean of within-cohort z-scored log2(x+1) over the available genes
sig_score <- function(d, genes) {
  g <- genes[genes %in% names(d)]
  g <- g[colSums(!is.na(d[, g, drop = FALSE])) > 0]
  m <- scale(log2(as.matrix(d[, g, drop = FALSE]) + 1))
  m <- m[, colSums(is.na(m)) < nrow(m), drop = FALSE]   # drop zero-variance genes
  list(score = rowMeans(m, na.rm = TRUE), n_used = ncol(m))
}

# ---------------------------------------------------------------------
# Association engine: partial Spearman per cohort + fixed-effect pooling
# ---------------------------------------------------------------------
assoc_one <- function(cohorts, x, y, Z) {
  per <- lapply(names(cohorts), function(k) {
    d <- cohorts[[k]]
    if (!all(c(x, y, Z) %in% names(d))) return(NULL)
    if (any(sapply(c(x, y, Z), function(v) sum(!is.na(d[[v]])) < 10))) return(NULL)
    v <- pspear(x, y, Z, d)
    data.frame(cohort = k, n = v[["n"]], r = v[["r"]], p = v[["p"]])
  })
  per <- do.call(rbind, per)
  if (is.null(per) || nrow(per) == 0) return(NULL)
  k <- length(Z)
  m <- meta_r(per$r, per$n, k)
  if (nrow(per) == 1) { m["p"] <- per$p; m["I2"] <- NA; m["Q_p"] <- NA }
  list(per = per, pooled = m, same_dir = length(unique(sign(per$r))) == 1)
}

assoc_table <- function(cohorts, x, outcomes, models, family) {
  rows <- list()
  for (y in outcomes) for (mn in names(models)) {
    Z <- models[[mn]]
    if (y %in% Z) next
    a <- assoc_one(cohorts, x, y, Z)
    if (is.null(a)) next
    row <- data.frame(family = family, outcome = y, model = mn, k = length(Z),
                      n_cohorts = nrow(a$per), n_total = sum(a$per$n),
                      r_pooled = a$pooled[["r"]], lo = a$pooled[["lo"]], hi = a$pooled[["hi"]],
                      p_pooled = a$pooled[["p"]], I2 = a$pooled[["I2"]], same_direction = a$same_dir)
    for (i in seq_len(nrow(a$per))) {
      row[[paste0("r_", a$per$cohort[i])]] <- a$per$r[i]
      row[[paste0("p_", a$per$cohort[i])]] <- a$per$p[i]
      row[[paste0("n_", a$per$cohort[i])]] <- a$per$n[i]
    }
    rows[[length(rows) + 1]] <- row
  }
  tab <- rbindlist(rows, fill = TRUE)
  # Benjamini-Hochberg within the hypothesis family, separately per model
  tab[, q_pooled := p.adjust(p_pooled, "BH"), by = model]
  tab[, robust := q_pooled < 0.05 & !is.na(I2) & I2 < 0.5 & same_direction & n_cohorts == 3]
  as.data.frame(tab)
}

print_assoc <- function(tab) {
  for (i in seq_len(nrow(tab))) {
    t <- tab[i, ]
    per <- paste(sapply(c("TCGA","CGGA_325","CGGA_693"), function(k) {
      rk <- t[[paste0("r_", k)]]; if (is.null(rk) || is.na(rk)) "   -  " else sprintf("%6.2f", rk)
    }), collapse = " ")
    say("%-16s %-22s | %s | pooled r = %5.2f [%5.2f, %5.2f]  q = %-7s I2 = %3s%%  %s",
        t$outcome, t$model, per, t$r_pooled, t$lo, t$hi, fmt_p(t$q_pooled),
        ifelse(is.na(t$I2), " NA", sprintf("%3.0f", 100 * t$I2)), ifelse(isTRUE(t$robust), "ROBUST", ""))
  }
}

# ---------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------
theme_pub <- theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10),
        strip.background = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "bottom")
COL_LOWHIGH <- c(Low = "#3B6FB6", High = "#C0392B")
COL_MODELS  <- c("Unadjusted" = "grey55", "+ proliferation" = "#3B6FB6",
                 "+ prolif + ESTIMATE" = "#C0392B", "+ prolif + ESTIMATE + subtype" = "#8E44AD",
                 "+ prolif + ABSOLUTE (TCGA)" = "#E67E22", "+ prolif + ESTIMATE + monocytic" = "#16A085",
                 "+ prolif + ESTIMATE + monocytic + endothelial" = "#2C3E50",
                 "+ prolif + ESTIMATE + monocytic + endothelial + fibroblast-like" = "#D35400")

PRETTY <- c(ESTIMATE_immune = "ESTIMATE immune", ESTIMATE_stromal = "ESTIMATE stromal",
            MCP_T_cell = "T cells", MCP_T_cell_CD8_ = "CD8 T cells", MCP_cytotoxicity_score = "Cytotoxic lymphocytes",
            MCP_NK_cell = "NK cells", MCP_B_cell = "B lineage", MCP_Monocyte = "Monocytic lineage (TAMs)",
            MCP_Myeloid_dendritic_cell = "Myeloid dendritic cells", MCP_Neutrophil = "Neutrophils",
            MCP_Endothelial_cell = "Endothelial cells", MCP_Cancer_associated_fibroblast = "Fibroblast-like",
            TAM_microglia = "Microglial TAMs (Mueller)", TAM_macrophage = "Blood-derived TAMs (Mueller)",
            He_ISG = "ISGs (He et al.)", MES_like = "MES-like (Neftel)", IFNa_hallmark = "Hallmark IFN-alpha",
            HAVCR2 = "HAVCR2 (TIM-3)", CD274 = "CD274 (PD-L1)", CD276 = "CD276 (B7-H3)", PDCD1 = "PDCD1 (PD-1)")
pretty_lab <- function(x) ifelse(x %in% names(PRETTY), PRETTY[x], x)

# Generic dodged forest plot (manual vertical offsets; robust across ggplot2 versions)
#   df: data with columns for category (cat), group (grp), estimate, lo, hi
forest_dodge <- function(df, cat, grp, est = "r", lo = "lo", hi = "hi", order = NULL, colours = NULL,
                         xlab = "", title = "", ref = 0, log_x = FALSE, shape_col = NULL, width = 0.7) {
  df <- as.data.frame(df)
  order <- if (is.null(order)) unique(as.character(df[[cat]])) else order
  order <- order[order %in% df[[cat]]]
  gl <- if (is.factor(df[[grp]])) levels(droplevels(df[[grp]])) else unique(as.character(df[[grp]]))
  ng <- length(gl)
  df$.g <- factor(as.character(df[[grp]]), levels = gl)
  # offsets computed within each category, so rows with fewer groups stay centred
  df$.rank <- ave(as.integer(df$.g), df[[cat]], FUN = function(v) rank(v))
  df$.n    <- ave(as.integer(df$.g), df[[cat]], FUN = length)
  df$.y <- (length(order) + 1 - match(as.character(df[[cat]]), order)) -
           (df$.rank - (df$.n + 1) / 2) * width / max(ng, 1)
  df$.e <- df[[est]]; df$.lo <- df[[lo]]; df$.hi <- df[[hi]]
  g <- ggplot(df, aes(x = .e, y = .y, colour = .g)) +
    geom_vline(xintercept = ref, linetype = 2, colour = "grey60") +
    geom_segment(aes(x = .lo, xend = .hi, yend = .y))
  g <- if (is.null(shape_col)) g + geom_point(size = 2) else
    g + geom_point(aes(shape = .data[[shape_col]]), size = 2) +
        scale_shape_manual(values = c("q < 0.05" = 16, "n.s." = 1), limits = c("q < 0.05", "n.s."), drop = FALSE, name = NULL)
  g <- g + scale_y_continuous(breaks = rev(seq_along(order)), labels = order, expand = expansion(add = 0.6)) +
    labs(x = xlab, y = NULL, title = title) + theme_pub + guides(colour = guide_legend(ncol = 2))
  g <- if (is.null(colours)) g + scale_colour_discrete(name = NULL) else g + scale_colour_manual(values = colours, name = NULL)
  if (log_x) g <- g + scale_x_log10(breaks = c(0.25, 0.5, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 2, 3, 4))
  g
}

forest_models <- function(tab, title, outcome_order = unique(tab$outcome)) {
  tab <- as.data.frame(tab)
  tab$model_lab <- factor(MODEL_LABELS[tab$model], levels = MODEL_LABELS)
  tab$sig <- ifelse(tab$q_pooled < 0.05, "q < 0.05", "n.s.")
  tab$outcome <- pretty_lab(tab$outcome)
  forest_dodge(tab, "outcome", "model_lab", "r_pooled", "lo", "hi", pretty_lab(outcome_order), COL_MODELS,
               "Pooled partial Spearman r with HMGA1 (95% CI)", title, shape_col = "sig")
}

save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}
