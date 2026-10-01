# =====================================================
# 00b_microenvironment_scores.R
# Computes ESTIMATE immune/stromal scores (input: log2[expression + 1])
# and MCP-counter populations (via immunedeconv) for TCGA, CGGA_325 and
# CGGA_693. Run ONCE after 00a_download_and_prepare_data.R, from the
# project folder (the one containing "share/").
# Output: share/scores_TCGA.csv, share/scores_CGGA325.csv, share/scores_CGGA693.csv
# Versions used: estimate 1.0.13, immunedeconv 2.1.4, MCPcounter 1.2.0
# =====================================================

library(estimate)
library(immunedeconv)

dir.create("share/tmp_estimate", recursive = TRUE, showWarnings = FALSE)

# ---------- helpers ----------
collapse_by_symbol <- function(mat, symbols) {
  keep <- !is.na(symbols) & symbols != ""
  mat <- mat[keep, , drop = FALSE]; symbols <- symbols[keep]
  mean_expr <- rowMeans(mat)
  idx <- tapply(seq_len(nrow(mat)), symbols, function(i) i[which.max(mean_expr[i])])
  out <- mat[unlist(idx), , drop = FALSE]
  rownames(out) <- names(idx)
  out
}

run_estimate <- function(mat, label) {
  in_f  <- file.path("share/tmp_estimate", paste0(label, "_input.txt"))
  cg_f  <- file.path("share/tmp_estimate", paste0(label, "_common.gct"))
  out_f <- file.path("share/tmp_estimate", paste0(label, "_estimate.gct"))
  df <- data.frame(NAME = rownames(mat), log2(mat + 1), check.names = FALSE)
  write.table(df, in_f, sep = "\t", quote = FALSE, row.names = FALSE)
  filterCommonGenes(input.f = in_f, output.f = cg_f, id = "GeneSymbol")
  estimateScore(input.ds = cg_f, output.ds = out_f, platform = "illumina")
  sc <- read.table(out_f, skip = 2, header = TRUE, sep = "\t", check.names = FALSE)
  rownames(sc) <- sc$NAME
  sc <- t(sc[, -(1:2)])
  data.frame(sample = colnames(mat),        # keep original IDs (read.table may alter them)
             ESTIMATE_stromal = as.numeric(sc[, "StromalScore"]),
             ESTIMATE_immune  = as.numeric(sc[, "ImmuneScore"]),
             ESTIMATE_score   = as.numeric(sc[, "ESTIMATEScore"]),
             check.names = FALSE)
}

run_mcp <- function(mat) {
  res <- immunedeconv::deconvolute(mat, "mcp_counter")
  m <- as.data.frame(t(as.matrix(res[, -1])))
  colnames(m) <- paste0("MCP_", gsub("[^A-Za-z0-9]+", "_", res$cell_type))
  m$sample <- rownames(m)
  m
}

process <- function(mat, label) {
  cat("\n=== ", label, ": ", nrow(mat), " genes x ", ncol(mat), " samples ===\n", sep = "")
  est <- run_estimate(mat, label)
  mcp <- run_mcp(mat)
  out <- merge(est, mcp, by = "sample")
  stopifnot(nrow(out) == ncol(mat))
  write.csv(out, paste0("share/scores_", label, ".csv"), row.names = FALSE)
  cat("✓ Saved share/scores_", label, ".csv\n", sep = "")
  print(summary(out$ESTIMATE_immune))
}

# ---------- TCGA ----------
tpm    <- readRDS("share/tcga_tpm.rds")
gi     <- readRDS("share/tcga_gene_info.rds")
cohort <- readRDS("share/tcga_cohort_clinical.rds")
sym    <- gi$gene_name[match(rownames(tpm), rownames(gi))]
tcga   <- collapse_by_symbol(tpm[, rownames(cohort)], sym)
process(tcga, "TCGA")

# ---------- CGGA ----------
for (id in c("325", "693")) {
  ex <- readRDS(paste0("share/cgga", id, "_expr_gbm.rds"))
  mat <- as.matrix(ex[, -1]); mode(mat) <- "numeric"
  mat <- collapse_by_symbol(mat, ex[[1]])
  process(mat, paste0("CGGA", id))
}

cat("\n✓ All done: share/scores_TCGA.csv, share/scores_CGGA325.csv, share/scores_CGGA693.csv\n")
