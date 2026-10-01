# =====================================================================
# 06_single_cell_TISCH2.R  (Supplementary Figure S1)
# TISCH2 cell-type mean expression, Neftel GSE131928 Smart-seq2 only.
# TISCH2 fills values for cell types ABSENT from a dataset, so only the
# confirmed Smart-seq2 cell types are kept (handover Sect. 7).
# =====================================================================
source("scripts/00_config.R", local = TRUE)

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
cat("OK Figure S1 saved\n")
