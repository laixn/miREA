# ---------------------------------------------------------------------------
# Data provenance: results/miRTarBase10_FunctionalMTI_unique.csv
#
# This file is treated as the RAW INPUT for this folder. 
# It was built once, offline, from miRTarBase release 10.0's human (hsa) miRNA-Target Interaction file 
# (downloaded from https://awi.cuhk.edu.cn/miRTarBase/downloads/files/10.0/miRTarBase_MTI.csv) by:
#   1. keeping only records with Support Type == "Functional MTI" (the strong,
#      reporter-assay/western-blot-equivalent evidence tier),
#   2. trimming whitespace from the miRNA and Target Gene fields,
#   3. deduplicating to unique (miRNA, gene) pairs,
#   4. formatting as a 2-column TERM2GENE table: ID = "miRTarBase" (constant),
#      MGI = "miRNA:gene".
# Result: 8903 unique Functional-MTI (miRNA, gene) pairs. This same pair set
# is used both as the GSEA pathway (TERM2GENE) below and as the tier-membership
# flag driving the ORA enrichment sweep.
# ---------------------------------------------------------------------------

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(fgsea)
  library(cowplot)
  library(enrichplot)
  library(enrichit)
})

# fgsea::fgseaMultilevel()'s enrichment-score sign is deterministic given a fixed ranked list
# and gene set -- this seed only controls the small Monte-Carlo noise in the p-value estimate
# (magnitude, not direction).
GSEA_SEED <- 1

BASE <- "C:/Users/rtd648/Desktop/project/miREA/submission/s3_ploscb/revise/NoteS5_miRTarBase_analysis"
RESULTS <- file.path(BASE, "results")
SCORE_DIR <- "C:/Users/rtd648/Desktop/project/miREA/00_revision/revise/major3/method13/results/mirtarbase"

cancers <- c("BLCA", "BRCA", "CESC", "COAD", "ESCA", "KICH", "KIRC",
             "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD", "READ",
             "STAD", "THCA", "UCEC")
THRESHOLDS <- c(0.01, 0.05, 0.10, 0.20, 0.30, 0.50, 0.75, 1.00)
THRESH_LABELS <- c("1", "5", "10", "20", "30", "50", "75", "100")
TIER_LABEL <- "miRTarBase: Functional MTI"

# ---------- shared raw input: the Functional-MTI pair set ----------
pathway <- fread(file.path(RESULTS, "miRTarBase10_FunctionalMTI_unique.csv"))
functional_mti_pairs <- pathway[, tstrsplit(MGI, ":", fixed = TRUE, names = c("miRNA", "gene"))]
cat("Functional-MTI pairs (raw input):", nrow(functional_mti_pairs), "\n")

# ---------- GSEA: called directly via fgsea::fgseaMultilevel(), NOT clusterProfiler::GSEA() ----------
# clusterProfiler::GSEA() in the installed version silently drops the scoreType argument (a known
# regression: clusterProfiler switched its GSEA backend to the 'enrichit' package in 4.19.3, and
# the GSEA() wrapper did not forward extra arguments -- including scoreType -- to
# enrichit::gsea_gson() until the fix landed in 4.21.1/4.21.2). Without scoreType, GSEA() silently
# falls back to the standard two-tailed test instead of the one-tailed "neg" test this analysis
# was designed around. Calling fgsea directly restores the intended one-tailed test, and
# reproduces (to the 3rd decimal) the original results computed on the original analysis
# machine's clusterProfiler 4.14.4, which predates the enrichit migration and forwarded
# scoreType correctly.
functional_mti_ids <- unique(pathway$MGI)
pathway_list <- list(miRTarBase = functional_mti_ids)

gsea_rows <- vector("list", length(cancers))
set.seed(GSEA_SEED)

for (i in seq_along(cancers)) {
  cancer <- cancers[i]
  score_data <- read.csv(file.path(SCORE_DIR, paste0(cancer, "_mirtarbase_merged.csv")),
                         header = TRUE)
  score_data$MGI <- paste0(score_data$miRNA, ":", score_data$gene)
  MGIList <- setNames(score_data$strength, score_data$MGI)
  MGIList <- sort(MGIList, decreasing = TRUE, na.last = NA)

  fg <- fgsea::fgseaMultilevel(pathways = pathway_list, stats = MGIList,
                                scoreType = "neg", minSize = 1, maxSize = length(MGIList),
                                eps = 0)

  if (cancer == "BLCA") {
    # Reconstruct the classic enrichplot::gseaplot(by="runningScore") layout by hand:
    # enrichplot's gseaplot.gseaResult() just calls enrichit::gseaScores() (the same
    # exponent-weighted running-sum statistic fgsea/GSEA have always used, independent of
    # scoreType) and draws the result -- no need to build a full gseaResult S4 object.
    gsdata <- enrichit::gseaScores(MGIList, pathway_list[["miRTarBase"]], exponent = 1, fortify = TRUE)
    gsdata$ymin <- 0; gsdata$ymax <- 0
    pos <- gsdata$position == 1
    h <- diff(range(gsdata$runningScore)) / 20
    gsdata$ymin[pos] <- -h
    gsdata$ymax[pos] <- h

    p1 <- ggplot(gsdata, aes(x = x)) +
      enrichplot::theme_dose() +
      xlab("Position in the Ranked List of Genes") +
      geom_linerange(aes(ymin = ymin, ymax = ymax), color = "black") +
      geom_line(aes(y = runningScore), color = "green", linewidth = 1) +
      geom_vline(data = data.frame(es = which.min(abs(gsdata$runningScore - fg$ES[1]))),
                 aes(xintercept = es), colour = "#FA5860", linetype = "dashed") +
      ylab("Running Enrichment Score") +
      geom_hline(yintercept = 0) +
      ggtitle("Rank-based Functional MTI Enrichment in BLCA") +
      theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1)
  }

  leading_edge_n <- length(fg$leadingEdge[[1]])
  gsea_rows[[i]] <- data.table(
    cancer = cancer,
    n_MGI = length(MGIList),
    NES = fg$NES,
    p_value = fg$pval,
    set_size = fg$size,
    leading_edge_n = leading_edge_n,
    leading_edge_ratio = leading_edge_n / fg$size
  )
}

gsea_table <- rbindlist(gsea_rows)
# As in the original code, apply BH correction to the 17 cancer-level p-values.
gsea_table[, p_adjusted := p.adjust(p_value, method = "BH")]
gsea_table <- gsea_table[, .(cancer, n_MGI, leading_edge_ratio, NES, p_value, p_adjusted)]
fwrite(gsea_table, file.path(BASE, "SupplementaryTable_GSEA.csv"))


# ---------- ORA: threshold-based fold enrichment, computed from raw input ----------
# For each cancer's full candidate (miRNA,gene) universe (SCORE_DIR, already used
# above for GSEA), flag membership in the Functional-MTI pair set, rank pairs by
# -strength, and test for enrichment of Functional-MTI pairs among the top-X%
# strongest pairs at each percentile threshold (one-sided Fisher exact test).
sig_stars <- function(p) fcase(p < 0.001, "***", p < 0.01, "**", p < 0.05, "*", default = "ns")

enrich_one <- function(d, thresh) {
  n <- nrow(d)
  n_top <- max(1L, round(n * thresh))
  in_top <- d$pct_rank <= thresh
  n_tier_total <- sum(d$is_functional_mti)
  n_tier_top <- sum(d$is_functional_mti[in_top])
  baseline_rate <- n_tier_total / n
  top_rate <- n_tier_top / n_top
  fold <- if (baseline_rate > 0) top_rate / baseline_rate else NA_real_
  ft <- matrix(c(n_tier_top, n_top - n_tier_top,
                 n_tier_total - n_tier_top, n - n_top - (n_tier_total - n_tier_top)), nrow = 2)
  p <- tryCatch(fisher.test(ft, alternative = "greater")$p.value, error = function(e) NA_real_)
  data.table(n = n, n_top = n_top, n_tier_total = n_tier_total, n_tier_top = n_tier_top,
             baseline_rate = baseline_rate, top_rate = top_rate, fold_enrichment = fold, p = p)
}

per_cancer_rows <- list()
for (cancer in cancers) {
  d <- fread(file.path(SCORE_DIR, paste0(cancer, "_mirtarbase_merged.csv")),
             select = c("miRNA", "gene", "strength"))
  d <- d[is.finite(strength)]
  d <- merge(d, functional_mti_pairs[, .(miRNA, gene, is_functional_mti = TRUE)],
             by = c("miRNA", "gene"), all.x = TRUE)
  d[is.na(is_functional_mti), is_functional_mti := FALSE]
  d[, pct_rank := rank(strength, ties.method = "average") / .N]

  for (th in THRESHOLDS) {
    res <- enrich_one(d, th)
    per_cancer_rows[[length(per_cancer_rows) + 1]] <- cbind(cancer = cancer, tier = TIER_LABEL, threshold = th, res)
  }
}

per_cancer <- rbindlist(per_cancer_rows)
per_cancer[, threshold_label := factor(as.character(threshold * 100), levels = THRESH_LABELS)]
per_cancer[, sig_label := sig_stars(p)]
fwrite(per_cancer, file.path(RESULTS, "mirtarbase_functional_mti_per_cancer.csv"))
cat("ORA: wrote", nrow(per_cancer), "rows (", length(cancers), "cancers x", length(THRESHOLDS), "thresholds ) to mirtarbase_functional_mti_per_cancer.csv\n")

## ---- SupplementaryTable_ORA: minimal per-cancer results only (no pooled stats) ----
ora_table <- per_cancer[, .(cancer, percentile = as.numeric(as.character(threshold_label)),
                             n_MGI = n, n_top, fold_enrichment, p, sig_label)]
setorder(ora_table, cancer, percentile)
if (nrow(ora_table) != length(cancers) * length(THRESH_LABELS) ||
    !setequal(ora_table$cancer, cancers) ||
    !setequal(ora_table$percentile, as.numeric(THRESH_LABELS))) {
  stop("ORA output is not the expected complete 17 x 8 table; check input CSVs.")
}
fwrite(ora_table, file.path(BASE, "SupplementaryTable_ORA.csv"))

# ---------- Figure S14: original panels A and C only ----------
per_cancer[, threshold_label := factor(threshold_label, levels = THRESH_LABELS)]
blca <- per_cancer[cancer == "BLCA"]
setorder(blca, threshold_label)
y_max_blca <- max(blca$fold_enrichment, na.rm = TRUE)

p2 <- ggplot(blca, aes(x = threshold_label, y = fold_enrichment,
                       fill = sig_label != "ns")) +
  geom_col(color = "black", alpha = 0.6, width = 0.6) +
  geom_text(aes(label = sig_label,
                y = fold_enrichment + 0.06 * y_max_blca), size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c(`TRUE` = "#F2918C", `FALSE` = "grey70"),
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Most Negative MGI Score Percentile (%)",
       y = "Fold Enrichment",
       title = "Threshold-based Functional MTI Enrichment in BLCA") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))

p_all <- cowplot::plot_grid(
  p1, p2,
  nrow = 1, ncol = 2,
  labels = c("A", "B"),
  align = "hv", axis = "tblr",
  rel_widths = c(1, 1)
)
ggsave(filename = file.path(BASE, "FigureS14.pdf"), plot = p_all,
       width = 12, height = 5, units = "in")

message("Saved FigureS14.pdf, SupplementaryTable_GSEA.csv, ",
        "and SupplementaryTable_ORA.csv to: ", BASE)
