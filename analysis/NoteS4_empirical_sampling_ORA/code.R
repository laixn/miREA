# This is the script to generate the self-contained summary figure for TG_ORA_emp
# (miRNet 2.0's empirical/permutation sampling ORA, a faithful re-implementation; see
# Enrichment_TG_based_emp.R for the method code), covering Note S4 of the revised manuscript.
# Generates 00_summary.pdf: 2x2 panels -- Positive benchmark (TPR), Negative benchmark (FPR),
# p-value distribution across TP/TN/hallmark pathway sets, and distinguish ability (relative
# rank of cancer-specific vs. non-specific pathways). No other method's data is shown.
#
# All panels are built from the minimal result tables already prepared under TP/, FP/,
# hallmark/, TP_all/ -- nothing is recomputed here (no permutation re-run).

setwd("/scratch/project_2011179/code/miREA_review_ploscb/miREA/") # change your own directory here
# .libPaths(c("/projappl/project_2011179/rpackages_440", .libPaths()))  ## Our special case, please remove this line when using.

library(dplyr)
library(readr)
library(ggplot2)
library(gghalves)
library(rstatix)
library(ggpubr)
library(patchwork)
library(scales)

result_dir <- "analysis/NoteS4_empirical_sampling_ORA/"

col_main <- "#1B98E0"
col_tp <- "#F4A7B9"      # light pink, the "TPR" color used in analysis/2.1_positive_benchmark/code.R (TPR_ht.pdf)
col_tp_dens <- "#99000d" # matches analysis/2.2_negative_benchmark/code.R's TP color (density panel only)
col_tn <- "#084594"      # matches the same reference's TN color
col_hallmark <- "#1B7837" # dark green, tone-matched to the TP/TN pair above

default.TP_cancer <- c("BLCA", "BRCA", "CESC", "COAD", "ESCA", "KIRC",
                        "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD",
                        "READ", "STAD", "THCA", "UCEC")

theme_panel <- theme_bw() + theme(
  axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
  axis.text.y = element_text(color = "black"),
  axis.title = element_text(color = "black"),
  panel.grid = element_blank(),
  panel.border = element_rect(color = "black", linewidth = 0.75, fill = NA),
  plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
  aspect.ratio = 1
)

## ===== top-left: Positive Benchmark (TPR) ===== ##
tp_summary <- read.csv(paste0(result_dir, "TP/TP_summary.csv"), stringsAsFactors = FALSE)
tpr_overall <- read.csv(paste0(result_dir, "TP/TPR_overall.csv"), stringsAsFactors = FALSE)
tp_summary$cancer <- factor(tp_summary$cancer, levels = default.TP_cancer)

p_TPR <- ggplot(tp_summary, aes(x = cancer, y = positive_rate)) +
  geom_col(fill = col_tp) +
  geom_hline(yintercept = tpr_overall$TPR, linetype = "dashed", color = "black") +
  annotate("text", x = 2, y = tpr_overall$TPR + 0.05,
           label = sprintf("Overall TPR = %.3f", tpr_overall$TPR), size = 3.2, hjust = 0) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Positive Benchmark", x = NULL, y = "True Positive Rate") +
  theme_panel

## ===== top-right: Negative Benchmark (FPR) ===== ##
tn_by_seed <- read.csv(paste0(result_dir, "FP/TN_summary_by_seed.csv"), stringsAsFactors = FALSE)
fpr_overall <- read.csv(paste0(result_dir, "FP/FPR_overall.csv"), stringsAsFactors = FALSE)
tn_by_seed$x <- "TG_ORA_emp"

p_FPR <- ggplot(tn_by_seed, aes(x = x, y = FPR)) +
  geom_violin(fill = col_main, alpha = 0.4, color = NA, width = 0.45) +
  geom_half_boxplot(fill = col_main, side = "r", outlier.shape = NA, width = 0.3) +
  geom_point(aes(x = 0.9), color = col_main, size = 1.5) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  annotate("text", x = 0.78, y = 0.05, label = "0.05", color = "red", vjust = -0.6, size = 3) +
  annotate("text", x = 1, y = 0.02,
           label = sprintf("Overall FPR = %.3f (n = 30 seeds)", fpr_overall$FPR), size = 3.2) +
  scale_x_discrete(expand = expansion(add = 0.9)) +
  scale_y_continuous(limits = c(0, max(tn_by_seed$FPR) * 1.15), expand = c(0, 0)) +
  labs(title = "Negative Benchmark", x = NULL, y = "False Positive Rate") +
  theme_panel +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## ===== bottom-left: p-value density across TP / TN / hallmark ===== ##
TP_pvals <- read_csv(paste0(result_dir, "TP/TP_pvalues.csv"), col_types = cols())$p_value
TN_pvals <- read_csv(paste0(result_dir, "FP/TN_pvalues.csv"), col_types = cols())$p_value
hallmark_pvals <- read_csv(paste0(result_dir, "hallmark/hallmark_pvalues.csv"), col_types = cols())$p_value

cat("Pooled p-values -- TP:", length(TP_pvals), " TN:", length(TN_pvals), " hallmark:", length(hallmark_pvals), "\n")

dens_df <- bind_rows(
  data.frame(pvalue = TP_pvals, type = "TP"),
  data.frame(pvalue = TN_pvals, type = "TN"),
  data.frame(pvalue = hallmark_pvals, type = "hallmark")
)
dens_df$type <- factor(dens_df$type, levels = c("TP", "TN", "hallmark"))

p_dens <- ggplot(dens_df, aes(x = pvalue)) +
  geom_density(aes(y = after_stat(density / max(density)), fill = type, group = type),
               alpha = 0.3, color = "black", linewidth = 0.5) +
  scale_fill_manual(name = NULL, values = c("TP" = col_tp_dens, "TN" = col_tn, "hallmark" = col_hallmark)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = percent_format(accuracy = 1)) +
  labs(title = "P-value Distribution", x = "P-values", y = "Density") +
  theme_panel +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.8), linewidth = 0.3),
        legend.key.size = unit(0.4, "cm"))

## ===== bottom-right: distinguish ability (single method) ===== ##
rank_df <- read.csv(paste0(result_dir, "TP_all/TP_all_full_is_cancer_pw.csv"), stringsAsFactors = FALSE) %>%
  mutate(method_short = "TG_ORA_emp")
rank_df$is_cancer_pw <- factor(rank_df$is_cancer_pw, levels = c("TRUE", "FALSE"))

stat_test <- rank_df %>% mutate(method = method_short) %>%
  group_by(method) %>%
  wilcox_test(relative_rank ~ is_cancer_pw) %>%
  adjust_pvalue(method = "BH") %>% add_significance()
# with a single x category, stat_pvalue_manual's automatic group-dodge bracket lookup
# collapses to a zero-width bracket (no visible line) -- set xmin/xmax explicitly to match
# geom_boxplot's position_dodge(0.7) geometry for 2 groups (TRUE/FALSE) instead.
stat_test$xmin <- 0.825
stat_test$xmax <- 1.175

p_rank <- ggplot(rank_df, aes(x = method_short, y = relative_rank, fill = is_cancer_pw)) +
  geom_boxplot(position = position_dodge(0.7), outlier.shape = NA, width = 0.4) +
  stat_pvalue_manual(
    stat_test, label = "p.adj.signif",
    y.position = max(rank_df$relative_rank, na.rm = TRUE) + 0.05,
    tip.length = 0.02
  ) +
  scale_fill_manual(name = "Specific cancer-related", values = c("TRUE" = "#F4A7B9", "FALSE" = "#9CC9E8")) +
  scale_color_manual(name = "Specific cancer-related", values = c("TRUE" = "#F4A7B9", "FALSE" = "#9CC9E8")) +
  theme_bw() +
  labs(title = "Relative Ranks Differences", y = "Relative Rank", x = NULL) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.75, fill = NA),
    legend.position = c(0.02, 0.95),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = alpha("white", 0.8), linewidth = 0.5),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    aspect.ratio = 1
  )

## ===== assemble: 2x2, square panels ===== ##
final <- (p_TPR | p_FPR) / (p_dens | p_rank) +
  plot_annotation(
    title = "TG_ORA_emp (miRNet 2.0 Empirical Sampling ORA)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15))
  )

out_pdf <- paste0(result_dir, "00_summary.pdf")
out_png <- paste0(result_dir, "00_summary.png")
ggsave(out_pdf, final, width = 11, height = 11)
ggsave(out_png, final, width = 11, height = 11, dpi = 200)

cat("Saved", out_pdf, "and", out_png, "\n")
