# This is the script to summary performance for different methods
# Generate:
# 1. Table S6: summary.csv
# 2. Figure 3: summary.pdf

setwd("/scratch/project_2011179/code/miREA/") # change your work directory here.
library(dplyr)
library(tidyr)
library(purrr)

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggforce)  # geom_arc_bar 用于绘制饼图

method_order <- c("TG_ORA", "TG_Score", "MiR_ORA", "MiR_Score",
                  "Edge_ORA", "Edge_Score", "Edge_2Ddist", "Edge_Topology", "Edge_Network")
metric_order <- c("TPR", "FPR", "median_rank", "sqrt_CMC", "sqrt_oncoKB", "sqrt_COSMIC", "hallmark_median_time", "Reactome_median_time") # "robustness",
metric_name <- c("TPR\n(Sensitivity)", "FPR\n(Specificity)", "Median Rank\n(Distinguish Ability)",
                 "Geometric Median of CMC\n(miRNA Biological Interpretability)", "Geometric Median of OncoKB\n(Gene Biological Interpretability)",
                 "Geometric Median of COSMIC\n(Gene Biological Interpretability)", "Median Hallmark Time (s)\n(Efficiency)",
                 "Median Reactome Time (s)\n(Efficiency)")
metric_range <- c(1, 1, 1, 1, 1, 1, 3000, 5000)
metric_map <- setNames(metric_name, metric_order)
metric_range <- setNames(metric_range, metric_name)

metric_pos <- c("TPR\n(Sensitivity)", "Geometric Median of CMC\n(miRNA Biological Interpretability)", "Geometric Median of OncoKB\n(Gene Biological Interpretability)",
                "Geometric Median of COSMIC\n(Gene Biological Interpretability)")
metric_neg <- setdiff(metric_name, metric_pos)

# 1. generate performance summary table ----
a_TP <- read.csv("analysis/2.1_positive_benchmark/TP_overall_performance.csv")
a_TP <- a_TP[, c("method", "TPR")]

b_TN_raw <- read.csv("analysis/2.2_negative_benchmark/TN_summary.csv")
b_TN <- b_TN_raw %>% group_by(method) %>%
  summarise(TN = sum(FP), n = sum(n_neg)) %>%
  ungroup() %>%
  mutate(FPR = TN / n) %>%
  select(method, FPR)

c_rank_raw <- read.csv("analysis/2.3_distinguish_ability/TP_all_result.csv")
c_rank <- c_rank_raw %>% group_by(method) %>%
  summarise(median_rank = median(median_rank)) %>%
  ungroup()

d_interpretability_raw <- read.csv("analysis/2.4_cancer_identify_ability/result.csv")
d_interpretability_raw <- na.omit(d_interpretability_raw)
d_interpretability <- d_interpretability_raw %>%
  group_by(method) %>%
  summarise(
    n_core_miR = sum(n_core_miR),
    n_CMC = sum(n_CMC),
    n_CMC_core = sum(n_CMC_core),
    n_core_gene = sum(n_core_gene),
    n_oncoKB = sum(n_oncoKB),
    n_oncoKB_core = sum(n_oncoKB_core),
    n_COSMIC = sum(n_COSMIC),
    n_COSMIC_core = sum(n_COSMIC_core)
  ) %>%
  ungroup() %>%
  mutate(
    prop_inter_core = n_CMC_core / n_core_miR,
    prop_inter_CMC = n_CMC_core / n_CMC,
    prop_oncoKB_inter_core = n_oncoKB_core / n_core_gene,
    prop_inter_oncoKB = n_oncoKB_core / n_oncoKB,
    prop_COSMIC_inter_core = n_COSMIC_core / n_core_gene,
    prop_inter_COSMIC = n_COSMIC_core / n_COSMIC,
    sqrt_CMC = sqrt(prop_inter_core * prop_inter_CMC),
    sqrt_oncoKB = sqrt(prop_oncoKB_inter_core * prop_inter_oncoKB),
    sqrt_COSMIC = sqrt(prop_COSMIC_inter_core * prop_inter_COSMIC)
  ) %>%
  select(method, sqrt_CMC, sqrt_oncoKB, sqrt_COSMIC)

e_robustness <- data.frame(
  method = c("Edge_ORA", "Edge_Score", "Edge_2Ddist", "Edge_Topology", "Edge_Network"),
  robustness = c(FALSE, TRUE, TRUE, TRUE, FALSE)
)

f_time_hallmark_raw <- read.csv("analysis/2.7_time_test/hallmark_time_summary.csv")
f_time_Reactome_raw <- read.csv("analysis/2.7_time_test/Reactome_time_summary.csv")

f_time_hallmark <-f_time_hallmark_raw %>%
  pivot_longer(
    cols = -cancer,
    names_to = "method",
    values_to = "time"
  ) %>%
  group_by(method) %>%
  summarise(
    #median_time = median(time, na.rm = TRUE),
    hallmark_median_time = median(time, na.rm = TRUE),
    .groups = "drop"
  )

f_time_Reactome <-f_time_Reactome_raw %>%
  pivot_longer(
    cols = -cancer,
    names_to = "method",
    values_to = "time"
  ) %>%
  group_by(method) %>%
  summarise(
    #median_time = median(time, na.rm = TRUE),
    Reactome_median_time = median(time, na.rm = TRUE),
    .groups = "drop"
  )




df_list <- list(a_TP, b_TN, c_rank, d_interpretability, f_time_hallmark, f_time_Reactome) # e_robustness,

summary_table <- reduce(
  df_list,
  full_join,
  by = "method"
)

summary_table$method <- factor(summary_table$method, levels = method_order)
summary_table <- summary_table %>% arrange(method)

write.csv(summary_table, file = "analysis/3_methods_performance_summary/summary.csv", row.names = FALSE)

# 2. visualzation ----


summary_table <- read.csv("analysis/0_methods_performance_summary/summary.csv")
plot_data <- summary_table %>%
  pivot_longer(cols = metric_order, names_to = "metric", values_to = "value") %>%
  mutate(
    label = ifelse(!is.na(value), sprintf("%.2f", round(as.numeric(value),2)), "")
  ) %>%
  filter(!is.na(value)) %>%
  rowwise()

plot_data$metric <- metric_map[plot_data$metric]
plot_data <- plot_data %>%
  mutate(start = 0.5*pi,
         end = 0.5*pi + 2*pi*value/metric_range[metric])

plot_data$metric <- factor(plot_data$metric, levels = rev(unique(metric_name)))

plot_data$method <- factor(plot_data$method, levels = unique(method_order))


plot_data <- plot_data %>%
  # mutate(start = 0,
  #        end = 2*pi*value/metric_range[metric]) %>%
  arrange(method, metric)

# num_metric <- ncol(summary_table) - 1
# num_method <- length(summary_table$method)
# ratio <- num_metric / num_method

# plot
p <- ggplot(plot_data) +
  geom_arc_bar(aes(
    y0 = as.numeric(factor(metric)),   # x
    x0 = as.numeric(factor(method)),   # y
    r0 = 0,
    r = 0.4,                            # radius
    start = start,
    end = end,
    #fill = ifelse(type=="quant","steelblue", ifelse(value==TRUE,"forestgreen","firebrick"))
    fill = ifelse(metric %in% metric_pos, "firebrick", "steelblue")
  ), color = "black") +
  geom_text(aes(
    y = as.numeric(factor(metric))+0.1,
    x = as.numeric(factor(method)),
    label = label
  ), size = 3) +
  scale_fill_identity() +
  scale_y_continuous(
    breaks = 1:length(unique(plot_data$metric)),
    labels = rev(unique(metric_name))
  ) +
  scale_x_continuous(
    limits = c(0.6, length(summary_table$method)+0.4),
    breaks = 1:length(summary_table$method),
    labels = unique(plot_data$method)
  ) +
  coord_fixed() + # ratio = ratio, clip = "off"
  theme_bw() +
  labs(x="",y="", title = "Methods Perfomance Metrics") +
  theme(
    axis.text.x = element_text(color = "black", angle=30, hjust=1, size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    panel.grid = element_blank(),
    plot.title = element_text(color = "black", face = "bold", hjust = 0.5, size = 14)
    # panel.background = element_rect(fill = "white"),   # 非透明背景
    # panel.border = element_rect(color = "black", linewidth = 0.75)
  )
ggsave("analysis/3_methods_performance_summary/summary.pdf", p, width = 10, height = 10)
