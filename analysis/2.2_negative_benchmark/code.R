# This is the script for summarize negative benchmark results, where the randomized gene/miR/MGI data were used for TP pathways
# FPR = number of pvalue < 0.05 / number of pathways (across all cancers)
# Outputs include:
# 1. FPR_ht.pdf: Figure 2B
# 2. TN_pval_hist.pdf: Figure 2C
# 3. TP_TN_pval_density: Figure S6

setwd("/scratch/project_2011179/code/miREA/") # change your own directory here

library(dplyr)
library(readr) # read_csv()
library(purrr)
library(tidyr) # pivot_wider()
library(ggplot2)
library(ComplexHeatmap)
library(tibble) # rownames_to_column()
# .libPaths(c("/projappl/project_2011179/rpackages_440", .libPaths()))  ## Our special case, please remove this line when using.
library(gghalves)
library(patchwork)
library(scales)
library(coin) # Fisher-Pitman permutation test



result_dir <- "analysis/2.2_negative_benchmark/"
if (!dir.exists(result_dir)){
  dir.create(result_dir)
}
fill_col = c("TG_ORA" = "#97D7F2", "TG_Score" = "#07AEE3", "MiR_ORA" = "#B3D49D", "MiR_Score" = "#35B257",
             "Edge_ORA" = "#EDC194", "Edge_Score" = "#F09137", "Edge_2Ddist" = "#9AA4D6",
             "Edge_Topology" = "#EAA5C2","Edge_Network" = "#FA8072")

methods <- c("TG_ORA", "TG_Score", "MiR_ORA", "MiR_Score", "Edge_ORA", "Edge_Score", "Edge_2Ddist", "Edge_Topology", "Edge_Network")
path_name = "TN"

# 1. generate summary table----
cat("Start processing", path_name, "...\n")
base_dir <- paste0("results/", path_name, "/")
seed_dirs <- list.dirs(base_dir, recursive = FALSE)
df_list <- list()

for (seed_dir in seed_dirs) {
  seed <- basename(seed_dir)
  cancer_dirs <- list.dirs(seed_dir, recursive = FALSE)
  for (cancer_dir in cancer_dirs){
    cancer <- basename(cancer_dir)
    file_path <- file.path(cancer_dir, paste0(cancer, "_", path_name, "_result.csv"))
    if (file.exists(file_path)) {
      df <- read_csv(file_path, col_types = cols())
      df <- df %>% mutate(cancer = cancer, seed = seed)
      df_list[[seed]][[cancer]] <- df
    } else {
      warning(paste("File not found:", file_path))
    }
  }
}

result_list <- df_list %>%
  list_transpose() %>%
  map(~ bind_rows(.x))


result_list <- lapply(result_list, function(x) {
    dplyr::bind_rows(x)
})

summary_df <- bind_rows(result_list) %>%
  select(cancer, seed, pathway, method, n_enrich, n_pathway, time, positive_rate)
method_order <- unique(summary_df$method)
summary_df$method <- factor(summary_df$method, levels = method_order)

summary_by_seed <- summary_df %>%
  group_by(seed, method) %>%
  summarize(FP = sum(n_enrich), n_neg = sum(n_pathway)) %>%
  ungroup() %>%
  mutate(FPR = FP / n_neg)

summary_by_cancer <- summary_df %>%
  group_by(cancer, method) %>%
  summarize(FP = sum(n_enrich), n_neg = sum(n_pathway)) %>%
  ungroup() %>%
  mutate(FPR = FP / n_neg)

summary_by_cancer <- summary_by_cancer %>%
  select(cancer, method, FPR) %>%
  pivot_wider(names_from = method, values_from = FPR)

write.csv(summary_df, file = paste0(result_dir, "TN_summary_df.csv"), row.names = FALSE)
write.csv(summary_by_seed, file = paste0(result_dir, "TN_summary_by_seed.csv"), row.names = FALSE)
write.csv(summary_by_cancer, file = paste0(result_dir, "TN_summary_by_cancer.csv"), row.names = FALSE)

# 2. TN boxplot: each dot is a seed ----
summary <- read.csv(paste0(result_dir, "TN_summary_by_seed.csv"))
# median_df <- summary %>% group_by(method) %>% summarise(med = median(FPR, na.rm = TRUE))
method_order <- unique(summary$method)
summary$method <- factor(summary$method, levels = method_order)
p <- ggplot(summary,
            aes(x = method, y = FPR, fill = method)) + #
  geom_violin(aes(fill = method), scale = "width", alpha = 0.5, color = NA) +
  geom_half_boxplot(position = position_dodge(width = 0.7),
                    side = "r",
               outlier.shape = NA,
               width = 0.6, alpha = 1) +
  geom_point(aes(color = method),
    position = position_nudge(x = -0.2),#position = position_dodge(width = 0.7),
             size = 1.5,
             alpha = 1) +
  # stat_boxplot(geom = "errorbar", width = 0.3, size = 0.75, color = "grey") +
  scale_fill_manual(values = fill_col) + # "#9CC9E870"
  scale_color_manual(values = fill_col)+ # fill_col
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed")+
  annotate("text", x = 0.05, y = 0.05, label = "0.05",
           color = "red", vjust = 1.5, hjust = -0.2, size = 3.5)+
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "grey", linewidth = 1) +
  annotate("text", x = 2,
           y = 0.13,
           label = "node-based", size = 4, fontface = "bold") +
  annotate("text", x = 7,
           y = 0.13,
           label = "edge-based", size = 4, fontface = "bold") +
  labs(
    title = "False Positive Rates Across 16 Cancer Types and 30 Seeds",
    y = "False Positive Rate",
    #x = "Method",
    x = NULL,
    fill = "Method", color = "Method"
  ) +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0,0.15,0.03), expand = c(0,0))+
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.75, fill = NA),
    legend.position = "none",
    legend.background = element_rect(fill = alpha("white", 0.8),
                                     color = "black",
                                     linewidth = 0.5),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    aspect.ratio = 1
  )

ggsave(paste0(result_dir, "FPR_ht.pdf"), p, width = 7, height = 7)

# 3. TN pvalue histogram ----
cat("Start processing", path_name, "...\n")
base_dir <- paste0("results/", path_name, "/")
seed_dirs <- list.dirs(base_dir, recursive = FALSE)
seeds <- basename(seed_dirs)
seeds <- seeds[seeds != "slurm_logs"]
cancer_dirs <- list.dirs(seed_dirs, recursive = FALSE)
cancers <- unique(basename(cancer_dirs))

TN_results <- list()
# seed - method
# build TN_results and TP_results list: wach sub-element is a dataframe, corresponds to one method, for results in all cancers
for (seed in seeds){
  seed <- as.character(seed)
  cat(seed,"\n")

  if (!seed %in% names(TN_results)) {
    TN_results[[seed]] <- list()
  }
  for (cancer in cancers) {
    cat("-", cancer,",")
    # TN
    load(paste0("results/TN/", seed, "/", cancer, "/RData/result.RData"))
    res <- result$result
    rm(result)
    for (method in names(res)) {
      df <- res[[method]]
      df$cancer <- cancer
      df$seed <- seed
      if (!method %in% names(TN_results[[seed]])) {
        TN_results[[seed]][[method]] <- df
      } else {
        TN_results[[seed]][[method]] <- rbind(TN_results[[seed]][[method]], df)
      }
    }
    rm(res)

  }
}

methods <- unique(unlist(lapply(TN_results, names)))
TN_results <- lapply(
  unique(unlist(lapply(TN_results, names))),
  function(m) {
    do.call(rbind,
            lapply(TN_results, function(seed_list) seed_list[[m]]))
  }
)
names(TN_results) <- methods

save(TN_results, file = paste0(result_dir, "result.RData"))

# plot
load(paste0(result_dir, "result.RData"))
methods <- names(TN_results)
fill_col <- fill_col[names(fill_col) %in% methods]

plot_list <- lapply(methods, function(method) {

  df <- data.frame(
    pvalue = TN_results[[method]]$p_value
  )

  ggplot(df, aes(x = pvalue)) +

    geom_histogram(
      bins = 30,
      fill = fill_col[method],
      # color = "black",
      color = "white",
      alpha = 1
    ) +

    theme_minimal(base_size = 12) +

    theme(
      panel.background = element_rect(
        fill = scales::alpha(fill_col[method], 0.1),
        color = NA
      ),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.title = element_blank(),
      axis.text = element_text(size = 10, color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75)
    ) +

    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      labels = label_number(accuracy = 0.01)
    ) +
    scale_y_continuous(
      limits = c(0, 1000),
      breaks = seq(0, 1000, 200)
    ) +

    ggtitle(method) +
    labs(x = "P-values", y = "Density")
})

names(plot_list) <- methods
p_all <- wrap_plots(plot_list, nrow = 3) +
  plot_annotation(
    title = "P-value Histogram Distribution for Negative Benchmark",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    )
  ) &
  labs(x = "p-value", y = "density")

ggsave(paste0(result_dir, "TN_pval_hist.pdf"), p_all, width = 10, height = 10)

# 4. pvalue density TN & TP----
# build TN_results and TP_results list: wach sub-element is a dataframe, corresponds to one method, for results in all cancers
## load data
cancer_dirs <- list.dirs("results/TP/", recursive = FALSE)
cancers <- basename(cancer_dirs)
TP_results <- list()
for (cancer in cancers) {
  cat(cancer,"\n")
  # TP
  load(paste0("results/TP/",cancer, "/RData/result.RData"))
  res <- result$result
  rm(result)
  for (method in names(res)) {
    df <- res[[method]]
    df$cancer <- cancer

    if (!method %in% names(TP_results)) {
      TP_results[[method]] <- df
    } else {
      TP_results[[method]] <- rbind(TP_results[[method]], df)
    }
  }
  rm(res)
}

load(paste0(result_dir, "result.RData")) # TN_results

## plot
format_pval <- function(x) {
  p <- if(length(x) > 1) x else NA_real_
  p_str <- if (p < 0.0001) "<0.0001" else sprintf("p = %.4f", p)
  return(p_str)
}

plot_list <- lapply(methods, function(method) {
  df <- data.frame(
    pvalue = c(TP_results[[method]]$p_value, TN_results[[method]]$p_value),
    type = rep(c("TP", "TN"), times = c(length(TP_results[[method]]$p_value),
                                        length(TN_results[[method]]$p_value)))
  )
  tp_pvalues <- TP_results[[method]]$p_value
  tn_pvalues <- TN_results[[method]]$p_value

  # K-S test
  ks_result <- ks.test(tp_pvalues, tn_pvalues)

  p_val <- ks_result$p.value
  p_text <- ifelse(p_val < 0.0001, "p < 0.0001", paste0("p = ", sprintf("%.4f", p_val)))

  # df$type <- factor(df$type, levels = c("TP", "TN"))
  # test <- coin::oneway_test(pvalue ~ type, data = df,
  #                           distribution = approximate(nresample = 10000)) # 置换 10000 次 approximate(nresample = 10000)
  # p_val <- pvalue(test)
  # p_text <- ifelse(p_val < 0.0001, "p < 0.0001", paste0("p = ", sprintf("%.4f", p_val)))

  ggplot(df, aes(x = pvalue)) + # , color = type, fill = type
    geom_density(aes(y = after_stat(density/max(density)), fill = type, group = type),
                 alpha = 0.3, color = "black", linewidth = 0.5)+
    scale_fill_manual(values = c("TP" = "#99000d", "TN" = "#084594")) + # "TP" = "#EDB2BF", "TN" = "#B2D0ED"
    theme_minimal(base_size = 12) +
    theme(
      #panel.background = element_rect(fill = scales::alpha(fill_col[method], 0.1), color = NA), # use method color as background fill
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.title = element_blank(),
      axis.text = element_text(size = 10, color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75)
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = label_number(accuracy = 0.01)) +
    scale_y_continuous(
      limits = c(0,1),
      breaks = seq(0, 1, 0.2),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggtitle(method) +
    labs(x = "P-values", y = "Density") +
    annotate("text", x = 0.5, y = 0.95, label = p_text, color = "black", size = 4) #fontface = "bold",
})

names(plot_list) <- methods

p_all <- wrap_plots(plot_list, nrow = 3, guides = "collect") +
  plot_annotation(
    title = "P-value Density Distribution for TP and TN",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    )
  ) &
  labs(x = "p-value", y = "density")


ggsave(paste0(result_dir, "TP_TN_pval_density.pdf"), width = 10, height = 10)

# # combined benchmark heatmap: TP + TN ----
# TP_input <- read.csv(paste0("analysis/2.1_positive_benchmark/", "TP", "_data_stat.csv"))
# cancer_order <- TP_input %>% arrange(desc(n_path)) %>% pull(cancer)
#
# TP_summary <- read.csv(paste0("analysis/2.1_positive_benchmark/", "TP", "_summary.csv"), header = TRUE)
# rownames(TP_summary) <- TP_summary$cancer
# TP_summary <- TP_summary[, -1]
# TP_matrix <- as.matrix(t(TP_summary))
# TP_matrix <- TP_matrix[,cancer_order]
#
# TN_summary <- read.csv(paste0(result_dir, "TN_summary_by_cancer.csv"), header = TRUE)
# rownames(TN_summary) <- TN_summary$cancer
# TN_summary <- TN_summary[, -1]
# TN_matrix <- as.matrix(t(TN_summary))
# TN_matrix <- TN_matrix[, cancer_order]
#
# TP_overall <- read.csv(paste0("analysis/2.1_positive_benchmark/", "TP", "_overall_performance.csv"))
# TN_overall <- read.csv(paste0(result_dir, "TN_summary_by_seed.csv"))
#
# benchmark_ht <- plot_benchmark_combined(TP_matrix, TN_matrix, TP_input, TP_overall, TN_overall, plot_path = paste0(result_dir, "/combined_benchmark.pdf"))
#

