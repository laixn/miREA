# This is the script for summarizing whether the methods could identify specific cancer-related pathways in cancer-related pathways for all cancer types.
# Generates Figure S7: cancer_specific_ranks.pdf
setwd("/scratch/project_2011179/code/miREA/") # change your own directory here

library(dplyr)
library(readr) # read_csv()
library(purrr)
library(ggplot2)
library(gghalves)
library(rstatix)
library(ggpubr)
library(patchwork)
#.libPaths(c("/projappl/project_2011179/rpackages_440", .libPaths()))  ## Our special case, please remove this line when using.


result_dir <- "analysis/2.3_distinguish_ability/"
if (!dir.exists(result_dir)){
  dir.create(result_dir)
}
fill_col = c("TG_ORA" = "#97D7F2", "TG_Score" = "#07AEE3", "MiR_ORA" = "#B3D49D", "MiR_Score" = "#35B257",
             "Edge_ORA" = "#EDC194", "Edge_Score" = "#F09137", "Edge_2Ddist" = "#9AA4D6",
             "Edge_Topology" = "#EAA5C2","Edge_Network" = "#FA8072")

method_order <- c("TG_ORA", "TG_Score", "MiR_ORA", "MiR_Score", "Edge_ORA", "Edge_Score", "Edge_2Ddist", "Edge_Topology", "Edge_Network")

# path_name = "TN"
# 1. TP pathways vs other cancers' TP pathways ----
## 1.0 data preparation
TP_all_pw <- read.csv("data/raw_data/pathway/TP/TP_pathway_all_cancers.csv")
TP_all_pw <- TP_all_pw %>% select(pathway, cancer) %>% distinct()

path_name = "TP_all"
cat("Start processing", path_name, "...\n")
base_dir <- paste0("results/", path_name, "/")
cancer_dirs <- list.dirs(base_dir, recursive = FALSE)
cancer_dirs <- cancer_dirs[!grepl("slurm_logs", cancer_dirs)]

df_list <- setNames(vector("list", length(method_order)), method_order)
for (cancer_dir in cancer_dirs) {
  cancer <- basename(cancer_dir)
  message("Processing cancer: ", cancer)
  cancer_pw <- TP_all_pw %>% filter(cancer == !!cancer) %>% pull(pathway)
  for (method in method_order){
    file_path <- file.path(cancer_dir, paste0("result/", cancer, "_", method, ".csv"))
    if (file.exists(file_path)){
      df <- read_csv(file_path, col_types = cols())
      df <- df %>% mutate(rank = rank(p_value, ties.method = "min", na.last = "keep"), # average
                          cancer = cancer,
                          method = method) %>%
        mutate(relative_rank = rank/nrow(df)) %>%
        filter(pathway %in% cancer_pw)
      df_list[[method]] <- bind_rows(df_list[[method]], df)
    } else {
      warning(paste("File not found:", file_path))
    }
  }
}


for (method in names(df_list)) {

  out_file <- file.path(result_dir, paste0(method, "_rank_result.csv"))

  write_csv(df_list[[method]], out_file)

  message("Saved: ", out_file)
}

result <- map_df(names(df_list), function(method){

  df_list[[method]] %>%
    group_by(cancer) %>%
    reframe(
      method = method,
      median_rank = median(relative_rank, na.rm = TRUE),
      min_rank    = min(relative_rank, na.rm = TRUE),
      max_rank    = max(relative_rank, na.rm = TRUE),
      mean_rank   = mean(relative_rank, na.rm = TRUE)
      # n_pw        = n()
    ) %>%
    select(method, cancer, everything()) %>% #
    distinct()
})

write.csv(result, file = paste0(result_dir, "TP_all_result.csv"), row.names = FALSE)

all_result <- list()

for (cancer_dir in cancer_dirs) {
  cancer <- basename(cancer_dir)
  message("Processing cancer: ", cancer)

  cancer_pw <- read.csv(paste0("data/raw_data/pathway/TP/cancer_specific/TP_", cancer, "_gene.csv"), header = TRUE)
  cancer_pw <- unique(cancer_pw$pathway)

  for (method in method_order){
    file_path <- file.path(cancer_dir, paste0("result/", cancer, "_", method, ".csv"))

    if (file.exists(file_path)){
      df <- read_csv(file_path, col_types = cols())

      df <- df %>%
        arrange(p_value) %>%   # 确保顺序
        mutate(
          rank = rank(p_value, ties.method = "min", na.last = "keep"),
          relative_rank = rank / nrow(df),
          cancer = cancer,
          method = method,
          is_cancer_pw = pathway %in% cancer_pw
        )

      all_result[[paste(cancer, method, sep = "_")]] <- df

    } else {
      warning(paste("File not found:", file_path))
    }
  }
}
all_result_df <- bind_rows(all_result)
method_order <- unique(all_result_df$method)
all_result_df$method <- factor(all_result_df$method, levels = method_order)
all_result_df$is_cancer_pw <- factor(all_result_df$is_cancer_pw,
                                     levels = c(TRUE, FALSE))

round_up_to <- function(x, step = 0.1) {
  if (x > 0){
    ceiling(x / step) * step
  } else {
    floor(x / step) * step
  }

}

## 1.1 median relative rank ----
result <- read.csv(paste0(result_dir, "TP_all_result.csv"))
# rank_range <- c(round_up_to(min(result$median_rank)), round_up_to(max(result$median_rank)))
# rank_breaks <- seq(from = ceiling(rank_range[1]), to = floor(rank_range[2]), length.out = 5)

pr_long <- result
pr_long$Methods <- factor(pr_long$method, levels = method_order)
pr_long$method_numeric <- as.numeric(pr_long$Methods) # convert method to numeric positions(same as heatmap rows)


method_range <- c(0.5, length(method_order)+0.5)
ratio <-diff(method_range)
p_median_rank <- ggplot(pr_long, aes(y = median_rank, x = method_numeric)) +
  geom_violin(aes(fill = Methods), scale = "width", alpha = 0.4, color = NA) +
  geom_half_boxplot(aes(fill = Methods), side = "r", outlier.shape = NA, alpha = 1) +
  geom_point(aes(x = method_numeric - 0.2, color = Methods), size = 1) +
  labs(y = "Median Relative Rank", x = NULL,
       title = "Median Relative Ranks for Specific Cancer-Related Pathways"
       ) + # subtitle = "(16 Cancer Types)"
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "grey", linewidth = 1) +
  scale_fill_manual(values = fill_col) +
  scale_color_manual(values = fill_col) +
  scale_x_continuous(
    breaks = 1:length(method_order),
    labels = method_order,
    limits = c(0.5, length(method_order)+0.5),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    limits = c(0,1),
    expand = c(0,0)
  ) +
  theme_minimal()+
  theme(
    legend.position = "none",
    # axis.text.x = element_text(color = "black", size = 8),
    axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12, margin = margin(t = 5)),
    axis.title.y = element_text(color = "black", size = 12, angle = 90, margin = margin(r = 5)),
    axis.ticks.length.x = unit(0.1, "cm"), # element_blank(),
    axis.ticks.length.y = unit(0.1, "cm"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10.5),
    plot.subtitle = element_text(hjust = 0.5),

    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid = element_blank()
  ) +
  coord_fixed(ratio = ratio, clip = "off")

# ggsave(paste0(result_dir, "TP_all_median_rank.pdf"), p_median_rank, width = 6, height = 6)

## 1.2 relative rank diff ----
pw_stat_test <- all_result_df %>%
  group_by(method) %>%
  wilcox_test(relative_rank ~ is_cancer_pw) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

p_diff <- ggplot(all_result_df,
       aes(x = method, y = relative_rank, fill = is_cancer_pw)) +
  geom_boxplot(position = position_dodge(0.7), outlier.shape = NA) +
  # geom_jitter(aes(color = is_cancer_pw),
  #             position = position_dodge(0.8)) +

  stat_pvalue_manual(
    pw_stat_test,
    label = "p.adj.signif",
    x = "method",
    group = "is_cancer_pw",
    y.position = max(all_result_df$relative_rank, na.rm = TRUE) + 0.05,
    tip.length = 0,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(name = "Specific cancer-related", values = c("TRUE" = "#F4A7B9", "FALSE" = "#9CC9E8")) +
  scale_color_manual(name = "Specific cancer-related",values = c("TRUE" = "#F4A7B9", "FALSE" = "#9CC9E8")) +
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "grey", linewidth = 1) +
  theme_bw() +
  labs(
    title = "Relative Ranks Differences",
    # subtitle = "(16 Cancer Types)",
    y = "Relative Rank",
    x = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.75, fill = NA),
    legend.position = c(0.02, 0.95),       # x, y ∈ [0, 1]
    legend.justification = c(0,1),         # 以右下角为锚点
    legend.background = element_rect(fill = alpha("white", 0.8),
                                     # color = "black",
                                     linewidth = 0.5),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10.5),
    plot.subtitle = element_text(hjust = 0.5),
    aspect.ratio = 1
  )
# ggsave(paste0(result_dir, "median_rank_diff.pdf"), p_diff, width = 6, height = 6)



## 1.3 merge plots ----
p <- (p_median_rank | p_diff) +
  plot_annotation(
    title = "Relative Ranks for Cancer-specific Pathways (16 Cancer Types)",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
    )
  )
ggsave(paste0(result_dir, "cancer_specific_ranks.pdf"), p, width = 12, height = 6)

