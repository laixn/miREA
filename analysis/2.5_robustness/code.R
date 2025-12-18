# This is the script to test methods' robustness to data, containing 2.5.1 for score-based methods, and 2.5.2 for DE-based methods
setwd("/scratch/project_2011179/code/miREA/") # change your own directory here
result_dir <- "analysis/2.5_robustness/"
if (!dir.exists(result_dir)){
  dir.create(result_dir)
}
thres_dir <- paste0(result_dir, "threshold_result/")

TN_summary_df <- read.csv("analysis/2.2_negative_benchmark/TN_summary_df.csv")

score_methods <- c("Edge_Score", "Edge_2Ddist", "Edge_Topology")
padj_thresholds <- c(0.1, 0.05, 0.2) # 0.05,
logFC_thresholds <- c(0, 0.584962500721156, 1)

library(dplyr)
library(ComplexHeatmap)
library(circlize) # colorRamp2()
library(ggplot2)
library(patchwork)


# 2.5.1. correlation among FPR for score-based methods and score statistics for cancers----
## generate data----
cancers <- c("BLCA", "BRCA", "CESC", "COAD", "ESCA", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD", "READ", "STAD", "THCA", "UCEC")
seeds <- c(1234,2345,3456,4567,5678,6789,7890,8901,9012,
           123,234,345,456,567,678,789,890,901,
           12,23,34,45,56,67,78,89,90,1,
           12345, 9999)
score_df <- data.frame(cancer = character(), MGI = character(), value = numeric(), stringsAsFactors = FALSE)
score_summary_df <- data.frame(cancer = character(), seed = integer(), sd = numeric(), range = numeric(), IQR = numeric(), MAD = numeric(), stringsAsFactors = FALSE)

for (cancer in cancers){
  cat(cancer,"\n")
  for (seed in seeds){
    if (seed == seeds[length(seeds)]){
      cat(seed,"\n")
    } else {
      cat(seed,"/")
    }
    load(paste0("data/input_data/TN/", cancer, "_TN_", seed, "_input_data.RData"))
    strength <- input_data_rand$data$Edge_Score
    data <- data.frame(MGI = names(strength), strength = strength, cancer = cancer, seed = seed, stringsAsFactors = FALSE)
    score_sum <- as.numeric(summary(data$strength))
    sd <- sd(data$strength)
    range <- diff(range(data$strength))
    IQR <- IQR(data$strength)
    MAD <- mad(data$strength)

    score_df <- rbind(score_df, data %>% select(cancer, seed, MGI, value = strength))

    score_row <- data.frame(cancer = cancer, seed = seed, sd = sd, range = range, IQR = IQR, MAD = MAD, stringsAsFactors = FALSE)
    score_summary_df <- rbind(score_summary_df, score_row)
    rm(input_data_rand)
  }
}

TN_summary <- TN_summary_df %>% filter(method %in% score_methods)
summary <- left_join(score_summary_df, TN_summary, by = c("cancer", "seed"))
# write.csv(score_summary_df, file = "analysis/2_benchmark/2.3.1_score_dist.csv", row.names = FALSE)
write.csv(summary, file = paste0(result_dir, "score_dist_vs_FPR.csv"), row.names = FALSE)

# ggplot(summary %>% filter(method == "Edge_2Ddist"), aes(x = sd, y = positive_rate)) +
#   geom_point(size = 2) +   # 散点
#   geom_smooth(method = "lm", se = TRUE) +   # 回归线 + 95% CI
#   theme_bw() +
#   labs(x = "SD", y = "Positive Rate", title = "SD vs Positive Rate")

## plot----
# pearson heatmap
summary <- read.csv(paste0(result_dir, "score_dist_vs_FPR.csv"))
x_cols <- c("sd", "range", "IQR", "MAD")
y_cols <- c("Edge_Score", "Edge_2Ddist", "Edge_Topology")

cor_pearson_mat <- matrix(NA, nrow = length(x_cols), ncol = length(y_cols),
                          dimnames = list(x_cols, y_cols))
p_pearson_mat <- matrix(NA, nrow = length(x_cols), ncol = length(y_cols),
                        dimnames = list(x_cols, y_cols))

for (i in seq_along(x_cols)){
  for (j in seq_along(y_cols)){
    df <- summary %>% filter(method == y_cols[j])
    x <- df[[x_cols[[i]]]]
    y <- df$positive_rate
    test <- suppressWarnings(cor.test(x, y, use = "complete.obs", method = "pearson"))
    cor_pearson_mat[i, j] <- test$estimate
    p_pearson_mat[i, j] <- test$p.value
  }
}

# star_mat <- matrix("", nrow(cor_pearson_mat), ncol(cor_pearson_mat),
#                    dimnames = list(rownames(cor_pearson_mat), colnames(cor_pearson_mat)))
#
# star_mat[p_pearson_mat < 0.001] <- "***"
# star_mat[p_pearson_mat >= 0.001 & p_pearson_mat < 0.01] <- "**"
# star_mat[p_pearson_mat >= 0.01 & p_pearson_mat < 0.05] <- "*"

pdf(paste0(result_dir, "1_score_dist_vs_FPR.pdf"), width = 5, height = 5)
ComplexHeatmap::Heatmap(
  cor_pearson_mat,
  name = "Pearson",
  width = ncol(cor_pearson_mat)*unit(2, "cm"),
  height = nrow(cor_pearson_mat)*unit(2, "cm"),
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_names_side = "bottom",
  column_names_rot = 30,
  column_names_gp = gpar(hjust = 1, fontsize = 10),
  row_names_gp = gpar(fontsize = 10),
  border_gp = gpar(col = "black", lwd = 1.5),
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill){
    grid.text(paste0(sprintf("%.4f", cor_pearson_mat[i, j]),"\n(p=", sprintf("%.4f", p_pearson_mat[i, j]), ")"), x, y, gp = gpar(fontsize = 9))
  }
)
dev.off()

# 2.5.2 positive rate change when adjusting thresholds for DE-based methods----
## generate data----
# get summary: results/thresholds_analysis/[method]_summary.csv
ora_summary <- data.frame(cancer = character(), padj_threshold = numeric(), logFC_threshold = numeric(), n_DEMGI = integer(),
                          n_enrich = integer(), n_total = integer(), PR = numeric(), stringsAsFactors = FALSE)
network_summary <- data.frame(cancer = character(), padj_threshold = numeric(), logFC_threshold = numeric(), n_DEMGI = integer(),
                              n_enrich = integer(), n_total = integer(), PR = numeric(), stringsAsFactors = FALSE)

for (cancer in cancers){
  cat(cancer,"\n")
  for (padj in padj_thresholds){
    cat(padj,"\n")
    for (logFC in logFC_thresholds){
      cat(logFC,"\n")
      hallmark_network <- read.csv(paste0(thres_dir, cancer, "/result/hallmark_padj", padj, "_logFC", logFC, "_Network.csv"))
      hallmark_ora <- read.csv(paste0(thres_dir, cancer, "/result/hallmark_padj", padj, "_logFC", logFC, "_ORA.csv"))
      reactome_network <- read.csv(paste0(thres_dir, cancer, "/result/Reactome_padj", padj, "_logFC", logFC, "_Network.csv"))
      reactome_ora <- read.csv(paste0(thres_dir, cancer, "/result/Reactome_padj", padj, "_logFC", logFC, "_ORA.csv"))
      ora_data <- rbind(hallmark_ora, reactome_ora)
      network_data <- rbind(hallmark_network, reactome_network)
      ora_enrich <- sum(ora_data$padj < 0.05)
      network_enrich <- sum(network_data$padj < 0.05)
      n_total = 2394+10
      n_DEMGI <- unique(network_data$n_DEMGI)
      ora_PR <- ora_enrich / n_total
      network_PR <- network_enrich / n_total

      ora_summary <- rbind(ora_summary, data.frame(
        cancer = cancer, padj_threshold = padj, logFC_threshold = logFC, n_DEMGI = n_DEMGI,
        n_enrich = ora_enrich, n_total = n_total, PR = ora_PR, stringsAsFactors = FALSE
      ))

      network_summary <- rbind(network_summary, data.frame(
        cancer = cancer, padj_threshold = padj, logFC_threshold = logFC, n_DEMGI = n_DEMGI,
        n_enrich = network_enrich, n_total = n_total, PR = network_PR, stringsAsFactors = FALSE
      ))

    }
  }
}
# network_summary$n_total <- 2394+10
# ora_summary$n_total <- 2394 + 10
# network_summary$PR <- network_summary$n_enrich / network_summary$n_total
# ora_summary$PR <- ora_summary$n_enrich / ora_summary$n_total

write.csv(network_summary, file = paste0(result_dir, "thres_network_summary.csv"), row.names = FALSE)
write.csv(ora_summary, file = paste0(result_dir, "thres_ora_summary.csv"), row.names = FALSE)

# start from this unless the first run
# network_summary <- read.csv(paste0(result_dir, "thres_network_summary.csv"))
# ora_summary <- read.csv(paste0(result_dir, "thres_ora_summary.csv"))

ora_total_summary <- ora_summary %>% group_by(padj_threshold, logFC_threshold) %>%
  summarise(n_enrich = sum(n_enrich),
            n_total = sum(n_total),
            mean_n_DEMGI = mean(n_DEMGI),
            .groups = "drop") %>%
  mutate(PR = n_enrich/n_total) %>%
  mutate(label = paste0(sprintf("%.4f", PR),"\n(",sprintf("%.0f", mean_n_DEMGI),")"))

network_total_summary <- network_summary %>% group_by(padj_threshold, logFC_threshold) %>%
  summarise(n_enrich = sum(n_enrich),
            n_total = sum(n_total),
            mean_n_DEMGI = mean(n_DEMGI),
            .groups = "drop") %>%
  mutate(PR = n_enrich/n_total) %>%
  mutate(label = paste0(sprintf("%.4f", PR),"\n(",sprintf("%.0f", mean_n_DEMGI),")"))



ora_total_summary$padj_label <- as.character(ora_total_summary$padj_threshold)
ora_total_summary <- ora_total_summary %>% mutate(
  logFC_label = case_when(
    logFC_threshold == 0 ~ "log2(1)",
    logFC_threshold == 1 ~ "log2(2)",
    TRUE ~ "log2(1.5)"
  )
)

network_total_summary$padj_label <- as.character(network_total_summary$padj_threshold)
network_total_summary <- network_total_summary %>% mutate(
  logFC_label = case_when(
    logFC_threshold == 0 ~ "log2(1)",
    logFC_threshold == 1 ~ "log2(2)",
    TRUE ~ "log2(1.5)"
  )
)
write.csv(ora_total_summary, file = paste0(result_dir, "thres_ora_total_summary.csv"), row.names = FALSE)
write.csv(network_total_summary, file = paste0(result_dir, "thres_network_total_summary.csv"), row.names = FALSE)

## plot ----
ora_total_summary <- read.csv(paste0(result_dir, "thres_ora_total_summary.csv"))
network_total_summary <- read.csv(paste0(result_dir, "thres_network_total_summary.csv"))

round_up_to <- function(x, step = 0.1) {
  ceiling(x / step) * step
}

# ora
v <- round_up_to(max(ora_total_summary$PR))
range_vals <- seq(0, v, length.out = 5)
padj_order <- sort(unique(ora_total_summary$padj_threshold))
logFC_order <- sort(unique(ora_total_summary$logFC_threshold), decreasing = TRUE)

p_ora <- ggplot(ora_total_summary, aes(
  x = factor(padj_threshold, levels = padj_order),
  y = factor(logFC_threshold, levels = logFC_order),
  fill = PR
)) +
  geom_tile(color = "black", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradientn(
    colours = c("white", "#EDC194"),
    limits = c(0, v),
    breaks = range_vals,
    labels = round(range_vals, 2),
    name = "Edge-ORA PR"
  ) +
  scale_x_discrete(
    labels = setNames(ora_total_summary$padj_label, ora_total_summary$padj_threshold)
  ) +
  scale_y_discrete(
    labels = setNames(ora_total_summary$logFC_label, ora_total_summary$logFC_threshold)
  ) +
  labs(
    x = "padj threshold",
    y = "logFC threshold",
    title = "Threshold Analysis for Edge-ORA Method"
  ) +
  coord_fixed() +
  theme_minimal() + # base_size = 13
  theme(
    # panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),  # 外黑框
    axis.text.x = element_text(angle = 0, color = "black", size = 9),
    axis.text.y = element_text(size = 9, color = "black"),
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5,margin = margin(b = 2)),
    panel.grid = element_blank(),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8),
    axis.title.x = element_text(color = "black", size = 9),
    axis.title.y = element_text(color = "black", size = 9),
  )
# network
v <- round_up_to(max(network_total_summary$PR))
range_vals <- seq(0, v, length.out = 5)

p_net <- ggplot(network_total_summary, aes(
  x = factor(padj_threshold, levels = padj_order),
  y = factor(logFC_threshold, levels = logFC_order),
  fill = PR
)) +
  geom_tile(color = "black", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradientn(
    colours = c("white", "#FA8072"),
    limits = c(0, v),
    breaks = range_vals,
    labels = round(range_vals, 2),
    name = "Edge-Network PR"
  ) +
  scale_x_discrete(
    labels = setNames(network_total_summary$padj_label, network_total_summary$padj_threshold)
  ) +
  scale_y_discrete(
    labels = setNames(network_total_summary$logFC_label, network_total_summary$logFC_threshold)
  ) +
  labs(
    x = "padj threshold",
    y = "logFC threshold",
    title = "Threshold Analysis for Edge-Network Method"
  ) +
  coord_fixed() +
  theme_minimal() + # base_size = 13
  theme(
    # panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),  # 外黑框
    axis.text.x = element_text(angle = 0, color = "black", size = 9),
    axis.text.y = element_text(size = 9, color = "black"),
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5,margin = margin(b = 2)),
    panel.grid = element_blank(),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8),
    axis.title.x = element_text(color = "black", size = 9),
    axis.title.y = element_text(color = "black", size = 9),
  )

combined <- (p_ora / p_net) +
  plot_layout(widths = c(1, 1), heights = c(1,1)) & #, guides = "collect"
  theme(plot.margin = margin(5,5,5,5))#&
# theme(plot.title = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)))

# 绘图
# print(p)
ggsave(paste0(result_dir, "2_threshold_heatmap.pdf"), combined, width = 5, height = 7)
