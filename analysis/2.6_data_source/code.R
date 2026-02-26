# This is the script to summarize results for methods based on miRNA data and methods based on both miRNA and gene data (only compare TPR)
# Generates Figure S10: TPR_heatmap.pdf

setwd("/scratch/project_2011179/code/miREA/") # change your own directory here
result_dir <- "analysis/2.6_data_source/"
data_dir <- paste0(result_dir, "miR_only_result/") # Here have results for Edge-ORA and Edge-Network, which use DEmiRs and all their target genes as the input MGIs (original is DEMGI)
# Therefore, it only use miRNA expression data, and the miRNA-gene interaction background
TP_dir <- "analysis/2.1_positive_benchmark/"

# Methods used for visualization
method_order2 <- c("TG_ORA", "MiR_ORA", "MiR_Score", "Edge_ORA (DEmiR-TG)", "Edge_Network (DEmiR-TG)",
                   "Edge_ORA", "Edge_Network")
mir_methods <- c("TG_ORA", "MiR_ORA", "MiR_Score", "Edge_ORA (DEmiR-TG)", "Edge_Network (DEmiR-TG)")
both_methods <- c("Edge_ORA", "Edge_Network")

fill_col = c("TG_ORA" = "#97D7F2", "MiR_ORA" = "#B3D49D", "MiR_Score" = "#35B257", "Edge_ORA (DEmiR-TG)" = "#EDC19470", "Edge_Network (DEmiR-TG)" = "#FA807270",
             "TG_Score" = "#07AEE3", "Edge_ORA" = "#EDC194", "Edge_Score" = "#F09137", "Edge_2Ddist" = "#9AA4D6",
             "Edge_Topology" = "#EAA5C2","Edge_Network" = "#FA8072")

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)


# generate data ----
## Edge-ORA and Edge-Network that based on miRNAs-TGs
path_name = "TP"
cat("Start processing", path_name, "...\n")

base_dir <- paste0(data_dir, path_name, "/")
cancer_dirs <- list.dirs(base_dir, recursive = FALSE)
df_list <- list()

for (cancer_dir in cancer_dirs) {
  cancer <- basename(cancer_dir)
  file_path <- file.path(cancer_dir, paste0(cancer, "_", path_name, "_result.csv"))

  if (file.exists(file_path)) {
    df <- read_csv(file_path, col_types = cols())
    df <- df %>% mutate(cancer = cancer)
    df_list[[cancer]] <- df
  } else {
    warning(paste("File not found:", file_path))
  }
}

summary_df <- bind_rows(df_list) %>%
  select(cancer, pathway, method, n_enrich, n_pathway, time, positive_rate)

method_order <- unique(summary_df$method)

summary <- summary_df %>%
  group_by(cancer, method) %>%
  summarize(positive_rate = mean(positive_rate, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = positive_rate) %>%
  select(cancer, all_of(method_order))

write.csv(summary, file = paste0(base_dir, "summary.csv"), row.names = FALSE)
write.csv(summary_df, file = paste0(base_dir, "summary_df.csv"), row.names = FALSE)


F1_TP <- read.csv(paste0(data_dir, "TP/summary_df.csv"), header = TRUE)

method_order = unique(F1_TP$method)
F1_TP <- F1_TP %>%
  group_by(method) %>%
  summarise(TP = sum(n_enrich, na.rm = TRUE),
            n_pos = sum(n_pathway, na.rm = TRUE)) %>%
  ungroup() %>%
  select(method, TP, n_pos) %>%
  mutate(TPR = TP/n_pos,
         method = factor(method, levels = method_order)) %>%
  arrange(method)


F1_TP$method <- paste0(F1_TP$method," (DEmiR-TG)")

df <- read.csv(paste0(TP_dir, "TP_overall_performance.csv"))
df <- rbind(df, F1_TP)

df <- df %>% filter(method %in% method_order2) %>% mutate(method = factor(method, levels = method_order2)) %>% arrange(method)

write.csv(df, file = paste0(result_dir, "TP_overall_performance.csv"), row.names = FALSE)



# plot ----
df <- read.csv(paste0(result_dir, "TP_overall_performance.csv"))
df <- df %>% filter(method %in% method_order2)


df_long <- df %>%
  pivot_longer(cols = TPR, names_to = "Metric", values_to = "Value")

df_long <- df_long %>%
  mutate(
    group = case_when(
      method %in% mir_methods ~ "only miR data",
      method %in% both_methods ~ "both miR and gene data",
      TRUE ~ "other"
    )
  )

df_long$method <- factor(df_long$method, levels = method_order2)

split_x <- length(mir_methods) + 0.5

fmt_label <- function(x){
  ifelse(x == 0, "0",
         ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
}

p <- ggplot(df_long, aes(x = method, y = Value, fill = method)) + # fill = Metric
  geom_col(position = position_dodge(width = 0.75),
           width = 0.7,
           color = "black",
           linewidth = 0.5) +
  geom_text(aes(label = fmt_label(Value), color = Metric),
            position = position_dodge(width = 0.75),
            vjust = -0.4,
            size = 3,
            angle = 0,
            show.legend = FALSE) +
  geom_vline(xintercept = split_x, color = "grey", linewidth = 0.5) +
  scale_fill_manual(values = fill_col) +
  scale_color_manual(values = c(TPR = "black")) +
  annotate("text", x = mean(seq_along(mir_methods)),
           y = max(df_long$Value, na.rm = TRUE) * 1.1,
           label = "only miR data", size = 4, fontface = "bold") +
  annotate("text", x = mean(seq_along(both_methods)) + length(mir_methods),
           y = max(df_long$Value, na.rm = TRUE) * 1.1,
           label = "miR and gene data", size = 4, fontface = "bold") +
  scale_y_continuous(
    limits = c(0, max(df_long$Value, na.rm = TRUE) * 1.2),
    breaks = seq(0, max(df_long$Value, na.rm = TRUE) * 1.2, by = 0.05),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(title = "Overall True Postive Rates for Methods",
       subtitle = "Divided by whether integrating gene data or not",
       x = NULL,
       y = "True Positive Rate",
       fill = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, color = "black", size = 9),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black", angle = 90, size = 12),
    axis.line = element_line(color = "black"),
    axis.ticks.length.x = unit(0.1, "cm"),
    axis.ticks.length.y = unit(0.1, "cm"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.75, fill = NA),
    legend.position = "none",
    legend.background = element_rect(fill = alpha("white", 0.8),
                                     color = "black",
                                     linewidth = 0.5),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    aspect.ratio = 1
  )

ggsave(paste0(result_dir, "TPR_heatmap.pdf"), p, width = 7, height = 7)
