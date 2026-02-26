# This is the script to summarize the time
# Outputs include:
# 1. Figure S11: 1_time_test.pdf
# 2. Figure 2E: 2_time_summary.pdf

setwd("/scratch/project_2011179/code/miREA/") # change your own directory here
result_dir <- "analysis/2.7_time_test/"
data_dir <- paste0(result_dir, "time_test_result/")

methods <- c("TG_ORA", "TG_Score", "MiR_ORA", "MiR_Score", "Edge_ORA", "Edge_Score", "Edge_2Ddist", "Edge_Topology", "Edge_Network")
fill_col = c("TG_ORA" = "#97D7F2", "TG_Score" = "#07AEE3", "MiR_ORA" = "#B3D49D", "MiR_Score" = "#35B257",
             "Edge_ORA" = "#EDC194", "Edge_Score" = "#F09137", "Edge_2Ddist" = "#626FB3",
             "Edge_Topology" = "#EAA5C2","Edge_Network" = "#FA8072")

library(dplyr)
library(readr)
library(tibble) # rownames_to_column
library(ggplot2)
library(tidyr)
library(scales) # trans_format
library(patchwork)

# .libPaths(c("/projappl/project_2011179/rpackages_440", .libPaths()))  ## Our special case, please remove this line when using.
library(gghalves)

# 1. test optimal number of cores on BLCA----
summary <- data.frame(cancer = character(), path_name = character(), method = character(), ncores = integer(), time = numeric(), stringsAsFactors = FALSE)
for (path_name in c("Reactome", "hallmark")){
  ncores <- c("1", "4", "8", "16", "32")
  methods <- c("Edge_2Ddist", "Edge_Topology", "Edge_Network")
  for (method in methods){
    base_dir <- paste0(data_dir, path_name, "/", method, "/BLCA/result/")
    for (ncore in ncores) {
      # cancer <- basename(cancer_dir)
      # cancers <- c(cancers, cancer)
      file_path <- file.path(base_dir, paste0("time_", ncore, "_ncores.csv"))

      if (file.exists(file_path)) {
        df <- read.csv(file_path)
        summary <- rbind(df, summary)
      } else {
        warning(paste("File not found:", file_path))
      }
    }
  }
}
write.csv(summary, file = paste0(result_dir, "time_test_BLCA.csv"), row.names = FALSE)

# plot
df <- summary %>%
  mutate(
    ncores = factor(ncores, levels = c(1,4,8,16,32),
                    labels = c("1 core (unparallel)", "4 cores", "8 cores", "16 cores", "32 cores")),
    method = factor(method, levels = c("Edge_2Ddist", "Edge_Topology", "Edge_Network"))
  )

df_full <- df %>%
  complete(
    path_name,
    method,
    ncores,
    fill = list(time = NA)
  )

df_full <- df_full %>%
  mutate(
    time_plot = ifelse(is.na(time), 10^0, time),
    is_failed = is.na(time),
    label = ifelse(is.na(time), "failed", floor(time))
  )

y_breaks_fun <- function(lims) {
  lo <- floor(log10(lims[1]))
  hi <- ceiling(log10(lims[2]))
  10^(lo:hi)
}

core_colors <- c(
  "1 core (unparallel)" = "#fdd0a2",
  "4 cores" = "#a1d99b",
  "8 cores" = "#c994c7",
  "16 cores" = "#ffffb2",
  "32 cores" = "#9ecae1"
)

plot_one <- function(df_sub, subtitle) {

  minor_breaks_fun <- function() {
    unlist(lapply(0:4, function(p) seq(10^p, 10^(p+1), length.out = 6)[-1]))
  }

  ggplot(df_sub, aes(x = method, y = time_plot, fill = ncores)) +
    geom_col(
      data = df_sub,
      position = position_dodge(width = 0.8),
      width = 0.8,
      color = "black"
    ) +
    geom_text(
      data = df_sub,
      aes(label = label),
      position = position_dodge(width = 0.8),
      vjust = -0.3,
      size = 3.5,
      color = "black"
    ) +
    scale_fill_manual(values = core_colors, name = NULL) +
    scale_y_log10(
      limits = c(1, 1e5),
      breaks = 10^(0:5),
      minor_breaks = minor_breaks_fun(),
      labels = trans_format("log10", math_format(10^.x)),
      expand = c(0,0)
    ) +
    annotation_logticks(sides = "l") +
    labs(
      x = NULL,
      y = "Calculation Time (s)",
      title = "Time under Different Number of Cores",
      subtitle = subtitle
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),

      panel.border = element_rect(color = "black", size = 1),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),

      plot.title = element_text(hjust = 0.5, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5, color = "black"),

      axis.text.x = element_text(color = "black"),

      legend.position = "none"
    )
}

p_hallmark <- df_full %>%
  filter(path_name == "hallmark") %>%
  plot_one("(10 Cancer Hallmarks for BLCA)")

p_reactome <- df_full %>%
  filter(path_name == "Reactome") %>%
  plot_one("(2394 Reactome Pathways for BLCA)")

combined <- (p_hallmark / p_reactome) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top",
        legend.justification = "center",
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black", hjust = 0.5))

ggsave(paste0(result_dir, "1_time_test.pdf"), width = 8, height = 10)

# 2. time summary ----
for (path_name in c("hallmark", "Reactome")){
  cat("Start processing", path_name, "...\n")
  base_dir <- paste0("results/", path_name, "/")
  cancer_dirs <- list.dirs(base_dir, recursive = FALSE)
  df_list <- list()

  for (cancer_dir in cancer_dirs) {
    cancer <- basename(cancer_dir)
    file_path <- file.path(cancer_dir, paste0(cancer, "_", path_name, "_result.csv"))

    if (file.exists(file_path)) {
      df <- read_csv(file_path, col_types = cols())
      df <- df %>% mutate(cancer = cancer) %>% filter(method %in% methods)
      df_list[[cancer]] <- df
    } else {
      warning(paste("File not found:", file_path))
    }
  }

  summary_df <- bind_rows(df_list) %>%
    select(cancer, pathway, method, n_enrich, n_pathway, time, positive_rate)

  method_order <- unique(summary_df$method)


  # ---- generate cancer × method time dataframe
  time <- summary_df %>%
    group_by(cancer, method) %>%
    summarize(time = mean(time, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = method, values_from = time) %>%
    select(cancer, all_of(method_order))

  write.csv(time, file = paste0(result_dir, path_name, "_time_summary.csv"), row.names = FALSE)

}

reactome_time <- read.csv(paste0(result_dir, "Reactome_time_summary.csv"))
hallmark_time <- read.csv(paste0(result_dir, "hallmark_time_summary.csv"))

rownames(reactome_time) <- reactome_time$cancer
reactome_time <- reactome_time[, -1]
reactome_mat <- as.matrix(t(reactome_time))

rownames(hallmark_time) <- hallmark_time$cancer
hallmark_time <- hallmark_time[, -1]
hallmark_mat <- as.matrix(t(hallmark_time))




# reactome_plot <-raincloud_plot(reactome_mat, object = "Reactome")
# plot(reactome_plot)

round_up_to <- function(x, step = 0.1) {
  if (x > 1){
    ceiling(x / step) * step
  } else {
    floor(x / step) * step
  }

}

# minor_breaks_fun <- function(time_range) {
#   unlist(lapply(time_range, function(p) seq(10^p, 10^(p+1), length.out = 6)[-1]))
# }

# minor_breaks_fun <- function(time_range) {
#   log_range <- log10(time_range)
#
#   unlist(lapply(seq(floor(log_range[1]), ceiling(log_range[2]) - 1),
#                 function(p)
#                   seq(10^p, 10^(p + 1), length.out = 10)[-1]
#   ))
# }

minor_breaks_fun <- function(time_range) {
  log_min <- floor(log10(time_range[1]))
  log_max <- ceiling(log10(time_range[2]))

  unlist(lapply(log_min:(log_max - 1), function(p)
    seq(10^p, 10^(p + 1), length.out = 10)[-1]
  ))
}

min_time <- min(c(min(reactome_mat), min(hallmark_mat)))
max_time <- max(c(max(reactome_mat), max(hallmark_mat)))
time_range <- c(round_up_to(min_time, 0.1), round_up_to(max_time, 1))

log_min <- floor(log10(time_range[1]))
log_max <- ceiling(log10(time_range[2]))
log_breaks <- 10^(log_min:log_max)
time_range <- c(10^log_min, 10^log_max)
# time_breaks <- seq(from = 0, to = time_range[2], by = 1000)
# log_breaks <- 10^(-5:10)
# log_breaks <- log_breaks[
#   log_breaks >= time_range[1] &
#     log_breaks <= time_range[2]
# ]
minor_breaks <- minor_breaks_fun(time_range)



matrix <- hallmark_mat
pr_long <- as.data.frame(matrix) %>%
  rownames_to_column("row") %>%
  pivot_longer(-row, names_to = "column", values_to = "value") %>%
  drop_na()
colnames(pr_long) <- c("Methods", "cancer", "time")
pr_long$Methods <- factor(pr_long$Methods, levels = rev(rownames(matrix)))
pr_long$method_numeric <- as.numeric(pr_long$Methods)# convert method to numeric positions(same as heatmap rows)
method_range <- c(0, nrow(matrix)+0.5)
ratio <- diff(method_range) /  diff(time_range)

# pr_long$time_raw <- pr_long$time
# pr_long$time <- log10(pr_long$time)

p_pr_hallmark <- ggplot(pr_long, aes(y = time, x = method_numeric)) +
  geom_half_violin(aes(fill = Methods), side = "l", scale = "width", alpha = 0.4, color = NA) +
  geom_half_boxplot(aes(fill = Methods), side = "r", outlier.shape = NA, alpha = 1) +
  geom_point(aes(x = method_numeric - 0.2, color = Methods), size = 1) +
  labs(y = NULL, x = NULL, title = "Cancer hallmark (10 pathways)") +
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "grey", linewidth = 1) +
  scale_fill_manual(values = fill_col) +
  scale_color_manual(values = fill_col) +
  scale_x_continuous(
    breaks = 1:nrow(matrix),
    labels = rev(rownames(matrix)),
    limits = method_range,
    expand = c(0,0)
  ) +
  scale_y_log10(
    limits = time_range,
    breaks = log_breaks,
    minor_breaks = minor_breaks,
    labels = trans_format("log10", math_format(10^.x)),
    expand = c(0,0)
  ) +
  annotation_logticks(sides = "b") +
  coord_flip(clip = "off") +
  theme_minimal()+
  theme(
    legend.position = "none",
    axis.text.x = element_text(color = "black", size = 7.5),
    axis.text.y = element_text(color = "black", size = 8),
    axis.title.x = element_text(color = "black", size = 9, margin = margin(t = 5)),
    # axis.title.y = element_text(color = "black", size = 12, angle = 90, margin = margin(r = 5)),
    axis.ticks.length.x = unit(0.1, "cm"),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.length.y = unit(0.1, "cm"),
    axis.ticks.y = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid = element_blank()
  ) +
  annotate("text",
           x = 1, y = minor_breaks[1],
           label = "Number of CPU cores:\nEdge_2Ddist: 1\nEdge_Topology: 10\nEdge_Network: 8",
           hjust = 0, vjust = 0,
           size = 2.5,
           color = "black")

matrix <- reactome_mat
# pr_long <- na.omit(melt(matrix))
pr_long <- as.data.frame(matrix) %>%
  rownames_to_column("row") %>%
  pivot_longer(-row, names_to = "column", values_to = "value") %>%
  drop_na()
colnames(pr_long) <- c("Methods", "cancer", "time")
pr_long$Methods <- factor(pr_long$Methods, levels = rev(rownames(matrix)))
pr_long$method_numeric <- as.numeric(pr_long$Methods) # convert method to numeric positions(same as heatmap rows)

method_range <- c(0, nrow(matrix)+0.5)
ratio <-diff(method_range) /  diff(time_range)
p_pr_reactome <- ggplot(pr_long, aes(y = time, x = method_numeric)) +

  geom_half_violin(aes(fill = Methods), scale = "width", side = "l", alpha = 0.4, color = NA) +
  geom_half_boxplot(aes(fill = Methods), side = "r", outlier.shape = NA, alpha = 1) +
  geom_point(aes(x = method_numeric - 0.2, color = Methods), size = 1) +
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "grey", linewidth = 1) +
  scale_fill_manual(values = fill_col, name = "Methods") +
  scale_color_manual(values = fill_col, name = "Methods") +
  scale_x_continuous(
    breaks = 1:nrow(matrix),
    labels = rev(rownames(matrix)),
    limits = method_range,
    expand = c(0,0)
  ) +
  scale_y_log10(
    limits = time_range,
    breaks = log_breaks,
    minor_breaks = minor_breaks,
    labels = trans_format("log10", math_format(10^.x)),
    expand = c(0,0)
  ) +
  annotation_logticks(sides = "b") +
  labs(y = "Calculation time (s)", x = NULL, title = "Reactome (2394 pathways)", fill = "Methods", color = "Methods", shape = "Methods") +
  guides(fill = guide_legend(title = "Methods"),
         color = guide_legend(title = "Methods")) +
  coord_flip() +
  theme_minimal()+
  theme(
    axis.text.x = element_text(color = "black", size = 7.5),
    axis.text.y = element_text(color = "black", size = 8),
    axis.title.x = element_text(color = "black", size = 9, margin = margin(t = 5)),
    axis.title.y = element_blank(),
    axis.ticks.length.x = unit(0.1, "cm"),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.length.y = unit(0.1, "cm"),
    axis.ticks.y = element_line(color = "black"),
    legend.position = "none",

    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),

    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid = element_blank()
  ) +
  annotate("text",
           x = 1, y = minor_breaks[1],
           label = "Number of CPU cores:\nEdge_2Ddist: 4\nEdge_Topology: 32\nEdge_Network: 8",
           hjust = 0, vjust = 0,
           size = 2.5,
           color = "black")

p_pr <- p_pr_hallmark / p_pr_reactome +
  plot_annotation(
    title = "Calculation Time Across 17 Cancer Types for Different Methods",
    theme = theme(
      plot.title = element_text(
        hjust = 0.8, face = "bold", size = 10.5
      )
    )
  )
ggsave(paste0(result_dir, "2_time_summary.pdf"), p_pr, width = 14, height = 16, units = "cm") # 2_

